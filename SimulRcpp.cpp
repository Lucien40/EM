#include <RcppArmadillo.h>

#include <armadillo>
#include <fstream>
#include <functional>
#include <iostream>
#include <string>
#include <tuple>
#include <vector>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

namespace {

// Physical constants
double eps0 = 8.854e-12;
double c = 299792458;
double pi = 3.14159265359;
double mu0 = 1.0 / (pow(eps0 * c, 2));
}  // namespace

std::tuple<double, double> getMaterialProps(int matCode, double R, double Ra) {
  // Takes as input a material code
  // and returns mu, eps, sigma for that material locally

  double sigma(1);
  double eps_r(0);
  double eps(0);
  switch (matCode) {
    case 0:  // vacuum
      sigma = 0;
      eps_r = 1;
      break;

    case 1:  // CARBON FIBRE; neglect polarisation

      // https://lib.dr.iastate.edu/cgi/viewcontent.cgi?article=2710&context=qnde
      sigma = 106;
      eps_r = sqrt(26.6 * 26.6 + 17.2 * 17.2);
      break;

    case 2:  // GLASS FIBRE
             // https://ieeexplore.ieee.org/document/7179196
      sigma = 1e-15;
      // for E glass;
      // https://lib.dr.iastate.edu/cgi/viewcontent.cgi?article=2710&context=qnde
      eps_r = sqrt(4.71 * 4.71 + 0.99 * 0.99);
      break;

    case 3:  // ALUMINUM
      sigma =
          3.8e7;  // https://sciencing.com/aluminum-vs-steel-conductivity-5997828.html
      eps = 82.3e-12;
      eps_r = eps / eps0;
      break;

    default:
      break;
  }

  double Ca = (1 - R * sigma / eps_r) / (1 + R * sigma / eps_r);
  double Cb = Ra / (eps_r + R * sigma);
  return std::make_tuple(Ca, Cb);
}

std::tuple<double, double> getMaterial(double x, double y, double z, double R,
                                       double Ra) {
  double R_rocket_outer = 5 * 0.0254;
  // inches to metres
  double R_rocket_inner = 4.9 * 0.0254;

  double GFRP_start_height = -0.3;
  double GFRP_end_height = 0.3;

  double Aluminum_start_height = 0;
  double Aluminum_end_height = .1;
  int matCode(0);

  double s = sqrt(x * x + y * y + z * z);
  if (R_rocket_inner < s && s < R_rocket_outer) {
    if (GFRP_start_height < z && z < GFRP_end_height) {
      matCode = 2;  // GFRP
    } else {
      matCode = 1;  // Carbon
    }
  } else if (s < R_rocket_inner && Aluminum_start_height < z &&
             z < Aluminum_end_height) {
    matCode = 3;  // ALUMINUM
  } else {
    matCode =
        0;  // VACUUM (also gets registered for PML region, but doesn't matter)
  }
  return getMaterialProps(matCode, R, Ra);
}

bool inFTDT(size_t N, size_t Depth, size_t i) {
  return (i >= Depth) && (i < N - Depth);
}

std::tuple<double, double> getCondx(double x, double depth) {
  double init(1);
  double final(100);
  double a(depth - final / depth);
  double sig(x * (x - a) + init);
  double sigs(sig * mu0 / eps0);
  return std::make_tuple(sig, sigs);
};
std::tuple<double, double> getCondy(double y, double depth) {
  double init(1);
  double final(100);
  double a(depth - final / depth);
  double sig(y * (y - a) + init);
  double sigs(sig * mu0 / eps0);
  return std::make_tuple(sig, sigs);
};
std::tuple<double, double> getCondz(double z, double depth) {
  double init(1);
  double final(100);
  double a(depth - final / depth);
  double sig(z * (z - a) + init);
  double sigs(sig * mu0 / eps0);
  return std::make_tuple(sig, sigs);
};

arma::vec getPMLR(
    size_t N, size_t PMLDepth, double d, double R,
    std::function<std::tuple<double, double>(double, double)> getCond) {
  arma::vec Rvec(N, arma::fill::zeros);
  for (size_t i = 0; i < N; i++) {
    if ((i < PMLDepth)) {
      Rvec(i) =
          exp(-std::get<0>(getCond((PMLDepth - i) * d, PMLDepth * d)) * 2 * R);
    } else if (i >= N - PMLDepth) {
      Rvec(i) = exp(
          -std::get<0>(getCond((N - PMLDepth + i) * d, PMLDepth * d)) * 2 * R);
    }
  }
  return Rvec;
}

arma::vec getPMLRS(
    size_t N, size_t PMLDepth, double d, double Rs,
    std::function<std::tuple<double, double>(double, double)> getCond) {
  arma::vec Rvec(N, arma::fill::zeros);
  for (size_t i = 0; i < N; i++) {
    if ((i < PMLDepth)) {
      Rvec(i) =
          exp(-std::get<1>(getCond((PMLDepth - i) * d, PMLDepth * d)) * Rs);
    } else if (i >= N - PMLDepth) {
      Rvec(i) =
          exp(-std::get<1>(getCond((N - PMLDepth + i) * d, PMLDepth * d)) * Rs);
    }
  }
  return Rvec;
}

arma::vec getPMLC(
    size_t N, size_t PMLDepth, double d, double R,
    std::function<std::tuple<double, double>(double, double)> getCond) {
  arma::vec C(N, arma::fill::zeros);
  for (size_t i = 0; i < N; i++) {
    if ((i < PMLDepth)) {
      double sig(std::get<0>(getCond((PMLDepth - i) * d, PMLDepth * d)));
      C(i) = (1 - exp(-sig * 2 * R)) / (sig * d);
    } else if (i >= N - PMLDepth) {
      double sig(std::get<0>(getCond((N - PMLDepth + i) * d, PMLDepth * d)));
      C(i) = (1 - exp(-sig * 2 * R)) / (sig * d);
    }
  }
  return C;
}

arma::vec getPMLCS(
    size_t N, size_t PMLDepth, double d, double Rs,
    std::function<std::tuple<double, double>(double, double)> getCond) {
  arma::vec C(N, arma::fill::zeros);
  for (size_t i = 0; i < N; i++) {
    if ((i < PMLDepth)) {
      double sigS(std::get<1>(getCond((PMLDepth - i) * d, PMLDepth * d)));
      C(i) = (1 - exp(-sigS * 2 * Rs)) / (sigS * d);
    } else if (i >= N - PMLDepth) {
      double sigS(std::get<1>(getCond((N - PMLDepth + i) * d, PMLDepth * d)));
      C(i) = (1 - exp(-sigS * 2 * Rs)) / (sigS * d);
    }
  }
  return C;
}

// [[Rcpp::export]]
List FTDT() {
  double xLen = 0.5;
  double yLen = 0.5;
  double zLen = 0.5;
  size_t NxInside = 10;
  size_t NyInside = 10;
  size_t NzInside = 10;
  size_t PML_depth = 3;

  size_t xN = NxInside + 2 * PML_depth;
  size_t yN = NyInside + 2 * PML_depth;
  size_t zN = NzInside + 2 * PML_depth;
  size_t totalN = xN * yN * zN;

  double dx = xLen / NxInside;
  double dy = yLen / NyInside;
  double dz = zLen / NzInside;

  // Time parameters
  double t = 0;
  double stab_factor = 1;
  double dt = dx / (c * sqrt(3) * stab_factor);
  size_t num_timesteps = totalN;
  double t_final = dt * num_timesteps;

  // std::cout << t_final;

  arma::field<arma::dcube> Ex(num_timesteps);
  arma::field<arma::dcube> Ey(num_timesteps);
  arma::field<arma::dcube> Ez(num_timesteps);
  arma::field<arma::dcube> Hx(num_timesteps);
  arma::field<arma::dcube> Hy(num_timesteps);
  arma::field<arma::dcube> Hz(num_timesteps);

  arma::field<arma::dcube> Hxy(num_timesteps);
  arma::field<arma::dcube> Hyx(num_timesteps);
  arma::field<arma::dcube> Hzy(num_timesteps);
  arma::field<arma::dcube> Hxz(num_timesteps);
  arma::field<arma::dcube> Hyz(num_timesteps);
  arma::field<arma::dcube> Hzx(num_timesteps);

  arma::field<arma::dcube> Exy(num_timesteps);
  arma::field<arma::dcube> Eyx(num_timesteps);
  arma::field<arma::dcube> Ezy(num_timesteps);
  arma::field<arma::dcube> Exz(num_timesteps);
  arma::field<arma::dcube> Eyz(num_timesteps);
  arma::field<arma::dcube> Ezx(num_timesteps);

  Ex(0) = arma::cube(xN, yN, zN, arma::fill::zeros);
  Ey(0) = arma::cube(xN, yN, zN, arma::fill::zeros);
  Ez(0) = arma::cube(xN, yN, zN, arma::fill::zeros);
  Hx(0) = arma::cube(xN, yN, zN, arma::fill::zeros);
  Hy(0) = arma::cube(xN, yN, zN, arma::fill::zeros);
  Hz(0) = arma::cube(xN, yN, zN, arma::fill::zeros);

  Hxy(0) = arma::cube(xN, yN, zN, arma::fill::zeros);
  Hyx(0) = arma::cube(xN, yN, zN, arma::fill::zeros);
  Hzy(0) = arma::cube(xN, yN, zN, arma::fill::zeros);
  Hxz(0) = arma::cube(xN, yN, zN, arma::fill::zeros);
  Hyz(0) = arma::cube(xN, yN, zN, arma::fill::zeros);
  Hzx(0) = arma::cube(xN, yN, zN, arma::fill::zeros);

  Exy(0) = arma::cube(xN, yN, zN, arma::fill::zeros);
  Eyx(0) = arma::cube(xN, yN, zN, arma::fill::zeros);
  Ezy(0) = arma::cube(xN, yN, zN, arma::fill::zeros);
  Exz(0) = arma::cube(xN, yN, zN, arma::fill::zeros);
  Eyz(0) = arma::cube(xN, yN, zN, arma::fill::zeros);
  Ezx(0) = arma::cube(xN, yN, zN, arma::fill::zeros);

  double R = dt / (2 * eps0);
  double Rs = dt / mu0;
  double Ra = pow(c * dt / dx, 2);
  double Rb = dt / (mu0 * dx);
  double invRb = 1.0 / Rb;
  arma::cube Ca(xN, yN, zN, arma::fill::zeros);
  arma::cube Cb(xN, yN, zN, arma::fill::zeros);

  arma::vec PMLRx = getPMLR(xN, PML_depth, dx, R, getCondx);
  arma::vec PMLRy = getPMLR(yN, PML_depth, dy, R, getCondy);
  arma::vec PMLRz = getPMLR(zN, PML_depth, dz, R, getCondz);

  arma::vec PMLCx = getPMLC(xN, PML_depth, dx, R, getCondx);
  arma::vec PMLCy = getPMLC(yN, PML_depth, dy, R, getCondy);
  arma::vec PMLCz = getPMLC(zN, PML_depth, dz, R, getCondz);

  arma::vec PMLRSx = getPMLRS(xN, PML_depth, dx, Rs, getCondx);
  arma::vec PMLRSy = getPMLRS(yN, PML_depth, dy, Rs, getCondy);
  arma::vec PMLRSz = getPMLRS(zN, PML_depth, dz, Rs, getCondz);

  arma::vec PMLCSx = getPMLCS(xN, PML_depth, dx, Rs, getCondx);
  arma::vec PMLCSy = getPMLCS(yN, PML_depth, dy, Rs, getCondy);
  arma::vec PMLCSz = getPMLCS(zN, PML_depth, dz, Rs, getCondz);

  for (size_t ix = 0; ix < xN; ix++) {
    for (size_t iy = 0; iy < yN; iy++) {
      for (size_t iz = 0; iz < zN; iz++) {
        double x_offset = xLen / 2 + dx * PML_depth;
        double y_offset = yLen / 2 + dy * PML_depth;
        double z_offset = zLen / 2 + dz * PML_depth;
        double x = ix * dx - x_offset;
        double y = iy * dy - y_offset;
        double z = iz * dz - z_offset;
        std::tie(Ca(ix, iy, iz), Cb(ix, iy, iz)) = getMaterial(x, y, z, R, Ra);
      }
    }
  }

  // DEFINE incident wave :
  double inc_frequency = 1e7;
  double inc_omega = 2 * pi * inc_frequency;
  double inc_amplitude = 1e5;
  // units are teslas
  double inc_harmonic = 1;

  for (size_t it = 1; it < num_timesteps; it++) {
    // std::cout << inc_amplitude * sin(inc_omega * inc_harmonic * it * dt)
    // << std::endl;
    Hx(it) = arma::cube(xN, yN, zN, arma::fill::zeros);
    Hy(it) = arma::cube(xN, yN, zN, arma::fill::zeros);
    Hz(it) = arma::cube(xN, yN, zN, arma::fill::zeros);
    Ex(it) = arma::cube(xN, yN, zN, arma::fill::zeros);
    Ey(it) = arma::cube(xN, yN, zN, arma::fill::zeros);
    Ez(it) = arma::cube(xN, yN, zN, arma::fill::zeros);

    Hxy(it) = arma::cube(xN, yN, zN, arma::fill::zeros);
    Hyx(it) = arma::cube(xN, yN, zN, arma::fill::zeros);
    Hzy(it) = arma::cube(xN, yN, zN, arma::fill::zeros);
    Hxz(it) = arma::cube(xN, yN, zN, arma::fill::zeros);
    Hyz(it) = arma::cube(xN, yN, zN, arma::fill::zeros);
    Hzx(it) = arma::cube(xN, yN, zN, arma::fill::zeros);

    Exy(it) = arma::cube(xN, yN, zN, arma::fill::zeros);
    Eyx(it) = arma::cube(xN, yN, zN, arma::fill::zeros);
    Ezy(it) = arma::cube(xN, yN, zN, arma::fill::zeros);
    Exz(it) = arma::cube(xN, yN, zN, arma::fill::zeros);
    Eyz(it) = arma::cube(xN, yN, zN, arma::fill::zeros);
    Ezx(it) = arma::cube(xN, yN, zN, arma::fill::zeros);

    for (size_t ix = 1; ix < xN - 1; ix++) {
      for (size_t iy = 1; iy < yN - 1; iy++) {
        for (size_t iz = 1; iz < zN - 1; iz++) {
          if (inFTDT(xN, PML_depth, ix) && inFTDT(yN, PML_depth, iy) &&
              inFTDT(zN, PML_depth, iz)) {
            Hx(it)(ix, iy, iz) =
                Hx(it - 1)(ix, iy, iz) +
                (Rb * Ey(it - 1)(ix, iy, iz + 1) - Rb * Ey(it - 1)(ix, iy, iz) +
                 Rb * Ez(it - 1)(ix, iy, iz) - Rb * Ez(it - 1)(ix, iy + 1, iz));

            Hy(it)(ix, iy, iz) =
                Hy(it - 1)(ix, iy, iz) +
                (Rb * Ex(it - 1)(ix, iy, iz) - Rb * Ex(it - 1)(ix, iy, iz + 1) +
                 Rb * Ez(it - 1)(ix + 1, iy, iz) - Rb * Ez(it - 1)(ix, iy, iz));

            Hz(it)(ix, iy, iz) =
                Hz(it - 1)(ix, iy, iz) +
                (Rb * Ey(it - 1)(ix, iy, iz) - Rb * Ey(it - 1)(ix + 1, iy, iz) +
                 Rb * Ex(it - 1)(ix, iy + 1, iz) - Rb * Ex(it - 1)(ix, iy, iz));

            Ex(it)(ix, iy, iz) =
                (Ca(ix, iy, iz) * Rb * Ex(it - 1)(ix, iy, iz) +
                 Cb(ix, iy, iz) *
                     (Hz(it)(ix, iy, iz) - Hz(it)(ix, iy - 1, iz) -
                      Hy(it)(ix, iy, iz) + Hy(it)(ix, iy, iz - 1))) *
                invRb;

            if (ix == PML_depth + 1) {
              Ex(it)(ix, iy, iz) +=
                  inc_amplitude * sin(inc_omega * inc_harmonic * it * dt);
            }
            // std::cout << inc_amplitude * sin(inc_omega * inc_harmonic * it *
            // dt)
            //<< std::endl;
            // std::cout << it * dt << std::endl;

            Ey(it)(ix, iy, iz) =
                (Ca(ix, iy, iz) * Rb * Ey(it - 1)(ix, iy, iz) +
                 Cb(ix, iy, iz) *
                     (Hx(it)(ix, iy, iz) - Hx(it)(ix, iy, iz - 1) -
                      Hz(it)(ix, iy, iz) + Hz(it)(ix - 1, iy, iz))) *
                invRb;
            Ez(it)(ix, iy, iz) =
                (Ca(ix, iy, iz) * Rb * Ez(it - 1)(ix, iy, iz) +
                 Cb(ix, iy, iz) *
                     (Hy(it)(ix, iy, iz) - Hy(it)(ix - 1, iy, iz) -
                      Hx(it)(ix, iy, iz) + Hx(it)(ix, iy - 1, iz))) *
                invRb;
          } else {
            Hxy(it)(ix, iy, iz) =
                (PMLRSy(iy)) * Hxy(it - 1)(ix, iy, iz) +
                PMLCSy(iy) *
                    (Ezx(it - 1)(ix, iy, iz) + Ezy(it - 1)(ix, iy, iz) -
                     Ezx(it - 1)(ix, iy - 1, iz) - Ezy(it - 1)(ix, iy - 1, iz));

            Hzy(it)(ix, iy, iz) =
                (PMLRSy(iy)) * Hzy(it - 1)(ix, iy, iz) -
                PMLCSy(iy) *
                    (Exy(it - 1)(ix, iy, iz) + Exz(it - 1)(ix, iy, iz) -
                     Exy(it - 1)(ix, iy - 1, iz) - Exz(it - 1)(ix, iy - 1, iz));
            inc_amplitude *sin(inc_omega * inc_harmonic * it * dt);
            Hzx(it)(ix, iy, iz) =
                (PMLRSx(ix)) * Hzx(it - 1)(ix, iy, iz) +
                PMLCSx(ix) *
                    (Eyx(it - 1)(ix, iy, iz) + Eyz(it - 1)(ix, iy, iz) -
                     Eyx(it - 1)(ix, iy - 1, iz) - Eyz(it - 1)(ix, iy - 1, iz));

            Hyx(it)(ix, iy, iz) =
                (PMLRSx(ix)) * Hyx(it - 1)(ix, iy, iz) -
                PMLCSx(ix) *
                    (Ezx(it - 1)(ix, iy, iz) + Ezy(it - 1)(ix, iy, iz) -
                     Ezx(it - 1)(ix, iy - 1, iz) - Ezy(it - 1)(ix, iy - 1, iz));

            Hyz(it)(ix, iy, iz) =
                (PMLRSz(iz)) * Hyz(it - 1)(ix, iy, iz) +
                PMLCSz(iz) *
                    (Exz(it - 1)(ix, iy, iz) + Exy(it - 1)(ix, iy, iz) -
                     Exz(it - 1)(ix, iy - 1, iz) - Exy(it - 1)(ix, iy - 1, iz));

            Hxz(it)(ix, iy, iz) =
                (PMLRSz(iz)) * Hxz(it - 1)(ix, iy, iz) -
                PMLCSz(iz) *
                    (Eyz(it - 1)(ix, iy, iz) + Eyx(it - 1)(ix, iy, iz) -
                     Eyz(it - 1)(ix, iy - 1, iz) - Eyx(it - 1)(ix, iy - 1, iz));
            Exy(it)(ix, iy, iz) =
                (PMLRy(iy)) * Exy(it - 1)(ix, iy, iz) +
                PMLCy(iy) *
                    (Hzx(it - 1)(ix, iy, iz) + Hzy(it - 1)(ix, iy, iz) -
                     Hzx(it - 1)(ix, iy - 1, iz) - Hzy(it - 1)(ix, iy - 1, iz));

            Ezy(it)(ix, iy, iz) =
                (PMLRy(iy)) * Ezy(it - 1)(ix, iy, iz) -
                PMLCy(iy) *
                    (Hxy(it - 1)(ix, iy, iz) + Hxz(it - 1)(ix, iy, iz) -
                     Hxy(it - 1)(ix, iy - 1, iz) - Hxz(it - 1)(ix, iy - 1, iz));

            Ezx(it)(ix, iy, iz) =
                (PMLRx(ix)) * Ezx(it - 1)(ix, iy, iz) +
                PMLCx(ix) *
                    (Hyx(it - 1)(ix, iy, iz) + Hyz(it - 1)(ix, iy, iz) -
                     Hyx(it - 1)(ix, iy - 1, iz) - Hyz(it - 1)(ix, iy - 1, iz));

            Eyx(it)(ix, iy, iz) =
                (PMLRx(ix)) * Eyx(it - 1)(ix, iy, iz) -
                PMLCx(ix) *
                    (Hzx(it - 1)(ix, iy, iz) + Hzy(it - 1)(ix, iy, iz) -
                     Hzx(it - 1)(ix, iy - 1, iz) - Hzy(it - 1)(ix, iy - 1, iz));

            Eyz(it)(ix, iy, iz) =
                (PMLRz(iz)) * Eyz(it - 1)(ix, iy, iz) +
                PMLCz(iz) *
                    (Hxz(it - 1)(ix, iy, iz) + Hxy(it - 1)(ix, iy, iz) -
                     Hxz(it - 1)(ix, iy - 1, iz) - Hxy(it - 1)(ix, iy - 1, iz));

            Exz(it)(ix, iy, iz) =
                (PMLRz(iz)) * Exz(it - 1)(ix, iy, iz) -
                PMLCz(iz) *
                    (Hyz(it - 1)(ix, iy, iz) + Hyx(it - 1)(ix, iy, iz) -
                     Hyz(it - 1)(ix, iy - 1, iz) - Hyx(it - 1)(ix, iy - 1, iz));

            Ex(it)(ix, iy, iz) = Exy(it)(ix, iy, iz) + Exz(it)(ix, iy, iz);
            Ey(it)(ix, iy, iz) = Eyx(it)(ix, iy, iz) + Eyz(it)(ix, iy, iz);
            Ez(it)(ix, iy, iz) = Ezy(it)(ix, iy, iz) + Ezx(it)(ix, iy, iz);
            Hx(it)(ix, iy, iz) = Hxy(it)(ix, iy, iz) + Hxz(it)(ix, iy, iz);
            Hy(it)(ix, iy, iz) = Hyx(it)(ix, iy, iz) + Hyz(it)(ix, iy, iz);
            Hz(it)(ix, iy, iz) = Hzy(it)(ix, iy, iz) + Hzx(it)(ix, iy, iz);
          }
        }
      }
    }
  }
  return List::create(Named("Ex") = Ex, Named("Ey") = Ey, Named("Ez") = Ez,
                      Named("Hx") = Hx, Named("Hy") = Hy, Named("Hz") = Hz);
}
