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
double mu0 = 1.0 / (eps0 * pow(c, 2));
}  // namespace

std::tuple<double, double, double> getMaterialProps(int matCode) {
  // Takes as input a material code
  // and returns mu, eps, sigma for that material locally

  double sigma(0);
  double eps_r(eps0);
  double eps(0);
  double mu(mu0);
  switch (matCode) {
    case 0:  // vacuum
      sigma = 0;
      eps_r = eps0;
      mu = mu0;
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

  return std::make_tuple(eps_r, sigma, mu);
}

int getMat(double x, double y, double z) {
  // Sets material geometry
  double R_rocket_outer = 5 * 0.0254 * 0.5;
  // inches to metres
  double R_rocket_inner = 4.9 * 0.0254 * 0.5;

  double GFRP_start_height = -0.3;
  double GFRP_end_height = 0.3;

  double Aluminum_start_height = 0;
  double Aluminum_end_height = .1;
  int matCode(0);

  double s = sqrt(x * x + y * y);
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
  return matCode;
}

std::tuple<double, double, double> getCondMaterial(double x, double y,
                                                   double z) {
  // get electromagnetic properties based on position
  return getMaterialProps(getMat(x, y, z));
}

bool inFTDT(size_t N, size_t Depth, size_t i) {
  // True if in FDTD region in 1D
  return (i >= Depth) && (i < N - Depth);
}

// PML layer coefficients, quadratically graded

std::tuple<double, double> getCondPMLx(size_t x, size_t depth, double d) {
  double init(0.001);
  double final(1e8);
  double a(depth * d - final / d / depth);
  double sig(x * d * (x * d - a) + init);
  double sigt((x + 0.5) * d * ((x + 0.5) * d - a) + init);
  double sigs(sigt * mu0 / eps0);
  if (sig == 0 || sigs == 0) {
    std::cout << x << std::endl;
  }
  return std::make_tuple(sig, sigs);
};
std::tuple<double, double> getCondPMLy(size_t y, size_t depth, double d) {
  double init(0.001);
  double final(1e8);
  double a(depth * d - final / d / depth);
  double sig(y * d * (y * d - a) + init);
  double sigt((y + 0.5) * d * ((y + 0.5) * d - a) + init);
  double sigs(sigt * mu0 / eps0);
  if (sig == 0 || sigs == 0) {
    std::cout << y << std::endl;
  }
  return std::make_tuple(sig, sigs);
};
std::tuple<double, double> getCondPMLz(size_t z, size_t depth, double d) {
  double init(0.001);
  double final(1e8);
  double a(depth * d - final / d / depth);
  double sig(z * d * (z * d - a) + init);
  double sigt((z + 0.5) * d * ((z + 0.5) * d - a) + init);

  double sigs(sigt * mu0 / eps0);
  if (sig == 0 || sigs == 0) {
    std::cout << z << std::endl;
  }
  return std::make_tuple(sig, sigs);
};

arma::cube getMatMat(size_t Nx, size_t Ny, size_t Nz, double d) {
  // get Material map
  arma::cube Rvec(Nx, Ny, Nz, arma::fill::zeros);
  double sig(0);
  double eps(0);
  for (size_t ix = 0; ix < Nx; ix++) {
    for (size_t iy = 0; iy < Ny; iy++) {
      for (size_t iz = 0; iz < Nz; iz++) {
        Rvec(ix, iy, iz) = getMat((ix - Nx * 0.5) * d, (iy - Ny * 0.5) * d,
                                  (iz - Nz * 0.5) * d);
      }
    }
  }
  return Rvec;
}

std::tuple<double, double, double, double> getCond(size_t Nx, size_t Ny,
                                                   size_t Nz, size_t PMLDepth,
                                                   double d, size_t ix,
                                                   size_t iy, size_t iz,
                                                   char dir) {
  // get electromagnetic properties based on lattice index, with PML layer

  std::tuple<double, double, double, double> cond;

  double sig(1);
  double sigS(1);
  double eps(eps0);
  double mu(mu0);
  if (inFTDT(Nx, PMLDepth, ix) && inFTDT(Ny, PMLDepth, iy) &&
      inFTDT(Nz, PMLDepth, iz)) {
    // In FTDT:
    std::tuple<double, double, double> condFTDT;
    switch (dir) {
      case 'x':
        mu = std::get<2>(getCondMaterial((ix - PMLDepth - Nx * 0.5) * d,
                                         (iy + 0.5 - PMLDepth - Ny * 0.5) * d,
                                         (iz + 0.5 - PMLDepth - Nz * 0.5) * d));
        condFTDT = getCondMaterial((ix + 0.5 - PMLDepth - Nx * 0.5) * d,
                                   (iy - PMLDepth - Ny * 0.5) * d,
                                   (iz - PMLDepth - Nz * 0.5) * d);
        break;
      case 'y':
        mu = std::get<2>(getCondMaterial((ix + 0.5 - PMLDepth - Nx * 0.5) * d,
                                         (iy - PMLDepth - Ny * 0.5) * d,
                                         (iz + 0.5 - PMLDepth - Nz * 0.5) * d));
        condFTDT = getCondMaterial((ix - PMLDepth - Nx * 0.5) * d,
                                   (iy + 0.5 - PMLDepth - Ny * 0.5) * d,
                                   (iz - PMLDepth - Nz * 0.5) * d);
        break;
      case 'z':
        mu = std::get<2>(getCondMaterial((ix + 0.5 - PMLDepth - Nx * 0.5) * d,
                                         (iy + 0.5 - PMLDepth - Ny * 0.5) * d,
                                         (iz - PMLDepth - Nz * 0.5) * d));
        condFTDT = getCondMaterial((ix - PMLDepth - Nx * 0.5) * d,
                                   (iy - PMLDepth - Ny * 0.5) * d,
                                   (iz + 0.5 - PMLDepth - Nz * 0.5) * d);
        break;

      default:
        break;
    }
    // Here the material is centered

    sig = std::get<1>(condFTDT);
    sigS = 0;
    eps = std::get<0>(condFTDT);

  } else {
    // In PML region
    std::tuple<double, double> condPML;

    switch (dir) {
      case 'x':
        if ((ix < PMLDepth)) {
          condPML = getCondPMLx((PMLDepth - ix) * d, PMLDepth, d);

        } else if (ix >= Nx - PMLDepth) {
          condPML = (getCondPMLx((Nx - PMLDepth + ix) * d, PMLDepth, d));
        }
        break;
      case 'y':
        if ((iy < PMLDepth)) {
          condPML = (getCondPMLy((PMLDepth - iy) * d, PMLDepth, d));
        } else if (iy >= Ny - PMLDepth) {
          condPML = (getCondPMLy((Ny - PMLDepth + iy) * d, PMLDepth, d));
        }
        break;
      case 'z':
        if ((iz < PMLDepth)) {
          condPML = (getCondPMLz((PMLDepth - iz) * d, PMLDepth, d));
        } else if (iz >= Nz - PMLDepth) {
          condPML = (getCondPMLz((Nz - PMLDepth + iy) * d, PMLDepth, d));
        }

      default:
        break;
    }

    sig = std::get<0>(condPML);
    sigS = std::get<1>(condPML);
  }
  return std::make_tuple(sig, sigS, eps, mu);
}

arma::cube getR(size_t Nx, size_t Ny, size_t Nz, size_t PMLDepth, double d,
                double dt, char dir) {
  arma::cube Rvec(Nx, Ny, Nz, arma::fill::zeros);
  double sig(0);
  double eps(0);
  for (size_t ix = 0; ix < Nx; ix++) {
    for (size_t iy = 0; iy < Ny; iy++) {
      for (size_t iz = 0; iz < Nz; iz++) {
        auto cond = getCond(Nx, Ny, Nz, PMLDepth, d, ix, iy, iz, dir);
        sig = std::get<0>(cond);
        eps = std::get<2>(cond);

        Rvec(ix, iy, iz) = exp(-sig * dt / eps);
      }
    }
  }
  return Rvec;
}

arma::cube getRs(size_t Nx, size_t Ny, size_t Nz, size_t PMLDepth, double d,
                 double dt, char dir) {
  arma::cube Rvec(Nx, Ny, Nz, arma::fill::zeros);
  double sigS(0);
  double mu(0);
  for (size_t ix = 0; ix < Nx; ix++) {
    for (size_t iy = 0; iy < Ny; iy++) {
      for (size_t iz = 0; iz < Nz; iz++) {
        auto cond = getCond(Nx, Ny, Nz, PMLDepth, d, ix, iy, iz, dir);
        sigS = std::get<1>(cond);
        mu = std::get<3>(cond);

        Rvec(ix, iy, iz) = exp(-sigS * dt / mu);
      }
    }
  }
  return Rvec;
}

arma::cube getC(size_t Nx, size_t Ny, size_t Nz, size_t PMLDepth, double d,
                double dt, char dir) {
  arma::cube C(Nx, Ny, Nz, arma::fill::zeros);
  double sig(0);
  double eps(0);
  for (size_t ix = 0; ix < Nx; ix++) {
    for (size_t iy = 0; iy < Ny; iy++) {
      for (size_t iz = 0; iz < Nz; iz++) {
        auto cond = getCond(Nx, Ny, Nz, PMLDepth, d, ix, iy, iz, dir);
        sig = std::get<0>(cond);
        eps = std::get<2>(cond);
        if ((inFTDT(Nx, PMLDepth, ix) && inFTDT(Ny, PMLDepth, iy) &&
             inFTDT(Nz, PMLDepth, iz)) ||
            sig == 0) {
          C(ix, iy, iz) = dt / (eps * d);

        } else {
          C(ix, iy, iz) = (1 - exp(-sig * dt / eps)) / (sig * d);
          // std::cout << sig << dir << std::endl;
        }
      }
    }
  }
  return C;
}

arma::cube getCs(size_t Nx, size_t Ny, size_t Nz, size_t PMLDepth, double d,
                 double dt, char dir) {
  arma::cube C(Nx, Ny, Nz, arma::fill::zeros);
  double sigS(0);
  double mu(0);
  for (size_t ix = 0; ix < Nx; ix++) {
    for (size_t iy = 0; iy < Ny; iy++) {
      for (size_t iz = 0; iz < Nz; iz++) {
        auto cond = getCond(Nx, Ny, Nz, PMLDepth, d, ix, iy, iz, dir);
        sigS = std::get<1>(cond);
        mu = std::get<3>(cond);

        if ((inFTDT(Nx, PMLDepth, ix) && inFTDT(Ny, PMLDepth, iy) &&
             inFTDT(Nz, PMLDepth, iz)) ||
            sigS == 0) {
          C(ix, iy, iz) = dt / (mu * d);
        } else {
          C(ix, iy, iz) = (1 - exp(-sigS * dt / mu)) / (sigS * d);
        }
      }
    }
  }
  return C;
}

// [[Rcpp::export]]
List FTDT(size_t NxIn, size_t NyIn, size_t NzIn) {
  std::cout << "Eps0: " << eps0 << "\n mu0: " << mu0 << std::endl;
  double xLen = 0.3;
  double yLen = 0.3;
  double zLen = 0.3;
  size_t NxInside = NxIn;
  size_t NyInside = NxIn;
  size_t NzInside = NxIn;
  size_t PML_depth = 5;

  size_t xN = NxInside + 2 * PML_depth;
  size_t yN = NyInside + 2 * PML_depth;
  size_t zN = NzInside + 2 * PML_depth;
  size_t totalN = xN * yN * zN;

  double dx = xLen / NxInside;
  double dy = yLen / NyInside;
  double dz = zLen / NzInside;

  // Time parameters
  double t = 0;
  double stab_factor = 2;
  double dt = dx / (c * sqrt(3) * stab_factor);
  size_t num_timesteps = 500;
  double t_final = dt * num_timesteps;

  // std::cout << t_final;

  arma::field<arma::dcube> Ex(num_timesteps);
  arma::field<arma::dcube> Ey(num_timesteps);
  arma::field<arma::dcube> Ez(num_timesteps);
  arma::field<arma::dcube> Hx(num_timesteps);
  arma::field<arma::dcube> Hy(num_timesteps);
  arma::field<arma::dcube> Hz(num_timesteps);

  arma::field<arma::dcube> Hxy(2);
  arma::field<arma::dcube> Hyx(2);
  arma::field<arma::dcube> Hzy(2);
  arma::field<arma::dcube> Hxz(2);
  arma::field<arma::dcube> Hyz(2);
  arma::field<arma::dcube> Hzx(2);

  arma::field<arma::dcube> Exy(2);
  arma::field<arma::dcube> Eyx(2);
  arma::field<arma::dcube> Ezy(2);
  arma::field<arma::dcube> Exz(2);
  arma::field<arma::dcube> Eyz(2);
  arma::field<arma::dcube> Ezx(2);

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

  // double R = dt / (2 * eps0);
  // double Rs = dt / mu0;
  // double Ra = pow(c * dt / dx, 2);
  // double Rb = dt / (mu0 * dx);
  // double invRb = 1.0 / Rb;
  // arma::cube Ca(xN, yN, zN, arma::fill::zeros);
  // arma::cube Cb(xN, yN, zN, arma::fill::zeros);

  // set integrating factors

  arma::cube Rx = getR(xN, yN, zN, PML_depth, dx, dt, 'y');
  arma::cube Ry = getR(xN, yN, zN, PML_depth, dx, dt, 'x');
  arma::cube Rz = getR(xN, yN, zN, PML_depth, dx, dt, 'z');

  arma::cube Cx = getC(xN, yN, zN, PML_depth, dx, dt, 'x');
  arma::cube Cy = getC(xN, yN, zN, PML_depth, dx, dt, 'y');
  arma::cube Cz = getC(xN, yN, zN, PML_depth, dx, dt, 'z');

  arma::cube RSx = getRs(xN, yN, zN, PML_depth, dx, dt, 'x');
  arma::cube RSy = getRs(xN, yN, zN, PML_depth, dx, dt, 'y');
  arma::cube RSz = getRs(xN, yN, zN, PML_depth, dx, dt, 'z');

  arma::cube CSx = getCs(xN, yN, zN, PML_depth, dx, dt, 'x');
  arma::cube CSy = getCs(xN, yN, zN, PML_depth, dx, dt, 'y');
  arma::cube CSz = getCs(xN, yN, zN, PML_depth, dx, dt, 'z');

  // DEFINE incident wave :
  double inc_frequency = c / (10 * dx);
  double inc_omega = 2 * pi * inc_frequency;
  double inc_amplitude = 1;
  // units are teslas
  double inc_harmonic = 1;

  std::cout << "f: " << inc_frequency << std::endl;
  std::cout << "dt: " << dt << std::endl;

  for (size_t it = 1; it < num_timesteps; it++) {
    // std::cout << inc_amplitude * sin(inc_omega * inc_harmonic * it * dt)
    //          << std::endl;
    Hx(it) = arma::cube(xN, yN, zN, arma::fill::zeros);
    Hy(it) = arma::cube(xN, yN, zN, arma::fill::zeros);
    Hz(it) = arma::cube(xN, yN, zN, arma::fill::zeros);
    Ex(it) = arma::cube(xN, yN, zN, arma::fill::zeros);
    Ey(it) = arma::cube(xN, yN, zN, arma::fill::zeros);
    Ez(it) = arma::cube(xN, yN, zN, arma::fill::zeros);

    Hxy(1) = arma::cube(xN, yN, zN, arma::fill::zeros);
    Hyx(1) = arma::cube(xN, yN, zN, arma::fill::zeros);
    Hzy(1) = arma::cube(xN, yN, zN, arma::fill::zeros);
    Hxz(1) = arma::cube(xN, yN, zN, arma::fill::zeros);
    Hyz(1) = arma::cube(xN, yN, zN, arma::fill::zeros);
    Hzx(1) = arma::cube(xN, yN, zN, arma::fill::zeros);

    Exy(1) = arma::cube(xN, yN, zN, arma::fill::zeros);
    Eyx(1) = arma::cube(xN, yN, zN, arma::fill::zeros);
    Ezy(1) = arma::cube(xN, yN, zN, arma::fill::zeros);
    Exz(1) = arma::cube(xN, yN, zN, arma::fill::zeros);
    Eyz(1) = arma::cube(xN, yN, zN, arma::fill::zeros);
    Ezx(1) = arma::cube(xN, yN, zN, arma::fill::zeros);

    for (size_t ix = 1; ix < xN - 1; ix++) {
      for (size_t iy = 1; iy < yN - 1; iy++) {
        for (size_t iz = 1; iz < zN - 1; iz++) {
          Hxy(1)(ix, iy, iz) =
              (RSy(ix, iy, iz)) * Hxy(0)(ix, iy, iz) +
              CSy(ix, iy, iz) *
                  (Ezx(0)(ix, iy, iz) + Ezy(0)(ix, iy, iz) -
                   Ezx(0)(ix, iy - 1, iz) - Ezy(0)(ix, iy - 1, iz));

          Hzy(1)(ix, iy, iz) =
              (RSy(ix, iy, iz)) * Hzy(0)(ix, iy, iz) -
              CSy(ix, iy, iz) *
                  (Exy(0)(ix, iy, iz) + Exz(0)(ix, iy, iz) -
                   Exy(0)(ix, iy - 1, iz) - Exz(0)(ix, iy - 1, iz));

          Hzx(1)(ix, iy, iz) =
              (RSx(ix, iy, iz)) * Hzx(0)(ix, iy, iz) +
              CSx(ix, iy, iz) *
                  (Eyx(0)(ix, iy, iz) + Eyz(0)(ix, iy, iz) -
                   Eyx(0)(ix - 1, iy, iz) - Eyz(0)(ix - 1, iy, iz));

          Hyx(1)(ix, iy, iz) =
              (RSx(ix, iy, iz)) * Hyx(0)(ix, iy, iz) -
              CSx(ix, iy, iz) *
                  (Ezx(0)(ix, iy, iz) + Ezy(0)(ix, iy, iz) -
                   Ezx(0)(ix - 1, iy, iz) - Ezy(0)(ix - 1, iy, iz));

          Hyz(1)(ix, iy, iz) =
              (RSz(ix, iy, iz)) * Hyz(0)(ix, iy, iz) +
              CSz(ix, iy, iz) *
                  (Exz(0)(ix, iy, iz) + Exy(0)(ix, iy, iz) -
                   Exz(0)(ix, iy, iz - 1) - Exy(0)(ix, iy, iz - 1));

          Hxz(1)(ix, iy, iz) =
              (RSz(ix, iy, iz)) * Hxz(0)(ix, iy, iz) -
              CSz(ix, iy, iz) *
                  (Eyz(0)(ix, iy, iz) + Eyx(0)(ix, iy, iz) -
                   Eyz(0)(ix, iy, iz - 1) - Eyx(0)(ix, iy, iz - 1));

          Exy(1)(ix, iy, iz) =
              (Ry(ix, iy, iz)) * Exy(0)(ix, iy, iz) +
              Cy(ix, iy, iz) *
                  (Hzx(0)(ix, iy, iz) + Hzy(0)(ix, iy, iz) -
                   Hzx(0)(ix, iy - 1, iz) - Hzy(0)(ix, iy - 1, iz));

          Ezy(1)(ix, iy, iz) =
              (Ry(ix, iy, iz)) * Ezy(0)(ix, iy, iz) -
              Cy(ix, iy, iz) *
                  (Hxy(0)(ix, iy, iz) + Hxz(0)(ix, iy, iz) -
                   Hxy(0)(ix, iy - 1, iz) - Hxz(0)(ix, iy - 1, iz));

          Ezx(1)(ix, iy, iz) =
              (Rx(ix, iy, iz)) * Ezx(0)(ix, iy, iz) +
              Cx(ix, iy, iz) *
                  (Hyx(0)(ix, iy, iz) + Hyz(0)(ix, iy, iz) -
                   Hyx(0)(ix - 1, iy, iz) - Hyz(0)(ix - 1, iy, iz));

          Eyx(1)(ix, iy, iz) =
              (Rx(ix, iy, iz)) * Eyx(0)(ix, iy, iz) -
              Cx(ix, iy, iz) *
                  (Hzx(0)(ix, iy, iz) + Hzy(0)(ix, iy, iz) -
                   Hzx(0)(ix - 1, iy, iz) - Hzy(0)(ix - 1, iy, iz));

          Eyz(1)(ix, iy, iz) =
              (Rz(ix, iy, iz)) * Eyz(0)(ix, iy, iz) +
              Cz(ix, iy, iz) *
                  (Hxz(0)(ix, iy, iz) + Hxy(0)(ix, iy, iz) -
                   Hxz(0)(ix, iy, iz - 1) - Hxy(0)(ix, iy, iz - 1));

          Exz(1)(ix, iy, iz) =
              (Rz(ix, iy, iz)) * Exz(0)(ix, iy, iz) -
              Cz(ix, iy, iz) *
                  (Hyz(0)(ix, iy, iz) + Hyx(0)(ix, iy, iz) -
                   Hyz(0)(ix, iy, iz - 1) - Hyx(0)(ix, iy, iz - 1));

          if (ix == PML_depth + 1) {
            Ezx(1)(ix, iy, iz) +=
                inc_amplitude * sin(inc_omega * inc_harmonic * it * dt);
            Ezy(1)(ix, iy, iz) +=
                inc_amplitude * sin(inc_omega * inc_harmonic * it * dt);
          }
          Ex(it)(ix, iy, iz) = Exy(1)(ix, iy, iz) + Exz(1)(ix, iy, iz);
          Ey(it)(ix, iy, iz) = Eyx(1)(ix, iy, iz) + Eyz(1)(ix, iy, iz);
          Ez(it)(ix, iy, iz) = Ezy(1)(ix, iy, iz) + Ezx(1)(ix, iy, iz);
          Hx(it)(ix, iy, iz) = Hxy(1)(ix, iy, iz) + Hxz(1)(ix, iy, iz);
          Hy(it)(ix, iy, iz) = Hyx(1)(ix, iy, iz) + Hyz(1)(ix, iy, iz);
          Hz(it)(ix, iy, iz) = Hzy(1)(ix, iy, iz) + Hzx(1)(ix, iy, iz);
        }
      }
    }
    Hxy(0) = Hxy(1);
    Hyx(0) = Hyx(1);
    Hzy(0) = Hzy(1);
    Hxz(0) = Hxz(1);
    Hyz(0) = Hyz(1);
    Hzx(0) = Hzx(1);

    Exy(0) = Exy(1);
    Eyx(0) = Eyx(1);
    Ezy(0) = Ezy(1);
    Exz(0) = Exz(1);
    Eyz(0) = Eyz(1);
    Ezx(0) = Ezx(1);
  }
  return List::create(Named("Ex") = wrap(Ex), Named("Ey") = wrap(Ey),
                      Named("Ez") = wrap(Ez), Named("Hx") = Hx,
                      Named("Hy") = Hy, Named("Hz") = Hz, Named("Cs") = CSx,
                      Named("C") = Cx, Named("Rs") = RSx, Named("R") = Rx,
                      Named("mat") = getMatMat(xN, yN, zN, dx));
}
