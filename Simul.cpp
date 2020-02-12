#include <RcppArmadillo.h>

#include <armadillo>
#include <iostream>
#include <tuple>
#include <vector>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

std::tuple<double, double> getMaterial(double x, double y, double z, double R,
                                       double Ra) {
  return std::make_tuple(x, y);
}
double getPMLRx(size_t ix) { return 0; }
double getPMLRy(size_t iy) { return 0; }
double getPMLRz(size_t iz) { return 0; }
double getPMLCx(size_t ix) { return 0; }
double getPMLCy(size_t iy) { return 0; }
double getPMLCz(size_t iz) { return 0; }

double getPMLRSx(size_t ix) { return 0; }
double getPMLRSy(size_t iy) { return 0; }
double getPMLRSz(size_t iz) { return 0; }
double getPMLCSx(size_t ix) { return 0; }
double getPMLCSy(size_t iy) { return 0; }
double getPMLCSz(size_t iz) { return 0; }

// [[Rcpp::export]]
List FTDT() {
  double xLen = 0.2;
  double yLen = 0.2;
  double zLen = 0.2;
  size_t NxInside = 5;
  size_t NyInside = 5;
  size_t NzInside = 5;
  size_t PML_depth = 2;

  size_t xN = NxInside + 2 * PML_depth;
  size_t yN = NyInside + 2 * PML_depth;
  size_t zN = NzInside + 2 * PML_depth;
  size_t totalN = xN * yN * zN;

  double dx = xLen / NxInside;
  double dy = yLen / NyInside;
  double dz = zLen / NzInside;

  // Physical constants
  double eps0 = 8.854e-12;
  double c = 299792458;
  double mu0 = 1.0 / (pow(eps0 * c, 2));

  // Time parameters
  double t = 0;
  double stab_factor = 1;
  double dt = dx / (c * sqrt(3) * stab_factor);
  size_t num_timesteps = totalN * 3;
  double t_final = dt * num_timesteps;

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

  double R = dt / (2 * eps0);
  double Ra = pow(c * dt / dx, 2);
  double Rb = dt / (mu0 * dx);
  double invRb = 1.0 / Rb;
  arma::cube Ca(xN, yN, zN, arma::fill::zeros);
  arma::cube Cb(xN, yN, zN, arma::fill::zeros);

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

  for (size_t it = 1; it < num_timesteps; it++) {
    Hx(it) = arma::cube(xN, yN, zN, arma::fill::zeros);
    Hy(it) = arma::cube(xN, yN, zN, arma::fill::zeros);
    Hz(it) = arma::cube(xN, yN, zN, arma::fill::zeros);
    Ex(it) = arma::cube(xN, yN, zN, arma::fill::zeros);
    Ey(it) = arma::cube(xN, yN, zN, arma::fill::zeros);
    Ez(it) = arma::cube(xN, yN, zN, arma::fill::zeros);

    for (size_t ix = 0; ix < xN; ix++) {
      for (size_t iy = 0; iy < yN; iy++) {
        for (size_t iz = 0; iz < zN; iz++) {
          if (ix >= PML_depth && ix < NxInside + PML_depth && iy >= PML_depth &&
              iy < NyInside + PML_depth && iz >= PML_depth &&
              iz < NzInside + PML_depth) {
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
            Ex(it)(ix, iy, iz) =
                (getPMLRy(iy) + getPMLRz(iz)) * Ex(it - 1)(ix, iy, iz) +
                getPMLCy(iy) *
                    (Hzx(it - 1)(ix, iy, iz) + Hzy(it - 1)(ix, iy, iz) -
                     Hzx(it - 1)(ix, iy - 1, iz) -
                     Hzy(it - 1)(ix, iy - 1, iz)) +
                getPMLCz(iz) *
                    (Hyx(it - 1)(ix, iy, iz) + Hyz(it - 1)(ix, iy, iz) -
                     Hyx(it - 1)(ix, iy - 1, iz) - Hyz(it - 1)(ix, iy - 1, iz));

            Ey(it)(ix, iy, iz) =
                (getPMLRx(iy) + getPMLRz(iz)) * Ey(it - 1)(ix, iy, iz) +
                getPMLCz(iy) *
                    (Hxy(it - 1)(ix, iy, iz) + Hxz(it - 1)(ix, iy, iz) -
                     Hxy(it - 1)(ix, iy - 1, iz) -
                     Hxz(it - 1)(ix, iy - 1, iz)) +
                getPMLCx(iz) *
                    (Hzx(it - 1)(ix, iy, iz) + Hzy(it - 1)(ix, iy, iz) -
                     Hzx(it - 1)(ix, iy - 1, iz) - Hzy(it - 1)(ix, iy - 1, iz));

            Ez(it)(ix, iy, iz) =
                (getPMLRy(iy) + getPMLRz(iz)) * Ex(it - 1)(ix, iy, iz) +
                getPMLCy(iy) *
                    (Hxz(it - 1)(ix, iy, iz) + Hxy(it - 1)(ix, iy, iz) -
                     Hxz(it - 1)(ix, iy - 1, iz) -
                     Hxy(it - 1)(ix, iy - 1, iz)) +
                getPMLCx(iz) *
                    (Hyx(it - 1)(ix, iy, iz) + Hyz(it - 1)(ix, iy, iz) -
                     Hyx(it - 1)(ix, iy - 1, iz) - Hyz(it - 1)(ix, iy - 1, iz));

            Hxy(it)(ix, iy, iz) =
                getPMLRSy(iy) * Hxy(it - 1)(ix, iy, iz) +
                getPMLCSy(iy) * (Ez(it)(ix, iy, iz) - Ez(it)(ix, iy + 1, iz));

            Hxz(it)(ix, iy, iz) =
                getPMLRSz(iz) * Hxz(it - 1)(ix, iy, iz) -
                getPMLCSz(iz) * (Ey(it)(ix, iy, iz) - Ey(it)(ix, iy, iz + 1));

            Hzy(it)(ix, iy, iz) =
                getPMLRSy(iy) * Hzy(it - 1)(ix, iy, iz) -
                getPMLCSy(iy) * (Ex(it)(ix, iy, iz) - Ex(it)(ix, iy + 1, iz));

            Hzx(it)(ix, iy, iz) =
                getPMLRSx(ix) * Hzx(it - 1)(ix, iy, iz) +
                getPMLCSx(ix) * (Ey(it)(ix, iy, iz) - Ey(it)(ix + 1, iy, iz));

            Hyx(it)(ix, iy, iz) =
                getPMLRSx(iy) * Hyx(it - 1)(ix, iy, iz) -
                getPMLCSx(iy) * (Ez(it)(ix, iy, iz) - Ez(it)(ix + 1, iy, iz));

            Hyz(it)(ix, iy, iz) =
                getPMLRSz(iy) * Hyz(it - 1)(ix, iy, iz) +
                getPMLCSz(iy) * (Ex(it)(ix, iy, iz) - Ex(it)(ix, iy, iz + 1));
          }
        }
      }
    }
  }
  return List::create(Named("Ex") = Ex, Named("Ey") = Ey, Named("Ez") = Ez,
                      Named("Hx") = Hx, Named("Hy") = Hy, Named("Hz") = Hz);
}
