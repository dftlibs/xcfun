/*
 * XCFun, an arbitrary order exchange-correlation library
 * Copyright (C) 2018 Ulf Ekström and contributors.
 *
 * This file is part of XCFun.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * For information on the complete list of contributors to the
 * XCFun library, see: <https://xcfun.readthedocs.io/>
 */

#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>

#include <XCFun/xcfun.h>

// This example contains calls to C++ interface routines that are needed to
// "talk to" the xcfun library and demonstrates how to use them.

// We will compute the kernel for an unpolarized system using total density and
// the gradient components as the variables. These are linear in the density
// matrix, which helps the code using the results from xcfun.

// we consider only one grid point
const int num_grid_points = 1;

// we will use XC_N_NX_NY_NZ
// N: density
// NX: x-gradient of the density
// NY: y-gradient of the density
// NZ: z-gradient of the density
const int num_density_variables = 4;

// forward declare function
double derivative(xcfun_t * fun,
                  int vector_length,
                  double density[][num_density_variables][num_grid_points]);

int main(int, char **) {
  // print some info and copyright about the library
  // please always include this info in your code
  std::cout << xcfun_splash() << std::endl;

  // create a new functional
  // we need this for interacting with the library
  auto fun = xcfun_new();

  {
    // in this example we use PBE
    std::cout << "Setting up PBE" << std::endl;
    auto ierr = xcfun_set(fun, "pbe", 1.0);
    if (ierr)
      std::cout << "functional name not recognized" << std::endl;
  }

  //----------------------------------------------------------------------------
  // in the first example we compute the XC energy ("order 0 derivative")
  {
    constexpr auto order = 0;
    // XC_CONTRACTED here has nothing to do with contracted basis sets
    // it means we will evaluated in the XC_CONTRACTED mode and internally
    // contract functional derivatives with the density taylor expansion
    // in other words: we will not have to explicitly assemble/contract partial
    // derivatives outside of XCFun
    auto ierr = xcfun_eval_setup(fun, XC_N_NX_NY_NZ, XC_CONTRACTED, order);
    if (ierr)
      std::cout << "xcfun_eval_setup failed" << std::endl;

    constexpr auto vector_length =
        1 << order; // bit shift to get power of two: 2**order
    double density[vector_length][num_density_variables][num_grid_points];

    for (auto i = 0; i < num_grid_points; i++) {
      // we use fantasy values here
      density[0][0][i] = 1.0; //         n
      density[0][1][i] = 2.0; // nabla_x n
      density[0][2][i] = 3.0; // nabla_y n
      density[0][3][i] = 4.0; // nabla_z n
    }

    auto res = derivative(fun, vector_length, density);
    std::cout << "The XC energy density is " << res << std::endl;

    // compare with reference
    auto diff = std::abs(-0.86494159400066051 - res);
    if (diff > 1.0e-6)
      std::cout << "derivatives do not match reference numbers" << std::endl;
  }

  //----------------------------------------------------------------------------
  // now we will compute the first derivatives ('potential')
  // and contract them with the first order densities
  {
    constexpr auto order = 1;
    auto ierr = xcfun_eval_setup(fun, XC_N_NX_NY_NZ, XC_CONTRACTED, order);
    if (ierr)
      std::cout << "xcfun_eval_setup failed" << std::endl;

    constexpr auto vector_length =
        1 << order; // bit shift to get power of two: 2**order
    double density[vector_length][num_density_variables][num_grid_points];

    for (auto i = 0; i < num_grid_points; i++) {
      // we use fantasy values here
      density[0][0][i] = 1.0; //         n zeroth order
      density[0][1][i] = 2.0; // nabla_x n zeroth order
      density[0][2][i] = 3.0; // nabla_y n zeroth order
      density[0][3][i] = 4.0; // nabla_z n zeroth order
      density[1][0][i] = 5.0; //         n first order
      density[1][1][i] = 6.0; // nabla_x n first order
      density[1][2][i] = 7.0; // nabla_y n first order
      density[1][3][i] = 8.0; // nabla_z n first order
    }

    auto res = derivative(fun, vector_length, density);

    // compare with reference
    auto diff = std::abs(-5.1509916226154067 - res);
    if (diff > 1.0e-6)
      std::cout << "derivatives do not match reference numbers" << std::endl;
  }

  //----------------------------------------------------------------------------
  // now we will compute a particular partial derivative
  // within order = 1
  // we do this with a trick: we set the perturbed density for
  // the density variable of interest to 1, and set other perturbed
  // densities to 0
  {
    constexpr auto order = 1;
    auto ierr = xcfun_eval_setup(fun, XC_N_NX_NY_NZ, XC_CONTRACTED, order);
    if (ierr)
      std::cout << "xcfun_eval_setup failed" << std::endl;

    constexpr auto vector_length = 1 << order; // bit shift to get 2**order
    double density[vector_length][num_density_variables][num_grid_points];

    for (auto i = 0; i < num_grid_points; i++) {
      // we use fantasy values here
      density[0][0][i] = 1.0; //         n zeroth order
      density[0][1][i] = 2.0; // nabla_x n zeroth order
      density[0][2][i] = 3.0; // nabla_y n zeroth order
      density[0][3][i] = 4.0; // nabla_z n zeroth order
      density[1][0][i] = 0.0;
      density[1][1][i] = 0.0;
      density[1][2][i] = 1.0; // we differentiate wrt this variable
      density[1][3][i] = 0.0;
    }

    auto res = derivative(fun, vector_length, density);

    // compare with reference
    auto diff = std::abs(-0.013470456737102541 - res);
    if (diff > 1.0e-6)
      std::cout << "derivatives do not match reference numbers" << std::endl;
  }

  //----------------------------------------------------------------------------
  // now we try 2nd order
  {
    constexpr auto order = 2;
    auto ierr = xcfun_eval_setup(fun, XC_N_NX_NY_NZ, XC_CONTRACTED, order);
    if (ierr)
      std::cout << "xcfun_eval_setup failed" << std::endl;

    constexpr auto vector_length = 1 << order; // bit shift to get 2**order
    double density[vector_length][num_density_variables][num_grid_points];

    for (auto i = 0; i < num_grid_points; i++) {
      // we use fantasy values here
      density[0][0][i] = 1.0; // zeroth order
      density[0][1][i] = 2.0; // zeroth order
      density[0][2][i] = 3.0; // zeroth order
      density[0][3][i] = 4.0; // zeroth order
      density[1][0][i] = 5.0; // first order
      density[1][1][i] = 6.0; // first order
      density[1][2][i] = 7.0; // first order
      density[1][3][i] = 8.0; // first order
      density[2][0][i] = 5.0; // first order
      density[2][1][i] = 6.0; // first order
      density[2][2][i] = 7.0; // first order
      density[2][3][i] = 8.0; // first order
      density[3][0][i] = 0.0; // second order
      density[3][1][i] = 0.0; // second order
      density[3][2][i] = 0.0; // second order
      density[3][3][i] = 0.0; // second order
    }

    auto res = derivative(fun, vector_length, density);

    // compare with reference
    auto diff = std::abs(-9.4927931153398468 - res);
    if (diff > 1.0e-6)
      std::cout << "derivatives do not match reference numbers" << std::endl;
  }

  //----------------------------------------------------------------------------
  // now we try 3nd order, contracted with perturbed densities
  {
    constexpr auto order = 3;
    auto ierr = xcfun_eval_setup(fun, XC_N_NX_NY_NZ, XC_CONTRACTED, order);
    if (ierr)
      std::cout << "xcfun_eval_setup failed" << std::endl;

    constexpr auto vector_length =
        1 << order; // bit shift to get power of two: 2**order
    double density[vector_length][num_density_variables][num_grid_points];

    for (auto i = 0; i < num_grid_points; i++) {
      // we use fantasy values here
      density[0][0][i] = 1.0;  // zeroth order
      density[0][1][i] = 2.0;  // zeroth order
      density[0][2][i] = 3.0;  // zeroth order
      density[0][3][i] = 4.0;  // zeroth order
      density[1][0][i] = 5.0;  // first order (1)
      density[1][1][i] = 6.0;  // first order (1)
      density[1][2][i] = 7.0;  // first order (1)
      density[1][3][i] = 8.0;  // first order (1)
      density[2][0][i] = 9.0;  // first order (2)
      density[2][1][i] = 10.0; // first order (2)
      density[2][2][i] = 11.0; // first order (2)
      density[2][3][i] = 12.0; // first order (2)
      density[3][0][i] = 5.0;  // second order (depending on (1) and (2))
      density[3][1][i] = 6.0;  // second order (depending on (1) and (2))
      density[3][2][i] = 7.0;  // second order (depending on (1) and (2))
      density[3][3][i] = 8.0;  // second order (depending on (1) and (2))
      density[4][0][i] = 9.0;  // first order (3)
      density[4][1][i] = 10.0; // first order (3)
      density[4][2][i] = 11.0; // first order (3)
      density[4][3][i] = 12.0; // first order (3)
      density[5][0][i] = 5.0;  // second order (depending on (1) and (3))
      density[5][1][i] = 6.0;  // second order (depending on (1) and (3))
      density[5][2][i] = 7.0;  // second order (depending on (1) and (3))
      density[5][3][i] = 8.0;  // second order (depending on (1) and (3))
      density[6][0][i] = 9.0;  // second order (depending on (2) and (3))
      density[6][1][i] = 10.0; // second order (depending on (2) and (3))
      density[6][2][i] = 11.0; // second order (depending on (2) and (3))
      density[6][3][i] = 12.0; // second order (depending on (2) and (3))
      density[7][0][i] = 0.0;  // third order (depending on (1), (2) and (3))
      density[7][1][i] = 0.0;  // third order (depending on (1), (2) and (3))
      density[7][2][i] = 0.0;  // third order (depending on (1), (2) and (3))
      density[7][3][i] = 0.0;  // third order (depending on (1), (2) and (3))
    }

    auto res = derivative(fun, vector_length, density);

    // compare with reference
    auto diff = std::abs(47.091223089835331 - res);
    if (diff > 1.0e-6)
      std::cout << "derivatives do not match reference numbers" << std::endl;
  }

  //----------------------------------------------------------------------------
  // we are done and can release the memory
  xcfun_delete(fun);
  std::cout << "Kernel test has ended properly!" << std::endl;
  return EXIT_SUCCESS;
}

// computes the derivative and takes care of offsetting
double derivative(xcfun_t * fun,
                  int vector_length,
                  double density[][num_density_variables][num_grid_points]) {
  std::vector<double> inp_array(vector_length * num_density_variables, 0.0);
  std::array<std::vector<double>, num_grid_points> out_array{
      std::vector<double>(vector_length, 0.0)};

  // put the densities into the right places
  // along the input array
  for (auto i = 0; i < num_grid_points; i++) {
    auto n = 0;
    for (auto j = 0; j < num_density_variables; j++) {
      for (auto k = 0; k < vector_length; k++) {
        inp_array[n++] = density[k][j][i];
      }
    }
    xcfun_eval(fun, inp_array.data(), out_array[i].data());
  }

  // The output_array holds a Taylor series expansion
  // and we pick here one particular element out of this array.
  return out_array[0][vector_length - 1];
}
