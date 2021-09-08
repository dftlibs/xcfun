/*
 * XCFun, an arbitrary order exchange-correlation library
 * Copyright (C) 2019 Ulf Ekström and contributors.
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

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include "XCFun/xcfun.h"
#include "XCFunctional.hpp"

// Not ideal, I think
#include "xcint.hpp"

namespace py = pybind11;
using namespace pybind11::literals;

namespace xcfun {
PYBIND11_MODULE(_xcfun, m) {
  m.doc() = "XCFun";

  py::enum_<xcfun_mode>(m, "xcfun_mode")
      .value("XC_MODE_UNSET", xcfun_mode::XC_MODE_UNSET)
      .value("XC_PARTIAL_DERIVATIVES", xcfun_mode::XC_PARTIAL_DERIVATIVES)
      .value("XC_POTENTIAL", xcfun_mode::XC_POTENTIAL)
      .value("XC_CONTRACTED", xcfun_mode::XC_CONTRACTED)
      .value("XC_NR_MODES", xcfun_mode::XC_NR_MODES)
      .export_values();

  py::enum_<xcfun_vars>(m, "xcfun_vars")
      .value("XC_VARS_UNSET", xcfun_vars::XC_VARS_UNSET)
      .value("XC_A", xcfun_vars::XC_A)
      .value("XC_N", xcfun_vars::XC_N)
      .value("XC_A_B", xcfun_vars::XC_A_B)
      .value("XC_N_S", xcfun_vars::XC_N_S)
      .value("XC_A_GAA", xcfun_vars::XC_A_GAA)
      .value("XC_N_S_GNN_GNS_GSS", xcfun_vars::XC_N_S_GNN_GNS_GSS)
      .value("XC_A_GAA_LAPA", xcfun_vars::XC_A_GAA_LAPA)
      .value("XC_A_GAA_TAUA", xcfun_vars::XC_A_GAA_TAUA)
      .value("XC_N_GNN_LAPN", xcfun_vars::XC_N_GNN_LAPN)
      .value("XC_N_GNN_TAUN", xcfun_vars::XC_N_GNN_TAUN)
      .value("XC_A_B_GAA_GAB_GBB_LAPA_LAPB",
             xcfun_vars::XC_A_B_GAA_GAB_GBB_LAPA_LAPB)
      .value("XC_A_B_GAA_GAB_GBB_TAUA_TAUB",
             xcfun_vars::XC_A_B_GAA_GAB_GBB_TAUA_TAUB)
      .value("XC_N_S_GNN_GNS_GSS_LAPN_LAPS",
             xcfun_vars::XC_N_S_GNN_GNS_GSS_LAPN_LAPS)
      .value("XC_N_S_GNN_GNS_GSS_TAUN_TAUS",
             xcfun_vars::XC_N_S_GNN_GNS_GSS_TAUN_TAUS)
      .value("XC_A_B_GAA_GAB_GBB_LAPA_LAPB_TAUA_TAUB",
             xcfun_vars::XC_A_B_GAA_GAB_GBB_LAPA_LAPB_TAUA_TAUB)
      .value("XC_A_B_GAA_GAB_GBB_LAPA_LAPB_TAUA_TAUB_JPAA_JPBB",
             xcfun_vars::XC_A_B_GAA_GAB_GBB_LAPA_LAPB_TAUA_TAUB_JPAA_JPBB)
      .value("XC_N_S_GNN_GNS_GSS_LAPN_LAPS_TAUN_TAUS",
             xcfun_vars::XC_N_S_GNN_GNS_GSS_LAPN_LAPS_TAUN_TAUS)
      .value("XC_A_AX_AY_AZ", xcfun_vars::XC_A_AX_AY_AZ)
      .value("XC_A_B_AX_AY_AZ_BX_BY_BZ", xcfun_vars::XC_A_B_AX_AY_AZ_BX_BY_BZ)
      .value("XC_N_NX_NY_NZ", xcfun_vars::XC_N_NX_NY_NZ)
      .value("XC_N_S_NX_NY_NZ_SX_SY_SZ", xcfun_vars::XC_N_S_NX_NY_NZ_SX_SY_SZ)
      .value("XC_A_AX_AY_AZ_TAUA", xcfun_vars::XC_A_AX_AY_AZ_TAUA)
      .value("XC_A_B_AX_AY_AZ_BX_BY_BZ_TAUA_TAUB",
             xcfun_vars::XC_A_B_AX_AY_AZ_BX_BY_BZ_TAUA_TAUB)
      .value("XC_N_NX_NY_NZ_TAUN", xcfun_vars::XC_N_NX_NY_NZ_TAUN)
      .value("XC_N_S_NX_NY_NZ_SX_SY_SZ_TAUN_TAUS",
             xcfun_vars::XC_N_S_NX_NY_NZ_SX_SY_SZ_TAUN_TAUS)
      .value("XC_A_2ND_TAYLOR", xcfun_vars::XC_A_2ND_TAYLOR)
      .value("XC_A_B_2ND_TAYLOR", xcfun_vars::XC_A_B_2ND_TAYLOR)
      .value("XC_N_2ND_TAYLOR", xcfun_vars::XC_N_2ND_TAYLOR)
      .value("XC_N_S_2ND_TAYLOR", xcfun_vars::XC_N_S_2ND_TAYLOR)
      .value("XC_NR_VARS", xcfun_vars::XC_NR_VARS)
      .export_values();

  m.def("xcfun_version",
        &xcfun_version,
        "XCFun version",
        py::return_value_policy::copy);
  m.def("xcfun_splash",
        &xcfun_splash,
        "XCFun splash screen",
        py::return_value_policy::copy);
  m.def("xcfun_test", &xcfun_test, "XCFun testing");

  py::class_<XCFunctional>(m, "xc_functional");

  m.def("xcfun_new",
        &xcfun::xcfun_new,
        "Create a new XC functional",
        py::return_value_policy::reference);
  m.def("xcfun_delete", &xcfun::xcfun_delete, "Free XC functional", "fun"_a);
  m.def("xcfun_set",
        &xcfun::xcfun_set,
        "Set a parameter in the XC functional",
        "fun"_a,
        "name"_a,
        "value"_a);
  m.def("xcfun_get",
        &xcfun::xcfun_get,
        "Get a parameter in the XC functional",
        "fun"_a,
        "name"_a,
        "value"_a);
  m.def("xcfun_is_gga",
        &xcfun::xcfun_is_gga,
        "Whether the functional is GGA",
        "fun"_a);
  m.def("xcfun_is_metagga",
        &xcfun::xcfun_is_metagga,
        "Whether the functional ia metaGGA",
        "fun"_a);

  m.def("xcfun_eval_setup",
        [](XCFunctional * fun, xcfun_vars vars, xcfun_mode mode, int order) {
          auto err_code = xcfun::xcfun_eval_setup(fun, vars, mode, order);
          if (err_code != 0)
            throw std::invalid_argument("Invalid options in xcfun_eval_setup " +
                                        std::to_string(err_code));
        },
        "Set up XC functional evaluation",
        "fun"_a,
        "vars"_a,
        "mode"_a,
        "order"_a);
  m.def("xcfun_eval",
        [](XCFunctional * fun,
           py::array_t<double, py::array::c_style | py::array::forcecast> density) {
          auto dens_len = xcfun::xcfun_input_length(fun);
          auto output_len = xcfun::xcfun_output_length(fun);

          auto dens_ndim = density.ndim();
          if (density.shape(dens_ndim - 1) != dens_len) {
            throw std::invalid_argument("Wrong dimension of density argument");
          }
          auto nr_points = density.shape(0);
          auto output =
              py::array_t<double, py::array::c_style | py::array::forcecast>({nr_points, output_len}, nullptr);

          if (dens_ndim == 1) {
            xcfun::xcfun_eval(fun, density.data(), output.mutable_data());
          } else if (dens_ndim == 2) {
            auto output_ndim = output.ndim();
            xcfun::xcfun_eval_vec(fun,
                                  nr_points,
                                  density.data(),
                                  density.shape(dens_ndim - 1),
                                  output.mutable_data(),
                                  output.shape(output_ndim - 1));
          } else {
            throw std::invalid_argument("Wrong shape of density argument");
          }

          return output;
        },
        "Evaluate XC functional",
        "fun"_a,
        "density"_a);
}
} // namespace xcfun
