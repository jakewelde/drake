# -*- python -*-
# This file contains rules for Bazel; see drake/doc/bazel.rst.

load("//tools:cpplint.bzl", "cpplint")
load("//tools:drake.bzl", "drake_cc_googletest", "drake_cc_library", "drake_cc_binary")

package(default_visibility = ["//visibility:public"])

drake_cc_library(
    name = "kalman_filter",
    srcs = ["kalman_filter.cc"],
    hdrs = ["kalman_filter.h"],
    deps = [
        ":luenberger_observer",
        "//drake/systems/controllers:linear_quadratic_regulator",
        "//drake/systems/framework",
        "//drake/systems/primitives:linear_system",
    ],
)

drake_cc_library(
    name = "luenberger_observer",
    srcs = ["luenberger_observer.cc"],
    hdrs = ["luenberger_observer.h"],
    deps = [
        "//drake/systems/framework",
    ],
)

# === test/ ===

drake_cc_googletest(
    name = "kalman_filter_test",
    deps = [
        ":kalman_filter",
        "//drake/common:eigen_matrix_compare",
    ],
)

drake_cc_googletest(
    name = "luenberger_observer_test",
    deps = [
        ":luenberger_observer",
        "//drake/common:eigen_matrix_compare",
        "//drake/systems/primitives:linear_system",
    ],
)

# === dev/ ===

drake_cc_binary(
    name = "pose_estimation_test",
    deps = [
        "//drake/common:drake_path",
        "//drake/multibody:rigid_body_tree",
        "//drake/multibody/parsers",
        "//drake/solvers:mathematical_program",
        "//drake/solvers:rotation_constraint",
        "@bot_core_lcmtypes//:lib",
    ],
    srcs = ["dev/pose_estimation_test.cc",
            "dev/rotation.h"],
)

cpplint()
