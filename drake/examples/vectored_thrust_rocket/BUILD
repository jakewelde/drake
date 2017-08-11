# -*- python -*-
# This file contains rules for Bazel; see drake/doc/bazel.rst.

load(
    "//tools:drake.bzl",
    "drake_cc_googletest",
    "drake_cc_library",
    "drake_cc_binary",
)
load("//tools:install_data.bzl", "install_data")
load("//tools:lint.bzl", "add_lint_tests")

drake_cc_binary(
    name = "rocket_run_passive",
    srcs = ["rocket_run_passive.cc"],
    add_test_rule = 1,
    data = [":models"],
    test_rule_args = ["-realtime_factor=0.0"],
    deps = [
        "//drake/common:find_resource",
        "//drake/lcm",
        "//drake/multibody:rigid_body_tree",
        "//drake/multibody/joints",
        "//drake/multibody/parsers",
        "//drake/multibody/rigid_body_plant:drake_visualizer",
        "//drake/multibody/rigid_body_plant:rigid_body_plant",
        "//drake/systems/analysis:simulator",
        "//drake/systems/framework:diagram",
        "@com_github_gflags_gflags//:gflags",
    ],
)

drake_cc_binary(
    name = "rocket_run_lqr",
    srcs = ["rocket_run_lqr.cc"],
    add_test_rule = 1,
    data = [":models"],
    test_rule_args = ["-realtime_factor=0.0"],
    deps = [
        "//drake/common:find_resource",
        "//drake/lcm",
        "//drake/multibody:rigid_body_tree",
        "//drake/multibody/parsers",
        "//drake/multibody/rigid_body_plant:drake_visualizer",
        "//drake/multibody/rigid_body_plant:rigid_body_plant",
        "//drake/systems/analysis:simulator",
        "//drake/systems/controllers:linear_quadratic_regulator",
        "@com_github_gflags_gflags//:gflags",
    ],
)

install_data()

add_lint_tests()