package(default_visibility = ["//visibility:public"])

SOURCE_FILES_NO_MAIN = glob(
    [
        "src/**/*.cc",
        "src/**/*.h",
        "src/**/*.inl",
    ],
    exclude = glob([
        "src/**/*.test.cc",
        "src/**/*.test.h",
        "src/**/*.perf.cc",
        "src/**/*.perf.h",
        "src/**/*.pybind.cc",
        "src/**/*.pybind.h",
        "src/**/main.cc",
    ]),
)

TEST_FILES = glob(
    [
        "src/**/*.test.cc",
        "src/**/*.test.h",
    ],
)

cc_binary(
    name = "pymatching",
    srcs = SOURCE_FILES_NO_MAIN + glob(["src/**/main.cc"]),
    copts = [
        "-std=c++20",
        "-O3",
        "-DNDEBUG",
    ],
    includes = ["src/"],
    deps = ["@stim//:stim_lib"],
)

cc_library(
    name = "libpymatching",
    srcs = SOURCE_FILES_NO_MAIN,
    copts = [
        "-std=c++20",
        "-O3",
        "-DNDEBUG",
    ],
    includes = ["src/"],
    deps = ["@stim//:stim_lib"],
)

cc_test(
    name = "pymatching_test",
    srcs = SOURCE_FILES_NO_MAIN + TEST_FILES,
    copts = [
        "-std=c++20",
    ],
    data = glob(["testdata/**"]),
    includes = ["src/"],
    deps = [
        "@gtest",
        "@gtest//:gtest_main",
        "@stim//:stim_lib",
    ],
)
