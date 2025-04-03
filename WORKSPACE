load("@bazel_tools//tools/build_defs/repo:http.bzl", "http_archive")
load("@bazel_tools//tools/build_defs/repo:git.bzl", "git_repository")

git_repository(
    name = "stim",
    commit = "1320ad7eac7de34d2e9c70daa44fbc6d84174450",
    remote = "https://github.com/quantumlib/stim.git",
)

http_archive(
    name = "gtest",
    sha256 = "353571c2440176ded91c2de6d6cd88ddd41401d14692ec1f99e35d013feda55a",
    strip_prefix = "googletest-release-1.11.0",
    urls = ["https://github.com/google/googletest/archive/refs/tags/release-1.11.0.zip"],
)

http_archive(
  name = "pybind11_bazel",
  sha256 = "9df284330336958c837fb70dc34c0a6254dac52a5c983b3373a8c2bbb79ac35e",
  strip_prefix = "pybind11_bazel-2.13.6",
  urls = ["https://github.com/pybind/pybind11_bazel/archive/v2.13.6.zip"],
)

http_archive(
  name = "pybind11",
  sha256 = "d0a116e91f64a4a2d8fb7590c34242df92258a61ec644b79127951e821b47be6",
  build_file = "@pybind11_bazel//:pybind11-BUILD.bazel",
  strip_prefix = "pybind11-2.13.6",
  urls = ["https://github.com/pybind/pybind11/archive/v2.13.6.zip"],
)

http_archive(
    name = "rules_python",
    sha256 = "2cc26bbd53854ceb76dd42a834b1002cd4ba7f8df35440cf03482e045affc244",
    strip_prefix = "rules_python-1.3.0",
    url = "https://github.com/bazel-contrib/rules_python/releases/download/1.3.0/rules_python-1.3.0.tar.gz",
)

load("@rules_python//python:repositories.bzl", "py_repositories")

py_repositories()

load("@pybind11_bazel//:build_defs.bzl", "pybind_extension")
