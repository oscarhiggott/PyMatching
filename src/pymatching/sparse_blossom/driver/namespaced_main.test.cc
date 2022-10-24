// Copyright 2022 PyMatching Contributors
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "pymatching/sparse_blossom/driver/namespaced_main.h"

#include "gtest/gtest.h"

struct RaiiTempNamedFile {
    int descriptor;
    std::string path;
    RaiiTempNamedFile();
    ~RaiiTempNamedFile();
};

RaiiTempNamedFile::RaiiTempNamedFile() {
    char tmp_stdin_filename[] = "/tmp/stim_test_named_file_XXXXXX";
    descriptor = mkstemp(tmp_stdin_filename);
    if (descriptor == -1) {
        throw std::runtime_error("Failed to create temporary file.");
    }
    path = tmp_stdin_filename;
}

RaiiTempNamedFile::~RaiiTempNamedFile() {
    if (!path.empty()) {
        remove(path.data());
        path = "";
    }
}

std::string result_of_running_main(const std::vector<std::string> args, const std::string &input) {
    std::vector<const char *> argv;
    argv.push_back("TEST_PROCESS");
    for (const auto &a : args) {
        argv.push_back(a.c_str());
    }
    RaiiTempNamedFile inp;
    RaiiTempNamedFile out;
    argv.push_back("--in");
    argv.push_back(inp.path.c_str());
    argv.push_back("--out");
    argv.push_back(out.path.c_str());
    FILE *f = fopen(inp.path.c_str(), "w");
    if (f == nullptr || fwrite(input.data(), 1, input.size(), f) != input.size()) {
        throw std::invalid_argument("Failed to write input.");
    }
    fclose(f);
    if (pm::main((int)argv.size(), argv.data()) != EXIT_SUCCESS) {
        throw std::invalid_argument("Returned not EXIT_SUCCESS");
    }
    f = fopen(out.path.c_str(), "r");
    if (f == nullptr) {
        throw std::invalid_argument("Failed to read output.");
    }
    std::string s;
    while (true) {
        int i = getc(f);
        if (i == EOF) {
            break;
        }
        s.push_back((char)i);
    }
    return s;
}

TEST(Main, predict) {
    RaiiTempNamedFile dem;
    FILE *f = fopen(dem.path.c_str(), "w");
    fprintf(f, "%s", R"DEM(
        error(0.1) D0 L0
        error(0.1) D0 D1 L1
        error(0.1) D1 L2
    )DEM");
    fclose(f);
    auto stdout = result_of_running_main(
        {"predict", "--dem", dem.path, "--in_format", "dets", "--out_format", "dets"},
        R"stdin(shot
shot D0
shot D1
shot D0 D1)stdin");
    ASSERT_EQ(stdout, R"stdout(shot
shot L0
shot L2
shot L1
)stdout");
}

TEST(Main, count_mistakes) {
    RaiiTempNamedFile dem;
    FILE *f = fopen(dem.path.c_str(), "w");
    fprintf(f, "%s", R"DEM(
        error(0.1) D0 L0
        error(0.1) D0 D1 L1
        error(0.1) D1 L2
    )DEM");
    fclose(f);
    auto stdout_text = result_of_running_main(
        {
            "count_mistakes",
            "--dem",
            dem.path,
            "--in_format",
            "dets",
            "--in_includes_appended_observables",
        },
        R"stdin(shot L0
shot D0 L0
shot D1 L2
shot D0 D1 L1)stdin");
    ASSERT_EQ(stdout_text, "1 / 4\n");
}
