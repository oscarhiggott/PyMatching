# Copyright 2020 Oscar Higgott

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#      http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import os
import re
import sys
import platform
import subprocess
import urllib.request
import shutil

from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext
from distutils.version import LooseVersion


def pip_install(package):
    return subprocess.check_call([sys.executable, '-m', 'pip', 'install', package])

# CMake is needed only if built from source (not for wheels)
pip_install("cmake")

root_dir = os.path.dirname(os.path.realpath(__file__))
lib_dir = os.path.join(root_dir, "lib")


def download_and_extract(pkg_url, pkg_fn, pkg_orig_dir=None, pkg_new_dir=None):
    if not os.path.isfile(pkg_fn):
        try:
            with urllib.request.urlopen(pkg_url) as r, open(pkg_fn, 'wb') as f:
                shutil.copyfileobj(r, f)
        except:
            print(f"Failed to download {pkg_url} dependency, aborting installation.")
            raise
    shutil.unpack_archive(pkg_fn, lib_dir)
    if pkg_orig_dir is not None and pkg_new_dir is not None:
        os.rename(pkg_orig_dir, pkg_new_dir)
    if os.path.isfile(pkg_fn):
        os.remove(pkg_fn)

lemon_url = "http://lemon.cs.elte.hu/hg/lemon-1.3/archive/e5af35e6c93f.tar.gz"
# lemon_url = "http://lemon.cs.elte.hu/pub/sources/lemon-1.3.1.tar.gz"
lemon_fn = os.path.join(root_dir, "lemon-1-3-e5af35e6c93f.tar.gz")
lemon_old_dir = os.path.join(lib_dir, "lemon-1-3-e5af35e6c93f")
lemon_new_dir = os.path.join(lib_dir, "lemon")

if not os.path.isfile(os.path.join(lemon_new_dir, "CMakeLists.txt")):
    download_and_extract(lemon_url, lemon_fn, lemon_old_dir, lemon_new_dir)


class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=''):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    def run(self):
        try:
            out = subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError("CMake must be installed to build the following extensions: " +
                               ", ".join(e.name for e in self.extensions))

        if platform.system() == "Windows":
            cmake_version = LooseVersion(re.search(r'version\s*([\d.]+)', out.decode()).group(1))
            if cmake_version < '3.1.0':
                raise RuntimeError("CMake >= 3.1.0 is required on Windows")

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        cmake_args = ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + extdir,
                      '-DPYTHON_EXECUTABLE=' + sys.executable]

        cfg = 'Debug' if self.debug else 'Release'
        build_args = ['--config', cfg]

        if platform.system() == "Windows":
            cmake_args += ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}'.format(cfg.upper(), extdir)]
            if sys.maxsize > 2**32:
                cmake_args += ['-A', 'x64']
            build_args += ['--', '/m']
        else:
            cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]
            build_args += ['--', '-j2']

        env = os.environ.copy()
        env['CXXFLAGS'] = '{} -DVERSION_INFO=\\"{}\\"'.format(env.get('CXXFLAGS', ''),
                                                              self.distribution.get_version())
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
        subprocess.check_call(['cmake', ext.sourcedir] + cmake_args, cwd=self.build_temp, env=env)
        subprocess.check_call(['cmake', '--build', '.'] + build_args, cwd=self.build_temp)


version = {}
here = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(here,"src/pymatching/_version.py")) as fp:
    exec(fp.read(), version)


with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()


setup(
    name='PyMatching',
    version=version['__version__'],
    author='Oscar Higgott',
    description='A package for decoding quantum error correcting codes using minimum-weight perfect matching.',
    url="https://github.com/oscarhiggott/PyMatching",
    ext_modules=[CMakeExtension('pymatching._cpp_mwpm')],
    packages=find_packages("src"),
    package_dir={'':'src'},
    cmdclass=dict(build_ext=CMakeBuild),
    install_requires=['scipy', 'numpy', 'networkx','matplotlib'],
    classifiers=[
        "License :: OSI Approved :: Apache Software License"
    ],
    long_description=long_description,
    long_description_content_type='text/markdown',
    python_requires='>=3.6',
    zip_safe=False,
)