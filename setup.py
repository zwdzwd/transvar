#!/usr/bin/env python

# The MIT License
# 
# Copyright (c) 2011 Seoul National University
# Copyright (c) 2014, 2015 The University of Texas MD Anderson Cancer Center
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:
#
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
# BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
# ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
# CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
# Contact: Wanding Zhou <zhouwanding@gmail.com>

import sys
from setuptools import setup, Extension

def main():
    if float(sys.version[:3])<2.6 or float(sys.version[:3])>=2.8:
        sys.stderr.write("CRITICAL: Python version must be 2.6 or 2.7!\n")
        sys.exit(1)
        
    ext_modules = [
        Extension("tabix",
                  sources = [
                      "pytabix/bgzf.c", "pytabix/bgzip.c", "pytabix/index.c",
                      "pytabix/knetfile.c", "pytabix/kstring.c", "pytabix/tabixmodule.c"
                  ],
                  include_dirs=["pytabix"],
                  libraries=["z"],
                  define_macros=[("_FILE_OFFSET_BITS", 64), ("_USE_KNETFILE", 1)],
                  extra_compile_args=["-w"],
              ),
        Extension("_sswlib",
                  sources = ['ssw/ssw.c', 'ssw/encode.c'],
                  include_dirs = ['ssw'],
                  extra_compile_args = ['-W', '-Wall', '-O2', '-finline-functions', '-fPIC', '-shared', '-Wl,-soname,sswlib'],
              ),
    ]


    setup(
        name = "TransVar",
        version = "2.0.2.20150407",
        description = "Transcript-based Variant annotator",
        url = "https://bitbucket.org/wanding/transvar",
        author = "Wanding Zhou",
        author_email = "zhouwanding@gmail.com",
        license = "MIT",
        keywords = ["bioinformatics", "genomics"],
        scripts = ['bin/transvar'],
        packages = ['transvar', 'ssw'],
        ext_modules = ext_modules,
        classifiers = [
            "Programming Language :: Python",
            "Development Status :: 4 - Beta",
            "Intended Audience :: Developers",
            "Intended Audience :: Science/Research",
            "License :: OSI Approved :: MIT License",
            "Operating System :: OS Independent",
            'Operating System :: POSIX',
            "Programming Language :: C",
            "Topic :: Scientific/Engineering :: Bio-Informatics"
        ],
        # long_description = """ """
        # install_requires=['numpy>=1.6']
    )

if __name__ == '__main__':
    main()
