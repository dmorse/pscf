#!/usr/bin/env python

# Copyright (c) 2002 Trent Mick
# 
# Permission is hereby granted, free of charge, to any person obtaining a
# copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:
# 
# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
# OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
# CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

"""Distutils setup script for 'preprocess'."""

import sys
import os
import shutil
from distutils.core import setup


#---- support routines

def _getVersion():
    import preprocess
    return '.'.join([str(n) for n in preprocess._version_])

def _getBinDir():
    """Return the current Python's bindir."""
    if sys.platform.startswith("win"):
        bindir = sys.prefix
    else:
        bindir = os.path.join(sys.prefix, "bin")
    return bindir

def _getSiteLibDir():
    if sys.platform.startswith("win"):
        sitelibdir = os.path.join(sys.prefix, "Lib", "site-packages")
    else:
        ver = '.'.join([str(n) for n in sys.version_info[:2]])
        sitelibdir = os.path.join(sys.prefix, "lib", "python%s" % ver,
                                  "site-packages")
    return sitelibdir


#---- setup mainline

if sys.platform.startswith('win'):
    scripts = []
    binFiles = ["preprocess.exe", "preprocess.py"]
else:
    if os.path.exists("preprocess"):
        os.remove("preprocess")
    shutil.copy2("preprocess.py", "preprocess")
    scripts = ["preprocess"]
    binFiles = []
siteFiles = ["content.types"]

setup(name="preprocess",
      version=_getVersion(),
      description="a multi-language preprocessor",
      author="Trent Mick",
      author_email="TrentM@ActiveState.com",
      url="http://starship.python.net/~tmick/",
      license="MIT License",
      platforms=["Windows", "Linux"],
      keywords=["preprocess", "preprocessor", "#if", "#else", "#endif"],

      py_modules=['preprocess'],
      scripts=scripts,
      data_files=[ (_getBinDir(), binFiles),
                   (_getSiteLibDir(), siteFiles) ],
     )

