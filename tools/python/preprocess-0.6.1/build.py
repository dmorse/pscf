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

"""
    preprocess build script

    Usage:
        python build.py [<options>...] [<targets>...]

    Options:
        --help, -h      Print this help and exit.
        --version       Print the version of this script and exit.
        --targets, -t   List all available targets.

    This is the primary build script for the preprocess project. It
    exists to assist in building, maintaining, and distributing this
    project.
    
    It is intended to have Makefile semantics. I.e. 'python build.py'
    will build execute the default target, 'python build.py foo' will
    build target foo, etc. However, there is no intelligent target
    interdependency tracking (I suppose I could do that with function
    attributes).
"""

import os
import sys
import getopt
import types
import getpass
import shutil
import glob


class BuildError(Exception):
    pass

#---- globals

# Site upload information.
gSiteUser = 'tmick'
gSite = 'starship.python.net'
gSiteDownloadsDir = '~tmick/public_html/downloads'


#---- internal logging facility

class _Logger:
    DEBUG, INFO, WARN, ERROR, CRITICAL = range(5)
    def __init__(self, threshold=None, streamOrFileName=sys.stderr):
        if threshold is None:
            self.threshold = self.WARN
        else:
            self.threshold = threshold
        if type(streamOrFileName) == types.StringType:
            self.stream = open(streamOrFileName, 'w')
            self._opennedStream = 1
        else:
            self.stream = streamOrFileName
            self._opennedStream = 0
    def __del__(self):
        if self._opennedStream:
            self.stream.close()
    def _getLevelName(self, level):
        levelNameMap = {
            self.DEBUG: "DEBUG",
            self.INFO: "INFO",
            self.WARN: "WARN",
            self.ERROR: "ERROR",
            self.CRITICAL: "CRITICAL",
        }
        return levelNameMap[level]
    def log(self, level, msg):
        if level < self.threshold:
            return
        message = "%s: " % self._getLevelName(level).lower()
        message = message + msg + "\n"
        self.stream.write(message)
        self.stream.flush()
    def debug(self, msg):
        self.log(self.DEBUG, msg)
    def info(self, msg):
        self.log(self.INFO, msg)
    def warn(self, msg):
        self.log(self.WARN, msg)
    def error(self, msg):
        self.log(self.ERROR, msg)
    def fatal(self, msg):
        self.log(self.CRITICAL, msg)

if 1:   # normal
    log = _Logger(_Logger.INFO)
else:   # debugging
    log = _Logger(_Logger.DEBUG)



#---- globals

_version_ = (0, 1, 0)



#---- support routines

def _rmtree_OnError(rmFunction, filePath, excInfo):
    if excInfo[0] == OSError:
        # presuming because file is read-only
        os.chmod(filePath, 0777)
        rmFunction(filePath)

def _rmtree(dirname):
    import shutil
    shutil.rmtree(dirname, 0, _rmtree_OnError)


def _getTargets():
    """Find all targets and return a dict of targetName:targetFunc items."""
    targets = {}
    for name, attr in sys.modules[__name__].__dict__.items():
        if name.startswith('target_'):
            targets[ name[len('target_'):] ] = attr
    return targets

def _listTargets(targets):
    """Pretty print a list of targets."""
    width = 77
    nameWidth = 15 # min width
    for name in targets.keys():
        nameWidth = max(nameWidth, len(name))
    nameWidth += 2  # space btwn name and doc
    format = "%%-%ds%%s" % nameWidth
    print format % ("TARGET", "DESCRIPTION")
    for name, func in targets.items():
        doc = func.__doc__ or ""
        if len(doc) > (width - nameWidth):
            doc = doc[:(width-nameWidth-3)] + "..."
        print format % (name, doc)


#---- the targets

def target_default():
    """all"""
    log.info("target: default")
    target_all()

def target_all():
    """Build the launcher and a source distro."""
    log.info("target: all")
    if sys.platform.startswith('win'):
        target_launcher()
    target_sdist()

def target_clean():
    log.info("target: clean")
    if sys.platform.startswith('win'):
        cmd = "nmake -f Makefile.win clean"
        retval = os.system(cmd)
        if retval:
            raise BuildError("Error cleaning launcher: cmd=%r" % cmd)

    import preprocess
    ver = '.'.join([str(part) for part in preprocess._version_])
    dirs = ["dist", "build", "preprocess-%s" % ver]
    for d in dirs:
        print "removing '%s'" % d
        if os.path.isdir(d): _rmtree(d)

    patterns = ["*.pyc", "*~", "MANIFEST", "preprocess",
                os.path.join("test", "*~"),
                os.path.join("test", "*.pyc"),
               ]
    for pattern in patterns:
        for file in glob.glob(pattern):
            print "removing '%s'" % file
            os.unlink(file)

def target_launcher():
    """Build the launcher executable."""
    log.info("target: launcher")
    if sys.platform.startswith("win"):
        cmd = "nmake -f Makefile.win"
        retval = os.system(cmd)
        if retval:
            raise BuildError("Error building launcher: cmd=%r" % cmd)

def target_sdist():
    """Build a source distribution for the current platform."""
    log.info("target: sdist")
    cmd = 'python setup.py sdist -f'
    retval = os.system(cmd)
    if retval:
        raise BuildError("Error building sdist: cmd=%r" % cmd)

def target_install():
    """Use the setup.py script to install."""
    log.info("target: install")
    cmd = 'python setup.py install'
    retval = os.system(cmd)
    if retval:
        raise BuildError("Error installing: cmd=%r" % cmd)

def target_upload():
    """Upload binary and source distribution to starship."""
    log.info("target: upload")
    import preprocess
    ver = '.'.join([str(i) for i in preprocess._version_])
    if sys.platform.startswith("win"):
        srcs = [os.path.join('dist', 'preprocess-%s.zip' % ver)]
        doNotExist = [s for s in srcs if not os.path.exists(s)]
        if doNotExist:
            raise BuildError("Cannot upload the following because they do "\
                             "not exist: %r" % doNotExist)
        passwd = getpass.getpass("%s@%s's password: " % (gSiteUser, gSite))
        dstDir = gSiteDownloadsDir + '/preprocess/%s/win32/' % ver
        # Ensure the dst dir is made.
        cmd = 'plink -pw %s %s@%s mkdir -p %s'\
              % (passwd, gSiteUser, gSite, dstDir)
        log.info("Create %s on %s." % (dstDir, gSite))
        retval = os.system(cmd)
        if retval:
            raise BuildError("Error creating %r on %r." % (dstDir, gSite))
        # Upload the src files.
        for src in srcs:
            log.info("Upload %s to %s." % (src, dstDir))
            cmd = 'pscp -pw %s %s %s@%s:%s'\
                  % (passwd, src, gSiteUser, gSite, dstDir)
            retval = os.system(cmd)
            if retval:
                raise BuildError("Error uploading to %r" % dstDir)
    else:
        srcs = [os.path.join('dist', 'preprocess-%s.tar.gz' % ver)]
        doNotExist = [s for s in srcs if not os.path.exists(s)]
        if doNotExist:
            raise BuildError("Cannot upload the following because they do "\
                             "not exist: %r" % doNotExist)
        dstDir = gSiteDownloadsDir + '/preprocess/%s/linux/' % ver
        # Ensure the dst dir is made.
        cmd = 'ssh %s@%s mkdir -p %s' % (gSiteUser, gSite, dstDir)
        log.info("Create %s on %s." % (dstDir, gSite))
        retval = os.system(cmd)
        if retval:
            raise BuildError("Error creating %r on %r." % (dstDir, gSite))
        # Upload the src files.
        for src in srcs:
            log.info("Upload %s to %s." % (src, dstDir))
            cmd = 'scp %s %s@%s:%s' % (src, gSiteUser, gSite, dstDir)
            retval = os.system(cmd)
            if retval:
                raise BuildError("Error uploading to %r" % dstDir)

def target_check_license():
    """List files whose license information is out of date."""
    license = """\
Copyright (c) 2002 Trent Mick

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""
    os.path.walk(os.curdir, _checkForLicenseInFiles, license)

def _checkForLicenseInFile(license, filename):
    licLines = license.split('\n')
    fin = open(filename, 'r')
    srcLines = fin.readlines()
    fin.close()
    # Try to find first license line in the first 30 lines.
    for i in range(30):
        if i >= len(srcLines):
            print filename
            return
        elif srcLines[i].find(licLines[0]) != -1:
            break
    else:
        print filename
        return
    # Ensure that the subsequent lines include the rest of the
    # license.
    for licLine in licLines:
        if i >= len(srcLines) or srcLines[i].find(licLine) == -1:
            print filename
            return
        i += 1

def _checkForLicenseInFiles(license, dirname, basenames): 
    """Check that the given license is near the top of the given source
    files. (Called via os.path.walk.)
    """
    exts = ['.cpp', '.py']
    for basename in basenames:
        base, ext = os.path.splitext(basename)
        ext = ext.lower()
        if ext in exts:
            filename = os.path.join(dirname, basename)
            _checkForLicenseInFile(license, filename)


#---- mainline

def build(targets=[]):
    log.info("build(targets=%r)" % targets)
    available = _getTargets()
    if not targets:
        if available.has_key('default'):
            return available['default']()
        else:   
            log.warn("No default target available. Doing nothing.")
    else:
        for target in targets:
            if available.has_key(target):
                retval = available[target]()
                if retval:
                    raise BuildError("Error running '%s' target: retval=%s"\
                                     % (target, retval))
            else:
                raise BuildError("Unknown target: '%s'" % target)

def main(argv):
    # Process options.
    optlist, targets = getopt.getopt(argv[1:], 'hvt',
        ['help', 'version', 'targets'])
    for opt, optarg in optlist:
        if opt in ('-h', '--help'):
            sys.stdout.write(__doc__ + '\n')
            return 0
        elif opt in ('--version',):
            ver = '.'.join([str(i) for i in _version_])
            sys.stdout.write("preprocess build script %s\n" % ver)
            return 0
        elif opt in ('-t', '--targets'):
            return _listTargets(_getTargets())

    return build(targets)

if __name__ == "__main__":
    sys.exit( main(sys.argv) )

