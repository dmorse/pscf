#!/bin/sh
#
# Author: Evan Bollig <boll0107@umn.edu>

PSCF_BUNDLE="`echo "$0" | sed -e 's/\/Contents\/MacOS\/pscf//'`"
PSCF_RESOURCES="$PSCF_BUNDLE/Contents/Resources"
PSCF_TEMP="/tmp/pscf/$UID"
PSCF_ETC="$PSCF_TEMP/etc"
PSCF_PANGO_RC_FILE="$PSCF_ETC/pango/pangorc"

echo "running $0"
echo "PSCF_BUNDLE: $PSCF_BUNDLE"

export "DYLD_LIBRARY_PATH=$PSCF_RESOURCES/lib:$DYLD_LIBRARY_PATH"
export "PATH=$PSCF_RESOURCES/bin:$PATH"

cat > $PSCF_TEMP/terminal <<EOM
export DYLD_LIBRARY_PATH=$PSCF_RESOURCES/lib
export PATH=$PSCF_RESOURCES/bin:$PATH
echo "Hello $USER! The Polymer Self-Consistent Field theory (PSCF) code is on your PATH and ready to run. Start with command \"pscf < [your input file]\""
EOM

# Start a new terminal window (shell) with the environment properly configured
# to run the PSCF command. No need to have a GUI :-D
chmod +x $PSCF_TEMP/terminal
osascript <<EOD
tell application "Terminal"
activate
set shell to do script ". $PSCF_TEMP/terminal"
end tell
EOD
