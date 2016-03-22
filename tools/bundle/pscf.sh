#!/bin/sh
#
# Author: Evan Bollig <boll0107@umn.edu>

PSCF_BUNDLE="`echo "$0" | sed -e 's/\/Contents\/MacOS\/pscf_terminal//'`"
echo $0
PSCF_RESOURCES="$PSCF_BUNDLE/Contents/Resources"
PSCF_TEMP="/tmp/pscf/$UID"
PSCF_ETC="$PSCF_TEMP/etc"
PSCF_PANGO_RC_FILE="$PSCF_ETC/pango/pangorc"

# TODO: write loop to iterate through each .app and append their paths
# TODO: make sure that we consider Resources dir for each fo the apps
MAIN_EXEC=$PSCF_RESOURCES/pscf.app/Contents/MacOS/
TEST_EXEC=$PSCF_RESOURCES/rpa_test.app/Contents/MacOS/
GROUP_EXEC=$PSCF_RESOURCES/new_2dgroups.app/Contents/MacOS/

echo "running $0"
echo "PSCF_BUNDLE: $PSCF_BUNDLE"

export "DYLD_LIBRARY_PATH=$PSCF_RESOURCES/lib:$MAIN_EXEC:$TEST_EXEC:$GROUP_EXEC:$DYLD_LIBRARY_PATH"
export "PATH=$PSCF_RESOURCES/bin:$MAIN_EXEC:$TEST_EXEC:$GROUP_EXEC:$PATH"

cat > $PSCF_TEMP/terminal <<EOM
export DYLD_LIBRARY_PATH=$PSCF_RESOURCES/lib
export PATH=$PSCF_RESOURCES/bin:$PATH
clear
echo -e "\n\nHello $USER! The Polymer Self-Consistent Field theory (PSCF) code is on your PATH and ready to run. Start with command \"pscf < [your input file]\"\n\n"
EOM

# Start a new terminal window (shell) with the environment properly configured
# to run the PSCF command. No need to have a GUI :-D
chmod +x $PSCF_TEMP/terminal
osascript <<EOD
tell application "Terminal"
    activate
    set newTab to do script ". $PSCF_TEMP/terminal"
    set current settings of newTab to settings set "Red Sands"
end tell
EOD
