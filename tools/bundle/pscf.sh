#!/bin/sh
#
# Author: Evan Bollig <boll0107@umn.edu>

# Original: 
# Author: Aaron Voisine <aaron@voisine.org>
# Inkscape Modifications: Michael Wybrow <mjwybrow@users.sourceforge.net>
# K-3D Modifications: Timothy M. Shead <tshead@k-3d.com>

PSCF_BUNDLE="`echo "$0" | sed -e 's/\/Contents\/MacOS\/pscf//'`"
PSCF_RESOURCES="$PSCF_BUNDLE/Contents/Resources"
PSCF_TEMP="/tmp/pscf/$UID"
PSCF_ETC="$PSCF_TEMP/etc"
PSCF_PANGO_RC_FILE="$PSCF_ETC/pango/pangorc"

echo "running $0"
echo "PSCF_BUNDLE: $PSCF_BUNDLE"
#
## Start X11 ...
#ps -wx -ocommand | grep -e '[X]11.app' > /dev/null
#if [ "$?" != "0" -a ! -f ~/.xinitrc ]; then
#    echo "rm -f ~/.xinitrc" > ~/.xinitrc
#    sed 's/xterm/# xterm/' /usr/X11R6/lib/X11/xinit/xinitrc >> ~/.xinitrc
#fi
#
#mkdir -p $PSCF_TEMP
#cat << __END_OF_GETDISPLAY_SCRIPT__ > "$PSCF_TEMP/getdisplay.sh"
##!/bin/sh
#mkdir -p "$PSCF_TEMP"
#
#if [ "\$DISPLAY"x == "x" ]; then
#    echo :0 > "$PSCF_TEMP/display"
#else
#    echo \$DISPLAY > "$PSCF_TEMP/display"
#fi
#__END_OF_GETDISPLAY_SCRIPT__
#chmod +x "$PSCF_TEMP/getdisplay.sh"
#rm -f $PSCF_TEMP/display
#open-x11 $PSCF_TEMP/getdisplay.sh || \
#open -a XDarwin $PSCF_TEMP/getdisplay.sh || \
#echo ":0" > $PSCF_TEMP/display
#
#while [ "$?" == "0" -a ! -f $PSCF_TEMP/display ];
#do
#  #echo "Waiting for display $PSCF_TEMP/display"
#  sleep 1;
#done
#export "DISPLAY=`cat $PSCF_TEMP/display`"
#
#ps -wx -ocommand | grep -e '[X]11' > /dev/null || exit 11
#
## Setup temporary runtime files
#rm -rf "$PSCF_TEMP"
#
## Because the bundle could be located anywhere at runtime, we have to
## create temporary copies of the Pango configuration files that
## reflect our current location
#mkdir -p "$PSCF_ETC/pango"
#sed -e 's|/opt/local/etc|'"$PSCF_ETC|g" "$PSCF_RESOURCES/etc/pango/pangorc" > "$PSCF_ETC/pango/pangorc"
#sed -e 's|/opt/local|\"'"$PSCF_RESOURCES|g" -e "s/\.so/.so\"/g" "$PSCF_RESOURCES/etc/pango/pango.modules" > "$PSCF_ETC/pango/pango.modules"
#cp -f "$PSCF_RESOURCES/etc/pango/pangox.aliases" "$PSCF_ETC/pango/pangox.aliases"

export "DYLD_LIBRARY_PATH=$PSCF_RESOURCES/lib"
export "FONTCONFIG_PATH=$PSCF_RESOURCES/etc/fonts"
export "PANGO_RC_FILE=$PSCF_PANGO_RC_FILE"
export "PATH=$PSCF_RESOURCES/bin:$PATH"

#export
#exec "$PSCF_RESOURCES/bin/pscf" u
#"--log-level=debug" "--plugins=$PSCF_RESOURCES/lib/pscf/plugins" "--share=$PSCF_RESOURCES/share/pscf" "--ui=$PSCF_RESOURCES/lib/pscf/uiplugins/pscf-ngui.module"

osascript -e 'display notification "Lorem ipsum dolor sit amet" with title "Title"'

osascript -e 'tell app "Terminal"
    do script "echo hello"
end tell'
