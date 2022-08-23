#!/bin/bash
set -e

find . -name "*.svg" | parallel --gnu -I% --max-args 1 inkscape % -w 100 -h 100 -o %.png

ffmpeg -framerate 30 -pattern_type glob -i '*.svg' -c:v libx264 -pix_fmt yuv420p out.mp4

#
#for filename in "${@:2}"; do
#  inkscape "${filename}" -w 100 -h 100 -o "${filename}.png"
#  echo "Converted ${filename}"
#done
#echo "Done"

