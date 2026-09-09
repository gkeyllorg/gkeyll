#!/bin/bash

source ./build-opts.sh

# Install prefix
PREFIX=$GKYLSOFT
ADAS_DIR="$GKYLSOFT/gkeyll/share"

mkdir -p $ADAS_DIR/adas
#cd -
python3 ../gyrokinetic/data/adas/process_adas.py
rm *.dat
cp *.npy $ADAS_DIR/adas/.
mv *.npy ../gyrokinetic/data/adas/.
echo "ADAS data downloaded to $ADAS_DIR/adas"
