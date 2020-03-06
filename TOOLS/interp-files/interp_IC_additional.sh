#!usr/bin/bash

###############################################
# interp_IC_additional.sh
# This script will perform single interpolation for a variable 
# with a single time record, for variables where the source 
# grid has already had one variable interpolated
###############################################

var=$1       #Input variable name
sfile=$2     #Source file
sourceid=$3  #Source ID Tag

# Fill land mask with zeros
ncatted -a _FillValue,$var,m,f,0 $sfile

#Fill land values
$SOSIEDIR/sosie3.x -f 1_initcd_${sourceid}_to_${sourceid}_${var}.namelist 

# Create weights
$SCRIPDIR/scripinterp.exe 2_${sourceid}_weights_${var}.namelist

# Fill values
$SOSIEDIR/sosie3.x -f 3_initcd_${sourceid}_to_nemo_${var}.namelist 

