#!usr/bin/bash

###############################################
# interp_IC_initial.sh
# This script will perform single interpolation for a variable 
# with a single time record, for the first variable from a source
###############################################

var=$1       #Input variable name
sfile=$2     #Source file
sourceid=$3  #Source ID Tag
stime=$4     #Source time variable

# Fill land mask with zeros
ncatted -a _FillValue,$var,m,f,0 $sfile

#Create mask file
ncks -d ${stime},0,0,1 -v $var ${sfile} ${sourceid}_mask.nc
ncrename -v $var,mask ${sourceid}_mask.nc
ncatted -a _FillValue,,d,, ${sourceid}_mask.nc
ncap2 -O -s 'where(mask>0) mask=1' ${sourceid}_mask.nc ${sourceid}_mask.nc

#Fill land values
$SOSIEDIR/sosie3.x -f 1_initcd_${sourceid}_to_${sourceid}_${var}.namelist 

# Create weights
$SCRIPDIR/scripgrid.exe 2_${sourceid}_weights_${var}.namelist # creates datagrid_file and nemogrid_file
$SCRIPDIR/scrip.exe 2_${sourceid}_weights_${var}.namelist
$SCRIPDIR/scripinterp.exe 2_${sourceid}_weights_${var}.namelist

#Create mask
ncks -d time_counter,0,0,1 -v $var initcd_${var}.nc sosie_initcd_mask.nc
ncrename -v $var,mask sosie_initcd_mask.nc
ncap2 -O -s 'where(mask>=0) mask=1' sosie_initcd_mask.nc sosie_initcd_mask.nc

# Fill values
$SOSIEDIR/sosie3.x -f 3_initcd_${sourceid}_to_nemo_${var}.namelist 

