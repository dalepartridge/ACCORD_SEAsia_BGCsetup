#!usr/bin/bash

###############################################
# interp_OBC_initial.sh
# This script will perform single interpolation for a variable 
# with a multiple time records, for the first variable from a source
###############################################

var=$1       #Input variable name
sfile=$2     #Source file
sourceid=$3  #Source ID Tag
stime=$4     #Source time variable

#Fnd mask with zeros
ncatted -a _FillValue,$var,m,f,0 $sfile

#Create mask file
ncks -d ${stime},0,0,1 -v $var ${sfile} ${sourceid}_mask.nc
ncrename -v $var,mask ${sourceid}_mask.nc
ncatted -a _FillValue,,d,, ${sourceid}_mask.nc
ncap2 -O -s 'where(mask>0) mask=1' ${sourceid}_mask.nc ${sourceid}_mask.nc

#Fill land values
$SOSIEDIR/sosie3.x -f 1_${sourceid}_to_${sourceid}_${var}.namelist 

# Create weights
$SCRIPDIR/scripgrid.exe 2_${sourceid}_weights_${var}.namelist # creates datagrid_file and nemogrid_file
$SCRIPDIR/scrip.exe 2_${sourceid}_weights_${var}.namelist
cdo splitsel,1 ${var}_${sourceid}-${sourceid}_OBC.nc split_
for f in split* 
do
    sed -i "64 c\ \ \ \ input_file = \"$f\"" 2_${sourceid}_weights_${var}.namelist
    sed -i "74 c\ \ \ \ output_file = \"init_$f\"" 2_${sourceid}_weights_${var}.namelist 
    $SCRIPDIR/scripinterp.exe 2_${sourceid}_weights_${var}.namelist
done
ncrcat init* initcd_${var}.nc
rm -rf split* init_*

#Create mask
ncks -d time_counter,0,0,1 -v $var initcd_${var}.nc sosie_initcd_mask.nc
ncrename -v $var,mask sosie_initcd_mask.nc
ncap2 -O -s 'where(mask>=0) mask=1' sosie_initcd_mask.nc sosie_initcd_mask.nc

# Fill values
$SOSIEDIR/sosie3.x -f 3_${sourceid}_to_nemo_${var}.namelist 

