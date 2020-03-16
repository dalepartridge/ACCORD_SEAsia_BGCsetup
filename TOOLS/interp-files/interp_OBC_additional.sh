#!usr/bin/bash

###############################################
# interp_OBC_additional.sh
# This script will perform interpolation for a variable 
# with a multiple time records, for variables where the source 
# grid has already had one variable interpolated
###############################################

var=$1       #Input variable name
sfile=$2     #Source file
sourceid=$3  #Source ID Tag

# Fill land mask with zeros
#ncatted -a _FillValue,$var,m,f,0 $sfile
python fill_land_mask.py ALK.nc ALK

#Fill land values
$SOSIEDIR/sosie3.x -f 1_initcd_${sourceid}_to_${sourceid}_${var}.namelist 

# Split file into individual files for each time record
cdo splitsel,1 ${var}_${sourceid}-${sourceid}_OBC.nc split_
for f in split* 
do
    sed -i "64 c\ \ \ \ input_file = \"$f\"" 2_woa18_weights_${var}.namelist
    sed -i "74 c\ \ \ \ output_file = \"init_$f\"" 2_woa18_weights_${var}.namelist 
    $SCRIPDIR/scripinterp.exe 2_${sourceid}_weights_${var}.namelist
done   
ncrcat init_split* initcd_${var}.nc
rm -rf split* init_split*

# Fill values
sed -i "88 ccf_z_src   = \'bdy_gdept.nc\'" 3_initcd_${sourceid}_to_nemo_${var}.namelist
sed -i "89 ccv_z_src   = \'gdept\'" 3_initcd_${sourceid}_to_nemo_${var}.namelist
$SOSIEDIR/sosie3.x -f 3_initcd_${sourceid}_to_nemo_${var}.namelist 

