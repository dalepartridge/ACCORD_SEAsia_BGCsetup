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

python fill_mask.py $sfile $var ${sourceid}_mask.nc

#Fill land values
$SOSIEDIR/sosie3.x -f 1_${sourceid}_to_${sourceid}_${var}.namelist 

python fix_first_record.py ${var}_${sourceid}-${sourceid}_OBC.nc $var

# Split file into individual files for each time record
cdo splitsel,1 ${var}_${sourceid}-${sourceid}_OBC.nc split_
for f in split* 
do
    sed -i "64 c\ \ \ \ input_file = \"$f\"" 2_${sourceid}_weights_${var}.namelist
    sed -i "74 c\ \ \ \ output_file = \"init_$f\"" 2_${sourceid}_weights_${var}.namelist 
    $SCRIPDIR/scripinterp.exe 2_${sourceid}_weights_${var}.namelist
done   
ncrcat init_split* initcd_${var}.nc
rm -rf split* init_split*

# Fill values
sed -i "88 ccf_z_src   = \'bdy_gdept.nc\'" 3_${sourceid}_to_nemo_${var}.namelist
sed -i "89 ccv_z_src   = \'gdept\'" 3_${sourceid}_to_nemo_${var}.namelist
$SOSIEDIR/sosie3.x -f 3_${sourceid}_to_nemo_${var}.namelist 

