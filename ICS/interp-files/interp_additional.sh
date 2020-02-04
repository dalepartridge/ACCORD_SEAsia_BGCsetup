module unload cray-netcdf-hdf5parallel cray-hdf5-parallel
module load cray-netcdf cray-hdf5
module load nco/4.5.0

var=$1
sfile=$2
sourceid=$3

SOSIEDIR=/home/n01/n01/dapa/code/git/sosie/bin/
SCRIPDIR=/work/n01/n01/annkat/SEAsia_ERSEM_R12/NEMO_ERSEM_CMEMS/trunk_NEMOGCM_r8395/TOOLS/WEIGHTS/BLD/bin/

# Fill land mask with zeros
ncatted -a _FillValue,$var,m,f,0 $sfile

#Fill land values
$SOSIEDIR/sosie3.x -f 1_initcd_${sourceid}_to_${sourceid}_${var}.namelist 

# Create weights
$SCRIPDIR/scripinterp.exe 2_${sourceid}_weights_${var}.namelist

# Fill values
/home/n01/n01/dapa/code/git/sosie/bin/sosie3.x -f 3_initcd_${sourceid}_to_nemo_${var}.namelist 

