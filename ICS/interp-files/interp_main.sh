module unload cray-netcdf-hdf5parallel cray-hdf5-parallel
module load cray-netcdf cray-hdf5
module load nco/4.5.0

var=$1
sfile=$2
sourceid=$3
stime=$4

SOSIEDIR=/home/n01/n01/dapa/code/git/sosie/bin/
SCRIPDIR=/work/n01/n01/annkat/SEAsia_ERSEM_R12/NEMO_ERSEM_CMEMS/trunk_NEMOGCM_r8395/TOOLS/WEIGHTS/BLD/bin/

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
/home/n01/n01/dapa/code/git/sosie/bin/sosie3.x -f 3_initcd_${sourceid}_to_nemo_${var}.namelist 

