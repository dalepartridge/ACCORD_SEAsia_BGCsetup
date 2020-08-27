
#Load modules
module unload cray-netcdf-hdf5parallel cray-hdf5-parallel
module load cray-netcdf cray-hdf5
module load nco/4.5.0
module load anaconda/python3

cd $WDIR/OBC/
. 1_interp.sh

mkdir bdyfiles
ln -s $RAWDATA/DOMAIN-HADGEM/coordinates.bdy.nc .
python 2_create_OBC_file.py

ln -s $RAWDATA/DOMAIN-HADGEM/bdy_gdept.nc .
python 3_extract_OBC.py

python 4_fill_phosphate.py

. 5_add_gdep.sh

cd bdyfiles
cp accord_bdytrc_y1980m01.nc accord_bdytrc_y1979m12.nc
cp accord_bdytrc_y2099m12.nc accord_bdytrc_y2100m01.nc
cd $WDIR

