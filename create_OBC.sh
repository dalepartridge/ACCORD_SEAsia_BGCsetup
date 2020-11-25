
#Load modules
module unload cray-netcdf-hdf5parallel cray-hdf5-parallel
module load cray-netcdf cray-hdf5
module load nco/4.5.0
module load anaconda/python3

cd $WDIR/OBC/
ln -s $DOMAINFILE domain_cfg.nc
python 0_make_mask.py
. 1_interp.sh

mkdir bdyfiles
ln -s $RAWDATA/DOMAIN/coordinates.bdy.nc .
python 2_create_OBC_file.py

ln -s $RAWDATA/DOMAIN/bdy_gdept.nc .
python 3_extract_OBC.py

python 4_clean_data.py

. 5_add_gdep.sh

cd bdyfiles
for i in *.nc; do;
cp $i ${a/2017/2018}
cp accord_bdytrc_y2017m01.nc accord_bdytrc_y2016m12.nc
cp accord_bdytrc_y2018m12.nc accord_bdytrc_y2019m01.nc
cd $WDIR

