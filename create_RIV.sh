#Load modules
module unload cray-netcdf-hdf5parallel cray-hdf5-parallel
module load cray-netcdf cray-hdf5
module load nco/4.5.0
module load anaconda/python3

cd $WDIR/RIV/
ln -s $DOMAINFILE domain_cfg.nc
ln -s $RAWDATA/RIVERS/*.nc .

#conda create -n xarr python=3.6
source activate xarr
pip install xarray netcdf4

python add_river_nutrients.py

source deactivate
#conda remove -n xarr --all

cd $WDIR

