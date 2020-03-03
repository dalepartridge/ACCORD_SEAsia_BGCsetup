module unload cray-netcdf-hdf5parallel cray-hdf5-parallel
module load cray-netcdf cray-hdf5
module load nco/4.5.0

ncrename -v TRNO3_bioalk,TRNO3_TA accord_bgc_ini.nc
ncrename -v TRBO3_bioalk,TRBO3_TA accord_bgc_ini.nc


