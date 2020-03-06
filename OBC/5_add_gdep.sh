
#module unload cray-netcdf-hdf5parallel cray-hdf5-parallel
#module load cray-netcdf cray-hdf5
#module load nco/4.5.0

for i in hadgem-data/accord*nc; do
  ncks -A hadgem-coordinates/bdy_gdept_hadgem.nc $i
done

