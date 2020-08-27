
#Load modules
module unload cray-netcdf-hdf5parallel cray-hdf5-parallel
module load cray-netcdf cray-hdf5
module load nco/4.5.0
module load anaconda/python3

cd $WDIR/ICS/
#. 1_interp.sh

ln -s $DOMAINFILE domain_cfg.nc
ln -s $TOOLS/python-scripts/cdl_parser.py .
ln -s $TOOLS/python-scripts/cdl/nemo_ini_trc.cdl  .
python 2_create_IC_file.py

ln -s $RAWDATA/DOMAIN-HADGEM/ACCORD_density.nc density.nc
python 3_extract_pelagic.py

python 4_extract_benthic.py

. 5_rename_bioalk_TA.sh

cd $WDIR

