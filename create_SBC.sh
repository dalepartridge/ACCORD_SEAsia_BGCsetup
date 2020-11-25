#Load modules
module unload cray-netcdf-hdf5parallel cray-hdf5-parallel
module load cray-netcdf cray-hdf5
module load nco/4.5.0
module load anaconda/python3
module load cdo

#####################
# NITROGEN DEPOSITION
#####################

# Create conda environment to use python package xarray
conda create -n ndep_env python=3.6
source activate ndep_env
pip install xarray pandas netcdf4 scipy

#Process Nitrogen Deposition
cd $WDIR/SBC/Ndep/
ln -s $DOMAINFILE domain_cfg.nc
ln -s $RAWDATA/SBC/Ndep/oxidized_reduced_Ndeposition.csv .
python process_Ndep.py

# Deactivate and remove environment
source deactivate
conda remove -n ndep_env --all

#####################
# LIGHT ABSORPTION
#####################

cd $WDIR/SBC/ady/
ln -s $RAWDATA/SBC/ady/adyBroadBandClimatology_ACCORD_extracted.nc ady_accord.nc
ln -s $DOMAINFILE domain_cfg.nc
ln -s $TOOLS/interp-files/interp_SBC_initial.sh .
python interp_ady.py $TOOLS/interp-files/namelist-templates/

####################
# Surface CO2
####################

cd $WDIR/SBC/pCO2
ln -s $RAWDATA/SBC/CO2/RCP85_MIDYR_CONC.DAT .
python process_pCO2.py
cdo -splitsel,1 pCO2a.nc pCO2_y
year=1989
for f in pCO2_y*; do
  mv $f ${f::6}${year}.nc
  year=$((year+1))
done

