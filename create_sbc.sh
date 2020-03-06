#Load modules
module load anaconda/python3

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


