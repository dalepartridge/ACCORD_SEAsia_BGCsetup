

################
# HadGEM
###############
cd hadgem
ln -s $RAWDATA/OBC/*.nc .
ln -s $DOMAINFILE domain_cfg.nc
ln -s $RAWDATA/DOMAIN-HADGEM/bdy_gdept.nc .
ln -s $TOOLS/interp-files/interp_OBC*.sh .
ln -s ../mesh_mask.nc .

python interp_hadgem.py $TOOLS/interp-files/namelist-templates/
cd ..




