module load cdo

cd woa18
ln -s $RAWDATA/OBC/ext*.nc .
ln -s $DOMAINFILE domain_cfg.nc
ln -s $RAWDATA/DOMAIN/bdy_gdept_reduced.nc bdy_gdept.nc
ln -s $TOOLS/interp-files/interp_OBC*.sh .
ln -s ../mesh_mask.nc .

python interp_woa18.py $TOOLS/interp-files/namelist-templates/
cd ..

cd woa18-oxy
ln -s $RAWDATA/OBC/ext_woa18_oxygen.nc .
ln -s $DOMAINFILE domain_cfg.nc
ln -s $RAWDATA/DOMAIN/bdy_gdept_reduced.nc bdy_gdept.nc
ln -s $TOOLS/interp-files/interp_OBC*.sh .
ln -s ../mesh_mask.nc .

python interp_woa18.py $TOOLS/interp-files/namelist-templates/
cd ..

