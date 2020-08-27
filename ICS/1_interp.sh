

################
# WOA18
###############
cd woa18
ln -s $RAWDATA/ICS/woa18_nitrate_extracted_mean.nc woa18_nitrate.nc
ln -s $RAWDATA/ICS/woa18_phosphate_extracted_mean.nc woa18_phosphate.nc
ln -s $RAWDATA/ICS/woa18_silicate_extracted_mean.nc woa18_silicate.nc
ln -s $DOMAINFILE domain_cfg.nc
ln -s $TOOLS/interp-files/interp_IC_*.sh .
python interp_woa.py $TOOLS/interp-files/namelist-templates/
cd ..

################
# WOA18-OXYGEN
###############
cd woa18-oxy
ln -s $RAWDATA/ICS/woa18_oxy_extracted_mean.nc woa18_oxygen.nc
ln -s $DOMAINFILE domain_cfg.nc
ln -s $TOOLS/interp-files/interp_IC_initial.sh .
python interp_woaoxy.py $TOOLS/interp-files/namelist-templates/
cd ..

################
# GLODAP
###############
cd glodap
ln -s $RAWDATA/ICS/GLODAP_TAlk_extracted.nc GLODAP_TAlk.nc 
ln -s $RAWDATA/ICS/GLODAP_TCO2_extracted.nc GLODAP_TCO2.nc
ln -s $DOMAINFILE domain_cfg.nc
ln -s $TOOLS/interp-files/interp_IC_*.sh .
python interp_glodap.py $TOOLS/interp-files/namelist-templates/
cd ..

################
# OCCCI
###############
cd occci
ln -s $RAWDATA/ICS/CCI-OC_chla_extracted_jan_mean.nc CCI-OC_chla.nc
ln -s $DOMAINFILE domain_cfg.nc
ln -s $TOOLS/interp-files/interp_IC_initial.sh .
python interp_occci.py $TOOLS/interp-files/namelist-templates/
cd ..

################
# ADY
###############
cd ady
ln -s $RAWDATA/ICS/adyBroadBandClimatology_ACCORD_extracted.nc adyBroadBandClimatology.nc
ln -s $DOMAINFILE domain_cfg.nc
ln -s $TOOLS/interp-files/interp_IC_initial.sh .
python interp_ady.py $TOOLS/interp-files/namelist-templates/
cd ..

################
# iMarNet
###############
cd imarnet
ln -s $RAWDATA/ICS/iMarNet_data.nc .
ln -s $DOMAINFILE domain_cfg.nc
ln -s $TOOLS/interp-files/interp_IC_*.sh .
python interp_imarnet.py $TOOLS/interp-files/namelist-templates/
cd ..




