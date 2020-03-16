

################
# HadGEM
###############
ln -s $RAWDATA/OBC/ALK.nc hadgem/
ln -s $RAWDATA/OBC/DIC.nc hadgem/
ln -s $RAWDATA/OBC/DIN.nc hadgem/
ln -s $RAWDATA/OBC/SIL.nc hadgem/
ln -s $RAWDATA/OBC/OXY.nc hadgem/
ln -s $DOMAINFILE domain_cfg.nc

python 0_make_mask.py
python 1_interp.py 1_interp_SIL.yaml
python 1_interp.py 1_interp_DIN.yaml
python 1_interp.py 1_interp_DIC.yaml
python 1_interp.py 1_interp_OXY.yaml
python 1_interp.py 1_interp_ALK.yaml





