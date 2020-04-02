import netCDF4
from cdl_parser import cdl_parser
import datetime as dt
import os 

nc = netCDF4.Dataset('domain_cfg.nc')
lat = nc.variables['nav_lat'][:]
lon = nc.variables['nav_lon'][:]
z = nc.variables['nav_lev'][:]
nc.close()

_,vars,_ = cdl_parser('nemo_ini_trc.cdl')

dims = {'y': lon.shape[0], 'x': lon.shape[1], 'z':len(z), 't': None}

nco = netCDF4.Dataset('bgc_ini.nc','w')
for d in dims:
    nco.createDimension(d,dims[d])
for v in vars:
    dat=nco.createVariable(v['name'],v['type'],v['dims'])
    dat.units = v['attr']['units']
    dat.long_name = v['attr']['long_name']

nco.createVariable('kt','double',{})

nco.author = os.getenv('USER')
nco.history = dt.datetime.now().strftime("Created on %a, %B %d, %Y at %H:%M")
nco.title = 'Nemo initial conditions'

# Set variables
nco.variables['time_counter'][:] = 1
nco.variables['nav_lat'][:] = lat
nco.variables['nav_lon'][:] = lon
nco.variables['nav_lev'][:] = z
nco.close()
