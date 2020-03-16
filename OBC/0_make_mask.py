import netCDF4
import numpy as np
import os
import datetime as dt

nc = netCDF4.Dataset('domain_cfg.nc')

mask = np.zeros((nc.variables['e3t_0'].shape))
hlev = np.squeeze(nc.variables['bottom_level'][:])
for i in range(hlev.shape[0]):
    for j in range(hlev.shape[1]):
        mask[0,hlev[i,j]:,i,j] = 1

dims = {x:nc.dimensions[x].size for x in nc.dimensions}
vars =[{'name':'nav_lon','type':'float32', 'dims':nc.variables['nav_lon'].dimensions},
        {'name':'nav_lat','type':'float32', 'dims':nc.variables['nav_lat'].dimensions},
        {'name':'tmask','type':'float32', 'dims':nc.variables['e3t_0'].dimensions}]


nco = netCDF4.Dataset('mesh_mask.nc','w')
for d in dims:
    nco.createDimension(d,dims[d])
for v in vars:
    dat=nco.createVariable(v['name'],v['type'],v['dims'])

nco.author = os.getenv('USER')
nco.history = dt.datetime.now().strftime("Created on %a, %B %d, %Y at %H:%M")
nco.title = 'Nemo mesh mask'

nco.variables['nav_lon'][:]=nc.variables['nav_lon'][:]
nco.variables['nav_lat'][:]=nc.variables['nav_lat'][:]
nco.variables['tmask'][:]=mask
nco.close()
nc.close()
