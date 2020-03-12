'''
Script to add benthic initial conditions from GOTM runs to ACCORD domain
'''
import netCDF4
import numpy as np

##########################
# Create/Open File
##########################
outfile = 'bgc_ini.nc'
nco = netCDF4.Dataset(outfile,'a')

pf = 'fabm_st2Dn'

dat = nco.variables[pf+'Y2_c'].shape

ben_dict = {
    'Y2':{'c':5669.74}, 'Y3':{'c':5104.38}, 'Y4':{'c':59.5165}, \
    'H1':{'c':39.0568}, 'H2':{'c':689.821}, \
    'K3':{'n':1.05327}, 'K4':{'n':0.0214783}, 'K1':{'p':11.8849}, 'K5':{'s':2.52971}, \
    'G2':{'o':7.02669, 'o_deep':0}, 'G3':{'c':134.886}, \
    'Q1':{'c':6.10334,'n':0.0624711,'p':0.00368191}, \
    'Q6':{'c':665.254,'n':3.79317,'p':0.24895,'s':5.14074, \
        'pen_depth_c':0.0143457,'pen_depth_n':0.00904731,'pen_depth_p':0.00903652,'pen_depth_s':0.00493341}, \
    'Q7':{'c':8795.6,'n':98.3962,'p':5.2213, \
        'pen_depth_c':0.180185,'pen_depth_n':0.166656,'pen_depth_p':0.226404}, \
    'Q17':{'c':0,'n':0,'p':0}, \
    'bL2':{'c':1.41338}, 'ben_nit':{'G4n':0}}

print('Filling fields')
for p in ben_dict:
    for q in ben_dict[p]:
        nco.variables[pf+p+'_'+q][:] = ben_dict[p][q]*np.ones(dat)

print('Calculating horizons')
ncbath = netCDF4.Dataset('/work/n01/n01/annkat/SEAsia_ERSEM_R12/DOMAIN/bathy_meter_ORCA12_adj.nc')
h = ncbath.variables['Bathymetry'][:]
ncbath.close()

D1 = 0.0026 * h**0.3286 
D1[D1<0.01] = 0.01 #Enforce minimum depth
nco.variables[pf+'ben_col_D1m'][:] = D1[np.newaxis,:,:]

D2 = 0.011 * h**0.4018
D2[D2<0.01] = 0.01 #Enforce minimum depth 1cm
tot = D1 + D2
D2[tot>0.3] = 0.3 - D1[tot>0.3] # Enforce total thickness less than 30cm
nco.variables[pf+'ben_col_D2m'][:] = D2[np.newaxis,:,:]



#############################################################################

print('Filling Before Fields')
fields = list(filter(lambda x: x.startswith('fabm_st2Dn'), nco.variables.keys()))
for i in fields:
    nco.variables[i.replace('st2Dn','st2Db')][:] = nco.variables[i][:]
nco.sync()

print('Complete')
nco.close()
