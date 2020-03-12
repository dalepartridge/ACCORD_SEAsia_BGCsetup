'''
Script to compile pelagic initial conditions from variety of sources into a single file, with some quality checks, for the ACCORD domain
'''

import netCDF4
import numpy as np

#Load Density
nc = netCDF4.Dataset('density.nc')
dens = np.squeeze(nc.variables['rhop'][:])
nc.close()

outfile = 'bgc_ini.nc'
nco = netCDF4.Dataset(outfile,'a')

##########################
# Set Nutrients + Oxygen
# Source - WOA18
##########################
print('Compiling pelagic variables into '+outfile)

# Nitrate
print('Nitrate')
nci = netCDF4.Dataset('woa18/TRNN3_n_woa18-SEAsia_IC.nc')
dat = nci.variables['TRNN3_n'][:]
dat[dat<0] = 0
nco.variables['TRNN3_n'][:] = dat
nco.sync()
nci.close()
_,z,y,x = dat.shape

# Phosphate
print('Phosphate')
nci = netCDF4.Dataset('woa18/TRNN1_p_woa18-SEAsia_IC.nc')
dat = nci.variables['TRNN1_p'][:]
dat[dat<0] = 0
nco.variables['TRNN1_p'][:] = dat
nco.sync()
nci.close()

# Silicate
print('Silicate')
nci = netCDF4.Dataset('woa18/TRNN5_s_woa18-SEAsia_IC.nc')
dat = nci.variables['TRNN5_s'][:]
dat[dat<0] = 0
nco.variables['TRNN5_s'][:] = dat
nco.sync()
nci.close()

# Oxygen
print('Dissolved Oxygen')
nci = netCDF4.Dataset('woa18-oxy/TRNO2_o_woa18-SEAsia_IC.nc')
dat = nci.variables['TRNO2_o'][:]
dat[dat<0] = 0
nco.variables['TRNO2_o'][:] = dat
nco.sync()
nci.close()

################################################
# Set Alkalinity and Dissolved inorganic carbon
# Source - GLODAP
################################################

# Total Alkalinity
print('Total Alkalinity')
nci = netCDF4.Dataset('glodap/TRNO3_bioalk_glodap-SEAsia_IC.nc')
dat = nci.variables['TRNO3_bioalk'][:]
dat[dat<0] = 0
nco.variables['TRNO3_bioalk'][:] = dat
nco.sync()
nci.close()

print('DIC')
nci = netCDF4.Dataset('glodap/TRNO3_c_glodap-SEAsia_IC.nc')
dat = nci.variables['TRNO3_c'][:]
dat[dat<0] = 0
nco.variables['TRNO3_c'][:] = dat
nco.sync()
nci.close()

##########################
# Set Light Attenuation
# Source - PML
##########################
print('ADY') 
nci = netCDF4.Dataset('ady/TRNlight_ADY_ady-SEAsia_IC.nc')
dat = nci.variables['TRNlight_ADY'][0,:]
dat[dat<0] = 0
nco.variables['TRNlight_ADY'][:] = np.tile(dat,(z,1,1))[np.newaxis,:,:,:]
nco.sync()
nci.close()

##########################################################
# Set homogeneous fields - Ammonium & DOC
# Source - Yuri Artioli
########################################################

dat = nco.variables['TRNN3_n'][:] 
unit = np.ones((dat.shape))

print('Ammonium')
nco.variables['TRNN4_n'][:] = 0.25 * nco.variables['TRNN3_n'][:] 

print('Dissolved Organic Carbon')
nco.variables['TRNR1_c'][:] = 12 * unit
nco.variables['TRNR1_n'][:] = 0.151 * unit
nco.variables['TRNR1_p'][:] = 0.0094 * unit
print('Semi-labile Organic Carbon')
nco.variables['TRNR2_c'][:] = 60 * unit
print('Semi-refractory Organic Carbon')
nco.variables['TRNR3_c'][:] = 60 * unit

print('Calcite')
nco.variables['TRNL2_c'][:] = 0.1 * unit

print('Bacteria')
nco.variables['TRNB1_c'][:] =  3.7523 * unit # Check with Luca if this is a sensible value
nco.variables['TRNB1_n'][:] =  0.0167 * nco.variables['TRNB1_c'][:] # B1: qnc in fabm.yaml
nco.variables['TRNB1_p'][:] =  0.0019 * nco.variables['TRNB1_c'][:] # B1: qpc in fabm.yaml

nco.sync()

########################################################
# Set phytoplankton values
# Source - OC-CCI monthly image for Chloraphyll
#          Chloraphyll conversion from NEMO-ERSEM 
#          Redfield ratio for nitrogen, phosphate, silicate
########################################################


nci = netCDF4.Dataset('occci/TRNChl_Tot_CCI-OC-SEAsia_IC.nc')
dat = nci.variables['TRNChl_Tot'][:]
dat[dat<0] = 0
# Set uniform down to pycnocline, exponentially step down below
Chl_tot = np.tile(dat,(z,1,1))
pycno = np.argmax(np.diff(dens,axis=0),axis=0)
x = nco.dimensions['x'].size
for i in np.arange(x):
    for j in np.arange(y):
        k = pycno[j,i]
        p1 = Chl_tot[k,j,i]
        p2 = 0.0001
        Chl_tot[pycno[j,i]:,j,i] = p1*np.exp((np.log(p2/p1)/(z-k))*np.arange(1,z-k+1)) 
Chl_tot = Chl_tot[np.newaxis,:,:,:]

# Split Chl using Brewin2010 Eq13-16, parameterised using Brewin2012 Table 1
Chl_pn = 0.937*(1-np.exp(-1.033*Chl_tot))
Chl_dia = 2.0/3.0 * (Chl_tot - Chl_pn) #Split micro 2:1 into diatom and micro
Chl_mic = 1.0/3.0 * (Chl_tot - Chl_pn)
Chl_pic = 0.170*(1-np.exp(-4.804*Chl_tot))
Chl_nan = Chl_pn - Chl_pic 

print('Diatom Phytoplankton')
nco.variables['TRNP1_Chl'][:] = Chl_dia
carbon = 25*Chl_dia # Chl:C = 0.04
nco.variables['TRNP1_c'][:] = carbon
nco.variables['TRNP1_n'][:] = (carbon/12.0)*(16.0/106.0)
nco.variables['TRNP1_p'][:] = (carbon/12.0)*(1.0/106.0)
nco.variables['TRNP1_s'][:] = (carbon/12.0)*(15.0/106.0)

print('Nano Phytoplankton')
nco.variables['TRNP2_Chl'][:] = Chl_nan
carbon = 50*Chl_nan # Chl:C = 0.02
nco.variables['TRNP2_c'][:] = carbon
nco.variables['TRNP2_n'][:] = (carbon/12.0)*(16.0/106.0)
nco.variables['TRNP2_p'][:] = (carbon/12.0)*(1.0/106.0)

print('Pico Phytoplankton')
nco.variables['TRNP3_Chl'][:] = Chl_pic
carbon = 80*Chl_pic # Chl:C = 0.0125
nco.variables['TRNP3_c'][:] = carbon
nco.variables['TRNP3_n'][:] = (carbon/12.0)*(16.0/106.0)
nco.variables['TRNP3_p'][:] = (carbon/12.0)*(1.0/106.0)

print('Micro Phytoplankton')
nco.variables['TRNP4_Chl'][:] = Chl_mic
carbon = 30*Chl_mic # Chl:C = 0.03'
nco.variables['TRNP4_c'][:] = carbon
nco.variables['TRNP4_n'][:] = (carbon/12.0)*(16.0/106.0)
nco.variables['TRNP4_p'][:] = (carbon/12.0)*(1.0/106.0)

nco.sync()

########################################################
# Set zooplankton and POM values
# Source - iMarNet Model (Provided by LEDM)
########################################################

meso_zoo = ['TRNZ4_c']
micro_zoo = ['TRNZ5_c','TRNZ5_n','TRNZ5_p']
heteroflagellates = ['TRNZ6_c','TRNZ6_n','TRNZ6_p']
small_pom = ['TRNR4_c','TRNR4_n','TRNR4_p']
med_pom = ['TRNR6_c','TRNR6_n','TRNR6_p','TRNR6_s']
large_pom = ['TRNR8_c','TRNR8_n','TRNR8_p','TRNR8_s']
vars = meso_zoo + micro_zoo + heteroflagellates + \
        small_pom + med_pom + large_pom

print('Zooplankton and POM')
for var in vars:
    nci = netCDF4.Dataset('imarnet/{}_imarnet-SEAsia_IC.nc'.format(var))
    dat = nci.variables[var][:]
    dat[dat<0] = 0
    nco.variables[var][:] = dat
    nci.close()
nco.sync()

##########################################################

print('Filling Before Fields')
fields = list(filter(lambda x: x.startswith('TRN'), nco.variables.keys()))
for i in fields:
    nco.variables[i.replace('TRN','TRB')][:] = nco.variables[i][:]
nco.sync()

print('Complete')
nco.close()
