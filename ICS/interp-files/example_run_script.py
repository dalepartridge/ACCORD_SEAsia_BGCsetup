import fileinput
import os
import re

def create_namelists(dat):
    templates = ['1_initcd_source_to_source_var.namelist.template',
                 '2_source_weights_var.namelist.template',
                 '3_initcd_source_to_nemo_var.namelist.template']
    namelists = ['1_initcd_{}_to_{}_{}.namelist'.format(dat['SOURCEID'],dat['SOURCEID'],dat['VAR']),
                 '2_{}_weights_{}.namelist'.format(dat['SOURCEID'],dat['VAR']),
                 '3_initcd_{}_to_nemo_{}.namelist'.format(dat['SOURCEID'],dat['VAR'])]
    for t,n in zip(templates,namelists):
        with open(t) as fin, \
             open(n,'w') as fout:
            for l in fin.readlines():
                a = l 
                for d in dat:
                    a = re.sub('__'+d+'__',dat[d],a)
                fout.write(a)
            fin.close()
            fout.close()
    return

data = {
        'SOURCEID':'glodap',
        'STIMEVAR':'missing_rec',
        'SLONVAR':'lon',
        'SLATVAR':'lat',
        'SZVAR':'Depth',
        'TARGETID': 'SEAsia',
        'TAG': 'IC',
        'DOMAIN': 'domain_cfg_ORCA12_adj.nc'
       }

vars = {
        'TAlk': {'SFILE': 'GLODAP_TAlk_extracted.nc',
                 'VARL': 'Total Alkalinity',
                 'OVAR': 'TRNO3_bioalk',
                 'SCALE': '1.0'
                },
        'TCO2': {'SFILE': 'GLODAP_TCO2_extracted.nc',
                 'VARL': 'DIC',
                 'OVAR': 'TRNO3_c',
                 'SCALE': '1.0'
                }
        }

t = 'snr' if data['STIMEVAR']=='missing_rec' else data['STIMEVAR']

for i,v in enumerate(vars):
    create_namelists({'VAR':v, **data, **vars[v]})
    if i == 0:
        os.system('sh ./interp_main.sh {} {} {} {}'.format(v, vars[v]['SFILE'], data['SOURCEID'], t))
    else:
        os.system('sh ./interp_additional.sh {} {} {}'.format(v, vars[v]['SFILE'], data['SOURCEID']))
