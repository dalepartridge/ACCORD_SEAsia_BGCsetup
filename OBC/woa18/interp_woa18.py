import fileinput
import os
import re
import sys

def create_namelists(template_dir,dat):
    templates = ['1_initcd_source_to_source_var.namelist.template',
                 '2_source_weights_var.namelist.template',
                 '3_initcd_source_to_nemo_var.namelist.template']
    namelists = ['1_{}_to_{}_{}.namelist'.format(dat['SOURCEID'],dat['SOURCEID'],dat['VAR']),
                 '2_{}_weights_{}.namelist'.format(dat['SOURCEID'],dat['VAR']),
                 '3_{}_to_nemo_{}.namelist'.format(dat['SOURCEID'],dat['VAR'])]
    for t,n in zip(templates,namelists):
        with open(template_dir+t) as fin, \
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
        'SOURCEID':'woa18',
        'STIMEVAR':'time',
        'SLONVAR':'lon',
        'SLATVAR':'lat',
        'SZVAR':'depth',
        'TARGETID': 'SEAsia',
        'TAG': 'OBC',
        'DOMAIN': 'domain_cfg.nc'
       }

vars = {'i_an': {'SFILE': 'ext_woa18_silicate.nc',
               'VARL': 'Silicate',
               'OVAR': 'SIL',
               'SCALE': '1.0'
               },   
        'n_an': {'SFILE': 'ext_woa18_nitrate.nc',
               'VARL': 'Nitrate',
               'OVAR': 'NO3',
               'SCALE': '1.0'
               },   
        'p_an': {'SFILE': 'ext_woa18_phosphate.nc',
               'VARL': 'Phosphate',
               'OVAR': 'PO4',
               'SCALE': '1.0'
               }}

for i,v in enumerate(vars):
    print('Processing {}'.format(v))
    create_namelists(sys.argv[1],{'VAR':v, **data, **vars[v]})
    if i == 0:
        os.system('sh ./interp_OBC_initial.sh {} {} {} {}'.format(v, vars[v]['SFILE'], data['SOURCEID'], data['STIMEVAR']))
    else:
        os.system('sh ./interp_OBC_additional.sh {} {} {}'.format(v, vars[v]['SFILE'], data['SOURCEID']))
