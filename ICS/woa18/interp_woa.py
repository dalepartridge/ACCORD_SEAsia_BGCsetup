import fileinput
import os
import re
import sys

def create_namelists(template_dir,dat):
    templates = ['1_initcd_source_to_source_var.namelist.template',
                 '2_source_weights_var.namelist.template',
                 '3_initcd_source_to_nemo_var.namelist.template']
    namelists = ['1_initcd_{}_to_{}_{}.namelist'.format(dat['SOURCEID'],dat['SOURCEID'],dat['VAR']),
                 '2_{}_weights_{}.namelist'.format(dat['SOURCEID'],dat['VAR']),
                 '3_initcd_{}_to_nemo_{}.namelist'.format(dat['SOURCEID'],dat['VAR'])]
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
        'TAG': 'IC',
        'DOMAIN': 'domain_cfg.nc'
       }

vars = {'n_an': {'SFILE': 'woa18_nitrate_extracted_mean.nc',
               'VARL': 'Nitrate',
               'OVAR': 'TRNN3_n',
               'SCALE': '1.0'
               },   
        'i_an': {'SFILE': 'woa18_silicate_extracted_mean.nc',
               'VARL': 'Silicate',
               'OVAR': 'TRNN5_s',
               'SCALE': '1.0'
               },   
        'p_an': {'SFILE': 'woa18_phosphate_extracted_mean.nc',
               'VARL': 'Phosphate',
               'OVAR': 'TRNN1_p',
               'SCALE': '1.0'
               }}

for i,v in enumerate(vars):
    create_namelists(sys.argv[1],{'VAR':v, **data, **vars[v]})
    
    if i == 0:
        os.system('sh ./interp_IC_initial.sh {} {} {} {}'.format(v, vars[v]['SFILE'], data['SOURCEID'], data['STIMEVAR']))
    else:
        os.system('sh ./interp_IC_additional.sh {} {} {}'.format(v, vars[v]['SFILE'], data['SOURCEID']))
