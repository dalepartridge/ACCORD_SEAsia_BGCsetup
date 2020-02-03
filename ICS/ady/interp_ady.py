import fileinput
import os
import re

def create_namelists(template_dir,dat,TWOD=False):
    templates = ['1_initcd_source_to_source_var.namelist.template',
                 '2_source_weights_var.namelist.template',
                 '3_initcd_source_to_nemo_var.namelist.template']
    namelists = ['1_initcd_{}_to_{}_{}.namelist'.format(dat['SOURCEID'],dat['SOURCEID'],dat['VAR']),
                 '2_{}_weights_{}.namelist'.format(dat['SOURCEID'],dat['VAR']),
                 '3_initcd_{}_to_nemo_{}.namelist'.format(dat['SOURCEID'],dat['VAR'])]
    TWOD_str = '2D' if TWOD else ''
    for t,n in zip(templates,namelists):
        with open(template_dir+t+TWOD_str) as fin, \
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
        'SOURCEID':'ady',
        'STIMEVAR':'day',
        'SLONVAR':'lon',
        'SLATVAR':'lat',
        'SZVAR':'',
        'TARGETID': 'SEAsia',
        'TAG': 'IC',
        'DOMAIN': 'domain_cfg_ORCA12_adj.nc'
       }

vars = {
        'ady': {'SFILE': 'adyBroadBandClimatology_ACCORD_extracted.nc',
                 'VARL': 'Light Attenuation',
                 'OVAR': 'TRNlight_ADY',
                 'SCALE': '1.0'
                },
        }

t = 'snr' if data['STIMEVAR']=='missing_rec' else data['STIMEVAR']

interp_dir='/work/n01/n01/dapa/INTERP/BGC/interp-files'
domain_file='/work/n01/n01/dapa/INTERP/BGC/'+data['DOMAIN']

#LINK TO FILES
os.system('ln -s {}/interp*.sh .'.format(interp_dir))
os.system('ln -s {} .'.format(domain_file))

for i,v in enumerate(vars):
    create_namelists(interp_dir+'/namelist-templates/',{'VAR':v, **data, **vars[v]},TWOD=True)
    
    if i == 0:
        os.system('sh ./interp_main.sh {} {} {} {}'.format(v, vars[v]['SFILE'], data['SOURCEID'], t))
    else:
        os.system('sh ./interp_additional.sh {} {} {}'.format(v, vars[v]['SFILE'], data['SOURCEID']))
