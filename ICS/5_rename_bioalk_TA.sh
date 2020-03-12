#!/usr/bin/bash

# Script to rename bioalk fields to TA

ncrename -v TRNO3_bioalk,TRNO3_TA bgc_ini.nc
ncrename -v TRBO3_bioalk,TRBO3_TA bgc_ini.nc


