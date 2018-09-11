import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import os
os.chdir(os.getcwd())

import os
for file in os.listdir("/mydir"):
    if file.endswith(".txt"):
        print(os.path.join("/mydir", file))



file_pop = nc.Dataset(r'/home/glofris1/IWC/Risk/example_exp_april/exposure_base_2010_pop.h5', 'r', format='h5'); #old Hyde land uses
pop = file_pop.variables['data'];
pop.set_auto_maskandscale(False);

geogunits=nc.Dataset(r'/home/glofris1/IWC/Risk/example_exp_april/geogunit_2.h5', 'r', format='h5');
data=geogunits.variables['data']
data.set_auto_maskandscale(False)
# SELECTING THE CELLS OF germany code: DEU number 85


Table = np.array([['FID','sum_pop','mean_pop_by_cell']])


for i in range(0,253):
	print(i)
	#cell_condition = np.where(data[:][:]==i)
	cell_condition = (data[:][:]==i)
	pop_by_country = np.extract(cell_condition,pop[:][:])
	sumpop=np.nansum(pop_by_country)
	meanpop=np.nanmean(pop_by_country)
	Table=np.vstack((Table, np.array([i,sumpop,meanpop])))
np.savetxt(r'/home/glofris1/IWC/Risk/example_exp_april/Population.txt',(Table),delimiter='\t',fmt='%s')

        #or only netherlands
	cell_condition = (data[:][:]==157)
	pop_by_country = np.extract(cell_condition,pop[:][:])
	sumpop=np.nansum(pop_by_country)
	meanpop=np.nanmean(pop_by_country)
        print(sumpop,meanpop)