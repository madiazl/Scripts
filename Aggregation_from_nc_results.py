import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import os
from xlrd import open_workbook

os.chdir(os.getcwd())
# os.chdir(r'C:\Glofris')


def func_aggregation (out_list,geogunit_list,geogunit_idnums):
        #agregation dependig on th condition of the agregation
        impact_list_index=np.where(out_list>0); # hacerlo como condicion revisar los tamanyos de los arreglos
        impact_list_index=np.array(impact_list_index);
        impact_list=out_list[impact_list_index];
        impact_list_geogunit = geogunit_list[impact_list_index];
        out=np.zeros((len(geogunit_idnums),2));
        out[:,0]=np.array(geogunit_idnums.astype(int));
        for m in range (0,len(geogunit_idnums)):
                auxiliar=np.where(impact_list_geogunit==geogunit_idnums[m]);
                out[m,1]=np.nansum(impact_list[auxiliar]);
        # Carry out analysis...
        return out;

table = []
print(os.getcwd(), table)

for file in os.listdir(os.getcwd()):
    if file.endswith(".h5"):
        print(os.path.join(os.getcwd(), file))
        table.append(file)
        # table=np.vstack((table,np.array([file])))

print table, type(table), table[0]

# file_pop = nc.Dataset(r'/home/glofris1/IWC/Risk/example_exp_april/exposure_base_2010_pop.h5', 'r', format='h5'); #old Hyde land uses
# pop = file_pop.variables['data'];
# pop.set_auto_maskandscale(False);
#
geogunits=nc.Dataset(r'C:\Glofris\WRI\geogunit_14.nc', 'r', format='netCDF4');
data=geogunits.variables['data']
data.set_auto_maskandscale(False)
print data
list_geog=(data[:][:]>=0)
print list_geog
data_list=np.extract(list_geog,data[:][:])
print(data_list)



book = open_workbook(r'C:\Glofris\WRI\geogunit_14_list.xlsx', on_demand=True)
book_sheet = book.sheet_by_index(0)
geogunitxls = [[book_sheet.cell_value(r,col)
for col in range(book_sheet.ncols)]
	for r in range(book_sheet.nrows)]
geogunitxls= np.array(geogunitxls)
geogunit_idnums=np.array(geogunitxls[:,0],dtype=np.float64)
geogunit_idnames=np.array(geogunitxls[:,1],dtype=np.str)
print geogunit_idnums, geogunit_idnames
print (len(geogunit_idnums))

Table = np.array([['','']])
for i in range(0,len(table)):
        print table[i]
	np.savetxt('exa.txt',np.asarray([i]),fmt="%s")
        a = nc.Dataset(table[i], 'r',format='h5')
	print a
        var = a.variables['data']
        #var.set_auto_maskandscale(False)
        var_list=np.extract(list_geog,var[:][:])
        a.sync()
        a.close()
        del(var,a)
        out=func_aggregation(var_list,data_list,geogunit_idnums)
        name=table[i]
        if name.endswith("gdp.h5"):
                np.savetxt(name.replace('GDP_Exposed_geogunit_14.txt','GDP_Exposed_geogunit_14.txt'),(out[:]),delimiter='\t',fmt='%f')
        if name.endswith("pop.h5"):
                np.savetxt(name.replace('POP_Exposed_geogunit_14.txt','POP_Exposed_geogunit_14.txt'),(out[:]),delimiter='\t',fmt='%f')

#       for j in range(np.int(geogunit_idnums[0]),len(geogunit_idnums)):
#               print j
#               condition=(data[:][:]==j)
#               cells_a=np.extract(condition,var[:][:]);
#               impact_list_index=np.where(cells_a>0); # hacerlo como condicion revisar los tamanyos de los arreglos
#               impact_list_index=np.array(impact_list_index);
#               impact_list=cells_a[impact_list_index];
#               value=np.nansum(impact_list)
#               Table=np.vstack((Table, np.array([j,value])))
#               del(condition,cells_a,impact_list,impact_list_index,value)
#       name=table[i]
#       print Table
#       np.savetxt(name.replace('.nc','Urban_Damage_geogunit_14.txt'),(Table),delimiter='\t',fmt='%f')
#       #np.savetxt(name.replace,(Table),delimiter='\t',fmt='%s')
#       a.sync()
#       a.close()

# SELECTING THE CELLS OF germany code: DEU number 85
# Table = np.array([['FID','sum_pop','mean_pop_by_cell']])


        # for i in range(0,253):
        #       print(i)
        #       #cell_condition = np.where(data[:][:]==i)
        #       cell_condition = (data[:][:]==i)
        #       pop_by_country = np.extract(cell_condition,pop[:][:])
        #       sumpop=np.nansum(pop_by_country)
        #       meanpop=np.nanmean(pop_by_country)
        #       Table=np.vstack((Table, np.array([i,sumpop,meanpop])))
        # np.savetxt(r'/home/glofris1/IWC/Risk/example_exp_april/Population.txt',(Table),delimiter='\t',fmt='%s')
        #
        #         #or only netherlands
        #       cell_condition = (data[:][:]==157)
        #       pop_by_country = np.extract(cell_condition,pop[:][:])
        #       sumpop=np.nansum(pop_by_country)
        #       meanpop=np.nanmean(pop_by_country)
        #         print(sumpop,meanpop)
