import netCDF4 as nc
import numpy as np
import os
from xlrd import open_workbook



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

os.chdir(os.getcwd())
# os.chdir(r'C:\Glofris')

table = []
print(os.getcwd(), table)

for file in os.listdir(os.getcwd()):
    if file.endswith(".nc"):
        print(os.path.join(os.getcwd(), file))
        table.append(file)
        # table=np.vstack((table,np.array([file])))

print table, type(table), table[0]

geogunits=nc.Dataset(r'/home/adiaz/IWC/geogunit_21.nc', 'r', format='netCDF4');
data=geogunits.variables['data']
data.set_auto_maskandscale(False)
print data
list_geog=(data[:][:]>=0)
print list_geog
data_list=np.extract(list_geog,data[:][:])
print(data_list)

book = open_workbook(r'/home/adiaz/IWC/geogunit_21_list.xlsx', on_demand=True)
book_sheet = book.sheet_by_index(0)
geogunitxls = [[book_sheet.cell_value(r,col)
for col in range(book_sheet.ncols)]
        for r in range(book_sheet.nrows)]
geogunitxls= np.array(geogunitxls)
geogunit_idnums=np.array(geogunitxls[:,0],dtype=np.float64)
geogunit_idnames=np.array(geogunitxls[:,1],dtype=np.str)
print geogunit_idnums, geogunit_idnames
print (len(geogunit_idnums))

for i in range(0,len(table)):
        print table[i]
        a = nc.Dataset(table[i], 'a', format='NETCDF4')
        var = a.variables['Risk_Results']
        var.set_auto_maskandscale(False)
        var_list=np.extract(list_geog,var[:][:])
        a.sync()
        a.close()
        del(var,a)
        out=func_aggregation(var_list,data_list,geogunit_idnums)
        name=table[i]
        if name.endswith("GDPexp.nc"):
                np.savetxt(name.replace('GDPexp.nc','GDP_Exposed_geogunit_21.txt'),(out[:]),delimiter='\t',fmt='%f')
        if name.endswith("POPexp.nc"):
                np.savetxt(name.replace('POPexp.nc','POP_Exposed_geogunit_21.txt'),(out[:]),delimiter='\t',fmt='%f')
        if name.endswith("Urban_Damage.nc"):
                np.savetxt(name.replace('Urban_Damage.nc','Urban_Damage_geogunit_21.txt'),(out[:]),delimiter='\t',fmt='%f')

