from netCDF4 import Dataset

a = Dataset("sis_ice_grid@sis.nc")
for key in a.variables.keys():
	if(key.find("lon")!=-1 and key.find("corner")==-1):
		lon = a[key][:]
	if(key.find("lat")!=-1 and key.find("corner")==-1):
		lat = a[key][:]

print(max(lon))
print(min(lon))
print(max(lat))
print(min(lat))
