import os
import csv
files=os.listdir('./')

checksumlist=[]
for file in files:
	if(file.find("all_grid_checksum")!=-1):
		checksumlist.append(file)
result={}
for i in checksumlist:
	parameter = str(i[17:]).split("_")
	for line in open(i,"r"):
		line_split=line.split(" ")
		line_split[0] = line_split[0].strip("\n")
		if(line_split[0] not in result.keys()):
			result[line_split[0]]=[]
		if(len(line_split)>1):
			result[line_split[0]].append((i,line_split[1]))
		else:
			result[line_split[0]].append((i,"-1"))

csv_write = csv.writer(open('differ.csv','w'),dialect="excel")
for i in result.keys():
	for j in result[i]:
		parameter = str(j[0][17:]).split("_")
		csv_write.writerow([i,"O"+parameter[0],parameter[1],parameter[2],j[0],j[1].strip("\n")])
# get the result
# next get the md5 count

count_dict={}
for i in result.keys():
	if( i not in count_dict.keys()):
		count_dict[i] = [ {},{},{},{} ]
		# for four optimization
	for j in result[i]:
		str_temp = j[1].strip("\n")
		if(str_temp not in count_dict[i][int(j[0][17:18])].keys() ):
			count_dict[i][int(j[0][17:18])][str_temp] = 0
		count_dict[i][int(j[0][17:18])][str_temp] = count_dict[i][int(j[0][17:18])][str_temp] + 1


csv_write = csv.writer(open('differ_with_optimize.csv','w'),dialect="excel")
csv_write.writerow(["grid_name","md5","optimize","number","total"])
for i in count_dict.keys():
	temp = 0
	for k in range(0,4):
		for j in count_dict[i][k].keys():
			temp+=count_dict[i][k][j]
	for k in range(0,4):
		for j in count_dict[i][k].keys():
			csv_write.writerow([i,j,k,count_dict[i][k][j],temp])

