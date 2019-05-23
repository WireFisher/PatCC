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
	for line in  open(i,"r"):
		line_split=line.split(" ")
		if(len(line_split)>1):
			if(line_split[0] not in result.keys()):
				result[line_split[0]]=[]
			result[line_split[0]].append((i,line_split[1]))
csv_write = csv.writer(open('differ.csv','w'),dialect="excel")
for i in result.keys():
	for j in result[i]:
		parameter = str(j[0][17:]).split("_")
		csv_write.writerow([i,"\-O"+parameter[0],parameter[1],parameter[2],j[0],j[1].strip("\n")])
# get the result
# next get the md5 count

count_dict={}
for i in result.keys():
	if( i not in count_dict.keys()):
		count_dict[i] = { }
	for j in result[i]:
		if(j[1].strip("\n") not in count_dict[i].keys()):
			count_dict[i][j[1].strip("\n")] = 0
		count_dict[i][j[1].strip("\n")] = count_dict[i][j[1].strip("\n")] + 1

csv_write.writerow([])
csv_write.writerow(["below is difference"])
for i in count_dict.keys():
	for j in count_dict[i].keys():
		csv_write.writerow([i,j,count_dict[i][j]])

