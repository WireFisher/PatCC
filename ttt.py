import os
import re
import csv
from operator import itemgetter
result=[]
#result.append(['compiler','optimize','process','thread','status','test case','time','failed grid'])
list = os.listdir('.')
for file_name in list:
	split_name=file_name.split('.')
	if(len(split_name)==4 and split_name[3]=='out'):
		complie_name=split_name[0]+'.'+split_name[1]+'.'+split_name[3]
		complier=""
		with open(complie_name) as f:
			lines = f.readline()
			line = lines.split(" ")
			complier=line[0]+line[1]+line[2]

		optimize_degree=int(split_name[0][6])
		process_num = int(split_name[1])
		thread_num = int(split_name[2])
		with open(file_name, 'r') as f:
			once_result=[]
			last_result=0
			lines = f.readlines()
			for line in lines:
				if(line.find('[')!=-1 and (line.find('OK')!=-1 or line.find('FAILED')!=-1)):
					if(line.find('(')!=-1):
						result1 = line.replace(' ','')
						result2 = re.split(r'[\[\]\(\)]',result1)
						once_result.append(complier)
						once_result.append('\-O'+str(optimize_degree))
						once_result.append(process_num)
						once_result.append(thread_num)
						once_result.append(result2[1])
						once_result.append(result2[2])
						once_result.append(result2[3])
						if(line.find('FAILED')!=-1):
							error_message=""
							for ii in range(last_result,lines.index(line)):
								if(lines[ii].find('Failure')!=-1):
									if(lines[ii-1].find(':')!=-1):
										error_message=error_message+','+str.strip(lines[ii-1].split(':')[1])
							once_result.append(error_message)
						else:
							once_result.append('all pass')
						result.append(once_result)
						once_result=[]
						last_result=lines.index(line)
table_sorted = sorted(result, key=itemgetter(2, 5))
csv_writer=csv.writer(open('all_result.csv','w'),dialect='excel')
csv_writer.writerows(table_sorted)

# write the result to csv
dict={ }
for line in result:
	if(line[5] not in dict.keys()):
		dict[line[5]]={}
	if(line[1] not in dict[line[5]].keys()):
		dict[line[5]][line[1]]=[]	
	dict[line[5]][line[1]].append([line[2],line[3],line[4]])
late_result=[]
for key in dict.keys():
	for key2 in dict[key].keys():
		dict[key][key2].sort()
		temp_result=[key,key2]
		for i in range(0,len(dict[key][key2])):
			temp_result.append(dict[key][key2][i][2])
		late_result.append(temp_result)

csv_writer=csv.writer(open('late_result.csv','w'),dialect='excel')
csv_writer.writerows(late_result)
