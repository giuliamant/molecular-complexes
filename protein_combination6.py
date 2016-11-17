from datetime import datetime
startTime = datetime.now()
import csv
import os
import sys
import pandas as pd
import plotly.plotly as py
import plotly.graph_objs as go
from textwrap import wrap
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns



os.chdir(sys.path[0])
csv.register_dialect(
    'mydialect',
    delimiter = '\t',
    quotechar = '"',
    doublequote = True,
    skipinitialspace = True,
    lineterminator = '\r\n',
    quoting = csv.QUOTE_MINIMAL)



#################################### Prepare the BIG MATRIX ####################################


oneh=open('onehour.csv','wb')		
import pandas as pd
df = pd.read_csv(sys.path[0]+"complex_protein_matrixS5.csv",delimiter='\t')

tag_to_pro={'TTTATT':'CTCF','CTCCTC':'RAD21','TGTACC':'SMC1','CTCGCT':'SMC3','ACCGAC':'Nipbl','GCAAAT':'E2F','AGTCCT':'Max','GAGGGC':'Myc','CCCTTA':'new_CTCF','GGCAAC':'new_Myc','GTACTT':'new_SMC1'}  


header_column=list(df.columns.values)


protein_values=[]
for i in range (1,len(header_column)-1):
	protein_values.append((df.iat[0,i],i))

protein_values.sort(reverse=True)

sel_protein=protein_values[0:11]


index=[x[1] for x in sel_protein]
index.sort()


col_names=[]
for i in index:
	col_names.append(header_column[i])


df1 = df[col_names]
df2 = df1.ix[1:]
#print df2
row_n=df2.shape[0]
col_n=df2.shape[1]

comb_dict={}
for r in range (row_n):
 	cou=0
 	comb=[]
 	for c in range (col_n):
		
		value=df2.iloc[r,c]
		if value != 0:
			
			cou+=1
			protein=col_names[c]
			prop=tag_to_pro[protein]
			comb.append(prop)
			
	
	tt=tuple(comb)
	if tt in comb_dict.keys():
		comb_dict[tt]+=1
	else:
		comb_dict[tt]=1


for w in sorted(comb_dict, key=comb_dict.get, reverse=True):
	oneh.write(','.join(''.join(s) for s in w))
	oneh.write('\t')
	oneh.write(str(comb_dict[w]))
	oneh.write('\n')
oneh.close()


############################################### PROTEINS COMBINATION PLOTS ######################################

data_chart={}
for i in range (1,12):
	data_chart[i]=[]
with open(sys.path[0]+'onehour.csv', 'rb') as mycsvfile:
	thedata = csv.reader(mycsvfile, dialect='mydialect')
	for row in thedata:
		a=row[0].split(',')
		
		#data_chart[len(a)].append(row)
		
		if len(data_chart[len(a)]) <= 50:
			#print len(data_chart[len(a)])
			data_chart[len(a)].append(row)
			#pass
	
for x in data_chart.keys():
	plot_dict={}
	for i in data_chart[x]:
		plot_dict[i[0]]=int(i[1])
	

	
	if x == 2:
		
		new_dict={}
		for i in data_chart[x]:
			
			k=i[0].split(',')
			
			new_dict[k[0],k[1]]=float(i[1])
				
		ser = pd.Series(list(new_dict.values()),
                  index=pd.MultiIndex.from_tuples(new_dict.keys()))
		df = ser.unstack().fillna(0)
		df.shape	
		
		
		

		
		sns.heatmap(df, annot=True, fmt="g", cmap='viridis')
		plt.show()
	
	

	
	
	
	
print datetime.now() - startTime




























