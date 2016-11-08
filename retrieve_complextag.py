# coding=utf-8
from datetime import datetime
startTime = datetime.now()
import random
import re
import numpy as np
import os
import subprocess
os.system("cat *R1* > read1.fastq.gz")
os.system("cat *R2* > read2.fastq.gz")

print 'number of initial sequences:'
print os.system("echo $(zcat <read1.fastq.gz | wc -l ) / 4 | bc")

os.system("gunzip read1.fastq.gz")
os.system("gunzip read2.fastq.gz")

################### Fastx_Toolkit commands #####################


 #remove redundancy
os.system("fastx_collapser -i read1.fastq -o read1collapsed.fasta -Q33")

#from fasta to txt

os.system("grep -v '^>' read1collapsed.fasta > read1collapsed.txt")
print 'sequences with no redundancy:'
print os.system("wc -l read1collapsed.txt")
#Preparing input file for analysis

os.system("grep CATTACTAGGAATCACACGC read1collapsed.txt > proteins.txt")
 
os.system("grep ACCTGAGACATCATAATAGCA proteins.txt > proteins_with_extension.txt")


file=open('proteins_with_extension.txt')

diz_compl={}

biglist=[]
ll=0.0
difpath=0.0

ptaglist=['CCTAAA','TGGAAT','TTGAGC','GGGATG','GTCCAA','ACACAG','ATGTCC','ACTTTG','CCGCTT','TTTATT','CTCCTC','TGTACC','CTCGCT','ACCGAC','GCAAAT','AGTCCT','GAGGGC','TTACGA','GTCTGG','TACATA','TATTGT','AAGCAT','TAAACT','CTTGTG','GTTGAT','CGACTA','TCTCTA','CGGGAC','GAAGTA','AACTTT','CCCTTA','GGCAAC','GTACTT','AAGATC','AGTTAA','ATTCGC','TAGCGG','GCTTCA','CTGTAT','CAGAAG','GATCTC','CGTCGG','ATGGAG','CGTTTT','ATGCTA','ACCCGA','TACCAC','TCTAGG']

molecule={}

complextag={}
for line in file:
	ll+=1
	ext=line[15:36]
	primer=re.search(r"ACCTGAGACATCATAATAGCA",line)
	
	
	UniPro=re.search(r"CATTACTAGGAATCACACGC",line)
	
	
	
	innerpath=line[primer.end():UniPro.start()]
	compl_tag=line[primer.start()-15:primer.start()]
	#print len(compl_tag)
	
	if len(innerpath) != 16:
		difpath+=1
	else:
		
		if len(compl_tag) == 15:
			mol_tag=innerpath[0:10]
			
			pro_tag=innerpath[10::]
			if pro_tag in ptaglist:#print len(compl_tag)
				complextag[compl_tag]=''
				biglist.append(line)
				
				molecule[(compl_tag,pro_tag,mol_tag)]=1
				








######### creo la big_matrix ###############




big_matrix=[]
for i in complextag:
	inside=[]
	for j in range(49):
		inside.append(0)
	big_matrix.append(inside)


compllist=[]


for i in sorted(complextag):
	compllist.append(i)
print 'complex tag list length', len(compllist)

		
for j in range (0,len(big_matrix)):
	big_matrix[j][0]=compllist[j]


header=['complex_tag\protein_tag']
for i in ptaglist:
	header.append(i)

big_matrix.insert(0, header)  



A=np.array(big_matrix)


compllist.insert(0,'')

# print A.shape
# print A[:,0]

################ riempio la big_matrix cercando complex_tag(y) and protein_tag(x) ##########	

# biglist2=biglist[0:12]
# for line in biglist2:
# 	ll+=1
# 	ext=line[15:36]
# 	primer=re.search(r"ACCTGAGACATCATAATAGCA",line)
# 	
# 	
# 	UniPro=re.search(r"CATTACTAGGAATCACACGC",line)
# 	
# 	
# 	
# 	innerpath=line[primer.end():UniPro.start()]
# 	compl_tag=line[primer.start()-15:primer.start()]
# 	pro_ind=0
# 	compl_ind=0
# 	
# 	
# 		
# 	mol_tag=innerpath[0:10]
# 		
# 	pro_tag=innerpath[10::]
# 			
# 	
# 	
# 	pro_ind=big_matrix[0].index(pro_tag)
# 	compl_ind=compllist.index(compl_tag)
# 	value=molecule[(compl_tag,pro_tag,mol_tag)]
# 		
# 				
# 	A[compl_ind][pro_ind].fill(value)
# #  	
# print A
	
########################## checking ######################


# molpetit={}
# for k,v in molecule.items():
# 	if k[0]=='TCCGCAGTAACGCTC':
# 		if k[1] == 'GGCAAC':
# 			molpetit[k]=v
#print molpetit	

##################################### set columns dataframe ###################
	

############################# big_matrix changes to fit in heatmap track ############################			
import pandas as pd
from pandas import DataFrame 
#pd.options.mode.chained_assignment = None
file_name='complex_protein_matrix.csv'
# A sparse matrix in dictionary form (can be a SQLite database). Tuples contains doc_id        and term_id. 
#doc_term_dict={('d1','t1'):12, ('d2','t3'):10, ('d3','t2'):5}

#extract all unique documents and terms ids and intialize a empty dataframe.

rows = set([c for (c,p,m) in molecule.keys()])  
cols = set([p for (c,p,m) in molecule.keys()])
df = DataFrame(index = rows, columns = cols )
df = df.fillna(0)
cols = df.columns.tolist()
index=[]
for i in cols:
	ind=ptaglist.index(i)
	index.append(ind)
propro_sorted=[x for (y,x) in sorted(zip(index,cols))]
#print cols
#print propro_sorted


#assign all nonzero values in dataframe
for key, value in molecule.items():
    df[key[1]][key[0]] += value

#sort the columns order by the ptaglist
#df1 = df[propro_sorted]
#df1.is_copy = False

# sum of molecules per complex tag
df['Total Molecule per Complex Tag'] = df.sum(axis=1)
#df1.loc[:,'Total Molecule per Complex Tag']*= 10 

################# sum of molecule per protein  ###########################

#Sum the columns:
sum_row = {col: df[col].sum() for col in df}
#Turn the sums into a DataFrame with one row with an index of 'Total':
sum_df = pd.DataFrame(sum_row, index=["Total Molecule per Protein"])
#Now append the row:
df = df.append(sum_df)

df2=df.sort_values(['Total Molecule per Complex Tag'], ascending=False)
propro_sorted.append('Total Molecule per Complex Tag')
df3 = df2[propro_sorted]
df3.to_csv(file_name, sep='\t')


	

	
	
	
	
	
	
	
print datetime.now() - startTime
	