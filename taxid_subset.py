#!/usr/bin/python
#script filters out all irrelevant species form centrifuge's results
#we only care about fungi under the kingdom of eukaryotes and we do not care about synthetic constructs 
from ete3 import NCBITaxa
#ncbi lib
new_taxid_list=[]
ncbi = NCBITaxa()
#open ip and op files
new_tax_id_to_be_aligned = open('new_tax_id_to_be_aligned', 'w')
with open("tax_id_to_be_aligned") as f:
  taxid_list = f.readlines()

#check superkingdom for each taxid
# we only care abt fungi under eukaryotes. we also do not want synthetic dna constructs
for i in range(0,len(taxid_list)):
 try:
  lineage=ncbi.get_lineage(taxid_list[i])
  superkingdom=lineage[2]
  if int(superkingdom) != 2759 and int(lineage[1])!=28384:
   new_taxid_list.append(taxid_list[i])   
   new_tax_id_to_be_aligned.write(taxid_list[i])
  elif int(lineage[4])== 4751 and int(lineage[1])!=28384:
   new_taxid_list.append(taxid_list[i])   
   new_tax_id_to_be_aligned.write(taxid_list[i])
 except:
  #new_tax_id_to_be_aligned.write(taxid_list[i])
  print("not found")
new_tax_id_to_be_aligned.close()

