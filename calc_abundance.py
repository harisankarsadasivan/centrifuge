#~/usr/bin/python

#code to analyze the aliIgnment results and to calculate the final abundance numbers.
#Input argument is the path to the NCBI downloaded database
import os
import sys
import sys

taxid_to_seqid={}
align_ctr={}
tax_ids=[]
rootdir = sys.argv[1]
sum_abd=0
#find all seqids that belong to a single taxid
for dirName, subdirList, fileList in os.walk(rootdir):
    #print('Found directory: %s' % dirName)
    for fname in fileList:
	if fname.endswith(".fa"):
		with open(dirName+"/"+fname) as f:
			for line in f:
				if line[0]=='>':
					tmp_var=line.split(None,1)[0][1:]
					#extracting seqid belonging to each taxid
					taxid_to_seqid[tmp_var]=fname[:-3]
					tax_ids.append(fname[:-3])
					align_ctr[fname[:-3]]=0
					print("Identified mapped reference ", tmp_var,"corresponding to taxid", taxid_to_seqid[tmp_var])


#analyze the alignment results
with open("/home/dna/hariss/minimap2/cumulative_alignment.sam") as f:
	for line in f:
		words=line.split('\t',11)
		if words[1]!=4 and words[2] in taxid_to_seqid.keys():
			align_ctr[taxid_to_seqid[words[2]]]=align_ctr[taxid_to_seqid[words[2]]]+1
			
#for key in align_ctr.keys():
#	print(align_ctr[key])

f=open("/home/dna/hariss/centrifuge_new/centrifuge/avg_rd_len.txt")
avg_rd_len=f.readline()
f.close()

f=open("abundance_results.txt",'w')
#add onto centrifug'e results
print(":::::::::::::::::::GENOMIC ABUNDANCE (in %)::::::::::::::::::")
print >>f, ":::::::::::::::::::GENOMIC ABUNDANCE (in %)::::::::::::::::::"
with open("/home/dna/hariss/centrifuge_new/centrifuge/centrifuge_report.tsv") as fh:
	for line in fh:
                 words=line.split('\t',9)
		 if words[1] in tax_ids:
			tmp=(float(words[5])*float(words[7])+ float(align_ctr[words[1]]*float(avg_rd_len) ))/ float(words[3])
			sum_abd=sum_abd+tmp
		 	print("Genomic abundance(in %) of ",words[0], "taxid:", words[1], ":", round(tmp*100,2))
			print >>f, "Genomic abundance(in %) of ",words[0], "taxid:", words[1], ":", round(tmp*100,2)


print(":::::::::::::::::::RELATIVE ABUNDANCE (in %)::::::::::::::::::")
print >>f, ":::::::::::::::::::RELATIVE ABUNDANCE (in %)::::::::::::::::::"
#finding relative abundance
with open("/home/dna/hariss/centrifuge_new/centrifuge/centrifuge_report.tsv") as fh:
	for line in fh:
                 words=line.split('\t',8)
		 if words[1] in tax_ids:
			tmp=(float(words[5])*float(words[7])+ float(align_ctr[words[1]]*float(avg_rd_len) ))/ float(words[3])
			tmp =float(tmp)/float(sum_abd)
			print("Relative abundance(in %) of ",words[0], "taxid:", words[1], ":",round(100*float(tmp),2))
			print >>f, "Relative abundance(in %) of ",words[0], "taxid:", words[1], ":",round(100* float(tmp),2)
f.close()
print(":::::::::::END of Software analysis pipeline::::::")
print("results saved in /home/dna/hariss/centrifuge_new/centrifuge/abundance_results.txt")
