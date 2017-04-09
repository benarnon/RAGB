import sys,os
from Bio import SeqIO

root = "/home/user/PycharmProjects/gene_block_identification/genomes_final/"
for path, subdirs, files in os.walk(root):
    for name in files:
        record_iterator = SeqIO.parse(path+"/"+name, "genbank")
	first_record = next(record_iterator)
	print(first_record.annotations["source"] + "\t" + first_record.name)