#!/usr/bin/python

from multiprocessing import Pool
import time
import os
from homolog import *
from Bio import SeqIO


# Copyright(C) 2014 David Ream
# Released under GPL version 3 licence. http://www.gnu.org/licenses/lgpl.html
# Do not remove this comment

# This exists to  make the main function easier to read. It contains code to run the argument parser, and does nothing else.

# this function will return all of the files that are in a directory. os.walk is recursive traversal.
class Gene(object):
    # The class "constructor" - It's actually an initializer
    def __init__(self, gene_id, start, stop):
        self.id = gene_id
        self.start = start
        self.stop = stop


def createGeneListFile(query_file, geneListFileName):
    handle2 = open(geneListFileName, "w")
    for line in [i.strip() for i in open(query_file).readlines()]:
        # print line
        if line[0] == '>':
            # Normal format accession, organism, locus, gene, product, note, str(start), str(stop), str(strand)
            print line

            if len(line.split('|')) == 11:
                a, b, c, gene_name, gene_id, gene_attribute, gene_start, gene_end, gene_strand, d, old_locus = line.split("|")
            else:
                a, b, c, gene_name, gene_id, gene_attribute, gene_start, gene_end, gene_strand, d = line.split("|")
            handle2.write(
                gene_name + "\t" + gene_id + "\t" + gene_start + "\t" + gene_end + "\t" + gene_attribute + "\t" + gene_strand + "\n")
    handle2.close()


def returnRecursiveDirFiles(root_dir):
    result = []
    for path, dir_name, flist in os.walk(root_dir):
        for f in flist:
            fname = os.path.join(path, f)
            if os.path.isfile(fname):
                result.append(fname)
    return result


# The filter file here i think should be either a vacant value (such as '') or a user defined
# value.  I do not think by default it should be given, like I have made the default behavior.    
def parallel_blast_parse_dict(in_folder, out_folder, filter_file):
    global accession, hlog
    operon_out_folder = out_folder
    if not os.path.isdir(operon_out_folder):
        os.makedirs(operon_out_folder)
    if filter_file != '':
        tmp = returnRecursiveDirFiles(in_folder)
        nc_list = [i.strip() for i in open(filter_file).readlines()]
        fname_list = [i for i in tmp if i.split('/')[-1].split('.')[0] in nc_list]
    else:
        fname_list = returnRecursiveDirFiles(in_folder)

    # print fname_list
    for fname in fname_list:
        homoList = []
        for line in [i.strip() for i in open(fname).readlines()]:
            try:
                if line != "":
                    # print line
                    hlog = Homolog.from_blast(line)
                    homoList.append(hlog)
            except:
                print "ERROR", line
            # hlog.Print()

            # this might have to be changed....
            # accession = hlog.accession()
            try:
                accession = hlog.accession()
            except:
                print ""


            # predicted_gene = hlog.blast_annotation()

            '''
            try: # faster implementation than "if predicted_gene in operon_dict.keys():"
                operon = operon_dict[predicted_gene]
                # Debugging the missing operon casABCDE12... no idea right now.
                if operon == 'casABCDE12':
                    print 'AFDFAFDSF'
                if operon in result.keys():    # check if the operon is in the result dict already, if not make a new entry in the else clause
                    if accession in result[operon].keys(): # Check if the organism has been added to the operon
                        result[operon][accession].append(hlog.ret_str())
                    else: # if the organims is not part of the operon, add it
                        result[operon].update({accession:[hlog.ret_str()]})
                else: # add the operon to the result
                    result.update({operon: {accession: [hlog.ret_str()]}})
            except:
                pass

            '''
        homoList.sort(key=lambda x: (x.percent_ident()), reverse=True)
        homoList.sort(key=lambda x: (int(x.start())), reverse=False)
        if len(homoList) > 0:
            handle = open(out_folder + "/" + accession + '.txt', 'w')
            for homo in homoList:
                # this code is omitted because it used to debugging purposes, and is currently unneeded
                '''
                outfile = intermediate_folder + unfilter_folder + operon + '.txt'
                #print "outfile", outfile
                handle = open(outfile, 'w')

                for accession in result[operon].keys():
                    handle.write('\n'.join(result[operon][accession]) + '\n')
                handle.close()
                '''

                # save results where i actually want them to go:
                # print "outfile", out_folder + operon + '.txt'
                start = str(homo.start())
                stop = str(homo.stop())
                e_val = str(homo.e_val())
                strand = str(homo.strand())

                # print start + " "   + stop +  " " + homo.blast_annotation() + " " + eval + " " + str(homo.aligned_length()) + " " + str(homo.percent_ident()) + "\n"
                # for ido's data
                # handle.write(start + " "   + stop +  " " + homo.blast_annotation() + " " + eval + " " + str(homo.aligned_length()) + " " + str(homo.percent_ident()) + "\n")
                # for island viewer
                if homo.query_gene_name() != "unknown" and False:
                    handle.write(
                        start + "|" + stop + "|" + strand + "|" + homo.query_gene_name() + "|" + homo.query_gene_id() + "|" + homo.query_description() + "|" + e_val + "|" + str(
                            homo.aligned_length()) + "|" + str(homo.percent_ident()) +  "|" + str(homo.gene_name()) + "|" + str(homo.description()) +  "\n")
                else:
                    handle.write(
                        start + "|" + stop + "|" + strand + "|" + homo.query_gene_name() + "|" + homo.query_gene_id() + "|" + homo.query_description() + "|" + e_val + "|" + str(
                            homo.aligned_length()) + "|" + str(homo.percent_ident()) +  "|" + str(homo.gene_name()) + "|" + str(homo.description()) + "\n")
            handle.close()


def parse_query(query_file):
    handle = open(query_file, "rU")
    genes = []
    for record in SeqIO.parse(handle, "fasta"):
        split = record.id.split("|")
        gene = Gene(split[3], split[4], split[5])
        genes.append(gene)
    handle.close()
    genes.sort(key=lambda x: (int(x.start)), reverse=False)


def parse_blast(infolder, outfolder, filter_file, num_proc, query_gene_list, query_fasta_file):
    start = time.time()
    folder = outfolder
    for the_file in os.listdir(folder):
        file_path = os.path.join(folder, the_file)
        try:
            if os.path.isfile(file_path):
                os.unlink(file_path)
                # elif os.path.isdir(file_path): shutil.rmtree(file_path)
        except Exception as e:
            print(e)
    print "Blast parse " + infolder, outfolder, filter_file, num_proc, query_gene_list, query_fasta_file
    createGeneListFile(query_fasta_file, query_gene_list)

    print infolder, outfolder, filter_file, num_proc
    parallel_blast_parse_dict(infolder, outfolder, filter_file)
    # parse_query(query_file,"./")


    print time.time() - start
