#!/usr/bin/python

from multiprocessing import Pool
import time
import os
import os, shutil
import sys
import argparse
from Bio.Blast.Applications import NcbiblastxCommandline

# Copyright(C) 2014 David Ream
# Released under GPL version 3 licence. http://www.gnu.org/licenses/lgpl.html
# Do not remove this comment

#this function will return all of the files that are in a directory. os.walk is recursive traversal.
def returnRecursiveDirFiles(root_dir):
    result = []
    for path, dir_name, flist in os.walk(root_dir):
        for f in flist:
            fname = os.path.join(path, f)
            if os.path.isfile(fname):
                result.append(fname)
    return result


# This code right now only deals with protein, but I will add functionality later for nucleotides. 
# Just moving the project along here, but this is a critical flaw moving forward.
def do_parallel_blast(arg_tuple):
    db, query_file, blast_result_folder, num_processors, eval_threshold = arg_tuple
    if (db.split('/')[-1] != query_file.split('/')[-1]):
        out_file = "%s%s.txt" % (blast_result_folder, db.split('/')[-1].split('.')[0])
        #print "db", db
        #print "got here"
        #cmd = "blastall -p tblastn -a %i -i %s -d %s -e %s -o %s -m 9" % (num_processors, query_file, db, eval_threshold, out_file)
        #cmd = "blastall -p tbla   stn -a %i -i %s -d %s -e %s -o %s -m 9" % (1, query_file, db, eval_threshold, out_file)
        #cmd = "blastall -p blastp -a %i -i %s -d %s -e %s -o %s -m 9" % (1, query_file, db, eval_threshold, out_file)
        print "before blast query file: " + query_file + ", db:" + db
        cmd = "blastall -p blastp -a %i -i %s -d %s -e %s -o %s -m 8" % (1, query_file, db, 0.01, out_file)
        #print cmd
        os.system( cmd )
        print "finish  blast on query file: " + query_file + ", db:" + db


#def parallel_blast(infile, query, folder, num_proc, e_val = '1e-10'):
def parallel_blast(database_folder, outfolder, filter_file, num_proc, query_file, e_val):
    #Delete all the current files in out folder
    import os, shutil

    print 'parallerl lbast'
    print query_file
    folder = outfolder
    print folder
    for the_file in os.listdir(folder):
        print "THE FILE"
        print the_file
        file_path = os.path.join(folder, the_file)
        try:
            if os.path.isfile(file_path):
                os.unlink(file_path)
            #elif os.path.isdir(file_path): shutil.rmtree(file_path)
        except Exception as e:
            print(e)

    # you kinda have to trust me here, but having blast run on as many threads per CPU as you have total processors is fastest
    # I have no idea why this is... ugh.
    
    unfiltered_db_list = [i for i in returnRecursiveDirFiles(database_folder) if i.split('/')[-1].split('.')[-1] == 'ffc']

    if filter_file == '':
        db_list = unfiltered_db_list
    else:
        filter_list = [i.strip() for i in open(filter_file).readlines()]
        db_list = [i for i in unfiltered_db_list if i.split('/')[-1].split('.')[0] in filter_list]
    
    #print len(unfiltered_db_list), len(db_list)
    
    #blast_arg_list = [(i, query_file, outfolder, num_proc, e_val) for i in db_list]
    blast_arg_list = [(i, query_file, outfolder, 1, e_val) for i in db_list]
    pool = Pool(processes = num_proc)
    pool.map(do_parallel_blast, blast_arg_list)


def blast(database_folder, outfolder, filter_file, num_proc, query_file, e_val):
    
    start = time.time()

    print database_folder, outfolder, filter_file, num_proc, query_file, e_val
    
    parallel_blast(database_folder, outfolder, filter_file, num_proc, query_file, e_val)

    print time.time() - start
    

