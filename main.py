from multiprocessing import Pool
import time
import os
import sys
import argparse
import shutil
import ParseGbk
import blast_script
import blast_parse
import computeBicliques
import high_throughput_tests
import json

__author__ = 'Arnon Benshahar'


def parser_code():
    parser = argparse.ArgumentParser(
        description='The purpose of this script is to run the full software suite that we have developed to study gene clusters.')

    parser.add_argument("-q", "--qfolder", dest="qfolder", metavar="DIRECTORY", default='./res/IslandViewerPAI/',
                        help="Folder containing the fasta and Island Viewer format files of the centroid query.")

    parser.add_argument("-g", "--dbfolder", dest="dbfolder", metavar="DIRECTORY", default='./res/genomes/',
                        help="Folder containing all genbank files for use by the program.")

    parser.add_argument("-o", "--outfolder", dest="outfolder", metavar="DIRECTORY", default='./OUT/',
                        help="Folder where the results of a run will be stored.")

    parser.add_argument("-d", "--window", dest="window_size", metavar="INT", default=15,
                        help="Folder where the results of a run will be stored.")

    parser.add_argument("-n", "--num_proc", dest="num_proc", metavar="INT", default=os.sysconf("SC_NPROCESSORS_CONF"),
                        type=int,
                        help="Number of processors that you want this script to run on. The default is every CPU that the system has.")

    parser.add_argument("-iv", "--island_viewer_format", dest="island_viewer_format", metavar="STRING", default='T',
                        help="IslandViewer queries format, T for islandviewer format and F for normal gbk file")

    parser.add_argument("-min_genomes", "--min_genomes_per_block", dest="min_genomes_per_block", metavar="INT", default=5,
                        help="Minimum genome in a gene-block.")

    parser.add_argument("-min_genes", "--min_genes_per_interval", dest="min_genes_per_interval", metavar="INT", default=5,
                        help="Minimum genes in a gene interval.")

    parser.add_argument("-rank", "--min_rank", dest="min_rank", metavar="INT", default=20,
                        help="Minimum ranking score that will be report")

    parser.add_argument("-parse", "--parse_input", dest="parse_input", metavar="STRING", default='T',
                        help="Parse the input files")

    return parser.parse_args()

def check_options(parsed_args):
    if os.path.isdir(parsed_args.dbfolder):
        dbfolder = parsed_args.dbfolder
    else:
        print "The folder %s does not exist." % parsed_args.dbfolder
        sys.exit()

    if os.path.isdir(parsed_args.qfolder):
        qfolder = parsed_args.qfolder
    else:
        print "The folder %s does not exist." % parsed_args.qfolder
        sys.exit()

    # if the directory that the user specifies does not exist, then the program makes it for them.
    if not os.path.isdir(parsed_args.outfolder):
        os.makedirs(parsed_args.outfolder)
    if parsed_args.outfolder[-1] != '/':
        outfolder = parsed_args.outfolder + '/'
    else:
        outfolder = parsed_args.outfolder

    # section of code that deals determining the number of CPU cores that will be used by the program
    if parsed_args.num_proc > os.sysconf("SC_NPROCESSORS_CONF"):
        num_proc = os.sysconf("SC_NPROCESSORS_CONF")
    elif parsed_args.num_proc < 1:
        num_proc = 1
    else:
        num_proc = int(parsed_args.num_proc)

    # validate the input for the window size
    try:
        window_size = int(parsed_args.window_size)
        if window_size <= 0:
            print "The window that you entered %s is a negative number, please enter a positive integer." % parsed_args.max_gap
            sys.exit()
        else:
            pass
    except:
        print "The window that you entered %s is not an integer, please enter a positive integer." % parsed_args.max_gap
        sys.exit()

    # validate the query input format (isalndviewer or gbk)
    if parsed_args.island_viewer_format == 'F' or parsed_args.island_viewer_format == 'T':
        island_viewer_format  = (parsed_args.island_viewer_format =='T')
    else:
        print "T for isalndviewer format and F for normal gbk file"
        sys.exit()

    # validate the input for the min_genomes_per_block
    try:
        min_genomes_per_block = int(parsed_args.min_genomes_per_block)
        if min_genomes_per_block <= 1:
            print "The min genomes per block that you entered %s is less than 2, please enter a positive integer greater than 2." % parsed_args.max_gap
            sys.exit()
        else:
            pass
    except:
        print "The min genomes per block you entered %s is not an integer, please enter a positive integer." % parsed_args.max_gap
        sys.exit()

    # validate the input for the min_genomes_per_block
    try:
        min_genes_per_interval = int(parsed_args.min_genes_per_interval)
        if min_genes_per_interval <= 1:
            print "The min_genes_per_interval you entered %s is less than 2, please enter a positive integer greater than 2." % parsed_args.min_genes_per_interval
            sys.exit()
        else:
            pass
    except:
        print "The min genomes per block you entered %s is not an integer, please enter a positive integer." % parsed_args.min_genes_per_interval
        sys.exit()

    # validate the input for the min_genomes_per_block
    try:
        min_rank = int(parsed_args.min_rank)
        if min_rank <= 0:
            print "The min rank you entered %s is not an integer, please enter a positive integer." % parsed_args.min_rank
            sys.exit()
        else:
            pass
    except:
        print "The min rank you entered %s is not an integer, please enter a positive integer." % parsed_args.min_rank
        sys.exit()

    # validate the query input format (isalndviewer or gbk)
    if parsed_args.parse_input == 'F' or parsed_args.parse_input == 'T':
        parse_input= (parsed_args.parse_input == 'T')
    else:
        print "T for isalndviewer format and F for normal gbk file"
        sys.exit()

    return dbfolder, qfolder, outfolder, num_proc, window_size, island_viewer_format, island_viewer_format, min_genomes_per_block, parse_input



def parse_cmd():
    return "./res/High-Throughput-Plamids/" , "./res/genomes/", "1e-50"
    # return "./res/High-Throughput-Plamids/", "./res/genomes/", "1e-50"


def biclustering(tuple_list, high_throughput):
    print 'Start Biclustering'
    print str(tuple_list)
    query_file, refernce_folder, ref_fasta, query_fasta, blast_results, e_val, blast_parse_dir, query_gene_list_dir, bicluster_results, max_genome_size = tuple_list

    print "Stage 3 parse blast results"
    list_file_name = query_gene_list_dir + query_file.split("/")[-1].split(".")[0] + ".txt"
    blast_parse.parse_blast(blast_results, blast_parse_dir, "", 10, list_file_name, query_file)

    print "Stage 4 biclustering"
    if not high_throughput:
        return computeBicliques.compute_bicluster(list_file_name, blast_parse_dir, bicluster_results, refernce_folder,
                                                  max_genome_size)


def parallel_high_throughput_test(tuple_list_array):
    num_proc = 1
    pool = Pool(processes=num_proc)
    pool.map(high_throughput_tests.compute_bicluster, tuple_list_array)


def main():
    print 'Start main'
    parsed_args = parser_code()

    db_folder, q_folder, outfolder, num_proc, window_size, island_viewer_format, island_viewer_format, min_genomes_per_block, parse_input = check_options(parsed_args)

    #only for high_throughput
    high_throughput = False
    start = time.time()

    # parsed_args = parser_code()

    # this json will contain all the stats about this run
    stats_json = {'query_list': []}

    tmp_dir = './TMP'
    ref_fasta_dir = tmp_dir + '/ffc/ref/'
    # query_fasta_dir = tmp_dir + '/ffc/qry/'/'
    query_fasta_dir = tmp_dir + '/ffc/qryheavy/'

    if parse_input:
        if os.path.exists(outfolder):
            shutil.rmtree(outfolder)
        os.makedirs(outfolder)

        if os.path.exists(tmp_dir):
            shutil.rmtree(tmp_dir)
        os.makedirs(tmp_dir)

        if os.path.exists(ref_fasta_dir):
            shutil.rmtree(ref_fasta_dir)
        os.makedirs(ref_fasta_dir)

        if os.path.exists(query_fasta_dir):
            shutil.rmtree(query_fasta_dir)
        os.makedirs(query_fasta_dir)

        # Stage 1: conver gbk files of the query into fasta format and the reference to ffc format.

        if island_viewer_format:
            query_json = ParseGbk.parseIslandViewer(q_folder, query_fasta_dir)

        else:
            query_json = ParseGbk.parse(q_folder, query_fasta_dir, "NONE", True,False)

        target_json = ParseGbk.parse(db_folder, ref_fasta_dir, "NONE", True, True)

        with open(outfolder + 'queries.json', 'w') as outfile1:
            json.dump(query_json, outfile1)

        with open(outfolder + 'targets.json', 'w') as outfile1:
            json.dump(target_json, outfile1)

    else:
        print 'Run biclustering'
        # create the queries file
        blast_results_dir = tmp_dir + '/blast_results-HEAVY/'
        if os.path.exists(blast_results_dir):
            shutil.rmtree(blast_results_dir)
        os.makedirs(blast_results_dir)

        blast_parse_dir = tmp_dir + '/blast_parse/'
        if os.path.exists(blast_parse_dir):
            shutil.rmtree(blast_parse_dir)
        os.makedirs(blast_parse_dir)

        query_gene_list_dir = tmp_dir + '/query_gene_list_dir/'
        if os.path.exists(query_gene_list_dir):
            shutil.rmtree(query_gene_list_dir)
        os.makedirs(query_gene_list_dir)

        query_fasta_list = []
        print 'Creating blast parse folders'
        for content in os.listdir(query_fasta_dir):  # "." means current directory
            if content.split(".")[-1] == "ffc":
                query_fasta_list.append(query_fasta_dir + content)
                if not os.path.exists(blast_results_dir + "" + content.split(".")[-2]):
                    os.makedirs(blast_results_dir + "" + content.split(".")[-2])
                if not os.path.exists(blast_parse_dir + "" + content.split(".")[-2]):
                    os.makedirs(blast_parse_dir + "" + content.split(".")[-2])

        print 'Create targets.json file'
        with open(outfolder + 'targets.json') as data_file:
            targets_json = json.load(data_file)

        sum = 0
        for target in targets_json:
            sum += target['number_of_genes']
        genome_size = sum / len(targets_json)

        print "Avg genome size " + str(genome_size)

        general_stats = []
        file_num = 1
        tuple_list_array = []
        e_val = 0.01
        for file in query_fasta_list:
            print 'File Number ' + str(file_num)
            file_num += 1

            # Stage 2: run blastall with the query fasta vs the ref fastas
            query_file_name = file.split("/")[-1].split(".")[-2]
            blast_output = blast_results_dir + query_file_name + "/"
            blast_script.blast(ref_fasta_dir, blast_output, "", os.sysconf("SC_NPROCESSORS_CONF"), file, e_val)

            blast_results_tmp = blast_results_dir + file.split("/")[-1].split(".")[-2] + "/"
            blast_parse_tmp = blast_parse_dir + file.split("/")[-1].split(".")[-2] + "/"
            bicluster_results_tmp = outfolder + file.split("/")[-1].split(".")[-2] + "/"

            print str(file) + "\n" + str(db_folder) + "\n" + str(ref_fasta_dir) + "\n" + str(
                query_fasta_dir) + "\n" + str(blast_results_tmp) + "\n" + str(blast_parse) + "\n" + str(
                query_gene_list_dir) + "\n" + outfolder
            tuple_list = (file, db_folder, ref_fasta_dir, query_fasta_dir, blast_results_tmp, e_val, blast_parse_tmp,
                          query_gene_list_dir, bicluster_results_tmp, genome_size)
            # file_stats = biclustering(tuple_list, high_throughput)
            # if not high_throughput:
            #     file_stats['accession'] = query_file_name
            #     with open(outfolder + 'queries.json') as data_file:
            #         query_json = json.load(data_file)
            #     for query in query_json:
            #         if query['accession'] == query_file_name:
            #             file_stats['length'] = query['length']
            #             file_stats['number_of_genes'] = query['number_of_genes']
            #     if file_stats['num_of_cliques'] > 0:
            #         general_stats.append(file_stats)
            #         with open(outfolder + 'resultStats.json', 'w') as outfile1:
            #             json.dump(general_stats, outfile1)
            # else:
            #     tuple_list_array.append((query_gene_list_dir + file.split("/")[-1].split(".")[0] + ".txt", blast_parse_tmp, './', db_folder))

        if island_viewer_format:
            print "high_throughput_test"
            parallel_high_throughput_test(tuple_list_array)


if __name__ == '__main__':
    main()
