import os
import math
import logging
from Bio import SeqIO

__author__ = 'Arnon Benshahar'


class Utilities(object):

    @staticmethod
    def ffc_files_in_dir(root_dir):
        results = []
        for file in os.listdir(root_dir):
            if file.endswith(".ffc"):
                results.append(file)
        return results

    # This function will return all the files that are in a directory. os.walk is recursive traversal.
    @staticmethod
    def return_recursive_dir_files(root_dir):
        result = []
        for path, dir_name, flist in os.walk(root_dir):
            for f in flist:
                fname = os.path.join(path, f)
                if os.path.isfile(fname):
                    result.append(fname)
        return result

    @staticmethod
    def calculate_lengths(root_dir, results):
        genomes_length = {}
        for result in results:
            path = root_dir + "/" + result
            with open(path, 'r') as f:
                count = 0
                for line in f:
                    if line[0] == ">":
                        count += 1
                genomes_length[result.split(".")[0]] = count
        return genomes_length

    @staticmethod
    def ncr(n, r):
        f = math.factorial
        return f(n) / f(r) / f(n - r)
