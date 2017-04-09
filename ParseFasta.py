__author__ = 'user'

def filesInDir(root_dir):
    results = []
    return results;

def calculate_lengths(root_dir,results):
    for result in results:
        path = root_dir + "/" + result
        result_file = open("./genomes_length.txt",'w')
        with open(path, 'r') as f:
            count = 0;
            for line in f:
                if line[0] == ">":
                    count = count + 1
            result_file.write(result + ": " + count)
    result_file.close()

def compute_length(root_dir):
    results = filesInDir(root_dir);
    calculate_lengths(root_dir,results);