__author__ = 'user'



def set_partition(query_file,out_path,d):
    print "Divide " + query_file + " into " + d + "-mers"
    island_name = query_file.split("/")[-1].split(".")[-2]
    lines = [];
    line_index = 0
    file_index = 0
    for line in [i.strip() for i in open(query_file).readlines()]:
        if line_index < d:
            lines.append(line)
        else:
             handle = open(out_path + "/" + island_name + file_index + '.ffc', 'w')
             for line2 in lines:
                 handle.write(line2)
             lines.pop(0)
             lines.append(line)