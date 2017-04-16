# RAGBI Program, 09/04/2017
**Version 1.0.0** 
By Arnon Benshahar

---

---


## How to execute RAGBI?
In order to run RAGBI you need:

- Python 2.7 or later
- **ncbi-blast-2.6.0+** program (please go to https://www.ncbi.nlm.nih.gov/books/NBK52640/ and follow the **downloading, configuration and execution steps**)

The easiest way to run the project is to execute the script named 'main.py'. The defaults that are provided are sufficient to run the project with the inputs provided.

```
./main.py -h usage: main.py [-q FOLDER] [-g FOLDER] [-o FOLDER] [-d INT] [-n INT] [-iv STRING] [-min_genomes INT] [-min_genes INT] [-rank INT] [-parse STRING] [-e FLOAT]
```

### Optional Arguments:
- **-q** : folder containing the gbk or *IslandViewer* format files (see **Input Format** section) of the centroid query (or queries in case of multiple run).
- **-g** : folder containing all reference genbank files for use by the program.
- **-o** : folder where the results of a run will be stored.
- **-d** : size of the window.
- **-n** : number of processors that you want this script to run on. The default is every CPU that the system has.
- **-iv** : IslandViewer queries format, T for islandviewer format and F for normal gbk files.
- **-min_genomes** : minimum genome in a gene-block.
- **-min_genes** : minimum genes in a gene interval.
- **-rank** : minimum ranking score that will be report.
- **-parse** : parse the input files.
- **--e** : eval for the BLAST search.
---

### Input Formt
Our program requires two input folders, the reference folder and the query folder:
1. **Reference folder:** this folder contains subfolder for each specie. Each subfolder contain the gbk file for that specie.
Example: inside ```/db``` there is ```/db/specie1``` folder which contains ```specie1.gbk```, this file includes all the genes in *specie1*. 
2. **Query folder:** for this folder we have two options:
    1. ***IslandViewer Format***:  this is for centroid queries that were predicted by the [IslandViewer 3 tool](http://www.pathogenomics.sfu.ca/islandviewer/browse/) as [genomic islands](https://en.wikipedia.org/wiki/Genomic_island). In this case the folder contains subfolder for each specie that was analysed by the islandviewer tool. Each folder contains the gbk file of the specie and a csv file which includes the details of all the islands that were discovered by the *IslandViewer* tool. 
**Example**:  inside ```/query/specie1``` there are two files: ```speice1.gbk``` and ```specie1.csv```
        - ```speice1.gbk``` includes all the genes in specie1 
        - ```specie1.csv```contains the information of all the islands that were predicted in *specie1*.
    2.  **Genebank format (gbk)** : the folder contains the gbk files of all the centroid queries.

## Output
Our program output a directory with the following files:
1. ```targets_file.csv```: csv file which contains information about the target/refernece genomes that were given as an input.
2. ```centroid_genome_file.csv```: csv file which contains information about the centroid genome that was given as an input.
3. ```general_results.csv```: csv file which contains a breif summary of the results. For each centroid genomes who got at least one **gene block** it displays the following information: it's top **gene block's** ranking score, number of cliques, number of gene blocks and an avarage number of gene block per clique.
4. Result's directory for each centroid genome, it contains two files:
    1.```***centroid_genome_name***_info_file.csv```: this file holds information about all the genes in the centroid genome.
    2.```***centroid_genome_name***_results_file.csv```: this file show an extensive report for the algorithm results for this centroid genome. The results **gene blocks** are divided into **cliques**. ***Example***:

    
```
Clique Number	1								
 									
Block Number	1	
Ranking Score	175.6785691449	
Number of Genes from the Query genome	4	
Number of Intervals from the Target genomes	30		
 									
Group A (genes)							
Gene Number|Start|End|Gene Id|Gene Name|Attribute|Strand|			
1|001|200|g1|ID_001|att1|-1|	
2|201|500|g2|ID_002|att2|1|
3|505|800|g3|ID_003|att3|-1|
 									
Group B (Genes Intervals)								

Interval 1
Specie              genome1	
Strain	            NC_00001	
Number of Genes	    3

Gene Number|Start|End|Gene Id|Gene Name|Attribute|Strand|Target Gene Id|Target Gene Attribute|Blast E-value|
1|500|700|g1|ID_001|att1|-1|g57|att80|0.001|
2|701|1000|g2|ID_002|att2|1|g58|att73|1e-10|
3|1005|1300|g3|ID_003|att3|1|g59|att22|7e-100|

Interval 2
Specie              genome2
Strain	            NC_00002	
Number of Genes	    4

Gene Number|Start|End|Gene Id|Gene Name|Attribute|Strand|Target Gene Id|Target Gene Attribute|Blast E-value|
3|2200|2495|g3|ID_003|att3|1|g120|att80|0.001|
2|2496|2795|g2|ID_002|att2|1|g121|att73|1e-10|
1|2796|2996|g1|ID_001|att1|1|g122|att22|7e-100|
1|2997|3117|g1|ID_001|att1|1|g123|att22|7e-100|
```
## Contact

- Voice: 00972-524541543
- Email: arnon.benshahar@gmail.com

---

## License & Copyright 