# RAGBI Program, 09/04/2017

**Version 1.0.0** 

By Arnon Benshahar

---

---


## Running Details
In order to run RAGBI you need:

- Python 2.7 or later
- Blast+ program (please go to https://www.ncbi.nlm.nih.gov/books/NBK52640/ and follow the **downloading, configuration and execution steps**. install the **ncbi-blast-2.6.0+**)

The easiest way to run the project is to execute the script named 'main.py'. The defaults that are provided are sufficient to run the project with the inputs provided.

```
./main.py -h usage: main.py [-q FOLDER] [-g FOLDER] [-o FOLDER] [-d INT] [-n INT] [-iv STRING] [-min_genomes INT] [-min_genes INT] [-rank INT] [-parse STRING] [-e FLOAT]
```

### Optional Arguments:
- **-q** : folder containing the fasta and Island Viewer format files of the centroid query.
- **-g** : folder containing all reference genbank files for use by the program.
- **-o** : folder where the results of a run will be stored.
- **-d** : size of the window.
- **-n** : number of processors that you want this script to run on. The default is every CPU that the system has.
- **-iv** : IslandViewer queries format, T for islandviewer format and F for normal gbk file.
- **-min_genomes** : minimum genome in a gene-block.
- **-min_genes** : minimum genes in a gene interval.
- **-rank** : minimum ranking score that will be report.
- **-parse** : parse the input files.
- **--e** : eval for the BLAST search.
---

### Input formt
Our program requires two input folders, the reference folder and the query folder:

1. **Reference folder:** this folder contains subfolder for each specie. Each subfolder contain the gbk file of that specie.
Example: inside ```/db``` there is ```/db/specie1``` folder which contains ```specie1.gbk```, this file includes all the genes in *specie1*. 
2. **Query folder:** for this folder we have two options:
    1. ***IslandViewer Format***:  this is for islands that were predicted by the [IslandViewer 3 tool](http://www.pathogenomics.sfu.ca/islandviewer/browse/). In this case the folder contains subfolder for each specie that was analysed by the islandviewer tool. Each folder contains the gbk file of the specie and a csv file which includes all the islands that were discovered by the *IslandViewer* tool. 
Example:  inside ```/query/specie1``` there are two files: ```speice1.gbk``` and ```specie1.csv```. ```speice1.gbk``` includes all the genes in specie1 and ```specie1.csv```contains the prediction of all the islands in specie1.
    2. ***gbk file***: the folder contains the gbk files of all the queries.


## Results Format
---
## Contact

- Voice: 00972-524541543
- Email: arnon.benshahar@gmail.com

---

## License & Copyright 