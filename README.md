A tool for the optimal design of sgRNAs for simultaneous cleavage of multiple genes.

Given a set of genomic sequences as input, CRISPys first clusters all the potential CRISPR-Cas9 targets that are located within them, and then designs the most promising sgRNAs to each cluster.

This is the standalone version of CRISPys. An online tool can be found at http://multicrispr.tau.ac.il/.
Running instructions:
Dependencies: Python 3.5; Biopython module; For the ‘considering genes homology’ service install MAFFT (https://mafft.cbrc.jp/alignment/software/) and keep the Protdist (http://evolution.genetics.washington.edu/phylip/doc/protdist.html) execution file (supplied here for convenience) in the source directory.

Running command:

python stage0.py path_of_input_file execution_directory_path

or:

python stage0.py –i path_of_input_file -d file execution_directory_path [-a A –w 1 –t 0 –o 0.43 –s cfd –p outfile –l 20 –m 20 –s 0 –n 20 –c 12]

alg: A for the basic CRISPys algorithm or E for considering genes homology when desining the sgRNAs. The considering homology option can run only on Unix. Default: A.

where_in_gene : floating point number between 0 to 1. Search for targets only in this primary fraction of the gene, in order to increase the likelihood that cleavage in the chosen target will be upstream to the active site and results in a loss of function mutation. Default: 1.

t: 0 for using sgRNA to gain maximal gaining score among all of the input genes or 1 for the maximal cleavage likelihood only among genes with score higher than the average. Default: 0.

v: the value of the threshold. Default: 0.43

s: the scoring function of the targets. Optional scoring systems are: cfd (default), …

p: protDist output file name. Default: "outfile".

l: minimal length of the target site. Default:20

m: maximal length of the target site, Default:20

g: 1 if the target sites are obligated to start with a G codon or 0 otherwise. Default: 0.

i: when choosing the consider homology option, this is the number of sgRNAs designed for each homology sub-group. Default: 10

ps the maximal number of possible polymorphic sites in a target. Default: 12

For any questions or commercial use, please contect the author at galhyams@gmail.com 
