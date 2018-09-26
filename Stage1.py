__author__ = 'GH'

##using UPGMA insted of NJ

import TreeConstruction_changed as TreeConstruction
import Stage2
import mafft_and_phylip
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SubsMat import MatrixInfo as matlist
from Bio import pairwise2
import copy
import os
import CasSites
import pickle




#from Bio.Phylo.TreeConstruction import DistanceCalculator


def print_res_to_csvV3(res, input_sg_genes_dict, genesList, genesNames, path, homology = False):
	def do_it(f,candidate_lst):
		for candidate in candidate_lst:
			f.write(str(candidate))
			f.write("\n")
	CSV_file_name = path+ "/output_simple.txt"
	f = open(CSV_file_name, "w")
	if homology:
		for subgroup in res:
			do_it(f,subgroup.candidate_lst)
	else:
		do_it(f, res)

	f.close()


def make_prot_lst(genes_lst):
    prot_lst = []
    for i in range(len(genes_lst)):
        prot_lst.append(DNA_to_prot(genes_lst[i]))  ##have to have the str  only for the names list
        #prot_lst.append(str(DNA_to_prot(genes_lst[i])))  ##have to have the str  only for the names list
    return prot_lst

def DNA_to_prot(dna_seq):
    '''the outpath is the file off al the protein seq to send to nadav
     the outpath should be an open file'''
    coding_dna = Seq(dna_seq, generic_dna)
    return coding_dna.translate()


def make_UPGMA(dm):
    '''use by the doc in http://biopython.org/DIST/docs/api/Bio.Phylo.TreeConstruction.DistanceTreeConstructor-class.html'''
    constructor = TreeConstruction.DistanceTreeConstructor()
    tree = constructor.upgma(dm)
    return tree

def df_genes(seq1, seq2):
	'''pairwise alingment'''
	#matrix = matlist.blosum80
	matrix = matlist.pam120
	#uncomment this part:
	#if(('*' in seq1) or ('*' in seq2)):
	#	print("There is stop codon in the input!! please use only coding DNA as input.")
	#	exit(1)
	seq1 = str(seq1).replace('*','A')
	seq2 = str(seq2).replace('*','A')
	alignments = pairwise2.align.globaldx(seq1, seq2, matrix)  ##scoring funcion by the blosum 62 matrix
	score = len(seq1) - alignments[0][2]  #temporal. To change to something more accurate.
	return score

def make_initiale_matrix(df, seqList):
    '''input: df: distance_function. seqList: list of sequences, in the order coresponding to the names of sequences
        output: the distance according to the given distance function, arranged as list of lists: a lower triangle matrix
    '''
    res = [] # an array of arrays: a matrix
    for i in range(len(seqList)):
        j = 0
        row = []
        while j <= i:
            tempDistance = df(seqList[i],seqList[j])
            row += [tempDistance]
            j += 1
        res += [row]
    #print(res)
    return res

def make_initiale_matrix_try(df, seqList):
    '''input: df: distance_function. seqList: list of sequences, in the order coresponding to the names of sequences
        output: the distance according to the given distance function, arranged as list of lists: a lower triangle matrix
    '''
    res = [] # an array of arrays: a matrix
    for i in range(len(seqList)):
        j = 0
        row = []
        while j <= i:
            tempDistance = df(seqList[i],seqList[j])
            row += [tempDistance]
            j += 1
        res += [row]
    #print(res)
    return res



def make_initiale_matrix_from_protdist(protdist_file):
	'''triangle'''
	import re
	os.rename(protdist_file, protdist_file + ".txt")
	res = []
	in_file = open(protdist_file + ".txt", 'r')
	#next(in_file)
	p = re.compile("[a-zA-Z][a-zA-Z]")
	#p = re.compile("a_gene_name_prefix_")
	temp_res = re.split(p,in_file.read())[1:] #without the number of genes
	res = []
	i = 1
	for line in temp_res:
		#print(line)
		l = re.split(" +",line)[1:] #without the name, I guess..
		#l = list(map(lambda x: float(x.rstrip()),l[:i]))
		#l = list(map(lambda x: max(res[-1]) if x == -1.0 else x, res[-1]))
		#res.append(l)
		to_append = list(map(lambda x: float(x.rstrip()),l[:i]))
		if to_append:
			res.append(to_append)
			
		#res.append(list(map(lambda x: float(x.rstrip()),l[:i]))) #the 1 in the begining of the line - for neglectig the gene name
		#res[-1] = list(map(lambda x: max(res[-1]) if x == -1.0 else x, res[-1]))
		i += 1
	#print(res)
	#print(len(res))
	return res

def make_initiale_matrix_from_protdist_original(protdist_file):
	'''triangle'''
	import re
	os.rename(protdist_file, protdist_file + ".txt")
	res = []
	in_file = open(protdist_file + ".txt", 'r')
	#next(in_file)
	p = re.compile("[a-zA-Z][a-zA-Z]")

	temp_res = re.split(p,in_file.read())[1:] #without the number of genes
	res = []
	i = 1
	for line in temp_res:
		#print(line)
		l = re.split(" +",line)[1:]
		#l = list(map(lambda x: float(x.rstrip()),l[:i]))
		#l = list(map(lambda x: max(res[-1]) if x == -1.0 else x, res[-1]))
		#res.append(l)
		res.append(list(map(lambda x: float(x.rstrip()),l[:i]))) #the 1 in the begining of the line - for neglectig the gene name
		#res[-1] = list(map(lambda x: max(res[-1]) if x == -1.0 else x, res[-1]))
		i += 1
	return res

def handle_minus1(res):
	#handle the -1:
	#find max_val:
	max_val = -1
	for i in range(len(res)):
		for j in range(len(res[i])):
			if res[i][j] > max_val:
				max_val = res[i][j]
	for i in range(len(res)):
		for j in range(len(res[i])):
			if res[i][j] == -1.0:
				res[i][j] = max_val * 1.1




def make_initiale_matrix_from_protdist_old0(protdist_file):
	'''triangle'''

	#first, make a matrix out of the file. then, make it a triangle matrix.

	temp_res, seq_names = [], []
	in_file = open(protdist_file, 'r')
	next(in_file)
	i=0  #row number
	New_line, to_append = True, False
	for line in in_file:  #the rows
		i+=1
		line_as_array = line.split()
		if not line_as_array[0].strip().split(".")[0].isdigit():
			if (to_append):
				temp_res.append(list(map(float, matrix_row)))
			seq_names.append(line_as_array[0].strip())
			matrix_row = []
			matrix_row += line_as_array[1:]

			to_append = True
		else:
			matrix_row += line_as_array
	if temp_res[-1] != list(map(float, matrix_row)):
		temp_res.append(list(map(float, matrix_row))) #last row
	#now, change to lower
	res = []
	for i in range(len(temp_res)):
		res.append(temp_res[i][:i+1])
	return res


def make_initiale_matrix_from_protdist_old01(protdist_file):
	'''triangle'''
	res, seq_names = [], []
	in_file = open(protdist_file, 'r')
	next(in_file)
	f = in_file.read().split(" ")
	matrix_row = []
	for word in f:  #the rows
		word = word.replace(" ","")
		if not word.replace(".","").isdigit(): #a new line
			if len(matrix_row)>0:
				res.append(matrix_row)
			matrix_row = []
		else:
			matrix_row.append(float(word))
	res.append(matrix_row)#appending the last row
	#os.remove(protdist_file)
	return res



def make_initiale_matrix_from_protdist_old1(protdist_file):
	'''triangle'''
	res, seq_names = [], []
	in_file = open(protdist_file, 'r')
	next(in_file)
	i=0  #row number
	for line in in_file:  #the rows
		i+=1
		line_as_array = line.split(" ")
		seq_names.append(line_as_array[0].strip())
		matrix_row = []
		num_of_white_spaces = 0
		for j in range(1, len(line_as_array)):
			if line_as_array[j] != '': ##this line was added at 12/12 - later then the rest of the code
				matrix_row.append(float(line_as_array[j].strip()))
			else:
				num_of_white_spaces += 1
			if j >= i + num_of_white_spaces:
				break
		res.append(matrix_row)
	#os.remove(protdist_file)
	return res


def make_distance_matrix(names, initiale_matrix):
    '''input: list of names of the sequences, and the output of 'make_initiale_matrix'
    output: a distance matrix, in a format adequate to the UPGMA function'''
    m = TreeConstruction._DistanceMatrix(names, initiale_matrix)
    return m

def make_initiale_matrix_from_protdist_old2(protdist_file):
	'''not triangle'''
	res, seq_names = [], []
	in_file = open(protdist_file, 'r')
	next(in_file)
	for line in in_file:
		line_as_array = line.split(" ")
		seq_names.append(line_as_array[0].strip())
		matrix_row = []
		for j in range(1, len(line_as_array)):
			if line_as_array[j] != '': ##this line was added at 12/12 - later then the rest of the code
				matrix_row.append(float(line_as_array[j]))
		res.append(matrix_row)
	os.remove(protdist_file)
	return res


def return_UPGMA_old(seq_list, names_list):
	'''input:  a list of names and a list of sequences, calibrated
	output: an nj tree instance.
	without protdist
	'''
	seq_list = make_prot_lst(seq_list)  #from DNA to protain
	#names_list = make_prot_lst(names_list)
	matrix = make_initiale_matrix(df,seq_list)
	distance_matrix = make_distance_matrix(names_list, matrix)
	upgma1 = make_UPGMA(distance_matrix)
	return upgma1

def return_UPGMA_working(seq_list, names_list, protdist_outfile, pylip_temps_path):
	'''
	with protdist
	:param seq_list: list of sequence for the UPGMA
	:param names_list: the name of the sequences at seq_list, at the same order
	:return: UPGMA tree
	'''
	#protdist_outfile = "protdist_outfile.ph"
	mafft_and_phylip.runIt(names_list,seq_list,protdist_outfile, pylip_temps_path) #to uncomment when really ranning the tool on the server
	matrix = make_initiale_matrix_from_protdist(protdist_outfile)
	#print(matrix)
	distance_matrix = make_distance_matrix(names_list, matrix)
	return make_UPGMA(distance_matrix), distance_matrix
	

def return_UPGMA(seq_list, names_list, protdist_outfile, pylip_temps_path):
	'''
	with protdist
	:param seq_list: list of sequence for the UPGMA
	:param names_list: the name of the sequences at seq_list, at the same order
	:return: UPGMA tree
	'''
	#protdist_outfile = "protdist_outfile.ph"
	new_names_list = list()
	for i in range(len(names_list)):
		new_names_list.append("GG"+str(i)) #for the regex
	
	#print(new_names_list)
	mafft_and_phylip.runIt(new_names_list, seq_list,protdist_outfile, pylip_temps_path) #to uncomment when really ranning the tool on the server
	matrix = make_initiale_matrix_from_protdist(protdist_outfile)
	#print(matrix)
	distance_matrix = make_distance_matrix(names_list, matrix)
	return make_UPGMA(distance_matrix), distance_matrix




def test_DNA_to_prot():
    path = "test_DNA_to_prot.txt"
    print(DNA_to_prot("AGGAGGAGAGGAGGAA"))
    DNA_to_prot("AGGAGGAGAGGAGGAA")


def test_UPGMA():
    lst = ["AGGAGGAGAGGAGGAA", "AGGAGGAGAGGATTAA","AGGAGGAGCCCAGGAA"]
    return_UPGMA(lst)

def test_df():
    seq1 = "GGCA"
    seq2 = "GGTC"
    print(df(seq1, seq2))


###until here: the nj. from here: the rapper algorithm###

#best_permutations_DS = []  #list of tuples: (max_seq, max_fraction, max_cut_prob, genes_list)
sg_genes_dict = {}#key: sgRNA. value: list of genes. mostly, this list will be of size 1.
genes_sg_dict = {}

def print_res_to_file(res, input_sg_genes_dict, path='', homology = False):
	''' input(res) format:	array of permutations_DS'''
	##first, get genes set, for the file name:
	genes_set = set()
	for sg, genes in input_sg_genes_dict.items():  #genes  is a list of genes
		for gene in genes:
			genes_set.add(gene.split("RNA")[0])
	if path != '':
		file_name = path+ "/output.txt"
	else:
		file_name = str(genes_set)[1:-1] + ".txt"
	f = open(file_name, "w")
	for permDS in res:
		f.write("SG:\n")
		f.write(permDS[0])
		f.write("\n")
		f.write("score:\n")
		f.write(str(permDS[2]))
		f.write("\n")
		f.write("Cleaved genes:\n")
		for gene_tup in permDS[3]:
			f.write(gene_tup[0] + "  ")
		f.write("\n")
		f.write("cleavage site for each gene:\n")
		next_line = str(permDS[4])
		next_line = next_line.replace("RNAfile","")
		f.write(next_line)
		f.write("\n\n")
	f.close()


def print_res_to_json(res, input_sg_genes_dict, genesList, genesNames, path):
	res = pickle.load(open("res.p", "rb"))
	genesList = pickle.load(open("genesList.p", "rb"))
	genesNames = pickle.load(open("genesNames.p", "rb"))
	json_file_name = path + "/resultTree.json"
	f = open(json_file_name, "w")
	f.write("[\t{\n")
	title = " ".join(genesNames).replace("'","")
	title = '"input genes: ' + title + 'Predicted sgRNAs:'  #change to 2 saperated titles
	f.write('"text": ' + title + '"\n')
	f.write("}\n")
	for candidate in res:
		f.write(',{\n"text": "sgRNA:' + candidate.seq + '"')
		f.write(',\n"nodes":[\n')#genes this sgRNA is expected to cleave
		number_of_genes_of_candidate = len(candidate.genes_score_dict)
		i = 0
		for gene_name in candidate.genes_score_dict.keys():
			i += 1
			geneIndex = genesNames.index(gene_name)
			f.write('{\n"text": "' + gene_name + '   Cleavage probability: ' + str(candidate.genes_score_dict[gene_name]) + '"')
			f.write(',\n"nodes":[\n')
			number_of_sites = len(candidate.match_sites_dict[gene_name]) #the number of sites of this candidate in this gene
			j = 0
			for site in candidate.match_sites_dict[gene_name]:
				j += 1
				#find the position in the gene
				site_position = str(genesList[geneIndex].find(site[0])) #item[0] is the target site seq
				if site_position == "-1":
					site_position = str(CasSites.give_complementary(genesList[geneIndex]).find(site[0])) + "R"
				to_write = '{\n"text": "' + site[0] + ' Position in gene: ' + site_position +'"\n}'
				if j != number_of_sites:
					to_write += ',\n'
				else:
					to_write += '\n'
				f.write(to_write)
			f.write(']\n}\n') #] ] is the end of the nodes. is the end of the gene.
			if i != number_of_genes_of_candidate:
				f.write(',\n')
			else:
				f.write('\n')
		f.write(']\n}\n') # ] is the end of the nodes. } is the end of the sgRNA
	f.write("]")
	f.close()

def print_res_to_csvV2(res, input_sg_genes_dict, genesList, genesNames, path, homology = False):
	def do_it(res, input_sg_genes_dict, genesList, genesNames, f):
		def str_of_mm(mm_lst):
			copy_mm_list = list(mm_lst)
			copy_mm_list.sort()
			res = ''
			for i in range(len(copy_mm_list)):
				res += str(copy_mm_list[i]+1) + " " #first pos will be considered as 1
			return res
		f.write("sgRNA, Sum of genes score, Targeted genes (score; mm; postions)\n")
		for candidate in res:
			#make the output in the wanted way
			genes_for_Thefile = ""
			genes = sorted(list(candidate.genes_score_dict.keys()), key = lambda gene: candidate.genes_score_dict[gene], reverse = True)
			sites = candidate.targets_dict
			for gene in genes:
				geneIndex = genesNames.index(gene)
				score = str(candidate.genes_score_dict[gene])
				if len(score)>5: #can ve >= as well
					score = score[:5]
				genes_for_Thefile += gene + " (" + str(candidate.genes_score_dict[gene])
				if len(sites) > 0:
					site_index = 1
					for site in sites.keys():
						if site_index >1:
							genes_for_Thefile += "; target_" + str(site_index) + ": "
						if site == gene:
							genes_for_Thefile += "; "
							mm = ''
							for item in sites[site]:  #item here is a single match site
								if len(item[1].keys()) > 0:
									genes_for_Thefile += str_of_mm(item[1].keys())#" " + str(list(map(lambda x: x+1, list(item[1].keys()))).sort())[1:-1].replace(",", "")
								#find pos:
								sg_position = str(genesList[geneIndex].find(item[0])) #item[0] is the target site seq
								if sg_position == "-1":
									sg_position = str(CasSites.give_complementary(genesList[geneIndex]).find(item[0])) + "R"

								genes_for_Thefile += "; pos: " + str(sg_position)
					genes_for_Thefile += ") "
				else:
					sg_position = str(genesList[geneIndex].find(candidate.seq))
					if sg_position == "-1":
						sg_position = str(CasSites.give_complementary(genesList[geneIndex]).find(candidate.seq)) + "R"
					genes_for_Thefile += "1; perfect match; " + sg_position + ") "
			row = candidate.seq + ", " + str(candidate.cut_expectation)+ ","  + genes_for_Thefile + "\n"
			f.write(row)
			#end of do_it()
	CSV_file_name = path+ "/output.csv"
	f = open(CSV_file_name, "w")
	if not homology:
		do_it(res, input_sg_genes_dict, genesList, genesNames, f)
	else:
		for subgroup in res:
			f.write("genes of subgroup: " + str(subgroup.genes_lst)[1:-1].replace("'","") + "\n")#genes_of_subgroup(subgroup) + "\n")
			do_it(subgroup.candidate_lst, input_sg_genes_dict, genesList, genesNames, f)
	f.close()

def genes_of_subgroup(subgroup):
	res = set()
	for c in subgroup:
		for g in c.genes_score_dict.keys():
			res.add(g)
	return str(list(res))[1:-1]


def print_res_to_file_CSV(res, input_sg_genes_dict, genesList, genesNames, path='', homology = False):
	''' input(res) format:	array of permutations_DS
	including finding the position on the gene'''
	##first, get genes set, for the file name:
	genes_set = set()
	for sg, genes in input_sg_genes_dict.items():  #genes is a list of genes
		for gene in genes:
			genes_set.add(gene.split("RNA")[0])
	if path != '':
		file_name = path+ "/output.csv"
	else:
		file_name = str(genes_set)[1:-1] + ".csv"
	f = open(file_name, "w")
	f.write("sgRNA,score, genes (score; mm; postions)\n")
	for sg in res:
		#make the output in the wanted way
		sites = sg[4]
		genes_for_Thefile = ""
		genes = sg[3]
		sites = sg[4]
		for gene in genes:
			geneIndex = genesNames.index(gene[0])
			genes_for_Thefile += gene[0] + " (" + str(gene[1])
			if len(sites) > 0:
				for site in sites:
					if site[0] == gene[0]:
						genes_for_Thefile += "; "
						mm = ''
						for item in site[1]:
							genes_for_Thefile += " " + str(list(map(lambda x: x+1, list(item[1].keys()))))[1:-1]
							##
							#find pos:
							sg_position = str(genesList[geneIndex].find(item[0])) #item[0] is the target site seq
							if sg_position == -1:
								sg_position = str(CasSites.give_complementary(genesList[geneIndex]).find(item[0])) + "R) "
							genes_for_Thefile += "; pos: " + str(sg_position)
			else:
				sg_position = str(genesList[geneIndex].find(sg[0]))
				if sg_position == -1:
					sg_position = str(CasSites.give_complementary(genesList[geneIndex]).find(sg[0])) + "R) "
				genes_for_Thefile += "1; perfect match; " + sg_position
		row = sg[0] + ", " + str(sg[2])+ ","  + genes_for_Thefile + "\n"
		f.write(row)
	f.close()

def stopping_condition(current_best_perm):
	'''rapper version: return True if need to stop'''
	if not current_best_perm:
		return False
	return current_best_perm[1] < 0.75  #if at least 0.75 of the genes in the subtree are cut, continue

def bottem_up(node, current_sg_genes_dict, current_genes_sg_dict, sgList, sgNames, Omega):
	'''caling the buttoms up algorithem with a sg genes dict sutable for the subtree'''
	if node.colour == 'b':
		return
	##making the genes_sg dict for this subtree and the sg_genes_dict to send to the intermadiate algorithm
	if not(current_sg_genes_dict):
		current_sg_genes_dict = {}
	if not(current_genes_sg_dict):
		current_genes_sg_dict = {}
	if not(sgList):
		sgList = []
	if not(sgNames):
		sgNames = []
	for leaf in node.leaves_DS: ##leaf here is a gene. taking only the relevant genes
		current_genes_sg_dict[leaf.name] = genes_sg_dict[leaf.name]
		##filling the sg genes dict
		for sg in current_genes_sg_dict[leaf.name]:
			current_sg_genes_dict[sg] = sg_genes_dict[sg] ##the checking if this sg is already in the dict just be more expensive overall
			if sg not in sgList:
				sgList.append(sg)
				sgNames.append(sg)
	##second, find the key sequence##
	current_res = None
	if len(current_genes_sg_dict) < 2 :  #only one gene
		current_best_perm, lowest_of_widest = Stage2.find_best_sg_for_single_gene(leaf.name, sgList )#lowest_of_widest is not in use in this function
	else:
		#get the set cover from the bottem up algorithm
		current_res = Stage2.call_it_all(sgList, sgNames, current_sg_genes_dict, Omega)##call_it_all(sgList, sgNames, input_sg_genes_dict, Omega)## Naive.find_Uno_sgRNA(current_genes_sg_dict, Omega) #current best perm is a tuple with the perm and metedata of this perm
		current_best_perm = current_res[0]  #the best sg at the current set cover
		if current_res == None:
			return
	global best_permutations_DS
	if (current_res):
		#remove unnided candidates from the current_res
		remove_unrelevant_candidates(current_res)
		#continue
		best_permutations_DS += current_res
	else:
		best_permutations_DS += [current_best_perm]
	node.set_colour('b')
	##continue up##
	if (node.parent) and not(stopping_condition(current_best_perm)):
		bottem_up(node.parent, current_sg_genes_dict, current_genes_sg_dict, sgList, sgNames, Omega)  ##this line is adopted to the rapper algorithm

def remove_unrelevant_candidates(res):
	'''
	:param res:
	:return: without only 1 genes to target (remove the adding of those, if this path will be chosen), and without candidates that there is another canidate which is equal of better then them on each gene.
			#first, implementation in quradratic time, and see if it work
	'''
	if len(res) < 1:
		return
	candidates_to_del = []
	for i in range(len(res) -1, -1,-1):
		if i < len(res) - 1:
			for j in range(i -1, -1, -1):
				#print(i, j, len(res))
				if j >= len(res):
					continue
				#check if res[j] is better or equale than res[i] on all genes
				if is_better_candidate(res[j], res[i]):
					del res[i]

def is_better_candidate(candidateA, candidateB):
	'''
	:param candidateA:
	:param candidateB:
	:return: return true if candidateA has a highr or equale score then candidateB on each of the genes in candidateB.genes_score_dict
	'''
	for gene, cut_exp in candidateB.genes_score_dict.items():
		if not gene in candidateA.genes_score_dict:
			return False
		if candidateA.genes_score_dict[gene] < cut_exp: #not to delete res[i]
			return  False
	return  True

def bottem_up_tree(upgmaTree, Omega):
	for i in range(len(upgmaTree.leaves_DS)):
		bottem_up(upgmaTree.leaves_DS[i], None,None,None,None, Omega)

def call_it_all_wighted(genesList, genesNames, input_sg_genes_dict, input_genes_sg_dict, Omega, protodist_outfile, pylip_temps_path):
	upgmaTree, distance_matrix = return_UPGMA(genesList, genesNames, protodist_outfile, pylip_temps_path) #to uncomment when using wighted
	Stage2.fill_leaves_sets(upgmaTree)  # as apposed to the intermediate algorithem, here leaves are genes
	fill_sg_genes_dict(input_sg_genes_dict)
	fill_genes_sg_dict(input_genes_sg_dict)
	#making the sgList for Algorithm B:
	sgList = list(input_sg_genes_dict.keys())
	sgNames = copy.deepcopy(sgList)
	best_permutations_DS = Stage2.call_it_all(sgList, sgNames, input_sg_genes_dict, Omega, df_targets)##call_it_all(sgList, sgNames, input_sg_genes_dict, Omega)## Naive.find_Uno_sgRNA(current
	best_permutations_DS.sort(key = lambda item: len(item[3]), reverse=True)# or (len(item[2]) and item[1]))  . sort for the print
	res =  find_w_set_cover(best_permutations_DS, distance_matrix)  ##if the output of the intermadiante is wanted
	return res

def call_it_all(genesList, genesNames, input_sg_genes_dict, input_genes_sg_dict, Omega, protodist_outfile, pylip_temps_path, df_targets, cfd_dict = None, PS_number = 12):
	fill_sg_genes_dict(input_sg_genes_dict)
	fill_genes_sg_dict(input_genes_sg_dict)
	sgList = list(input_sg_genes_dict.keys())
	sgNames = copy.deepcopy(sgList)
	best_permutations_DS = Stage2.call_it_all(sgList, sgNames, input_sg_genes_dict, Omega, df_targets, cfd_dict, PS_number)##call_it_all(sgList, sgNames, input_sg_genes_dict, Omega)## Naive.find_Uno_sgRNA(current
	best_permutations_DS.sort(key = lambda item: item.cut_expectation, reverse=True)# or (len(item[2]) and item[1]))  . sort for the print
	return best_permutations_DS


def find_set_cover(candidates_lst, genes_lst):
	'''	 the standard greedy aproximation algorithm
	:param candidates_lst ; a toy example:  ([['ACGCACCC', 0.6, 3.7961971025899606e-06, [('gene1', 0.01608054522924407), ('gene2', 0.01424411400247827), ('gene4', 0.01657343550446999)], [['gene1', [('ACGTACGT', {3: {'C', 'T'}, 6: {'C', 'G'}, 7: {'C', 'T'}})]], ['gene2', [('ACGTAGCT', {3: {'C', 'T'}, 5: {'C', 'G'}, 7: {'C', 'T'}})]], ['gene3', [('ACGTATTG', {3: {'C', 'T'}, 5: {'C', 'T'}, 6: {'C', 'T'}, 7: {'C', 'G'}})]], ['gene4', [('ATGCACGT', {1: {'C', 'T'}, 6: {'C', 'G'}, 7: {'C', 'T'}})]], ['gene5', [('ATGCATGC', {1: {'C', 'T'}, 5: {'C', 'T'}, 6: {'C', 'G'}})]]]], ['ACGCACCT', 0.6, 4.393668496780193e-05, [('gene1', 0.036452247191011256), ('gene2', 0.03157967032967035), ('gene4', 0.03816764705882347)], [['gene1', [('ACGTACGT', {3: {'C', 'T'}, 6: {'C', 'G'}})]], ['gene2', [('ACGTAGCT', {3: {'C', 'T'}, 5: {'C', 'G'}})]], ['gene3', [('ACGTATTG', {3: {'C', 'T'}, 5: {'C', 'T'}, 6: {'C', 'T'}, 7: {'T', 'G'}})]], ['gene4', [('ATGCACGT', {1: {'C', 'T'}, 6: {'C', 'G'}})]], ['gene5', [('ATGCATGC', {1: {'C', 'T'}, 5: {'C', 'T'}, 6: {'C', 'G'}, 7: {'C', 'T'}})]]]] ...
	:param distance_matrix:
	:return:
	'''
	res = []
	names_lst = genes_lst #names of genes by their  order in the distance matrix
	covered_genes = [0 for i in range(len(names_lst))] # 0 if not covered, 1 if covered
	used_candidates = [0 for i in range(len(candidates_lst))]  # the same here
	prices_lst = [0 for i in range(len(candidates_lst))]
	while(0 in covered_genes):
		for c in range(len(candidates_lst)): #going over all of the candidates
			num_of_coveres_genes = 0
			if used_candidates[c] == 1:#this candidate was allready chosen
				continue
			for i in range(len(candidates_lst[c][3])): #the covered genes lst by this candidate
				if covered_genes[names_lst.index(candidates_lst[c][3][i][0])] == 0:# #for knowing which to include in take_into_agv list
					num_of_coveres_genes += 1
			if num_of_coveres_genes == 0:  #this candidate is irelevant for now
				prices_lst[c] = 100000 #just a big number
			else:
				if num_of_coveres_genes > 1:
					print("num_of_coveres_genes > 1")
				prices_lst[c] = 1/num_of_coveres_genes #unwighted set cover
		#stage 2: find the minimal priced set
		min_priced = [-1, max(prices_lst)] #index and price
		for c in range(len(prices_lst)):
			if used_candidates[c] == 0 and prices_lst[c] < min_priced[1]:
				min_priced = [c, prices_lst[c]]
		#add the result to the cover and update the covered genes:
		print("min priced",min_priced)
		res.append(candidates_lst[min_priced[0]])
		used_candidates[min_priced[0]] = 1
			#update the covered genes:
		for i in range(len(candidates_lst[min_priced[0]][3])): #the covered genes lst
			gene_index = names_lst.index(candidates_lst[min_priced[0]][3][i][0])
			print("gene index", gene_index)
			print(covered_genes[gene_index])
			covered_genes[gene_index] = 1
	return res

def find_w_set_cover(candidates_lst ,distance_matrix):
	annealing_coef = 0
	return find_w_set_cover_heated(candidates_lst ,distance_matrix, annealing_coef)

def find_w_set_cover_sevral(candidates_lst ,distance_matrix):
	res = [] #list of lists f lists. each level 2 list represent the size of the set cover, and contains the best candidate with the lowest considering in homology and with highest consideting in homology, correspondingly
	annealing_coef = 0
	min_set_cover = find_w_set_cover_heated(candidates_lst ,distance_matrix, annealing_coef)
	res.append([0,min_set_cover])
	current_set_cover = min_set_cover
	for i in range(1,10):
		prev_set_cover = current_set_cover
		current_set_cover = find_w_set_cover_heated(candidates_lst ,distance_matrix, i/10)
		if len(current_set_cover)>len(prev_set_cover):
			res[len(res)-1].append(prev_set_cover)
			res.append([i/10,current_set_cover])
	return res[0]

def find_w_set_cover_heated(candidates_lst ,distance_matrix, annealing_coef):
	'''	 the standard greedy aproximation algorithm
	:param candidates_lst ; a toy example:  ([['ACGCACCC', 0.6, 3.7961971025899606e-06, [('gene1', 0.01608054522924407), ('gene2', 0.01424411400247827), ('gene4', 0.01657343550446999)], [['gene1', [('ACGTACGT', {3: {'C', 'T'}, 6: {'C', 'G'}, 7: {'C', 'T'}})]], ['gene2', [('ACGTAGCT', {3: {'C', 'T'}, 5: {'C', 'G'}, 7: {'C', 'T'}})]], ['gene3', [('ACGTATTG', {3: {'C', 'T'}, 5: {'C', 'T'}, 6: {'C', 'T'}, 7: {'C', 'G'}})]], ['gene4', [('ATGCACGT', {1: {'C', 'T'}, 6: {'C', 'G'}, 7: {'C', 'T'}})]], ['gene5', [('ATGCATGC', {1: {'C', 'T'}, 5: {'C', 'T'}, 6: {'C', 'G'}})]]]], ['ACGCACCT', 0.6, 4.393668496780193e-05, [('gene1', 0.036452247191011256), ('gene2', 0.03157967032967035), ('gene4', 0.03816764705882347)], [['gene1', [('ACGTACGT', {3: {'C', 'T'}, 6: {'C', 'G'}})]], ['gene2', [('ACGTAGCT', {3: {'C', 'T'}, 5: {'C', 'G'}})]], ['gene3', [('ACGTATTG', {3: {'C', 'T'}, 5: {'C', 'T'}, 6: {'C', 'T'}, 7: {'T', 'G'}})]], ['gene4', [('ATGCACGT', {1: {'C', 'T'}, 6: {'C', 'G'}})]], ['gene5', [('ATGCATGC', {1: {'C', 'T'}, 5: {'C', 'T'}, 6: {'C', 'G'}, 7: {'C', 'T'}})]]]] ...
	:param distance_matrix:
	:return:
	'''
	res = []
	names_lst = distance_matrix.names #names of genes by their  order in the distance matrix
	covered_genes = [0 for i in range(len(names_lst))] # 0 if not covered, 1 if covered
	used_candidates = [0 for i in range(len(candidates_lst))]  # the same here
	prices_lst = [0 for i in range(len(candidates_lst))]
	while(0 in covered_genes):
		print(covered_genes)
		for c in range(len(candidates_lst)): #going over all of the candidates
			single_gene_by_the_set = -1
			weight = 0 # a list of tuples to consider to the avg that makes the weight
			num_of_pairs = 0 #to be later davide the weight
			num_of_coveres_genes = 0
			if used_candidates[c] == 1:#this candidate was allready chosen
				continue
			for i in range(len(candidates_lst[c][3])):
				if covered_genes[names_lst.index(candidates_lst[c][3][i][0])] == 1:# #for knowing which to include in take_into_agv list
					continue
				if num_of_coveres_genes == 0:
					single_gene_by_the_set = i
				num_of_coveres_genes += 1
				for j in range(i+1, len(candidates_lst[c][3])):
					if covered_genes[names_lst.index(candidates_lst[c][3][j][0])] == 1:# #for knowing which to include in take_into_agv list
						continue
					num_of_pairs += 1 # also = num_of_coveres_genes(num_of_coveres_genes-1)
					single_gene_by_the_set = -1
					if annealing_coef == 0:
						weight += (candidates_lst[c][3][i][1] + candidates_lst[c][3][j][1])/2
					else: #dosen't have to have the 'else' here
						weight += ((1-distance_matrix[i,j])**annealing_coef) *(candidates_lst[c][3][i][1] + candidates_lst[c][3][j][1])/2 #is it correct? # (1/cleaving_probobility_i+ 1/cleaving_proboblity_j)/2 # the distance times the avereged covering probobility is the weight
						print(1-distance_matrix[i,j])
			if num_of_coveres_genes == 0:  #this candidate is irelevant for now
				prices_lst[c] = 100000 #just a big number
			else:
				if num_of_coveres_genes == 1: #assuming taking an sgRNA with perfect mach to the target site
					weight = 2 - candidates_lst[c][3][single_gene_by_the_set][1]  #just for testing - have to find what to do with those perfect mached sg.
				else:
					weight = 2 - (weight/num_of_pairs)
				prices_lst[c] = 1/num_of_coveres_genes #unwighted set cover

		#stage 2: find the minimal priced set
		min_priced = [-1, max(prices_lst)] #index and price
		for c in range(len(prices_lst)):
			#print(prices_lst[c])
			if used_candidates[c] == 0 and prices_lst[c] < min_priced[1]:
				min_priced = [c, prices_lst[c]]
		#add the result to the cover and update the covered genes:
		print("min priced",min_priced)
		res.append(candidates_lst[min_priced[0]])
		used_candidates[min_priced[0]] = 1
			#update the covered genes:
		for i in range(len(candidates_lst[min_priced[0]][3])): #the covered genes lst
			gene_index = names_lst.index(candidates_lst[min_priced[0]][3][i][0])
			covered_genes[gene_index] = 1
	return res

def find_set_cover_old(best_permutations_DS):
	res = []
	temp_best_perm_DS = best_permutations_DS
	uncovered_genes = set()
	for sg, genesLst in sg_genes_dict.items():
		for gene in genesLst:
			uncovered_genes.add(gene)
	while(len(uncovered_genes)) > 0 and len(temp_best_perm_DS) > 0:
	##going over all the permutations, and return the permutation that cover the maximal amount of genes haven't been covered yet, in the highest probability among the maximal covered permutations
		best_current_perm, best_num_of_coverd, best_prob_of_covered = None, 0,0  #best_current_perm is the whole tuple
		i = 0
		while i < (len(temp_best_perm_DS)):
			num_of_coverd = 0
			for gene_name in temp_best_perm_DS[i].genes_score_dict.keys():
				if gene_name in uncovered_genes:
					num_of_coverd += 1
			if num_of_coverd > best_num_of_coverd:## and temp_best_perm_DS[i][2] > best_prob_of_covered:
				best_num_of_coverd, best_prob_of_covered = num_of_coverd, temp_best_perm_DS[i].cut_prob
				best_current_perm = temp_best_perm_DS[i]

			elif (num_of_coverd >= best_num_of_coverd and temp_best_perm_DS[i].fraction_of_cut >= best_prob_of_covered) :## and temp_best_perm_DS[i][2] > best_prob_of_covered:  ##need to check if 2 is the right index, and not 1.
				best_num_of_coverd, best_prob_of_covered = num_of_coverd, temp_best_perm_DS[i].fraction_of_cut
				best_current_perm = temp_best_perm_DS[i]
			i+=1
		if(best_current_perm):
			res.append(best_current_perm)
			for gene_name in best_current_perm.genes_score_dict.keys():
				if gene_name in uncovered_genes: #there is a probability that this gene had already been covered bya prevuis sgRNA
					uncovered_genes.remove(gene_name)
	return res

def fill_sg_genes_dict(input_sg_genes_dict):
	global sg_genes_dict
	sg_genes_dict = input_sg_genes_dict

def fill_genes_sg_dict(input_genes_sg_dict):
	global genes_sg_dict
	genes_sg_dict = input_genes_sg_dict

def test_tree_leavesDS(tree):
	for l in tree.leaves_DS:
		print(l)
