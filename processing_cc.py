import pickle, numpy as np
import csv

     # Import SPRING helper functions
from preprocessing_python import *
#http://localhost:8000/springViewer.html?datasets/rbcd

     # Import expression matrix; rows are cells and columns are genes
     ### ****** Make sure E.npy is unzipped *************
print 'Loading expression matrix'
#E = np.load('datasets/rbc/K3_1.csv')
#E = np.loadtxt('datasets/rbc/K3_2.csv')
#E = np.load('example_inputs/python_E.npy')
E = np.genfromtxt('datasets/flag_fli1/FLI1_Flag.csv')

#print E

     # Filter out cells with fewer than 1000 UMIs
print 'Filtering cells'
E,cell_filter = filter_cells(E,0)

     # Normalize gene expression data
     # Only use genes that make up <
print 'Row-normalizing'
E = row_normalize(E)

     # Filter genes with mean expression < 0.1 and fano factor < 3
print 'Filtering genes'
_,gene_filter = filter_genes(E,0.0,0)

     # Z-score the gene-filtered expression matrix and do PCA with 20 pcs
print 'Zscoring and PCA'
Epca = get_PCA(Zscore(E[:,gene_filter]), 25)

     # get euclidean distances in the PC space
print 'Getting distance matrix'
D = get_distance_matrix(Epca)

     # load additional data (gene_list, cell_groupings, custom_colors)
     # gene_list is a list of genes with length E.shape[1]
     # cell_groupings is a dict of the form: { <grouping_name> : [<cell1_label>, <cell2_label>,...] }
     # a "grouping" could be the sample id, cluster label, or any other categorical variable
     # custom_colors is a dict of the form { <color_track_name> : [<cell1_value>, <cell2_value>,...] }
     # a "custom color" is any continuous variable that you would like to use for coloring cels. 

#cell_group = np.loadtxt('datasets/rbc/K3_days.csv')
#gene_list, cell_groupings, custom_colors = {NULL, NULL, NULL}
#gene_list, cell_groupings, custom_colors = pickle.load(open('datasets/rbc/K3_days.csv'))
#print pickle.load(open('example_inputs/python_data.p'))
#gene_list, cell_groupings, custom_colors = pickle.load(open('example_inputs/python_data.p'))

gene_list = ["ATRX", "CD123", "CD235ab", "CD33", "CD34", "CD41", "CD44", "CD45RA", "CD49f", "CD71", "CD90", "CEBPa", "c_Jun", "c_Myc", "CXCR4", "GATA1", "GATA2", "HBA", "IKZF1", "KAT3B", "KLF1", "MAFG", "PU1", "RUNX1", "TAL1"]  

print len(gene_list)

reader = csv.reader(open('datasets/flag_fli1/FLI1_Flag_day.csv', 'rb'))
cell_groupings = {}

for row in reader:
	key = row[0]
	if key in cell_groupings:
		pass
	cell_groupings[key] = row[1:]

#print cell_groupings

     # save a SPRING plots with k=5 edges per node in the directory "datasets/frog_python/"
print 'Saving SPRING plot'
save_spring_dir(E, D, 5, gene_list, 'datasets/flag_fli1', cell_groupings=cell_groupings)

print E
#save_spring_dir(E,D,5,gene_list,'datasets/frog_python', cell_groupings=cell_groupings, custom_colors=custom_colors)
