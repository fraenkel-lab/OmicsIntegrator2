#!/usr/bin/env python3

import argparse
import sys

def createOutFile(proteinf, garnetf, rnaf, metf):
	#Output is a dict keyed by geneSymbol/ID, like {name:{prize:float, proteinChange:float, geneChange:float, terminalType:string}}
	prizes = {}

	#protein file
	if proteinf:
		for line in proteinf:
			words = line.strip().split('\t')
			if len(words) != 2:
				sys.exit('Protein file does not have two columns')
			name = words[0]
			try:
				logfc = float(words[1])
			except TypeError:
				sys.exit('Second column of protein file should be a number')
			if name in prizes:
				if prizes[name]['prize'] < abs(logfc):
					prizes[name]['prize'] = abs(logfc)
					prizes[name]['proteinChange'] = logfc
					prizes[name]['terminalType'] = 'Proteomic'
			else:
				prizes[name] = {'prize': abs(logfc),'proteinChange':logfc,'geneChange':0.0,'terminalType':'Proteomic'}	

	#metabolite file
	if metf:
		for line in metf:
			words = line.strip().split('\t')
			if len(words) != 2:
				sys.exit('Metabolite file does not have two columns')
			name = words[0]
			try:
				logfc = float(words[1])
			except TypeError:
				sys.exit('Second column of metabolite file should be a number')
			if name in prizes:
				if prizes[name]['prize'] < abs(logfc):
					prizes[name]['prize'] = abs(logfc)
					prizes[name]['proteinChange'] = logfc
					prizes[name]['terminalType'] = 'Metabolite'
			else:
				prizes[name] = {'prize': abs(logfc),'proteinChange':logfc,'geneChange':0.0,'terminalType':'Metabolite'}

	#garnet file
	if garnetf:
		for line in garnetf:
			words = line.strip().split('\t')
			if len(words) != 2:
				sys.exit('Garnet file should have 2 columns')
			name = words[0]
			try:
				prize = float(words[1])
			except TypeError:
				sys.exit('Second column of garnet file should be a number')
			if name in prizes:
				if prizes[name]['prize'] < prize:
					prizes[name]['prize'] = prize
					prizes[name]['terminalType'] = 'TF'
			else:
				prizes[name] = {'prize': prize,'proteinChange':0.0,'geneChange':0.0,'terminalType':'TF'}
	
	#RNA file - This needs to go last, since we only want RNA information on nodes that are already prizes
	if rnaf:
		for line in rnaf:
			words = line.strip().split('\t')
			if len(words) != 2:
				sys.exit('RNA file does not have two columns')
			name = words[0]
			try:
				logfc = float(words[1])
			except TypeError:
				sys.exit('Second column of RNA file should be a number')
			if name in prizes:
				if abs(prizes[name]['geneChange']) < abs(logfc):
					prizes[name]['geneChange'] = logfc
			else:
				continue

	return prizes




def main():
	#Parsing arguments (run python PCSF.py -h to see all these decriptions)
	parser = argparse.ArgumentParser(description='Create one prize file as expected by graph.py from multiple user input data types')
	parser.add_argument("-o" "--outFile", dest="outfile", type=argparse.FileType('w'), required=True,
			help='Name of the output prize file with node attributes as columns, as expected by graph.py')
	parser.add_argument("-p", "--proteinFile", dest='proteinFile', type=argparse.FileType('r'), required=False, 
			help='Proteomic data file with two columns: GeneSymbol\tlogFC of Protein')
	parser.add_argument("-g", "--garnetFile", dest='garnetFile', type=argparse.FileType('r'), required=False, 
			help='Output file from Garnet, showing relevant TFs, with two columns: GeneSymbol\tprize')
	parser.add_argument("-r", "--rnaFile", dest='rnaFile', type=argparse.FileType('r'), required=False, 
			help='mRNA data file with two columns: GeneSymbol\tlogFC of gene. Note that this file will only be used in the visualization, not in the choice of pathways in the protein interaction network. To use mRNA data to assign protein prizes, use Garnet.')
	parser.add_argument("-m", "--metaboliteFile", dest='metFile', type=argparse.FileType('r'), required=False, 
			help='Metabolomics data file with two columns: MetaboliteID\tlogFC of metabolite. Make sure that your interactome includes metabolite interactions, the default does not.')

	args = parser.parse_args()

	prizes = createOutFile(args.proteinFile, args.garnetFile, args.rnaFile, args.metFile)
	
	args.outfile.write('name\tprize\tProteinChange\tGeneChange\tTerminalType\n')
	for p in prizes:
		args.outfile.write('%s\t%s\t%s\t%s\t%s\n'%(p,prizes[p]['prize'],prizes[p]['proteinChange'],prizes[p]['geneChange'],prizes[p]['terminalType']))
	

if __name__ == '__main__':
	main()
