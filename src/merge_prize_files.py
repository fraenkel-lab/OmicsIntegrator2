

### this needs to be merged into __main__.py

def main():

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


if __name__ == '__main__':
	main()


