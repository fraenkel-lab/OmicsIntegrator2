

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



def merge_prize_files(prize_files, prize_types):
	"""
	Arguments:
		prize_files (list of str or FILE): a filepath or FILE object with a tsv of name(\t)prize(\t)more...
		prize_types (list of str): a node type name to associate with the nodes from each prize_file

	Returns:
		pandas.DataFrame: a DataFrame of prizes with duplicates removed (first entry kept)
	"""

	dataframes = []

	for prize_file, prize_type in zip(prize_files, prize_types):

		prize_df = pd.read_csv(prize_file, sep='\t')
		prize_df.columns = ['name', 'prize'] + prize_df.columns[2:].tolist()
		prize_df['type'] = prize_type
		dataframes.append(prize_df)

	return merge_prize_dataframes(dataframes)


def merge_prize_dataframes(prize_dataframes):
	"""
	Arguments:
		prize_dataframes (list of pandas.DataFrame): a list of dataframes, each of which must at least have columns 'name' and 'prize'

	Returns:
		pandas.DataFrame: a DataFrame of prizes with duplicates removed (first entry kept)
	"""

	prizes_dataframe = pd.concat(prizes_dataframes)
	prizes_dataframe.drop_duplicates(subset=['name'], inplace=True) # Unclear if we should do this?

	return prizes_dataframe
