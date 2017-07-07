#!/usr/bin/env python3

# Core python modules
import sys, os

# Peripheral python modules
import argparse
import logging

# python external libraries
import numpy as np
import pandas as pd

# list of classes and methods we'd like to export:
__all__ = [ "merge_prize_files", "merge_prize_dataframes" ]


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
handler.setLevel(logging.INFO)
handler.setFormatter(logging.Formatter('%(asctime)s - Prizes: %(levelname)s - %(message)s', "%I:%M:%S"))
logger.addHandler(handler)


parser = argparse.ArgumentParser(description="""
	Create one prize file as expected by graph.py from multiple user input data types""")

class FullPaths(argparse.Action):
	"""Expand user- and relative-paths"""
	def __call__(self,parser, namespace, values, option_string=None):
		setattr(namespace, self.dest, os.path.abspath(os.path.expanduser(values)))

def directory(dirname):
	if not os.path.isdir(dirname): raise argparse.ArgumentTypeError(dirname + " is not a directory")
	else: return dirname

parser.add_argument("-o" "--out", dest="output_file", type=argparse.FileType('w'), required=True,
	help='Name of the output prize file with node attributes as columns, as expected by graph.py')

parser.add_argument("-p", "--proteins", dest='protein_file', type=argparse.FileType('r'), required=False,
	help='Proteomic data file with two tab-spaced columns: GeneSymbol\tprize')

parser.add_argument("-g", "--garnet", dest='garnet_file', type=argparse.FileType('r'), required=False,
	help='Output file from Garnet, showing relevant TFs')

parser.add_argument("-r", "--rna", dest='rna_file', type=argparse.FileType('r'), required=False,
	help='mRNA data file with two tab-spaced columns: GeneSymbol\tprize. Note that this file will only be used in the visualization, not in the choice of pathways in the protein interaction network. To use mRNA data to assign protein prizes, use Garnet.')

parser.add_argument("-m", "--metabolites", dest='metabolite_file', type=argparse.FileType('r'), required=False,
	help='Metabolomics data file with two columns: MetaboliteID\tlogFC of metabolite. Make sure that your interactome includes metabolite interactions, the default does not.')


def output_dataframe_to_tsv(dataframe, output_file):
	"""
	Output the dataframe to a csv
	"""
	dataframe.to_csv(output_file, sep='\t', header=True, index=False)


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

	prizes_dataframe = pd.concat(prize_dataframes)
	prizes_dataframe.drop_duplicates(subset=['name'], inplace=True) # Unclear if we should do this?

	return prizes_dataframe


def main():

	args = parser.parse_args()

	prizes_dataframe = merge_prize_files(args.protein_file, args.garnet_file, args.metabolite_file, ['protein', 'garnet', 'metabolites'])

	output_dataframe_to_tsv(prizes_dataframe, args.output_file)



if __name__ == '__main__':
	main()
