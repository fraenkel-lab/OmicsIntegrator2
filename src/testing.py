#!/usr/bin/env python3

# Core python modules
import sys
import os

# Peripheral python modules
import argparse
import pickle
import logging
import random
from collections import Counter
from copy import copy

# Core python external libraries
import numpy as np
import pandas as pd

# Peripheral python external libraries
import networkx as nx

# Lab modules
from pcst_fast import pcst_fast

import matplotlib.pyplot as plt
from graph import *



logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
handler.setLevel(logging.INFO)
handler.setFormatter(logging.Formatter('%(asctime)s - Forest: %(levelname)s - %(message)s', "%I:%M:%S"))
logger.addHandler(handler)



INTERACTOME_FILE = "/nfs/latdata/iamjli/ALS/data/interactome/iRefIndex_v13_MIScore_interactome.txt"
PRIZE_FILE = "/nfs/latdata/iamjli/ALS/data/iMNs/proteomics/20170323_ALS_CTR_iMNs_protein_log2FC_NO_FILTER.tsv"
OUTPUT_DIR = "../output/"


def output_dataframe_to_tsv(dataframe, output_dir, filename):
	"""
	Output the dataframe to a csv
	"""
	path = os.path.join(os.path.abspath(output_dir), filename)
	dataframe.to_csv(path, sep='\t', header=True)


def main():
	
	g = Graph(INTERACTOME_FILE, {})
	prizes, terminals, terminal_attributes = g.prepare_prizes(PRIZE_FILE)

	vertices, edges = g.pcsf(prizes)
	node_attributes = g.output_forest_as_node_attributes_tsv(vertices, edges, terminal_attributes)

	# forest, augmented_forest = graph.randomizations(prizes, terminals, 5, 0)

	# output_networkx_graph_as_gml_for_cytoscape(forest, os.path.join(OUTPUT_DIR, 'output.gml'))

	output_dataframe_to_tsv(node_attributes, OUTPUT_DIR, "node_attributes.tsv")


main()