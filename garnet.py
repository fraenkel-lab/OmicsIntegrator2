#!/usr/bin/env python3

from garnet import *

options = {}
motifs_and_genes = map_known_genes_and_motifs_to_peaks(known_genes_file, motifs_file, peaks_file, options)

motif_prizes = motif_regression(motifs_and_genes, expression_file)


