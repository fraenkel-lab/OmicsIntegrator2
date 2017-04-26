#!/usr/bin/env python3

import numpy as np
import pandas as pd


from graph import Graph


def main(Ws, Bs, As):

	graph = Graph(interactome_file, {})

	defaults = graph.params


	Ws = np.linspace(start, stop, num=50)
	Bs = np.linspace(start, stop, num=50)
	As = np.linspace(start, stop, num=50)

	# generate a grid over w, b, a with https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.optimize.brute.html

	def func(w, b, a):

		graph.params.w = w
		graph.params.b = b
		graph.params.a = a   # this won't work since this is computed during startup.

		prizes, terminals, terminal_attributes = graph.prepare_prizes(prize_file)
		vertex_indices, edge_indices = graph.pcsf(prizes)


	graph._aggregate_pcsf(results, 'frequency')


if __name__ == '__main__':
	main()
