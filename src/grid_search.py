import numpy as np
import pandas as pd


from graph import Graph


class Options(object):
	def __init__(self, options):
		self.__dict__.update(options)


graph = Graph()

defaults = {"w": 6, "b": 1, "mu": 0, "a": 20, "noise": 0.1, "mu_squared": False, "exclude_terminals": False, "dummy_mode": "terminals", "seed": None}
params = Options(defaults)
graph.params = params



Ws = np.linspace(start, stop, num=50)
Bs = np.linspace(start, stop, num=50)
As = np.linspace(start, stop, num=50) # this won't work since this is computed during startup.

# generate a grid over w, b, a with https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.optimize.brute.html

# inside loop:
	# prepare_prizes
	# pcsf
