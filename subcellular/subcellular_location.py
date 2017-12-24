import csv
import argparse

def sc_parser(filename):
	sc_dict = {}
	with open(filename, "r") as infile:
		reader = csv.reader(infile, delimiter = "\t")
		for row in reader:
			proteinName = row[1]
			sc_compartment = row[3]
			confidence = row[6]
			try:
				sc_dict[proteinName].append((sc_compartment, confidence))
			except:
				sc_dict[proteinName] = [(sc_compartment, confidence)]

	return sc_dict

def sc_combine_dict(exp_dict, know_dict):
	for key in exp_dict.keys():
		try:
			know_dict[key].extend(exp_dict[key])
		except:
			know_dict[key] = exp_dict[key]
	return know_dict

def get_broad_location(location):
	location = location.lower()
	if "mitochond" in location or "respir" in location:
		return "mitochondria"
	elif "golgi" in location:
		return "golgi"
	elif "endosom" in location:
		return "endosome"
	elif "perox" in location or "microbod" in location:
		return "peroxisome"
	elif "endoplasmic reticulum" in location or "sarcoplas" in location or "signal recog" in location:
		return "ER"
	elif "lytic" in location or "lyso" in location:
		return "lysosome"
	elif "nucleolus" in location or "fibrillar" in location:
		return "nucleolus"
	elif  "snrnp" in location or \
		"repli" in location or "splice" in location or "replic" in location or \
		"kineto" in location or "nucl" in location or "polymerase" in location or \
		"chrom" in location or "dna" in location or "centrosome" in location:
		return "nucleus"
	elif "endomemb" in location or "copi" in location or "granule" in location:
		return "endomembrane"
	elif "intracell" in location:
		return "intracellular"
	elif "riboso" in location:
		return "ribosome"
	elif "aggresome" in location or "proteasome" in location:
		return "proteasome"
	elif "exocyt" in location or "exosome" in location or "extracellular" in location or "cell-cell" in location:
		return "extracellular"
	elif "leading edge" in location or "filopo" in location or "cortex" in location or "desmo" in location or "centrio" in location or \
		"myosin" in location or "lamell" in location or "microtub" in location or "ruffle" in location or "cytoske" in location or \
		"cili" in location or "actin" in location or "focal ad" in location or "projection" in location or "adhere" in location or \
		"junction" in location or "contractile" in location or "kinesin" in location or "filament" in location or \
		"spindle" in location:
		return "cytoskeleton"
	elif "phago" in location or "pino" in location:
		return "phagocytic"
	elif "organelle" in location:
		return "organelle"
	elif "vesic" in location or "vacuol" in location:
		return "vesicle"
	elif "apical" in location or "membr" in location or "periph" in location or "furrow" in location or "surface" in location:
		return "membrane"
	elif "macromol" in location or "complex" in location or "supramol" in location:
		return "protein complex"
	elif "cytoplasm" in location or "cytosol" in location:
		return "cytoplasm"
	elif "cell body" in location or "ranvier" in location or "synap" in location or "cajal" in location or \
		"bouton" in location or "neuron" in location or "cone" in location or "axon" in location or \
		"dendrit" in location or "perik" in location: 
		return "neuron"
	elif "tubule" in location or "sarco" in location or "zone" in location or "fibril" in location or "myelin" in location or "disc" in location or "band" in location:
		return "muscle"
	elif "lipopro" in location or "chylo" in location:
		return "lipoprotein"
	else:
		return location



def get_final_loc(locations):

	if len(locations) == 1:
		return get_broad_location(locations[0])
	else:
		cleaned_locs = [get_broad_location(x) for x in locations]
		counts = {"endomembrane": 0,
					"endosome" : 0,
					"golgi" : 0,
					"ER" : 0,
					"mitochondria" : 0,
					"nucleolus" : 0,
					"nucleus" : 0,
					"proteasome" : 0,
					"cytoskeleton": 0,
					"lysosome" : 0,
					"proteasome" : 0,
					"phagocytic" : 0,
					"extracellular" :0,
					"vesicle" : 0,
					"intracellular" : 0,
					"protein complex" : 0,
					"lipoprotein" :0,
					"organelle" : 0,
					"membrane" : 0,
					"cytoplasm" : 0,
					"ribosome" : 0,
					"neuron" : 0,
					"peroxisome" : 0,
					"muscle" : 0}
		endomembrane = ["golgi", "endosome", "ER", "endomembrane"]
		for location in cleaned_locs:
			try:
				if location in endomembrane:
					counts["endomembrane"] += 1
					if location != "endomembrane":
						counts[location] += 1
				else:
					counts[location] += 1
			except:
				if "go" not in location:
					pass
					#print(location)
		if counts["mitochondria"] >= 1:
			return "mitochondria"
		elif counts["lysosome"] >= 2:
			return "lysosome"
		elif counts["peroxisome"] >= 1:
			return "peroxisome"
		elif counts["golgi"] >= 2 and counts["ER"] >= 2:
			return "endomembrane"
		elif counts["endomembrane"] >= 2:
			if counts["golgi"] >= 1 and counts["golgi"] > counts["ER"]:
				return "golgi"
			elif counts["ER"] >= 1 and counts["ER"] > counts["golgi"]:
				return "ER"
			return "endomembrane"
		elif counts["nucleolus"] >= 1:
			return "nucleolus"
		elif counts["ribosome"] >= 2:
			return "ribosome"
		elif counts["nucleus"] >= 1:
			return "nucleus"
		elif counts["cytoskeleton"] >= 2:
			return "cytoskeleton"
		elif counts["proteasome"] >= 1:
			return "proteasome"
		elif counts["phagocytic"] >= 1:
			return "phagocytic"
		elif counts["protein complex"] >= 1:
			return "protein complex"
		elif counts["lipoprotein"] >= 1:
			return "lipoprotein"
		else:
			max_so_far = ["other", 0]
			for loc in counts.keys():
				if counts[loc] > max_so_far[1]:
					max_so_far = [loc, counts[loc]]

			return max_so_far[0]

def buildLocationDict(sc_exp_fn, sc_know_fn, outfile):
	exp_dict = sc_parser(sc_exp_fn)
	know_dict = sc_parser(sc_know_fn)
	combined_dict = sc_combine_dict(exp_dict, know_dict)

	final_dict = {}
	badwords = ["cellular_component",
				"Cell",
				"Cell part",
				"Intracellular"]

	for proteinName in combined_dict.keys():
		possible_locs = combined_dict[proteinName] 
		sorted(possible_locs, key = lambda loc: loc[1], reverse = True)
		filt_badwords = [x for x in possible_locs if x[0] not in badwords]
		try:
			max_top = int(filt_badwords[0][1])
			filt = [x for x in filt_badwords if int(x[1]) == max_top]
		except:
			filt = [("nonspecific", 0)]
		#print(possible_locs[0:min(len(possible_locs), 4)])
		final_loc = get_final_loc([x[0] for x in filt])
		final_dict[proteinName] = [final_loc, [x[0] for x in filt[0:min(len(possible_locs), 10)]]]

	with open(outfile, "w", newline = '') as ofile:
		writer = csv.writer(ofile, delimiter = "\t")
		for proteinName in final_dict.keys():
			row = [proteinName, final_dict[proteinName][0], final_dict[proteinName][1]]
			writer.writerow(row)




def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("subcellular_experiments_info", help = "Filename for subcellular experiments information from COMPARTMENT.")
	parser.add_argument("subcellular_knowledge_info", help = "Filename for subcellular knowledge information from COMPARTMENT.")
	parser.add_argument("-o", "--outfile", help = "Filename for output file")
	args = parser.parse_args()

	sc_exp_fn = args.subcellular_experiments_info
	sc_know_fn = args.subcellular_knowledge_info
	isOutfile = args.outfile is not None
	if isOutfile:
		outfile = args.outfile
	else:
		outfile = "./subcellular_location_by_protein.tsv"
	buildLocationDict(sc_exp_fn, sc_know_fn, outfile)

if __name__ == "__main__":
	main()

