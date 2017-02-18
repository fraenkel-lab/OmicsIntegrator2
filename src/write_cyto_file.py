#Write Simple Interaction Format files to store edges in a format supported by  Cytoscape 3.0
optSif = open('%s/%s_optimalForest.sif'%(outputpath,outputlabel), 'wb')
augSif = open('%s/%s_augmentedForest.sif'%(outputpath,outputlabel), 'wb')
dumSif = open('%s/%s_dummyForest.sif' %(outputpath, outputlabel), 'wb')

#Write attribute files in a format supported by Cytoscape
#The first lines of these files contains the variable names
noa = open('%s/%s_nodeattributes.tsv'%(outputpath, outputlabel), 'wb')
noa.write('Protein\tPrize\tBetweennessCentrality\tFractionOfOptimalForestsContaining\tTerminalType\n')

eda = open('%s/%s_edgeattributes.tsv'%(outputpath, outputlabel), 'wb')
eda.write('Edge\tWeight\tFractionOfOptimalForestsContaining\n')

undirEdgesAdded = {}

edgesSorted = self.augForest.edges(data=True)
edgesSorted.sort(key = itemgetter(0, 1))

#iterate through edges to record edge types and edge attributes
for (node1,node2,data) in edgesSorted:
	#Check if interaction between node1 and node2 is directed
	try:
		w = self.inputObj.dirEdges[node1][node2]
		augSif.write(node1+'\tpd\t'+node2+'\n')
		eda.write(node1+' (pd) '+node2+'\t'+str(data['weight'])+'\t'+
				  str(data['fracOptContaining'])+'\n')
		#Check if this edge is in optForest
		if data['fracOptContaining'] > 0.0:
			optSif.write(node1+'\tpd\t'+node2+'\n')
			dumSif.write(node1+'\tpd\t'+node2+'\n')
	#Else undirected
	except KeyError:
		#Don't want to write undirected interactions twice (A pp B, B pp A).
		try:
			undirEdgesAdded[node2][node1] = 2
		except KeyError:
			augSif.write(node1+'\tpp\t'+node2+'\n')
			eda.write(node1+' (pp) '+node2+'\t'+str(data['weight'])+'\t'+
					  str(data['fracOptContaining'])+'\n')
			if data['fracOptContaining'] > 0.0:
				optSif.write(node1+'\tpp\t'+node2+'\n')
				dumSif.write(node1+'\tpp\t'+node2+'\n')
			if node1 in undirEdgesAdded:
				undirEdgesAdded[node1][node2] = 1
			else:
				undirEdgesAdded[node1] = {node2:1}

nodesSorted = self.augForest.nodes(data=True)
nodesSorted.sort(key = itemgetter(0, 1))
#iterate through nodes to record node attributes
for (node,data) in nodesSorted:
	noa.write(node+'\t'+str(data['prize'])+'\t'+str(data['betweenness'])+'\t'+str(data['fracOptContaining'])+'\t'+data['TerminalType']+'\n')

dumSorted = self.dumForest.edges()
dumSorted.sort(key = itemgetter(0, 1))
#Record dummy edges
for (node1,node2) in dumSorted:
	if node1 == 'DUMMY':
		dumSif.write(node1+'\tpd\t'+node2+'\n')

optSif.close()
augSif.close()
dumSif.close()
noa.close()
eda.close()
print 'Wrote output files for Cytoscape, in directory %s, with names starting with '\
	  '"%s".\n' %(outputpath, outputlabel)
