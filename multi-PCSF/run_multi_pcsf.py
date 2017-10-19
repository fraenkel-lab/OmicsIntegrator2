#!/usr/bin/env python3

#Run dendrogram-based multi-sample PCSF
#constraining similar samples to result in similar networks

#Inputs: 
#   -prize file per sample
#   -One interactome
#   -Dendrogram representing heirarchical clustering of samples

import argparse
import sys, os
import pickle
from graph import Graph, output_networkx_graph_as_edgelist

def run_single_PCSF(prizeFile, edgeFile, paramDict, outdir):
    #One straightforward application of PCSF
    graph = Graph(edgeFile, paramDict)
    graph.prepare_prizes(prizeFile)
    vertex_indices, edge_indices = graph.pcsf()
    forest, augmented_forest = graph.output_forest_as_networkx(vertex_indices, edge_indices)
    output_networkx_graph_as_edgelist(augmented_forest, outdir)
    return forest

def run_multi_PCSF(dendrogram, prizefileslist, edgeFile, paramDict, alpha, lbda, outdir):
    #Iterate through dendrogram, and at each clade, run forest for each sample, adding artificial prizes
    
    dendrogram = pickle.load(open(dendrogram,'rb'))
    with open(prizefileslist, 'r') as pf:
        prizefiles = pf.readlines()
    N = len(prizefiles)
    lastF = [[] for _ in range(N)] #list of N lists of nodes in latest version of each sample
    last_iteration_for_samples = [-1 for _ in range(N)] #keep track of last iteration where this sample appeared
    origP = [[] for _ in range(N)] #keep track of original prizes per sample so we know which are Steiners

    #Run the first iteration with unaltered prizes and note which nodes are terminals
    os.makedirs(outdir + '/initial', exist_ok=True)
    names = []
    for i,p in enumerate(prizefiles):
        p = p.strip()
        p_name = os.path.basename(p)
        names.append(p_name)
        i_outdir = outdir + '/initial/%s'%p_name
        os.makedirs(i_outdir, exist_ok=True)
        unadjusted_forest = run_single_PCSF(p, edgeFile, paramDict, i_outdir)
        lastF[i] = unadjusted_forest.nodes()
        with open(p,'r') as pf:
            for line in pf:
                prot = line.split('\t')[0]
                origP[i].append(prot)
      # Change the disimilarity to similarity score: 
    Height=[]               
    for i, c in enumerate(dendrogram):
         Height.append(float(c[2]))
       #print(i)
       
    Similarity=[]             

    for i in range(0,len(Height)):
        Similarity.append((1-(Height[i]-min(Height))/(max(Height)-min(Height)))+ 0.00001)


    #now interate over dendrogram, and at each merge, re-run PCSF for samples in that merge
    for i,c in enumerate(dendrogram):
        os.makedirs(outdir + '/iter%i'%i, exist_ok=True)
        num_samples_in_clade = int(c[3])
        s_c = float(Similarity[i])
        #s_c = float(c[2])
        d_c = 1
        sample1 = int(c[0])
        sample2 = int(c[1])
        samples_in_clade = calc_original_samples(sample1, N, dendrogram) + calc_original_samples(sample2, N, dendrogram)
        if len(samples_in_clade) != num_samples_in_clade: sys.exit('Calculation of samples in clade %i went wrong'%i)
        
        #get frequency of nodes in these networks
        forestFreq = nodeFrequency([lastF[k] for k in samples_in_clade])

        #Update forests for each sample
        for s in samples_in_clade:
            s_outdir = outdir + '/iter%i/%s'%(i,names[s])
            os.makedirs(s_outdir, exist_ok=True)
            #read in artificial prizes we already assigned to this sample in previous iterations
            last_iter = last_iteration_for_samples[s]
            if last_iter == -1:
                artificial_prizes = {}
            else:
                s_lastdir = outdir + 'iter%i/%s'%(last_iter,names[s])
                with open('%s/artificial_prizes.txt'%(s_lastdir), 'r') as a:
                    artificial_prizes = {}
                    for line in a:
                        node, prize = line.strip().split('\t')
                        artificial_prizes[node] = float(prize)
            #Calculate new artificial prizes for this iteration
            for node in forestFreq:
                if node not in origP[s]:
                    if node in artificial_prizes:
                        artificial_prizes[node] = artificial_prizes[node] + lbda*((s_c*forestFreq[node])**(alpha-d_c)) 
                    else:
                        artificial_prizes[node] = lbda*((s_c*forestFreq[node])**(alpha-d_c))
            #write new prizes + original prizes to a file
            with open('%s/updated_prizes.txt'%(s_outdir),"w") as f:
                with open('%s/artificial_prizes.txt'%(s_outdir),'w') as a:
                    for item in artificial_prizes:
                        f.write("%s\t%s\n" % (str(item), str(artificial_prizes[item])))
                        a.write("%s\t%s\n" % (str(item), str(artificial_prizes[item])))
                with open(prizefiles[s].strip(), "r") as p:
                    f.writelines(p.readlines())

            #submit new artificial prizes + orig prize list to run_single_pcsf
            new_forest = run_single_PCSF('%s/updated_prizes.txt'%(s_outdir), edgeFile, paramDict, s_outdir)

            #update lastF
            lastF[s] = new_forest.nodes()
            
def calc_original_samples(sample, N, Z):
    #recursively find all original samples in this merge
    #if the sample number is greater than N, this is a clade of several other samples
    if sample < N:
        return [sample]
    else:
        newIdx = sample-N
        return calc_original_samples(int(Z[newIdx][0]), N, Z) + calc_original_samples(int(Z[newIdx][1]), N, Z)
        
def nodeFrequency(list_of_node_lists):
    #return dict with node:freq of all nodes in these lists
    n = float(len(list_of_node_lists)) #want floating point division
    countDict = {}
    for forest in list_of_node_lists:
        for node in forest:
            if node in countDict:
                countDict[node] += 1
            else:
                countDict[node] = 1
    freqDict = {}
    for node in countDict:
        freqDict[node] = countDict[node] / n
    return freqDict
        

def main():
    parser = argparse.ArgumentParser(description="""
	Run the PCSF algorithm several times on hierarchically clustered samples, refining the tree for each sample so that similar
    samples result in similar trees.""")

    parser.add_argument("-e", "--edge", dest='edge_file', required=True,
	    help ='(Required) Path to the text file containing the edges. Should be a tab delimited file with 3 columns: "nodeA\tnodeB\tcost"')
    parser.add_argument("-p", "--prizefiles", dest='prize_files', required=True,
	    help='(Required, one or more) A text file containing paths to the files containing the prizes. The list should be in the same order as was provided to create the dendrogram. Each path should lead to a tab delimited files with lines: "nodeName(tab)prize"')
    parser.add_argument("-d", "--dendrogram", dest="dendrogram", required=True,
        help='(Required) pickled object denoting hierarchical clustering of samples, of the type returned by scipy.cluster.heirarchy\'s linkage().. Should be an array of length n-1, where dendrogram[i] indicates which clusters are merged at the i-th iteration.')
    parser.add_argument('-o', '--output', dest='output_dir', default='.',
	    help='Output directory path. Default current directory')

    parser.add_argument("-w", dest="w", default=5, type=float, help="Omega: the weight of the edges connecting the dummy node to the nodes selected by dummyMode [default: 5]")
    parser.add_argument("-b", dest="b", default=1, type=float, help="Beta: scaling factor of prizes [default: 1]")
    parser.add_argument("-g", dest="g", default=0, type=float, help="Gamma: Edge penalty on hubs [default: 0]")
    parser.add_argument("-l", "--lambda",dest="lbda", default=1, type=float, help="Lambda: scaling factor on artificial prizes [default: 1]")
    parser.add_argument("-a", "--alpha",dest="alpha", default=1, type=float, help="Alpha: non-linear scaling on artificial prizes [default: 1]")
    
    args = parser.parse_args()
    paramDict = {'w':args.w, 'b':args.b, 'g':args.g}

    run_multi_PCSF(args.dendrogram, args.prize_files, args.edge_file, paramDict, args.alpha, args.lbda, args.output_dir)


if __name__ == '__main__':
	main()
