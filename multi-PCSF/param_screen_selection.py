import json
import os
import sys
import argparse
from graph import *

#Run parameter screen and then select the "best" parameters
#In this case, meaning the average degree of hidden nodes should be less than, or no more than 25% greater than, the average degree of steiner nodes. After that is satisfied, choose the one with the highest number of terminals

def run_param_screen(prize_file, edge_file, w_list, b_list, g_list):
    goodparams = []
    for w in w_list:
        for a in a_list:
            for b in b_list:
                #run pcsf for these parameter combinations
                outdir = '%s_w%s_b%s_g%s'%(prize_file.rsplit('.',1)[0],w,b,g)
                if not os.path.exists(outdir): os.mkdir(outdir)
                params = {"w":int(w), "b":int(b), "g":int(g)}
                graph = Graph(edge_file, params)
                graph.prepare_prizes(prize_file)
                vertex_indices, edge_indices = graph.pcsf()
                forest, augmented_forest = graph.output_forest_as_networkx(vertex_indices, edge_indices)
                output_networkx_graph_as_json_for_cytoscapejs(augmented_forest, outdir)


                #Evaluate parameter set
                terminals=0
                hiddens=0
                totalterm=0
                totalhidden=0
                graphf = open('%s_w%s_b%s_g%s/graph_json.json'%(prize_file.rsplit('.',1)[0],w,b,g),'r')
                nodes = json.load(graphf)['elements']['nodes']
                for node in nodes:
                    if 'prize' in node['data']:
                        terminals = terminals+1
                        totalterm = totalterm + int(node['data']['degree'])
                    else:
                        hiddens = hiddens+1
                        totalhidden = totalhidden + int(node['data']['degree'])
                if terminals>0:
                    avgterm = totalterm/terminals
                    avghidden = totalhidden/hiddens
                    diff = avghidden/float(avgterm)
                    #print('For w = %s, b = %s, g = %s, there are %i nodes, and the average degree of hidden nodes is %.2f, while the average degree of terminals is %.2f'%(w,b,g,totalterm+totalhidden, avghidden, avgterm))
                    if diff>0 and diff<1.25: goodparams.append([w,b,g,totalterm])
                #else:
                    #print('For w = %s, b = %s, g = %s, there are 0 terminals'%(w,b,g))

    #Determine optimal set (if it exists)
    print (prize_file)
    print ('There are %i parameter sets with acceptable difference between average degrees.'%len(goodparams))
    if len(goodparams)>0:
        maxtermsindex = -1
        maxterms = 0
        for i, params in enumerate(goodparams):
            if params[3] > maxterms:
                maxtermindex = i
                maxterms = params[3]
        bestw, bestb, bestg, num = goodparams[maxtermsindex]
        print ('Of these, the parameter set with maximum terminals is w = %s, b = %s, g = %s.\n'%(bestw, bestb, bestg))

        
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""Run OmicsIntegrator2 on several parameter sets and suggest the best one""")
    parser. add_argument("-ws", "--w_list", dest="w_list", nargs="+", default=[5,10], help="A list of integers for the parameter w (number of trees). default='5 10'")
    parser. add_argument("-bs", "--b_list", dest="b_list", nargs="+", default=[5,10], help="A list of integers for the parameter b (size of network). default='5 10'")
    parser. add_argument("-gs", "--g_list", dest="g_list", nargs="+", default=[0,100000,500000], help="A list of integers for the parameter a (negative prize on hubs). default='0 100000 500000'")
    parser.add_argument("-e", "--edge", dest='edge_file', required=True,
	help ='(Required) Path to the text file containing the edges. Should be a tab delimited file with 3 columns: "nodeA\tnodeB\tweight(between 0 and 1)"')
    parser.add_argument("-p", "--prize", dest='prize_file', required=True,
	help='(Required) Path to the text file containing the prizes. Should be a tab delimited file with lines: "nodeName(tab)prize"')

    args = parser.parse_args()

    run_param_screen(args.prize_file, args.edge_file, args.w_list, args.b_list, args.g_list)
