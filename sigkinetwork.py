
#Modules necessary to construct a kinase network from the user uploaded file
import pandas as pd
from pyvis.network import Network
import networkx as nx
from networkx import *
import plotly.express as px
import plotly.graph_objects as go 
import csv
import re
import itertools
from itertools import chain
import numpy as np
from scipy.stats import norm
import statistics
import operator
import functools
from networkx.readwrite import json_graph
import json, community
from networkx.algorithms import centrality as cn
from community import community_louvain
from statistics import median


networkfile = pd.read_csv('Desktop/PROJECT_DATASET.csv')
networkfile = networkfile.dropna()  

global headers
headers = list(networkfile.columns.values)

#Convert first column (edges) of the document to a list
edges_list= networkfile[networkfile.columns[0]].to_list()
#Convert the list into a list of tuples, each tuple being one edge (since this is the format needed to add edges onto pyvis.network)
edges_list2 = [tuple(x.split('.',1)) for x in edges_list]

#Separate edges into two lists
#List with positive edges and the value of each of the edges
positive_edges = [[(*edges_list2[n],+ networkfile[networkfile.columns[na]].to_list()[n]) for n in range(len(edges_list2)) if networkfile[networkfile.columns[na]].to_list()[n] > 0 ] for na in range(1,len(headers))]
#List with negative edges and the value of each of the edges
#Values are converted to positive values to calculate the centrality values
global negative_edges
negative_edges = [[(*edges_list2[n],+ networkfile[networkfile.columns[na]].to_list()[n] * -1) for n in range(len(edges_list2)) if networkfile[networkfile.columns[na]].to_list()[n] < 0 ] for na in range(1,len(headers))]

#Get list of nodes from positive edges
#Get list of nodes from negative edges

#Get list of allnodes
nodes = list(chain.from_iterable(edges_list2))
#Remove duplicates from the list
global nodes_all
nodes_all = list(dict.fromkeys(nodes))

networkgraph = Network()

##SELECTED

weighted_strength_graphs = 'n' * len (negative_edges)

global sp_edges
sp_edges = 'n' * len (negative_edges)

#Function to run community_louvain to return pair of kinases that are grouped together more than 80 times when
#running a community louvain partition algorithm a hundred times
def asso(network, n_runs):
    partition = community_louvain.best_partition(network,resolution=1, randomize=100)
    kinase_groups = {key:[k for k in partition.keys() if partition[k]==partition[key] and k!=key] for key in partition.keys()}
    for n in range(n_runs-1):
        partition2 = community_louvain.best_partition(network, randomize=1)
        kinase_groups2 = {key:[k for k in partition.keys() if partition2[k]==partition2[key] and k!=key] for key in partition.keys()}
        for k, v in kinase_groups2.items():
            kinase_groups[k] = kinase_groups[k]+ v if k in kinase_groups else v
    
    kg_counts = {combination:kinase_groups[combination[0]].count(combination[1])  for combination in nodes_comb}
    confidence = {key:value for key,value in kg_counts.items() if value>=80}
    #Now get groups
    confidence0 = list(set([key[0] for key in confidence.keys()]))
    confidence_lists = [ [x for x in confidence.keys() if x[0]==kin ] for kin in confidence0 ]
    confidence1 = [[ x[n][1] for n in range(len(x)) ] for x in confidence_lists]
    
    groups_dict = {}
    for n in range(len(confidence0)):
        groups_dict[confidence0[n]] = [[confidence1[n][na]] + [confidence1[n][ne] for ne in range(len(confidence1[n])) if (confidence1[n][na], confidence1[n][ne]) or (confidence1[n][ne], confidence1[n][na])  in confidence.keys()]  for na in range(len(confidence1[n]))]
    
    def IntersecOfSets(kinases_list): 
    # Converting the arrays into sets 
        s1 = set(kinases_list[0])  
        for n in range(len(kinases_list)):
            s2 = set(kinases_list[n]) 
            s1 = s1.intersection(s2)   
        return list(s1)

    groups_dict = {key: IntersecOfSets(groups_dict[key]) for key in groups_dict.keys()}
    
    global groups_dict2
    groups_dict2 = [ [key] + groups_dict[key]  for key in groups_dict.keys() if len([key] + groups_dict[key])>4]
    #groups_dict.sort(key= lambda x:len(x), reverse=True)
    groups_dict2 = { key:sorted([x for x in groups_dict2 if key in x], key= lambda x:len(x), reverse=True) for key in groups_dict.keys()}
    groups_dict2 = {key:groups_dict2[key] for key in groups_dict2.keys() if groups_dict2[key]}
    groups_dict2 = {key:groups_dict2[key][0] for key in groups_dict2.keys() }
    groups_dict2= list(groups_dict2.values())
    
    for x in groups_dict2:
        for y in groups_dict2:
            if set(x)==set(y):
                groups_dict2.remove(y)


    groups_dict2.sort()
    groups_dict2 = list(groups_dict2 for groups_dict2,_ in itertools.groupby(groups_dict2))    
    
    #Average WSC in the group
    global groups_score
    groups_score = [ median([ weighted_strength_graphs[kin] for kin in x]) for x in groups_dict2]
    
    sp = groups_dict2[groups_score.index(max(groups_score))]
    #sp = groups_dict2[-3]
    
    
    neg3 = ['AKT1_2' if x== 'AKT1' else x for x in neg2]
    sp_edges[n] = [edge for edge in nodes_bf if edge[0] in neg3 and edge[1] in neg3]
    sp_edges[n] = [('AKT1', edge[1]) if edge[0]=='AKT1_2' else (edge[0], 'AKT1') if edge[1]=='AKT1_2' else edge for edge in sp_edges ]
        
    return sp_edges

#Analysis of each of the negative edges networks in the dataset
#Repeat with positive edges networks
for n in len(negative_edges):
    selected_edges = negative_edges[n]
	
      #Get list of nodes
      #Keep only the kinases names (i.e. remove the z-scores) from the list of edges
    nodes_bf = [ (a,b) for a,b,c in selected_edges]
      #Flatten list of edges tuples to list of nodes
    nodes = list(itertools.chain(*nodes_bf))
      #Remove duplicates from the list
    nodes_list = list(dict.fromkeys(nodes))

        #Add nodes & edges to graph
    networkgraph.add_nodes(nodes_list, value=[10]*len(nodes_list), title=nodes_list, label=nodes_list)
    networkgraph.add_edges(selected_edges)

	#Show to which kinases each kinase is connected
      #Show to which kinases each kinase is connected
    kinases_map = networkgraph.get_adj_list()
    network_map = networkgraph.get_edges()
    for key,value in kinases_map.items():
        kinases_map[key] = ([ functools.reduce(operator.add, ((x,str([(y["width"]) for y in network_map if x in y["from"] and key in y["to"] or x in y["to"] and key in y["from"]])))) for x in value])



    nodes_numerical = {nodes_all[n]:n for n in range(len(nodes_all))} 
    nodes_names = dict(zip(list(nodes_numerical.values()), nodes_all))
    edges_numerical = [(nodes_numerical[selected_edges[n][0]],nodes_numerical[selected_edges[n][1]], abs(selected_edges[n][2])) for n in range(len(selected_edges))]  
    
    ###Centrality measures
    nxx = nx.Graph()
        #all_edges = [(*selected_edges[noo][0],+ selected_edges[noo][1]) for noo in range(len(selected_edges))]
    #nxx.add_weighted_edges_from(edges_numerical)
    nxx.add_weighted_edges_from(selected_edges)
 

    edges_weight = {x:sorted([ abs(selected_edges[n][2]) for n in range(len(selected_edges)) if x in selected_edges[n]], reverse=True) for x in nodes_list}

    strength_graphs = {kin:sum(edges_weight[kin]) for kin in edges_weight}
    degree_graphs = {kin:len(edges_weight[kin]) for kin in edges_weight}
    betweenness_graphs = betweenness_centrality(nxx,normalized=False, weight='weight')
  
    weighted_degree_graphs = {}
    weighted_strength_graphs = {}
    for kin in list(edges_weight.keys()):
        x_parameter = sum([ sum(edges_weight[kin][n:len(edges_weight[kin])]) for n in range(1,len(edges_weight[kin]))]) / strength_graphs[kin] + 0.5
        r_parameter = x_parameter / (degree_graphs[kin]/2) 
        weighted_degree_graphs[kin] = r_parameter * degree_graphs[kin]
        weighted_strength_graphs[kin] = np.sqrt( strength_graphs[kin] * weighted_degree_graphs[kin] )
    
    #if n < 7 :
    #weighted_strength_graphs[kin] = (weighted_strength_graphs[kin]) * -1

    #Get list of all unique possible node_node combinations
    #FOR EACH EDGE LIST

    nodes_comb = [ [(nodes_list[n],nodes_list[n+na]) for na in range(1,len(nodes_list)-n)] for n in range(len(nodes_list)) ]
    nodes_comb = list(itertools.chain(*nodes_comb))
    
    asso(nxx,100)
                
