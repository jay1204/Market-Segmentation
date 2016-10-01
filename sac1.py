import numpy as np
import math
import igraph as ig
import copy
import sys

if len(sys.argv) != 2:
	print "python sac1.py <alpha>"
	sys.exit(1)
alpha = float(sys.argv[1])

def main():
	# read and format the graph
	g = read_data()

	# pairwise modularity and similarity 
	newmanQ = [[0 for x in range(g.vcount())] for x in range(g.vcount())] 
	attrQ = [[0 for x in range(g.vcount())] for x in range(g.vcount())] 
	for i in range(g.vcount()):
		for j in range(i+1,g.vcount()):
			newmanQ[i][j] = g[i, j] - float(g.degree(i) * g.degree(j))/(2*g.ecount())
			attrQ[i][j] = cosineSimilarity(g.vs[i]['attr'], g.vs[j]['attr'])
			newmanQ[j][i] = newmanQ[i][j]
			attrQ[j][i] = attrQ[i][j]

	nodes = map(lambda x: [x], range(g.vcount()))
	communities = map(lambda x: [x], range(g.vcount()))
	nodeDict = dict((i,i) for i in range(len(nodes)))

	# Phase 1 and Phase 2
	count = len(communities)
	for i in range(15):
		phase1(nodes, communities, alpha, nodeDict, newmanQ, attrQ, count,g)
		nodes = list(x for x in communities if x != [])
		communities = copy.deepcopy(nodes)
		nodeDict = dict((i,i) for i in range(len(nodes)))

	outputfile = open("communities.txt",'w')

	for comm in communities:
		output = ','.join(map(str,comm))
		outputfile.write(output+'\n')

	outputfile.close()

def read_data():
	attrData = np.genfromtxt('data/fb_caltech_small_attrlist.csv', delimiter = ",", names = True, dtype = None)
	edgeData = np.genfromtxt('data/fb_caltech_small_edgelist.txt', delimiter = " ", dtype = None)

	g = ig.Graph()
	g.add_vertices(len(attrData))
	edgelists = map(lambda x: (x[0],x[1]), edgeData)
	g.add_edges(edgelists)
	g.vs['attr'] = attrData

	return g

'''
cosine similarity
'''
def dotproduct(v1, v2):
	return sum((a*b) for a, b in zip(v1, v2))

def length(v):
	return math.sqrt(dotproduct(v, v))

def cosineSimilarity(v1, v2):
	return dotproduct(v1, v2) / (length(v1) * length(v2))


"""
Phase 1: mergeCommunity
"""
def deltaQEachVertex(vertex, community, alpha, newmanQ, attrQ, count,g):
	newman = 0.0
	attrs = 0.0
	for index in community:
	    newman = newman + newmanQ[vertex][index]
	    attrs = attrs + attrQ[vertex][index]
	newman = newman / (2*g.ecount())
	attrs = attrs / (len(community) * len(community)*count)
	return newman * alpha + (1 - alpha) * attrs

def deltaQ(node, community, alpha, newmanQ, attrQ, count,g):
	value = 0.0
	#count = len(node)
	for nod in node:
		value = value + deltaQEachVertex(nod, community, alpha, newmanQ, attrQ, count,g)
	return value

def moveNode(node, communities,alpha, nodeDict, nodeIndex, newmanQ, attrQ, count,g):
	maxValue = 0.0
	maxIndex = -1
	change = False
	#for i in range(len(communities)):
	    #if communities[i] != []:
	for i in range(len(communities)):
		if communities[i] != []:
			value = deltaQ(node, communities[i], alpha, newmanQ, attrQ, count,g)
			if(value > maxValue):
				maxValue = value
				maxIndex = i
	if maxIndex != -1 and maxIndex != nodeDict[nodeIndex]:
		nodeDict[nodeIndex] = maxIndex
		communities[maxIndex].extend(node)
		change = True
	else:
		communities[nodeDict[nodeIndex]].extend(node)
	    
	return change

def phase1(nodes,communities, alpha, nodeDict, newmanQ, attrQ, count,g):
	gain = 1.0
	iteration  = 1
	change = True
	while change and iteration <= 15:
		change = False
		gain = 0.0
		iteration = iteration + 1
		for i in range(len(nodes)):
			for vertex in range(len(nodes[i])):
				communities[nodeDict[i]].remove(nodes[i][vertex])
			result = moveNode(nodes[i], communities, alpha, nodeDict, i, newmanQ, attrQ, count,g)
			if result:
				change = True


if __name__ == "__main__":
    main()


