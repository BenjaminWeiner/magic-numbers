#memory efficient version that does things per line
#adapted to work in systems with multiple lengths, by reading polyFile

#this version looks at number of monomers in a cluster
import sys
import numpy as np
import pandas as pd
import networkx 
from networkx.algorithms.components.connected import connected_components
import pickle

burnTime=0

testName=sys.argv[1]


#given an array of coords for many proteins at a single time, construct a dictionary 
#with the occupancy at each point, then return a list of connected proteins
def getOccupancies(timeStep):

	#split timeStep into individual proteins
	splitArray=np.zeros(numProteins+numRubisco-1)
	splitArray[0]=lengths[0]
	for j in range(1,numProteins+numRubisco-1):
		splitArray[j]=splitArray[j-1]+lengths[j]

	indivProteins=np.split(timeStep,splitArray.astype(int))    #list of arrays
	occupancyDict={}
	for i in range(numProteins+numRubisco):
		for j in range(lengths[i]):
			thisKey=np.array_str(indivProteins[i][j])   #key is coordinate
			if thisKey not in occupancyDict:
				occupancyDict[thisKey]=[i]
			else:
				occupancyDict[thisKey].append(i)
				
	return list(occupancyDict.values())


#using the networkX package, we convert list of connected nodes into a graph that can be broken into clusters
def to_graph(l):
	G = networkx.Graph()
	for part in l:
		# each sublist is a bunch of nodes
		G.add_nodes_from(part)
		# it also imlies a number of edges:
		G.add_edges_from(to_edges(part))
	return G

def to_edges(l):
	""" 
		treat `l` as a Graph and returns it's edges 
		to_edges(['a','b','c','d']) -> [(a,b), (b,c),(c,d)]
	"""
	it = iter(l)
	last = next(it)

	for current in it:
		yield last, current
		last = current  


#load data*****************************
path=''
filename='polyOutput_'+testName+'.dat'



lengths=[]
polyFile=open(path+'params_'+filename)
polyLines=polyFile.readlines()


proteinLength=int(polyLines[11])
numProteins=int(polyLines[13])
numRubisco=int(polyLines[14])

for i in range(numProteins):
	lengths.append(proteinLength)
for i in range(numRubisco):
	lengths.append(8)


polyOutputFileName='polyOutput_'+filename


#calculate************************************************************
sizeRecord=[]

linesPerStep=np.sum(np.array(lengths))
lineCounter=1
timeStepIndex=0
thisStep=[]
for line in open(path+polyOutputFileName):

	splitLine=line.split()
	intLine=[int(i) for i in splitLine]
	thisStep.append(intLine)
	if lineCounter==linesPerStep:   

		if timeStepIndex>=burnTime:
			stepArray=np.asarray(thisStep)
			thisOccupancy=getOccupancies(stepArray)
			thisGraph = to_graph(thisOccupancy)                   #convert site occupancies to graph
			clusterGen=connected_components(thisGraph)            #generator of sets of connected clusters

			#now go from clusters of protein IDs to monomer number
			sizeList=[]
			for i in clusterGen: #each i is a set of proteins
				thisClusterMonomers=0
				for k in i:               #each k is a protein in the cluster
					thisClusterMonomers+=lengths[k]
				sizeList.append(thisClusterMonomers)

			sizeRecord.append(sizeList)

		thisStep=[]
		lineCounter=1
		timeStepIndex+=1

	else:
		lineCounter+=1

#save your list
outFile='clusterSize_'+testName+'.dat'
pickle.dump(sizeRecord,open(outFile, "wb"))