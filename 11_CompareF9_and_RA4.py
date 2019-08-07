import matplotlib
matplotlib.use('Agg')
import numpy as np
import pandas as pd
import pickle as pk
import os

import matplotlib.pyplot as plt
from matplotlib_venn import venn2

# ********************** Compare total sites found in both F9 and RA4 ****************** # 
def CombineSites(samples, chrSeq):
	""" 
	Function to combine sites from all samples

	Output: 
		- Dictionary with individual chromosomes as keys. The values are combined sites across all samples. 
	"""

	sitesCombine = {str(k): pd.DataFrame(columns = ['Chr', 'Sites', 'CgarType', 'Count', 'SiteType', 'RPM']) for k in chrSeq}

	for i in range(len(samples)):
		print(samples[i])
		inputFile = str(samples[i]) + '_readsCatalogue/' + str(samples[i]) + '_allSites_noDups_final.tsv'
		dat = pd.read_csv(inputFile, delimiter = '\t', low_memory = False)

		for j in chrSeq:
			subcomb = sitesCombine[j]
			subdat  = dat.loc[dat.Chr == j]
			appsubs = pd.concat([subcomb, subdat]).reset_index(drop = True)

			sitesCombine[j] = appsubs

	return sitesCombine

def CompareExperiments(combine_RA4, combine_F9, chrSeq):
	"""
	Function to find similar sites at each chromosomes for two experiments

	Output:
		- Dictionary of results for overlapping sites, number that exist in sample 1, 
		and number in sample 2 only 

	"""

	numOverlap = {str(k): 0 for k in chrSeq}
	numDiffRA4 = {str(k): 0 for k in chrSeq}
	numDiffF9  = {str(k): 0 for k in chrSeq}

	for i in chrSeq:
		sub_RA4 = set(combine_RA4[i].Sites)
		sub_F9  = set(combine_F9[i].Sites)

		olapSites = sub_RA4.intersection(sub_F9)
		diffSites_RA4 = sub_RA4 - sub_F9
		diffSites_F9  = sub_F9  - sub_RA4

		numOverlap[i] = len(olapSites)
		numDiffRA4[i] = len(diffSites_RA4)
		numDiffF9[i]  = len(diffSites_F9)

	result = {'TotOverlap': sum(numOverlap.values()), 'TotDiffRA4': sum(numDiffRA4.values()), 
	'TotDiffF9': sum(numDiffF9.values())}

	return result


RA4_samples = ["UD", "RA4", "RA4_PG", "RA4_TCP"]
F9_samples  = ["F9_UD", "F9_D4", "F9_D4_PG", "F9_D4_TCP"]

chrseq = list(range(1, 20, 1))
chrSeq = [format(x, '01d') for x in chrseq]
chrSeq.extend(('X', 'Y', 'MT'))

combine_RA4 = CombineSites(RA4_samples, chrSeq)
combine_F9  = CombineSites(F9_samples, chrSeq)

resultCompareExp = CompareExperiments(combine_RA4, combine_F9, chrSeq)

# Plot Venn Diagram 
plt.figure(figsize=(6,6))
venn2(subsets = [resultCompareExp['TotDiffRA4'], resultCompareExp['TotDiffF9'], resultCompareExp['TotOverlap']], 
	set_labels = (str('RA4'), str('F9')))
plt.savefig("FiguresCompareF9vsRA4/VennOfSitesOverlap.pdf", bbox_inches = 'tight')


# ******************* Compare F9_D4 - UD samples to RA4 - UD samples ****************** # 

RA4_samples = ["UD", "RA4", "RA4_PG", "RA4_TCP"]
F9_samples  = ["F9_UD", "F9_D4", "F9_D4_PG", "F9_D4_TCP"]

chrseq = list(range(1, 20, 1))
chrSeq = [format(x, '01d') for x in chrseq]
chrSeq.extend(('X', 'Y', 'MT'))

# RA4_diff_samples = ["RA4_diff", "RA4_PG_diff", "RA4_TCP_diff"]
# F9_diff_samples = ["F9_D4_diff", "F9_D4_PG_diff", "F9_D4_TCP_diff"]

def CompareDiffSamples(RA4_samples, F9_samples, chrSeq):
	"""
		Function to compare two bed files 

	"""

	allResult = {}
	
	for i in range(1, len(RA4_samples)):
		print(RA4_samples[i])

		input1 = "Results_RA4/UD_vs_" + str(RA4_samples[i]) + "/" + str(RA4_samples[i]) + "_diff_sites.bed"
		input2 = "Results/F9_UD_vs_" + str(F9_samples[i]) + "/" + str(F9_samples[i]) + "_diff_sites.bed"

		dat1 = pd.read_csv(input1, header = None, delimiter = '\t', low_memory = False)
		dat1.columns = ["Chr", "Start", "End"]
		dat2 = pd.read_csv(input2, header = None, delimiter = '\t', low_memory = False)
		dat2.columns = ["Chr", "Start", "End"]

		allResult[str(RA4_samples[i]) + "_Min_UD_vs_" + str(F9_samples[i]) + "_Min_UD"] = {}
		subResult = {}

		olapVec  = []
		diff1Vec = []
		diff2Vec = []

		for j in chrSeq:
			newj = "chr" + str(j)

			subDat1 = set(dat1.loc[dat1.Chr == newj, 'Start'])
			subDat2 = set(dat2.loc[dat2.Chr == newj, 'Start'])

			olap  = subDat1 & subDat2
			diff1 = subDat1 - subDat2
			diff2 = subDat2 - subDat1 

			olapVec.append(len(olap))
			diff1Vec.append(len(diff1))
			diff2Vec.append(len(diff2))

		subResult = {'Intersection': sum(olapVec), str(RA4_samples[i]) + "_Min_UD": sum(diff1Vec), 
		str(F9_samples[i]) + "_Min_UD": sum(diff2Vec)}

		allResult[str(RA4_samples[i]) + "_Min_UD_vs_" + str(F9_samples[i]) + "_Min_UD"] = subResult

	return allResult

compTreat = CompareDiffSamples(RA4_samples, F9_samples, chrSeq)

for key, value in compTreat.items():
	print(key)
	
	figName = "FiguresCompareF9vsRA4/" + str(key) + ".pdf"

	valKeys = list(value.keys())
	
	# Plot Venn Diagram 
	plt.figure(figsize=(6,6))
	venn2(subsets = [value[valKeys[1]], value[valKeys[2]], value[valKeys[0]]], 
		set_labels = (valKeys[1], valKeys[2]))
	plt.savefig(figName, bbox_inches = 'tight')
	

# ******************* Compare F9_D4 samples to RA4 samples ****************** # 

RA4_samples = ["UD", "RA4", "RA4_PG", "RA4_TCP"]
F9_samples  = ["F9_UD", "F9_D4", "F9_D4_PG", "F9_D4_TCP"]

chrseq = list(range(1, 20, 1))
chrSeq = [format(x, '01d') for x in chrseq]
chrSeq.extend(('X', 'Y', 'MT'))

def CompareDiffEachSamples(RA4_samples, F9_samples, chrSeq):
	"""
		Function to compare two bed files 

	"""

	allResult = {}
	
	for i in range(0, len(RA4_samples)):
		print(RA4_samples[i])

		input1 = str(RA4_samples[i]) + "_readsCatalogue/" + str(RA4_samples[i]) + "_allSites_noDups_final.bed"
		input2 = str(F9_samples[i]) + "_readsCatalogue/" + str(F9_samples[i]) + "_allSites_noDups_final.bed"

		dat1 = pd.read_csv(input1, header = None, delimiter = '\t', low_memory = False)
		dat1.columns = ["Chr", "Start", "End", "Count", "RPM"]
		dat2 = pd.read_csv(input2, header = None, delimiter = '\t', low_memory = False)
		dat2.columns = ["Chr", "Start", "End", "Count", "RPM"]

		allResult[str(RA4_samples[i]) + "_vs_" + str(F9_samples[i])] = {}
		subResult = {}

		olapVec  = []
		diff1Vec = []
		diff2Vec = []

		for j in chrSeq:
			newj = "chr" + str(j)

			subDat1 = set(dat1.loc[dat1.Chr == newj, 'Start'])
			subDat2 = set(dat2.loc[dat2.Chr == newj, 'Start'])

			olap  = subDat1 & subDat2
			diff1 = subDat1 - subDat2
			diff2 = subDat2 - subDat1 

			olapVec.append(len(olap))
			diff1Vec.append(len(diff1))
			diff2Vec.append(len(diff2))

		subResult = {'Intersection': sum(olapVec), str(RA4_samples[i]): sum(diff1Vec), 
		str(F9_samples[i]): sum(diff2Vec)}

		allResult[str(RA4_samples[i]) + "_vs_" + str(F9_samples[i])] = subResult

	return allResult

compTreat = CompareDiffEachSamples(RA4_samples, F9_samples, chrSeq)

for key, value in compTreat.items():
	
	figName = "FiguresCompareF9vsRA4/" + str(key) + ".pdf"

	valKeys = list(value.keys())
	
	# Plot Venn Diagram 
	plt.figure(figsize=(6,6))
	venn2(subsets = [value[valKeys[1]], value[valKeys[2]], value[valKeys[0]]], 
		set_labels = (valKeys[1], valKeys[2]))
	plt.savefig(figName, bbox_inches = 'tight')




















