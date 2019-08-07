import matplotlib
matplotlib.use('Agg')
import numpy as np
import pandas as pd
import pickle as pk
import os
# import matplotlib.pyplot as plt
# from matplotlib_venn import venn2
# from matplotlib_venn import venn3, venn3_circles
from matplotlib_venn import venn3_unweighted 
	
# ************************* Find overlapping sites for every pair of samples ******************************* #
def CompareSamples(allSamples, sampleName, chrSeq):
	"""
		Function to find overlapping and differenting sites amongst all pairs of samples

		Output: 
			- SaveData/SampleSitesCompare.pkl contains a dictionary of a summary of the comparisons
			- SaveData/SampleSitesCompareList.pkl contains a dictionary of all the sites that are overlapping and 
				differentiating between all pair of samples 
	"""

	compareSite = {}
	compareSiteList = {}

	for i in range(len(allSamples)):
		
		tt = i + 1

		while tt < len(allSamples):
			smp1 = pd.read_csv(str(allSamples[i]) + "/" + sampleName[i] + "_allSites_noDups_final.tsv", delimiter = '\t', low_memory = False)
			smp2 = pd.read_csv(str(allSamples[tt]) + "/" + sampleName[tt] + "_allSites_noDups_final.tsv", delimiter = '\t', low_memory = False)

			vennStats = {str(sampleName[i]) + '_diff': np.zeros(len(chrSeq)), str(sampleName[tt]) + '_diff' : np.zeros(len(chrSeq)), 
			'Intersection': np.zeros(len(chrSeq))}
			siteStats = {str(sampleName[i]) + '_diff': {}, str(sampleName[tt]) + '_diff' : {}, 'Intersection': {}}

			for j in range(len(chrSeq)):
				print(j)
				dat1 = smp1.loc[smp1['Chr'] == chrSeq[j]]
				dat2 = smp2.loc[smp2['Chr'] == chrSeq[j]]

				site1 = set(dat1['Sites'].values)
				site2 = set(dat2['Sites'].values)

				olapSites = site1 & site2
				olapDiff1 = site1 - site2
				olapDiff2 = site2 - site1

				siteStats[str(sampleName[i]) + '_diff'][chrSeq[j]]  = olapDiff1
				siteStats[str(sampleName[tt]) + '_diff'][chrSeq[j]] = olapDiff2
				siteStats["Intersection"][chrSeq[j]] = olapSites

				vennStats[str(sampleName[i]) + '_diff'][j]  = len(olapDiff1)
				vennStats[str(sampleName[tt]) + '_diff'][j] = len(olapDiff2)
				vennStats["Intersection"][j] = len(olapSites)

			vennStats['Total'] = {str(sampleName[i]) + '_diff': np.sum(vennStats[str(sampleName[i]) + '_diff']), 
			str(sampleName[tt]) + '_diff': np.sum(vennStats[str(sampleName[tt]) + '_diff']), 
			"Intersection": np.sum(vennStats["Intersection"])}

			dictName = str(sampleName[i]) + '_vs_' + str(sampleName[tt])
			compareSite[dictName] = vennStats
			compareSiteList[dictName] = siteStats

			tt = tt + 1

	f = open("SaveData/SampleSitesCompare.pkl", "wb")
	pk.dump(compareSite, f)
	f.close()

	f = open("SaveData/SampleSitesCompareList.pkl", "wb")
	pk.dump(compareSiteList, f)
	f.close()

allSamples = ["F9_UD_readsCatalogue", "F9_D4_readsCatalogue", "F9_D4_PG_readsCatalogue", "F9_D4_TCP_readsCatalogue"]
sampleName = ['F9_UD', 'F9_D4', 'F9_D4_PG', 'F9_D4_TCP']
chrseq = list(range(1, 20, 1))
chrSeq = [format(x, '01d') for x in chrseq]
chrSeq.extend(('X', 'Y', 'MT'))

CompareSamples(allSamples, sampleName, chrSeq)

# ************************* Find overlapping sites TCP + PG - UD ******************************* #
def ComparePGplusTCPminUD(smp1, smp2, smp3):
	"""
		Function to compare sites between overlapping sites in TCP and PG to those of UD

		Output: 
			- SaveData/SampleSitesCompare.pkl contains a dictionary of a summary of the comparisons
			- SaveData/SampleSitesCompareList.pkl contains a dictionary of all the sites that are overlapping and 
				differentiating between all pair of samples 
	"""

	compareSite = pk.load(open("SaveData/SampleSitesCompare.pkl", "rb"))
	compareSiteList = pk.load(open("SaveData/SampleSitesCompareList.pkl", "rb"))

	uniSites = {}

	count = 0

	for j in range(len(chrSeq)):
		dat1 = smp1.loc[smp1['Chr'] == chrSeq[j]]
		dat2 = smp2.loc[smp2['Chr'] == chrSeq[j]]

		site1 = set(dat1['Sites'].values)
		site2 = set(dat2['Sites'].values)

		uniSites[chrSeq[j]] = site1 & site2
		# print(len(uniSites[chrSeq[j]]))

		count = count + len(uniSites[chrSeq[j]])

	siteStats = {'F9_D4_PG_TCP_Inter_diff': {}, 'F9_UD_diff' : {}, 'Intersection': {}}
	raud = {'F9_D4_PG_TCP_Inter_diff': np.zeros(len(chrSeq)), 'F9_UD_diff' : np.zeros(len(chrSeq)), 
	'Intersection': np.zeros(len(chrSeq))}

	for j in range(len(chrSeq)):
		dat3 = smp3.loc[smp3['Chr'] == chrSeq[j]]

		site1 = uniSites[chrSeq[j]]
		site2 = set(dat3['Sites'].values)

		olapSites = site1 & site2
		olapDiff1 = site1 - site2
		olapDiff2 = site2 - site1
	
		siteStats['F9_D4_PG_TCP_Inter_diff'][chrSeq[j]]  = olapDiff1
		siteStats['F9_UD_diff'][chrSeq[j]] = olapDiff2
		siteStats["Intersection"][chrSeq[j]] = olapSites

		raud['Intersection'][j] = len(olapSites)
		raud['F9_D4_PG_TCP_Inter_diff'][j] = len(olapDiff1)
		raud['F9_UD_diff'][j] = len(olapDiff2)

	raud['Total'] = {'F9_D4_PG_TCP_Inter_diff': np.sum(raud['F9_D4_PG_TCP_Inter_diff']), 'F9_UD_diff': np.sum(raud['F9_UD_diff']), 
		"Intersection": np.sum(raud["Intersection"])}

	compareSite['F9_D4_PG_TCP_vs_F9_UD'] = raud
	compareSiteList['F9_D4_PG_TCP_vs_F9_UD'] = siteStats

	f = open("SaveData/SampleSitesCompare.pkl", "wb")
	pk.dump(compareSite, f)
	f.close()

	f = open("SaveData/SampleSitesCompareList.pkl", "wb")
	pk.dump(compareSiteList, f)
	f.close()

######### Run function ########
smp1 = pd.read_csv("F9_D4_PG_readsCatalogue/F9_D4_PG_allSites_noDups_final.tsv", delimiter = '\t', low_memory = False)
smp2 = pd.read_csv("F9_D4_TCP_readsCatalogue/F9_D4_TCP_allSites_noDups_final.tsv", delimiter = '\t', low_memory = False)
smp3 = pd.read_csv("F9_UD_readsCatalogue/F9_UD_allSites_noDups_final.tsv", delimiter = '\t', low_memory = False)

ComparePGplusTCPminUD(smp1, smp2, smp3)

# ************************* Save Compare data into bed files ******************************* #
def createFolder(directory):
    ''' Create directory 
        Input: name of directory 
    '''
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print('Error: Creating directory. ' + directory)

def SaveCompare(compareSiteList):	
	""" 
		Function to save comparison files to a bed file

		Output:
			- Results/F9_UD_vs_F9_D4_sites.bed
	"""
	createFolder('Results')

	for keys in compareSiteList:
		filename = 'Results/' + str(keys)
		createFolder(filename)

		for kk in compareSiteList[keys]:
			tmp = compareSiteList[keys][kk]
			df  = pd.DataFrame(columns = ['Chr', 'Start', 'End'])
			print(keys)
			print(kk)
			for ii in chrSeq:
				cc = np.repeat("chr" + str(ii), len(tmp[ii]))
				tmparr = list(tmp[ii])
				tmparr = sorted(tmparr)
				tmpend = list((np.asarray(tmparr) + np.ones(len(tmparr))).astype('int'))
				df_tmp = pd.DataFrame({'Chr': cc, 'Start': tmparr, 'End': tmpend}, columns = ['Chr', 'Start', 'End'])
				df = pd.concat([df, df_tmp]).reset_index(drop = True)
			rsltName = str(filename) + '/' + str(kk) + '_sites.bed'
			df.to_csv(rsltName, sep = '\t', index = False, header = False)

compareSiteList = pk.load(open("SaveData/SampleSitesCompareList.pkl", "rb"))

SaveCompare(compareSiteList)

# ******************* Compare sites amongst all RA_diffs ******************************* #

## Create function for overlaps between two dictionaries
def FindOverlap(dict1, dict2, dictNames, seqChr):
	''' Function to compute overlap in sets 
	Input: 
		dict1, dict2 = dictionary containing sets in each seqChr
		names  = array containing strings for names of set1 and set2
		seqChr = array of chromosomes sequence
	Output: 
		dictionary with se1 differences, set2 differences, and intersection
	'''

	setNames = [str(i) + "_diff" for i in dictNames]
	siteList = {str(setNames[0]): {}, str(setNames[1]): {}, 'Intersection': {}}
	siteNum  = {str(setNames[0]): np.zeros(len(chrSeq)), str(setNames[1]) : np.zeros(len(chrSeq)), 
	'Intersection': np.zeros(len(chrSeq))}

	for j in range(len(seqChr)):

		site1 = dict1[str(seqChr[j])]
		site2 = dict2[str(seqChr[j])]

		olapSites = site1 & site2
		olapDiff1 = site1 - site2
		olapDiff2 = site2 - site1

		siteList[str(setNames[0])][chrSeq[j]] = olapDiff1
		siteList[str(setNames[1])][chrSeq[j]] = olapDiff2
		siteList["Intersection"][chrSeq[j]]  = olapSites

		siteNum[str(setNames[0])][j] = int(len(olapDiff1))
		siteNum[str(setNames[1])][j] = int(len(olapDiff2))
		siteNum["Intersection"][j]  = int(len(olapSites))

	siteNum['Total'] = {str(setNames[0]): int(np.sum(siteNum[setNames[0]])), 
	str(setNames[1]): int(np.sum(siteNum[setNames[1]])), 'Intersection': int(np.sum(siteNum['Intersection']))}

	return({'siteNum': siteNum, 'siteList': siteList})

compareSiteList = pk.load(open("SaveData/SampleSitesCompareList.pkl", "rb"))

allSamples = ["F9_UD_readsCatalogue", "F9_D4_readsCatalogue", "F9_D4_PG_readsCatalogue", "F9_D4_TCP_readsCatalogue"]
sampleName = ['F9_UD', 'F9_D4', 'F9_D4_PG', 'F9_D4_TCP']

chrseq = list(range(1, 20, 1))
chrSeq = [format(x, '01d') for x in chrseq]
chrSeq.extend(('X', 'Y', 'MT'))

# ******************************************************************** #

dict1 = compareSiteList['F9_UD_vs_F9_D4']['F9_D4_diff']
dict2 = compareSiteList['F9_UD_vs_F9_D4_PG']['F9_D4_PG_diff']
d4ud_pgud = FindOverlap(dict1, dict2, dictNames = ["F9_D4_Min_UD", "F9_D4_PG_Min_UD"], seqChr = chrSeq)

A_E_list = d4ud_pgud['siteList']['F9_D4_Min_UD_diff']
A_E_len  = np.sum([len(values) for keys, values in A_E_list.items()])

B_F_list = d4ud_pgud['siteList']['F9_D4_PG_Min_UD_diff']
B_F_len  = np.sum([len(values) for keys, values in B_F_list.items()])

D_G_list = d4ud_pgud['siteList']['Intersection']
D_G_len  = np.sum([len(values) for keys, values in D_G_list.items()])

# ******************************************************************** #

dict1 = compareSiteList['F9_UD_vs_F9_D4']['F9_D4_diff']
dict2 = compareSiteList['F9_UD_vs_F9_D4_TCP']['F9_D4_TCP_diff']
d4ud_tcpud = FindOverlap(dict1, dict2, dictNames = ["F9_D4_Min_UD", "F9_D4_TCP_Min_UD"], seqChr = chrSeq)

A_D_list = d4ud_tcpud['siteList']['F9_D4_Min_UD_diff']
A_D_len  = np.sum([len(values) for keys, values in A_D_list.items()])

C_F_list = d4ud_tcpud['siteList']['F9_D4_TCP_Min_UD_diff']
C_F_len  = np.sum([len(values) for keys, values in C_F_list.items()])

E_G_list = d4ud_tcpud['siteList']['Intersection']
E_G_len  = np.sum([len(values) for keys, values in E_G_list.items()])

# ******************************************************************** #

dict1 = compareSiteList['F9_UD_vs_F9_D4_PG']['F9_D4_PG_diff']
dict2 = compareSiteList['F9_UD_vs_F9_D4_TCP']['F9_D4_TCP_diff']
pgud_tcpud = FindOverlap(dict1, dict2, dictNames = ["F9_D4_PG_UD", "F9_D4_TCP_UD"], seqChr = chrSeq)

B_D_list = pgud_tcpud['siteList']['F9_D4_PG_UD_diff']
B_D_len  = np.sum([len(values) for keys, values in B_D_list.items()])

C_E_list = pgud_tcpud['siteList']['F9_D4_TCP_UD_diff']
C_E_len  = np.sum([len(values) for keys, values in C_E_list.items()])

G_F_list = pgud_tcpud['siteList']['Intersection']
G_F_len  = np.sum([len(values) for keys, values in G_F_list.items()])

# ******************************************************************** #

dict1 = D_G_list
dict2 = compareSiteList['F9_UD_vs_F9_D4_TCP']['F9_D4_TCP_diff']
findG = FindOverlap(dict1, dict2, dictNames = ['D', 'C'], seqChr = chrSeq)

G_len = np.sum([len(values) for keys, values in findG['siteList']['Intersection'].items()])

G = findG['siteNum']['Total']['Intersection']
E = E_G_len - G
A = A_E_len - E
F = G_F_len - G
B = B_F_len - F 
D = D_G_len - G
C = C_E_len - E

plt.figure(figsize = (6,6))
out = venn3_unweighted(subsets = (A, B, D, C, E, F, G), set_labels = ('F9_D4_Min_UD', 'F9_D4_PG_Min_UD', 'F9_D4_TCP_Min_UD'))
out.get_patch_by_id('100').set_alpha(1.0)
for text in out.set_labels:
	text.set_fontsize(6)
for text in out.subset_labels:
	text.set_fontsize(6)
plt.savefig("Figures/F9_D4_Min_UD_diffVenn.pdf", bbox_inches = 'tight')

#### Save sites into bed files ####

def dict2df(dicty, seqChr, filename):
	df = pd.DataFrame(columns = ['Chr', 'Start', 'End'])
	for ii in seqChr:
		cc = np.repeat("chr" + str(ii), len(dicty[ii]))
		tmparr = list(dicty[str(ii)])
		tmparr = sorted(tmparr)
		tmpend = list((np.asarray(tmparr) + np.ones(len(tmparr))).astype('int'))
		df_tmp = pd.DataFrame({'Chr': cc, 'Start': tmparr, 'End': tmpend}, columns = ['Chr', 'Start', 'End'])
		df = pd.concat([df, df_tmp]).reset_index(drop = True)
	rsltName = str(filename) + '_sites.bed'	
	df.to_csv(rsltName, sep = '\t', index = False, header = False)
		
dict1 = A_E_list
dict2 = C_E_list
Adiff = FindOverlap(dict1, dict2, dictNames = ['A', 'C'], seqChr = chrSeq)

# A
dict2df(dicty = Adiff['siteList']['A_diff'], seqChr = chrSeq, 
	filename = "Results/F9_D4_D_vs_PG_D_vs_TCP_D/F9_D4_Min_UD_only")
# E
dict2df(dicty = Adiff['siteList']['Intersection'], seqChr = chrSeq, 
	filename = "Results/F9_D4_D_vs_PG_D_vs_TCP_D/F9_D4_Min_UD_vs_TCP_Min_UD_intersection")

dict1 = B_F_list
dict2 = C_F_list
Bdiff = FindOverlap(dict1, dict2, dictNames = ['B', 'C'], seqChr = chrSeq)

# B
dict2df(dicty = Bdiff['siteList']['B_diff'], seqChr = chrSeq, 
	filename = "Results/F9_D4_D_vs_PG_D_vs_TCP_D/F9_D4_PG_Min_UD_only")
# F
dict2df(dicty = Bdiff['siteList']['Intersection'], seqChr = chrSeq, 
	filename = "Results/F9_D4_D_vs_PG_D_vs_TCP_D/F9_D4_PG_Min_UD_vs_TCP_Min_UD_intersection")
# C 
dict2df(dicty = Bdiff['siteList']['C_diff'], seqChr = chrSeq, 
	filename = "Results/F9_D4_D_vs_PG_D_vs_TCP_D/F9_D4_TCP_Min_UD_only")

dict1 = A_D_list
dict2 = B_D_list
Cdiff = FindOverlap(dict1, dict2, dictNames = ['A', 'B'], seqChr = chrSeq)

#D
dict2df(dicty = Cdiff['siteList']['Intersection'], seqChr = chrSeq, 
	filename = "Results/F9_D4_D_vs_PG_D_vs_TCP_D/F9_D4_Min_UD_vs_PG_Min_UD_intersection")
#G
dict2df(dicty = findG['siteList']['Intersection'], seqChr = chrSeq, 
	filename = "Results/F9_D4_D_vs_PG_D_vs_TCP_D/F9_D4_Min_UD_vs_PG_Min_UD_vs_TCP_Min_UD_intersection")




