import numpy as np
import pandas as pd
import pickle as pk

# *********************************** Find Union of all Enhancers ********************************************* #
def DiffCounts(fileName, allSamples, sampleName):

	smplEnh = {'F9_UD':{}, 'F9_D4': {}, 'F9_D4_PG': {}, 'F9_D4_TCP': {}}

	for i in range(len(sampleName)):
		inputFile = str(allSamples[i]) + "/" + str(sampleName[i]) + "_" + str(fileName)
		annoDf = pd.read_csv(inputFile, delimiter = '\t', header = None, low_memory = False, usecols = list(range(9)))
		annoDf.columns = ['ChrSite', 'SiteStart', 'SiteEnd', 'Counts', 'RPM', 'ChrEnh', 'EnhStart', 'EnhEnd', 'EnhName']
		smplEnh[sampleName[i]] = set(annoDf.EnhName)

	smplEnhLen = {k: len(v) for k, v in smplEnh.items()}
	listEnhAll = list(smplEnh.values())
	enhUnionSet = set().union(*listEnhAll)
	enhUnion = list(enhUnionSet)
	enhUnion.sort()

	# Create data frame for the union set
	initData  = {k: np.repeat(0, len(enhUnion)) for k, v in smplEnh.items()}
	smplEnhDf = pd.DataFrame.from_dict(initData, orient='index')
	smplEnhDf.columns = enhUnion

	# *********************************** Input values in DF Enhancers ********************************************* #
	for i in range(len(sampleName)):
		print(sampleName[i])
		inputFile = str(allSamples[i]) + "/" + str(sampleName[i]) + "_" + str(fileName)
		annoDf = pd.read_csv(inputFile, delimiter = '\t', header = None, low_memory = False, usecols = list(range(9)))
		annoDf.columns = ['ChrSite', 'SiteStart', 'SiteEnd', 'Counts', 'RPM', 'ChrEnh', 'EnhStart', 'EnhEnd', 'EnhName']

		setInUnion = smplEnh[sampleName[i]]

		# # Average values for RPM in regions
		# for j in setInUnion: 
		# 	tmp = annoDf.loc[annoDf.EnhName == str(j)]
		# 	avg = np.mean(tmp.RPM.values)
		# 	smplEnhDf.loc[sampleName[i], j] = avg

		# Sum values for counts in regions
		for j in setInUnion: 
			tmp = annoDf.loc[annoDf.EnhName == str(j)]
			sumval = sum(tmp.Counts.values)
			smplEnhDf.loc[sampleName[i], j] = sumval

	# # *********************************** Compute Log2FoldChange ********************************************* #
	# offSet = 0.0001
	# smplEnhDf.loc['RA4/UD'] = np.round(np.log2(smplEnhDf.loc['RA4'] + offSet)- np.log2(smplEnhDf.loc['UD'] + offSet), 4)
	# smplEnhDf.loc['RA4_PG/UD'] = np.round(np.log2(smplEnhDf.loc['RA4_PG'] + offSet)- np.log2(smplEnhDf.loc['UD'] + offSet), 4)
	# smplEnhDf.loc['RA4_TCP/UD'] = np.round(np.log2(smplEnhDf.loc['RA4_TCP'] + offSet)- np.log2(smplEnhDf.loc['UD'] + offSet), 4)

	# *********************************** Compute Difference Counts ********************************************* #
	offSet = 0.0001
	smplEnhDf = smplEnhDf + offSet
	smplEnhDf.loc['F9_D4_Min_UD'] = smplEnhDf.loc['F9_D4'] - smplEnhDf.loc['F9_UD']
	smplEnhDf.loc['F9_D4_PG_Min_UD'] = smplEnhDf.loc['F9_D4_PG'] - smplEnhDf.loc['F9_UD']
	smplEnhDf.loc['F9_D4_TCP_Min_UD'] = smplEnhDf.loc['F9_D4_TCP'] - smplEnhDf.loc['F9_UD']

	# Compute Percent Change
	smplEnhDf.loc['F9_D4_Min_UD_PerChange'] = np.round((smplEnhDf.loc['F9_D4'] - smplEnhDf.loc['F9_UD'])/smplEnhDf.loc['F9_UD']*100, 2)
	smplEnhDf.loc['F9_D4_PG_Min_UD_PerChange'] = np.round((smplEnhDf.loc['F9_D4_PG'] - smplEnhDf.loc['F9_UD'])/smplEnhDf.loc['F9_UD']*100, 2)
	smplEnhDf.loc['F9_D4_TCP_Min_UD_PerChange'] = np.round((smplEnhDf.loc['F9_D4_TCP'] - smplEnhDf.loc['F9_UD'])/smplEnhDf.loc['F9_UD']*100, 2)

	return(smplEnhDf)

allSamples = ["F9_UD_readsCatalogue", "F9_D4_readsCatalogue", "F9_D4_PG_readsCatalogue", "F9_D4_TCP_readsCatalogue"]
sampleName = ['F9_UD', 'F9_D4', 'F9_D4_PG', 'F9_D4_TCP']
chrseq = list(range(1, 20, 1))
chrSeq = [format(x, '01d') for x in chrseq]
chrSeq.extend(('X', 'Y', 'MT'))

# ************* Computing FC for LSD1 Enhancers ******************** # 
fileName = "allSites_noDups_LSD1_intersect.bed"
lsdFC = DiffCounts(fileName, allSamples, sampleName)
lsdFC = lsdFC.transpose()

fname  = "Results/FoldChanges/FC_allSamples_LSD1_counts.tsv"
lsdFC.to_csv(fname, sep = '\t', index = True)

# ************* Computing FC for ESC J1 Enhancers ******************* # 
fileName = "allSites_noDups_ESC_J1_enhancers_intersect.bed"
escFC = DiffCounts(fileName, allSamples, sampleName)
escFC = escFC.transpose()

fname = "Results/FoldChanges/FC_allSamples_ESC_J1_counts.tsv"
escFC.to_csv(fname, sep = '\t', index = True)

# ************* Computing FC for all other regions ******************* # 

def Log2FC_NoEnhance(fileName, allSamples, sampleName):

	smplEnh = {'F9_UD':{}, 'F9_D4': {}, 'F9_D4_PG': {}, 'F9_D4_TCP': {}}

	for i in range(len(sampleName)):
		inputFile = str(allSamples[i]) + "/" + str(sampleName[i]) + "_" + str(fileName)
		annoDf = pd.read_csv(inputFile, delimiter = '\t', header = 0, low_memory = False)
		annoDf = annoDf.rename(columns = {'V4': 'rpm'})
		smplEnh[sampleName[i]] = set(annoDf.ENSEMBL)

	smplEnhLen = {k: len(v) for k, v in smplEnh.items()}
	listEnhAll = list(smplEnh.values())
	enhUnionSet = set().union(*listEnhAll)
	enhUnion = np.asarray(list(enhUnionSet))
	enhUnion = enhUnion[enhUnion != 'nan']
	enhUnion.sort()

	# Create data frame for the union set
	gnRpmDict = {k: {} for k, v in smplEnh.items()}
	# smplEnhDf = pd.DataFrame.from_dict(initData, orient='index')
	# smplEnhDf.columns = enhUnion

	# *********************************** Input values in DF Enhancers ********************************************* #
	for i in range(len(sampleName)):
		print(sampleName[i])
		inputFile = str(allSamples[i]) + "/" + str(sampleName[i]) + "_" + str(fileName)
		annoDf = pd.read_csv(inputFile, delimiter = '\t', header = 0, low_memory = False)
		annoDf = annoDf.rename(columns = {'V4': 'rpm'})
		setInUnion = np.asarray(list(smplEnh[sampleName[i]]))
		setInUnion = setInUnion[setInUnion != 'nan']
		setInUnion.sort()

		smplGnRpm = np.repeat(0, len(enhUnion)).astype(np.float64)
		for j in setInUnion: 
			# print(j)
			tmp = annoDf.loc[annoDf.ENSEMBL == str(j)]
			avg = np.round(np.mean(tmp.rpm.values), 4)
			ind = np.where(enhUnion == j)
			smplGnRpm[ind] = avg

		gnRpmDict[sampleName[i]] = smplGnRpm

	# *********************************** Compute Log2FoldChange ********************************************* #
	offSet = 0.0001
	gnRpmDict['RA4/UD'] = np.round(np.log2(gnRpmDict['RA4'] + offSet)- np.log2(gnRpmDict['UD'] + offSet), 4)
	gnRpmDict['RA4_PG/UD'] = np.round(np.log2(gnRpmDict['RA4_PG'] + offSet)- np.log2(gnRpmDict['UD'] + offSet), 4)
	gnRpmDict['RA4_TCP/UD'] = np.round(np.log2(gnRpmDict['RA4_TCP'] + offSet)- np.log2(gnRpmDict['UD'] + offSet), 4)

	dictRes = {'gnRpmDict': gnRpmDict, 'enhUnion': enhUnion}

	return(dictRes)

allSamples = ["UD_readsCatalogue", "RA4_readsCatalogue", "RA4_PG_readsCatalogue", "RA4_TCP_readsCatalogue"]
sampleName = ['UD', 'RA4', 'RA4_PG', 'RA4_TCP']
chrseq = list(range(1, 20, 1))
chrSeq = [format(x, '01d') for x in chrseq]
chrSeq.extend(('X', 'Y', 'MT'))

fileName = "sitesAnnotation_noEnhancers.tsv"
regFC = Log2FC_NoEnhance(fileName, allSamples, sampleName)
regFCdf = pd.DataFrame.from_dict(regFC['gnRpmDict'], orient = 'index')
regFCdf.columns = list(regFC['enhUnion'])
regFCdf = regFCdf.transpose()

fname = "Results/FoldChanges/FC_allSamples_noEnhancers.tsv"
regFCdf.to_csv(fname, sep = '\t', index = True)

# ************* Computing FC for promoters ******************* # 

def Log2FC_Promoters(fileName, allSamples, sampleName):

	smplEnh = {'UD':{}, 'RA4': {}, 'RA4_PG': {}, 'RA4_TCP': {}}

	for i in range(len(sampleName)):
		inputFile = str(allSamples[i]) + "/" + str(sampleName[i]) + "_" + str(fileName)
		annoDf = pd.read_csv(inputFile, delimiter = '\t', header = None, low_memory = False, usecols = [0, 1, 2, 3, 10])
		annoDf.columns = ['Chr', 'Start', 'End', 'Rpm', 'Ens_id']
		smplEnh[sampleName[i]] = set(annoDf.Ens_id)

	smplEnhLen = {k: len(v) for k, v in smplEnh.items()}
	listEnhAll = list(smplEnh.values())
	enhUnionSet = set().union(*listEnhAll)
	enhUnion = np.asarray(list(enhUnionSet))
	enhUnion = enhUnion[enhUnion != 'nan']
	enhUnion.sort()

	# Create data frame for the union set
	gnRpmDict = {k: {} for k, v in smplEnh.items()}
	# smplEnhDf = pd.DataFrame.from_dict(initData, orient='index')
	# smplEnhDf.columns = enhUnion

	# *********************************** Input values in DF Enhancers ********************************************* #
	for i in range(len(sampleName)):
		print(sampleName[i])
		inputFile = str(allSamples[i]) + "/" + str(sampleName[i]) + "_" + str(fileName)
		annoDf = pd.read_csv(inputFile, delimiter = '\t', header = None, low_memory = False, usecols = [0, 1, 2, 3, 10])
		annoDf.columns = ['Chr', 'Start', 'End', 'Rpm', 'Ens_id']
		setInUnion = np.asarray(list(smplEnh[sampleName[i]]))
		setInUnion = setInUnion[setInUnion != 'nan']
		setInUnion.sort()

		smplGnRpm = np.repeat(0, len(enhUnion)).astype(np.float64)
		for j in setInUnion: 
			# print(j)
			tmp = annoDf.loc[annoDf.Ens_id == str(j)]
			avg = np.round(np.mean(tmp.Rpm.values), 4)
			ind = np.where(enhUnion == j)
			smplGnRpm[ind] = avg

		gnRpmDict[sampleName[i]] = smplGnRpm

	# *********************************** Compute Log2FoldChange ********************************************* #
	offSet = 0.0001
	gnRpmDict['RA4/UD'] = np.round(np.log2(gnRpmDict['RA4'] + offSet)- np.log2(gnRpmDict['UD'] + offSet), 4)
	gnRpmDict['RA4_PG/UD'] = np.round(np.log2(gnRpmDict['RA4_PG'] + offSet)- np.log2(gnRpmDict['UD'] + offSet), 4)
	gnRpmDict['RA4_TCP/UD'] = np.round(np.log2(gnRpmDict['RA4_TCP'] + offSet)- np.log2(gnRpmDict['UD'] + offSet), 4)

	dictRes = {'gnRpmDict': gnRpmDict, 'enhUnion': enhUnion}

	return(dictRes)

allSamples = ["UD_readsCatalogue", "RA4_readsCatalogue", "RA4_PG_readsCatalogue", "RA4_TCP_readsCatalogue"]
sampleName = ['UD', 'RA4', 'RA4_PG', 'RA4_TCP']
chrseq = list(range(1, 20, 1))
chrSeq = [format(x, '01d') for x in chrseq]
chrSeq.extend(('X', 'Y', 'MT'))

fileName = "allSites_noDups_noenhancer_promoters_intersect.bed"
regFC = Log2FC_Promoters(fileName, allSamples, sampleName)
regFCdf = pd.DataFrame.from_dict(regFC['gnRpmDict'], orient = 'index')
regFCdf.columns = list(regFC['enhUnion'])
regFCdf = regFCdf.transpose()

fname = "Results/FoldChanges/FC_allSamples_promoters.tsv"
regFCdf.to_csv(fname, sep = '\t', index = True)


