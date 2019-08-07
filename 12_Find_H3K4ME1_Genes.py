import numpy as np
import pandas as pd
import pickle as pk

def DiffCounts_Sites(fileName, allSamples, sampleName):

	smplEnh = {str(k): {} for k in sampleName}
	annoDfPos = {str(k): {} for k in sampleName}
	annoDfTot = {str(k): {} for k in sampleName}

	for i in range(len(sampleName)):
		inputFile = str(allSamples[i]) + "/" + str(sampleName[i]) + "_" + str(fileName)
		annoDf = pd.read_csv(inputFile, delimiter = '\t', header = None, low_memory = False, usecols = list(range(9)))
		annoDf.columns = ['ChrSite', 'SiteStart', 'SiteEnd', 'Counts', 'RPM', 'ChrEnh', 'EnhStart', 'EnhEnd', 'EnhName' , 'Ens', 'H3Change']
		annoDf = annoDf.drop(['SiteEnd', 'RPM', 'ChrEnh','EnhStart', 'EnhEnd', 'Ens'], axis = 1)
		smplEnh[sampleName[i]] = set(annoDf.EnhName)
		annoDfTot[sampleName[i]] = len(annoDf)
		annoDfPos[sampleName[i]] = annoDf.loc[annoDf.H3Change > -1, :]

	# Find length of H3K4ME1 genes found
	smplEnhLen = {k: len(v) for k, v in smplEnh.items()}

	# Find length of H3K4ME1 sites found
	annoDfPosLen = {k: len(v) for k, v in annoDfPos.items()}

	ud_sample = annoDfPos['F9_UD']

	joinname  = ['F9_D4_vs_UD', 'F9_D4_PG_vs_UD', 'F9_D4_TCP_vs_UD']
	joinlist  = {str(k): {} for k in joinname}
	joinlist  = {'F9_D4_vs_UD': {}, 'F9_D4_PG_vs_UD': {}, 'F9_D4_TCP_vs_UD': {}}

	for i in range(1, len(sampleName)):
		tmp = pd.merge(ud_sample, annoDfPos[sampleName[i]], how = "outer", on = ['ChrSite', 'SiteStart'])
		tmp[['Counts_y', 'Counts_x']] = tmp[['Counts_y', 'Counts_x']].fillna(value = 0.0)
		tmp.loc[:, str(sampleName[i]) + '_Min_UD'] = tmp.loc[:, 'Counts_y'] - tmp.loc[:, 'Counts_x']
		countList = [sum(tmp.iloc[:, 8] > 0), sum(tmp.iloc[:, 8] <= 0), len(tmp)]
		joinlist[joinname[i-1]] = countList

	return(smplEnhDf)

allSamples = ["F9_UD_readsCatalogue", "F9_D4_readsCatalogue", "F9_D4_PG_readsCatalogue", "F9_D4_TCP_readsCatalogue"]
sampleName = ['F9_UD', 'F9_D4', 'F9_D4_PG', 'F9_D4_TCP']
chrseq = list(range(1, 20, 1))
chrSeq = [format(x, '01d') for x in chrseq]
chrSeq.extend(('X', 'Y', 'MT'))

fileName = "allSites_noDups_H3K4ME1_intersect.bed"
lsdFC = DiffCounts_Sites(fileName, allSamples, sampleName)

def DiffCounts(fileName, allSamples, sampleName):

	smplEnh = {str(k): {} for k in sampleName}

	for i in range(len(sampleName)):
		inputFile = str(allSamples[i]) + "/" + str(sampleName[i]) + "_" + str(fileName)
		annoDf = pd.read_csv(inputFile, delimiter = '\t', header = None, low_memory = False)
		annoDf.columns = ['ChrSite', 'SiteStart', 'SiteEnd', 'Counts', 'RPM', 'ChrEnh', 'EnhStart', 'EnhEnd', 'EnhName' , 'Ens', 'H3Change']
		annoDf = annoDf.loc[annoDf.H3Change > -1, :]
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
		annoDf = pd.read_csv(inputFile, delimiter = '\t', header = None, low_memory = False)
		annoDf.columns = ['ChrSite', 'SiteStart', 'SiteEnd', 'Counts', 'RPM', 'ChrEnh', 'EnhStart', 'EnhEnd', 'EnhName' , 'Ens', 'H3Change']

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

fileName = "allSites_noDups_H3K4ME1_intersect.bed"
lsdFC = DiffCounts(fileName, allSamples, sampleName)

lsdFC = lsdFC.transpose()

fname  = "Results/FoldChanges/FC_allSamples_H3K4ME1_increased_counts.tsv"
lsdFC.to_csv(fname, sep = '\t', index = True)




