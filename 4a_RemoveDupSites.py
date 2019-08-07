import matplotlib
matplotlib.use('Agg')
import numpy as np
import pandas as pd
import pickle as pk
import os

def createFolder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print('Error: Creating directory. ' + directory)

# *************** Collect all sites from every samples into one dictionary & save ***********************
def CombinePatterns(allSamples, siteType, sampleName):
	""" Combine all patterns into one file 

		Output:
			- Save pickled sites into SaveData/SampleSites.pkl
			- Save each samples into (sample_name)/(sample_name)_allSites.tsv
	"""

	samples = {}

	# createFolder("SaveData")

	for i in range(len(allSamples)):
		sampleSite = pd.DataFrame(columns = ["Chr", "Sites", "CgarType", "Count", "SiteType"])

	# Reading in all sites from sample
		for j in range(len(siteType)):
			inputFile = str(allSamples[i]) + "/" + str(siteType[j]) + "_finalSites.tsv"
			countSite = pd.read_csv(inputFile, delimiter = '\t', header = 0, low_memory=False)
			siteRep   = pd.DataFrame({'SiteType': np.repeat(siteType[j], len(countSite))}).astype('str')

			frames  = [countSite, siteRep]
			siteDat = pd.concat(frames, axis = 1) 
			sampleSite = pd.concat([sampleSite, siteDat]).reset_index(drop=True)

			sampleSite = sampleSite.sort_values(by = ['Chr', 'Sites'], ascending=[True, True])

			samples[allSamples[i]] = sampleSite

	f = open("SaveData/SampleSites.pkl", "wb")
	pk.dump(samples, f)
	f.close()

	# Create tsv files for the sites for annotation 
	# SampleSites = pk.load(open('SaveData/SampleSites.pkl', 'rb'))
	SampleSites = samples
	i = 0
	for key, value in SampleSites.items():
		filename = str(key) + "/" + str(sampleName[i]) + "_allSites.tsv"
		value.to_csv(filename, sep = '\t', index = False)
		i = i + 1

# ******************************* Find duplicates & remove ********************************* #
def FindDupes(dat1):
	""" Function to find duplicate sites across all patterns 

		Output:
			- All duplicated sites and the number of duplications
	"""

	seen = {}
	dupes = []
	a = dat1.loc[:, "Sites"].values
	for x in a:
		if x not in seen:
			seen[x] = 1
		else:
			if seen[x] == 1:
				dupes.append(x)
			seen[x] +=1
	return(dupes)

def ConcDupes(smp1, chrSeq):
	""" Function to delete duplicate and combine counts 
		
		Output:
			- Data frame of sites without duplications

	"""

	nodupes = pd.DataFrame(columns = ['Chr', 'Sites', 'CgarType', 'Count', 'SiteType'])
	for j in range(len(chrSeq)):
		print(chrSeq[j])
		dat1 = smp1[smp1.Chr == chrSeq[j]]
		dupes = FindDupes(dat1)

		if (len(dupes) != 0):
			newdf = pd.DataFrame(columns = ['Chr', 'Sites', 'CgarType', 'Count', 'SiteType'])
			for k in dupes:
				dfdupes = dat1.loc[dat1.Sites == k]
				dfdupes.Count.values[0] = np.sum(dfdupes.Count.values)
				dfdupes = dfdupes[:1]
				newdf   = newdf.append(dfdupes)
				dat1 = dat1[dat1.Sites != k]

			dat2 = dat1.append(newdf).sort_values(by = 'Sites', ascending = True).reset_index(drop = True)
			nodupes = nodupes.append(dat2)
		else:
			nodupes = nodupes.append(dat1)

	return(nodupes)

def RPM(smp1, totReads):
	""" 
	Function to compute Count per Million 
	"""

	tmp = smp1.Count.values/totReads*1e6
	smp1['rpm'] = np.round(tmp, 4)
	return(smp1)

def RemoveDups(allSamples, sampleName, chrSeq):
	""" 
		Function to remove duplications and include RPM

		Output:
			- (sample_directory)/(sample_name)_allSites_noDups.tsv, files without rpm 
			- (sample_directory)/(sample_name)_allSites_noDups_final.tsv, files with rpm and counts
	"""

	for i in range(len(allSamples)):
		print(allSamples[i])
		smp1 = pd.read_csv(str(allSamples[i]) + "/" + sampleName[i] + "_allSites.tsv", delimiter = '\t', low_memory = False)
		editSmp1  = ConcDupes(smp1, chrSeq)
		filename  = str(allSamples[i]) + "/" + str(sampleName[i]) + "_allSites_noDups.tsv"
		editSmp1.to_csv(filename, sep = '\t', index = False)

	## Total number of mapped reads F9 
	## F9_UD = 2221852, F9_D4 = 2446243
	## F9_D4_PG = 2012215, F9_D4_TCP = 2336059

	totMapReads = {'F9_UD': 2221852, 'F9_D4': 2446243, 'F9_D4_PG': 2012215, 'F9_D4_TCP': 2336059}

	for i in range(len(allSamples)):
		editSmp1 = pd.read_csv(str(allSamples[i]) + "/" + sampleName[i] + "_allSites_noDups.tsv", delimiter = '\t', low_memory = False)
		editCount = RPM(editSmp1, totMapReads[sampleName[i]])
		editCount.columns = ['Chr', 'Sites', 'CgarType', 'Count', 'SiteType', 'RPM']
		filename  = str(allSamples[i]) + "/" + str(sampleName[i]) + "_allSites_noDups_final.tsv"
		editCount.to_csv(filename, sep = '\t', index = False)

def main():

	allSamples = ["F9_UD_readsCatalogue", "F9_D4_readsCatalogue", "F9_D4_PG_readsCatalogue", "F9_D4_TCP_readsCatalogue"]
	siteType   = ["CCGG", "CCAGG", "CCTGG"]
	sampleName = ['F9_UD', 'F9_D4', 'F9_D4_PG', 'F9_D4_TCP']

	chrseq = list(range(1, 20, 1))
	chrSeq = [format(x, '01d') for x in chrseq]
	chrSeq.extend(('X', 'Y', 'MT'))

	CombinePatterns(allSamples, siteType, sampleName)
	RemoveDups(allSamples, sampleName, chrSeq)

if __name__ == "__main__":
    main()






