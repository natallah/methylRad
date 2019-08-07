import numpy as np
import pandas as pd
import pickle as pk

# ********************* Read in Annotation File ******************** #
geneBed = []

with open("Mus_musculus.GRCm38.93.chr.bed") as f:
	for line in f:
		geneBed.append(line.strip().split())

geneBed = np.asarray(geneBed)

genePds = pd.DataFrame({'Chr': geneBed[:, 0], 'Start': geneBed[:, 1], 
	'End': geneBed[:, 2], 'GeneName': geneBed[0, 3], 'Score': geneBed[0, 4],
	'Strand': geneBed[:, 5]})

cols = ['Chr', 'Start', 'End', 'GeneName', 'Score', 'Strand']
genePds = genePds[cols]

# ********************* Read in Peak File *********************** #
chrseq = list(range(1, 20, 1))
chrSeq = [format(x, '01d') for x in chrseq]
chrSeq.extend(('X', 'Y', 'MT'))

allSamples = ["UD_readsCatalogue", "RA4_readsCatalogue", "RA4_PG_readsCatalogue", "RA4_TCP_readsCatalogue"]
siteType   = ["CCGG", "CCAGG", "CCTGG"]

sampleSites = pk.load(open("SaveData/SampleSites.pkl", "rb"))

for i in allSamples:
	geneAnno = np.zeros(len(sampleSites[i]), dtype=str)

	for j in chrSeq:
		tmp = sampleSites[i].loc[sampleSites[i]['Chr'] == j]
		loc = tmp['Sites']

		for k in range(len(loc)):
			



