import pandas as pd
# import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO
import sys, getopt
import re
import os

samplePath  = ["UD_readsCatalogue", "RA4_readsCatalogue", "RA4_TCP_readsCatalogue", "RA4_PG_readsCatalogue"]
inputGenome = ["genomeCatalogue/CCGG_genomeRefLoc.tsv", "genomeCatalogue/CCAGG_genomeRefLoc.tsv", "genomeCatalogue/CCTGG_genomeRefLoc.tsv"]

siteType = ["CCGG", "CCAGG", "CCTGG"]

for k in range(len(samplePath)):
    inputFile   = [str(samplePath[k]) + "/CCGG_sitesLoc.tsv", str(samplePath[k]) + "/CCAGG_sitesLoc.tsv", str(samplePath[k]) + "/CCTGG_sitesLoc.tsv"]

    chrseq = list(range(1, 20, 1))
    chrSeq = [format(x, '01d') for x in chrseq]
    chrSeq.extend(('X', 'Y', 'MT'))

    for j in range(len(inputFile)):
        print(inputFile[j])
        countSite = pd.read_csv(inputFile[j], delimiter = '\t', header = None, names = ["Chr", "Site", "Count"], low_memory=False)
        countSite['Count'] = countSite['Count'].astype('float')
        wholeGenome = pd.read_csv(inputGenome[j], delimiter = '\t', header = None, names = ["Chr", "Start", "End", "Site"],
                                  dtype = {'Chr': 'str', 'Start': 'int', 'End': 'int', 'Site': 'int'})
        lenGenome = len(wholeGenome)

        # Create data frame of whole genome and counts of reads found at site
        dfWholeGenomeCount = pd.DataFrame(wholeGenome, columns=['Chr','Site'])
        dfWholeGenomeCount['Count'] = np.zeros((lenGenome, 1), dtype='float')

        countSiteFil = countSite.loc[countSite['Count'] >= 5]
        print("Number of sites greater than 5", str(len(countSiteFil)))
        countDiff = {}
        countBoth = {}

        for i in chrSeq:
            # print("Chromosomes: " + str(i))
            tmp  = dfWholeGenomeCount.loc[(dfWholeGenomeCount['Chr'] == i)]
            tmp2 = tmp.loc[tmp['Site'].isin(countSiteFil['Site']), ]
            tmp3 = countSiteFil.loc[countSiteFil['Chr'] == i, ]
            setGenome = set(tmp2['Site'])
            setCounts = set(tmp3['Site'])
            inboth  = setGenome.intersection(setCounts)
            inCount = setCounts.difference(setGenome)
            countDiff[i] = len(inCount)
            countBoth[i] = len(inboth)

        print("Site type:" + str(siteType[j]))
        print("Number of Differences:" + str(sum(countDiff.values())))
        print("Number of Intersection:" + str(sum(countBoth.values())))






