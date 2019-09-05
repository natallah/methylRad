import pandas as pd
import re
import numpy as np
import sys, getopt

from collections import Counter

pd.options.mode.chained_assignment = None 

def main(argv):
    inputFile = ''
    try:
        opts, args = getopt.getopt(argv, "hi:", ["ifile="])
    except getopt.GetoptError:
        print("4_compareSites.py -i <inputfile>")
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print("4_compareSites.py -i <inputfile>")
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputFile = arg

    # inputFile = "UD_readsCatalogue"
    siteTypes = ["CCGG", "CCAGG", "CCTGG"]

    for tt in range(len(siteTypes)):
        siteFile = str(inputFile) + "/" + str(siteTypes[tt]) + "_sitesLoc.tsv"
        countSite = pd.read_csv(siteFile, delimiter = '\t', header = None, names = ["Chr", "Site", "Count"],
                                low_memory = False)
        countSite['Count'] = countSite['Count'].astype('float')

        inputGenome = "genomeCatalogue/" + str(siteTypes[tt]) + "_genomeRefLoc.tsv"
        wholeGenome = pd.read_csv(inputGenome, delimiter = '\t', header = None, names = ["Chr", "Start", "End", "Site"],
                                  dtype = {'Chr': 'str', 'Start': 'int', 'End': 'int', 'Site': 'int'})
        lenGenome = len(wholeGenome)

        cigarFile = str(inputFile) + "/" + str(siteTypes[tt]) + '_sitesLocRaw.tsv'
        cigarSite = pd.read_csv(cigarFile, delimiter = '\t', header = 0, low_memory = False)

        # Create data frame of whole genome and counts of reads found at site
        dfWholeGenome = pd.DataFrame(wholeGenome, columns = ['Chr', 'Site'])

        # List of all the chromosomes
        chrseq = list(range(1, 20, 1))
        chrSeq = [format(x, '01d') for x in chrseq]
        chrSeq.extend(('X', 'Y', 'MT'))

        # Only include sites that have greater than 5 counts
        countSiteFil = countSite.loc[countSite['Count'] >= 5]

        countDiff = {}
        countBoth = {}

        # Loop through every chromosomes to find matches and mismatches
        for i in chrSeq:
            print("Chromosomes: " + str(i))
            tmp  = dfWholeGenome.loc[dfWholeGenome['Chr'] == i]
            tmp1 = countSiteFil.loc[countSiteFil['Chr'] == i]

            setGenome = set(tmp['Site'])
            setCounts = set(tmp1['Site'])

            inboth = list(setGenome.intersection(setCounts))
            inCount = list(setCounts.difference(setGenome))

            countBoth[i] = inboth
            countDiff[i] = inCount
            
            print("Number of Sites in reference genome: " + str(len(inboth)))
            print("Number of Sites not in reference genome: " + str(len(inCount)))

        # Create data frame to determine substitution, deletion, and insertion for all sites that are different
        mcgarDict = {}

        for j in chrSeq:
            inCount = np.sort(countDiff[j])
            chrLoc = cigarSite.loc[cigarSite['Chr'] == j]
            mcgar = np.zeros((1, len(inCount)), dtype=bool)
            cgar = np.zeros((1, len(inCount)), dtype=str)

            # Loop over values that are not in reference genome
            for k in range(len(inCount)):
                
                # Find values not in reference genome in the set of Cigars with D/I
                val = chrLoc.loc[chrLoc['Loc'] == inCount[k]]

                # If not in the set of Cigars with D/I, it must be a substitution because then its Cigar is a Matched
                # Else, if it is in the set of Cigars with D/I, we must correct for the location accordingly
                if (val.empty == True):
                    mcgar[0, k] = True
                    cgar[0, k] = 'S'
                else:
                    cc = val['Cigar'].values[0]
                    newLoc = val['Start'].values[0]
                    mCount = 0
                    for num, let in re.findall('(\d+)([IDM])', cc):
                        if (let == 'M'):
                            newLoc = newLoc + int(num)
                            mCount = mCount + 1
                            if ((inCount[k] <= newLoc) & (mCount == 1)):
                                cgar[0, k] = 'S'
                                break
                        elif (let == 'I'):
                            if ((inCount[k] > newLoc)):
                                inCount[k] = inCount[k] - int(num)
                                cgar[0, k] = 'I'
                            newLoc = newLoc + int(num)
                            mCount = mCount + 1
                        elif (let == 'D'):
                            if ((inCount[k] > newLoc)):
                                inCount[k] = inCount[k] + int(num)
                                cgar[0, k] = 'D'
                            mCount = mCount + 1

            mcgardf = {}
            mcgardf["Chr"] = j
            mcgardf["Sites"] = inCount
            mcgardf["CgarType"] = cgar[0]
            mcgarDict[j] = pd.DataFrame(mcgardf)

            print(Counter(cgar[0]))

            tmp = countSiteFil.loc[(countSiteFil['Chr'] == j)]
            inVals = tmp.loc[tmp['Site'].isin(np.sort(countDiff[j]))]

            mcgarDict[j]['Count'] = inVals['Count'].values

        # Create data frame to determine substitution, deletion, and insertion for all sites that are the same
        mcgarDict2 = {}
        siteDict = {}
        df = pd.DataFrame(columns = ['Chr', 'Sites', 'CgarType', 'Count'])

        for j in chrSeq:
            print(j)
            bCount  = np.sort(countBoth[j])
            mcgardf = {}
            mcgardf["Chr"] = j
            mcgardf["Sites"] = bCount
            mcgardf["CgarType"] = np.repeat('M', len(bCount))
            mcgarDict2[j] = pd.DataFrame(mcgardf)

            # Include counts in data frame
            inVals = countSiteFil.loc[(countSiteFil['Chr'] == j) & (countSiteFil['Site'].isin(bCount))]
            newinVals = inVals
            newinVals['Total'] = newinVals.groupby(['Chr', 'Site'])['Count'].transform('sum')
            newinVals2 = newinVals.drop_duplicates(subset=['Chr', 'Site'])
            newinVals2['Count'] = newinVals2['Total']
            inVals = newinVals2.drop(['Total'], axis=1)
            mcgarDict2[j]['Count'] = inVals['Count'].values

            # Concatenate dataframe with sites marked as 'S/I/D' and 'M'
            result = pd.concat([mcgarDict[j], mcgarDict2[j]])
            siteDict[j] = result.sort_values(by=['Sites']).reset_index(drop=True)

            # Find duplicates in each chromosomes
            dupes = siteDict[j].duplicated(subset=['Chr', 'Sites'], keep = False)
            indupes = np.where(dupes == True)

            # Fix duplicates
            if (len(indupes[0]) > 0):
                dfDupes   = siteDict[j].loc[indupes[0], :]
                dupeSites = dfDupes['Sites'].values
                uniqDupes = np.unique(dupeSites)

                conSites  = []
                conCounts = []
                conCgar   = []

                for kk in range(len(uniqDupes)):
                    tmp = dfDupes.loc[dfDupes['Sites'] == uniqDupes[kk]]
                    combCount = np.sum(tmp['Count'].values)
                    conSites.append(uniqDupes[kk])
                    conCounts.append(combCount)

                    # Set Cgar as 'M' for every duplicates for now
                    conCgar.append('M')

                noDupesDf = pd.DataFrame({'Chr': np.repeat(str(j), len(uniqDupes)), 'Sites': conSites, 
                'CgarType': conCgar, 'Count': conCounts})

                # Remove duplicates and replace with combined values 
                siteDict[j] = siteDict[j].drop(indupes[0])
                siteDict[j] = pd.concat([siteDict[j], noDupesDf])

            siteDict[j] = siteDict[j].sort_values(by=['Sites']).reset_index(drop=True)

            frames = [df, siteDict[j]]
            df = pd.concat(frames)
            df = df[['Chr', 'Sites', 'CgarType', 'Count']]
            df['Sites'] = pd.to_numeric(df['Sites'], downcast = 'integer')

        df = df.reset_index(drop=True)
        
        filename = str(inputFile) + "/" + str(siteTypes[tt]) + "_finalSites.tsv"
        df.to_csv(filename, sep = "\t", header = True, index = False)

if __name__ == "__main__":
    main(sys.argv[1:])





