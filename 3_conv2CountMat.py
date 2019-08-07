import pandas as pd
import re
import os
import sys, getopt

def createFolder(directory):
    ''' Create directory 
        Input: name of directory 
    '''
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print('Error: Creating directory. ' + directory)

def main(argv):
    inputFile = ''
    try:
        opts, args = getopt.getopt(argv, "hi:o:", ["ifile=", "ofile="])
    except getopt.GetoptError:
        print("3_conv2CountMat.py -i <inputfile>")
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print("3_conv2CountMat.py -i <inputfile>")
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputFile = arg
        elif opt in ("-o", "--ofile"):
            outputFile = arg

    # inputFile = 'UD_MAPQ10_coord_sorted.tsv'
    udReads = pd.read_csv(inputFile, delimiter = '\t', header = None, names = ["Chr", "Start", "Cigar", "Sequence"])
    udReads.loc[:, "Start"] = udReads.loc[:, "Start"] - 1
    numrows = udReads.shape[0]

    # patterns = ["CCGG", "CCAGG", "CCTGG"]
    # methSites = [1, 1, 1]
    patterns  = ["CCTGG"]
    methSites = [1]

    chrseq = list(range(1, 20, 1))
    chrSeq = [format(x, '01d') for x in chrseq]
    chrSeq.extend(('X', 'Y', 'MT'))

    cigarList = ['I', 'D', 'S', 'H']
    cigarSet  = set(cigarList)

    # Test 
    # test = pd.DataFrame({'Chr': [1, 1], 'Start': [10, 20], 'Sequence': ['ATCCGGTG', 'ATGCGECCGGAT']})

    createFolder(outputFile)

    # Loop through all patterns of interest
    for patNum in range(len(patterns)):
        print("Finding sites: " + str(patterns[patNum]))
        pat = patterns[patNum]
        prog = re.compile(pat)
        outputFileSite = str(outputFile) + '/' + str(pat) + '_sitesLoc.tsv'
        outputFileSiteRaw =  str(outputFile) + '/' + str(pat) + '_sitesLocRaw.tsv'

        with open(outputFileSiteRaw, 'w') as f: 
            f.write('\t'.join(['Chr', 'Row', 'Start', 'Cigar', 'Loc' + '\n']))

        dictPos = {}
        # Loop through every reads in the bamfile
        for row in range(numrows):
            for match in prog.finditer(udReads.loc[row, 'Sequence']):
                pos = udReads.loc[row, 'Start'] + match.start() + methSites[patNum]
                chrm = udReads.loc[row, 'Chr']
                cgar = list(udReads.loc[row, 'Cigar'])
                cgar = set(cgar)

                # Determine if CGAR strings contain any I, D, S, H
                insect = cgar.intersection(cigarSet)

                # Save all reads with CGAR strings with either I, D, S, H
                if (bool(insect) == True):
                	with open(outputFileSiteRaw, 'a') as f:
                		f.write('\t'.join([str(chrm), str(row), str(udReads.loc[row, 'Start']), str(udReads.loc[row, 'Cigar']), 
                        str(pos) + '\n']))

                # Add up the counts for each sites
                # If chrm and sites exist
                if chrm in dictPos:
                    if pos in dictPos[chrm]:
                        dictPos[chrm][pos] += 1
                    else:
                        dictPos[chrm][pos] = 1
                
                # If chrm and sites don't exist       
                else:
                    dictPos[chrm] = {}
                    dictPos[chrm][pos] = 1

        # Save only counts and sites 
        with open(outputFileSite, 'w') as f:
            for chrmKey in dictPos.keys():
                for posKey in dictPos[chrmKey].keys():
                    f.write('\t'.join([str(chrmKey), str(posKey), str(dictPos[chrmKey][posKey])]) + '\n')

if __name__ == "__main__":
    main(sys.argv[1:])
