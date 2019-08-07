import sys, getopt
import re
import os
from Bio import SeqIO

def createFolder(directory):
	try: 
		if not os.path.exists(directory):
			os.makedirs(directory)
	except OSError:
		print('Error: Creating directory. ' + directory)

def main(argv):
	inputFile = ''
	try: 
		opts, args = getopt.getopt(argv, "hi:", "ifile=")
	except getopt.GetoptError:
		print("genomeCountSites.py -i <inputfile>")
		sys.exit(2)
	for opt, arg in opts: 
		if opt == '-h':
			print("genomeCountSites.py -i <inputfile>")
			sys.exit()
		elif opt in ("-i", "--ifile"):
			inputFile = arg

	# inputFile = 'genome_ref.fasta'

	createFolder('./genomeCatalogue/')

	seqDict = {}
	fastaSeq = SeqIO.parse(open(inputFile), 'fasta')

	for fasta in fastaSeq:
		name, sequence = fasta.id, str(fasta.seq)
		seqDict[name] = sequence
		print(name)
		print(len(sequence))

	patterns  = ["CCGG", "CCAGG", "CCTGG"]
	methSites = [1, 1, 1]

	chrseq = list(range(1, 20, 1))
	chrSeq = [format(x, '01d') for x in chrseq]
	chrSeq.extend(('X', 'Y', 'MT'))

	for patNum in range(len(patterns)):
		print("Finding sites: " + str(patterns[patNum]))
		pat  = patterns[patNum]
		prog = re.compile(pat)
		outputFile = './genomeCatalogue/' + str(pat) + '_genomeRefLoc.tsv'

		with open(outputFile, 'w') as f:
			for chrm in chrSeq:
				for match in prog.finditer(seqDict[str(chrm)]):
					f.write('\t'.join([str(chrm), str(match.start()), str(match.end()), str(match.start() + methSites[patNum])]) + '\n')

if __name__ == "__main__":
	main(sys.argv[1:])


