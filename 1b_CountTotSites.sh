wc -l CCGG_genomeRefLoc.tsv 
wc -l CCAGG_genomeRefLoc.tsv 
wc -l CCTGG_genomeRefLoc.tsv 

wc -l CCGG_sitesLoc.tsv >> sitesLocCount.out
wc -l CCAGG_sitesLoc.tsv >> sitesLocCount.out
wc -l CCTGG_sitesLoc.tsv >> sitesLocCount.out

# Number of Sites in Reference Genome:

# CCGG  = 1,630,178
# GGCC  = 6,884,821

# CCAGG = 3,934,781
# GGACC = 1,667,862

# CCTGG = 3,932,097
# GGTCC = 1,663,995  

# Number of mapped reads in each bamfiles: 
# UD  = 50,959,931
# RA4 = 73,712,891
# RA4_TCP = 74,553,237
# RA4_PG  = 67,235,311