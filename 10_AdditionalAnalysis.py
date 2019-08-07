import numpy as np
import pandas as pd
import pickle as pk
import os

# ************* Find Common Genes for RA4 only sites in promoters and enhancers *************** #
inputFile = "Results/RA4_D_vs_PG_D_vs_TCP_D/RA4_UD_only_sites_annoESC.bed"
escDat = pd.read_csv(inputFile, delimiter='\t', header=None, low_memory=False, usecols=[0, 1, 2, 3, 7, 8], 
	names=['Chr', 'Start', 'End', 'Rpm', 'Ens_id', 'Gene_name'])

inputFile = "Results/RA4_D_vs_PG_D_vs_TCP_D/RA4_UD_only_sites_annoPromoters.bed"
promDat = pd.read_csv(inputFile, delimiter='\t', header=None, low_memory=False, usecols=[0, 1, 2, 3, 10], 
	names=['Chr', 'Start', 'End', 'Rpm', 'Ens_id'])

inputFile = 
ensESC = set(escDat.Ens_id)
ensProm = set(promDat.Ens_id)

# Find overlap
intersectEns = list(ensESC & ensProm)
intersectEns_df = pd.DataFrame(intersectEns)

# Subset data frame
escDatInt  = escDat.loc[escDat.Ens_id.isin(intersectEns),:]
promDatInt = promDat.loc[promDat.Ens_id.isin(intersectEns),:] 

# Write to tsv files
rsltName = "Results/RA4_D_vs_PG_D_vs_TCP_D/RA4_UD_only_sites_geneList_in_ESC_and_promoters.bed"
intersectEns_df.to_csv(rsltName, sep = '\t', index = False, header = False)

# ***************** Combine sites for Analysis between RA4_PG-UD and RA4_TCP-UD ******************* #
# For sites in RA4_PG-UD only
inputFile = "Results/RA4_D_vs_PG_D_vs_TCP_D/RA4_UD_vs_PG_intersection_sites.bed"
dat1 = pd.read_csv(inputFile, delimiter='\t', header=None, low_memory=False, names=['Chr', 'Start', 'End'])

inputFile = "Results/RA4_D_vs_PG_D_vs_TCP_D/RA4_PG_UD_only_sites.bed"
dat2 = pd.read_csv(inputFile, delimiter='\t', header=None, low_memory=False, names=['Chr', 'Start', 'End'])

combSite = pd.concat([dat1, dat2]).reset_index(drop=True)

combSite = combSite.sort_values(by = ['Chr', 'Start'], ascending=[True, True]).reset_index(drop=True)

combName = 'Results/RA4_D_vs_PG_D_vs_TCP_D/RA4_PG_Min_UD_vs_RA4_TCP_Min_UD/RA4_PG_Min_UD_only_sites.bed'	
combSite.to_csv(combName, sep = '\t', index = False, header = False)

# For sites in RA4_TCP-UD only
inputFile = "Results/RA4_D_vs_PG_D_vs_TCP_D/RA4_UD_vs_TCP_intersection_sites.bed"
dat1 = pd.read_csv(inputFile, delimiter='\t', header=None, low_memory=False, names=['Chr', 'Start', 'End'])

inputFile = "Results/RA4_D_vs_PG_D_vs_TCP_D/RA4_TCP_UD_only_sites.bed"
dat2 = pd.read_csv(inputFile, delimiter='\t', header=None, low_memory=False, names=['Chr', 'Start', 'End'])

combSite = pd.concat([dat1, dat2]).reset_index(drop=True)

combSite = combSite.sort_values(by = ['Chr', 'Start'], ascending=[True, True]).reset_index(drop=True)

combName = 'Results/RA4_D_vs_PG_D_vs_TCP_D/RA4_PG_Min_UD_vs_RA4_TCP_Min_UD/RA4_TCP_Min_UD_only_sites.bed'	
combSite.to_csv(combName, sep = '\t', index = False, header = False)

# For sites in the intersection of RA4_TCP-UD and RA4_Pg-UD
inputFile = "Results/RA4_D_vs_PG_D_vs_TCP_D/RA4_vs_PG_vs_TCP_intersection_sites.bed"
dat1 = pd.read_csv(inputFile, delimiter='\t', header=None, low_memory=False, names=['Chr', 'Start', 'End'])

inputFile = "Results/RA4_D_vs_PG_D_vs_TCP_D/RA4_PG_vs_TCP_intersection_sites.bed"
dat2 = pd.read_csv(inputFile, delimiter='\t', header=None, low_memory=False, names=['Chr', 'Start', 'End'])

combSite = pd.concat([dat1, dat2]).reset_index(drop=True)

combSite = combSite.sort_values(by = ['Chr', 'Start'], ascending=[True, True]).reset_index(drop=True)

combName = 'Results/RA4_D_vs_PG_D_vs_TCP_D/RA4_PG_Min_UD_vs_RA4_TCP_Min_UD/RA4_TCP_minus_UD_intersect_RA4_PG_minus_UD_sites.bed'	
combSite.to_csv(combName, sep = '\t', index = False, header = False)

# *********************** Count the number of CCGG, CCAGG, CCTGG in each samples ************************** #
# UD
inputFile = "UD_readsCatalogue/UD_allSites_noDups_final.tsv"
dat = pd.read_csv(inputFile, delimiter='\t', header=0, low_memory=False)
dat['SiteType'].value_counts()

# RA4
inputFile = "RA4_readsCatalogue/RA4_allSites_noDups_final.tsv"
dat = pd.read_csv(inputFile, delimiter='\t', header=0, low_memory=False)
dat['SiteType'].value_counts()

# RA4_PG
inputFile = "RA4_PG_readsCatalogue/RA4_PG_allSites_noDups_final.tsv"
dat = pd.read_csv(inputFile, delimiter='\t', header=0, low_memory=False)
dat['SiteType'].value_counts()

# RA4_TCP
inputFile = "RA4_TCP_readsCatalogue/RA4_TCP_allSites_noDups_final.tsv"
dat = pd.read_csv(inputFile, delimiter='\t', header=0, low_memory=False)
dat['SiteType'].value_counts()













