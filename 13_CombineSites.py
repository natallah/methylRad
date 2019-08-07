import numpy as np
import pandas as pd
import pickle as pk

inputlsd1 = "FinalAnnotation/LSD1occupied_enhancers.mm10.use.edit.bed"
lsd1_data = pd.read_csv(inputlsd1, delimiter = '\t', header = None, low_memory = False)

lsd1_sub  = lsd1_data.loc[:, [0, 1, 2, 3]]
lsd1_sub.columns = ['Chr', 'SiteStart', 'SiteEnd', 'Enhancers']

inputfold = "Results/FoldChanges/FC_allSamples_LSD1_counts.tsv"
foldData  = pd.read_csv(inputfold, delimiter = '\t', header = 0, low_memory = False)
foldData  = foldData.rename(columns = {'Unnamed: 0': 'Enhancers'})

tmp = pd.merge(lsd1_sub, foldData, how = "inner", on = ["Enhancers"])

fname  = "Results/FoldChanges/FC_allSamples_LSD1_counts_locations.tsv"
tmp.to_csv(fname, sep = '\t', index = True)