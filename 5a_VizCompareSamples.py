import matplotlib
matplotlib.use('Agg')
import numpy as np
import pandas as pd
import pickle as pk
import matplotlib.pyplot as plt
import os

from matplotlib_venn import venn2

# ************************* Load comparisons and Venn Diagram ******************************** #
def createFolder(directory):
    ''' Create directory 
        Input: name of directory 
    '''
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print('Error: Creating directory. ' + directory)

def CreateVennPairs(allSamples, compareSite):
    """
        Create Venn diagram of paired comparisons of sites

        Output: 
            - creates folder named Figures
            - store indiviudal venn diagram in Figures 
    """
    # createFolder('Figures')

    dfCompareSite = {}

    for keys, values in compareSite.items():
        dfCompareSite[keys] = pd.DataFrame(values['Total'], index=[0])
        print(compareSite[keys])

    for keys, v in dfCompareSite.items():
        print(keys)
        print(v)

        vennLab = list(v.columns)

        vennLab[0] = vennLab[0].replace("_diff", "")
        vennLab[1] = vennLab[1].replace("_diff", "")

        plt.figure(figsize=(6,6))
        venn2(subsets = [int(v.iloc[0, 0]), int(v.iloc[0, 1]), int(v.iloc[0, 2])], set_labels = (vennLab[0], vennLab[1]))
        plt.savefig("Figures/" + str(keys) + ".pdf", bbox_inches = 'tight')

if __name__ == "__main__":

    allSamples  = ["F9_UD_readsCatalogue", "F9_D4_readsCatalogue", "F9_D4_PG_readsCatalogue", "F9_D4_TCP_readsCatalogue"]
    compareSite = pk.load(open('SaveData/SampleSitesCompare.pkl', 'rb'))

    CreateVennPairs(allSamples, compareSite)

