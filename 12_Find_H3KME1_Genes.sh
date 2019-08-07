# ************************ Find Genes in H3K4ME1 in each samples **************************** #
bedtools intersect -wa -wb -a F9_UD_readsCatalogue/F9_UD_allSites_noDups_final.bed -b H3K4ME1_Sub.bed  > F9_UD_readsCatalogue/F9_UD_allSites_noDups_H3K4ME1_intersect.bed
bedtools intersect -wa -wb -a F9_D4_readsCatalogue/F9_D4_allSites_noDups_final.bed -b H3K4ME1_Sub.bed  > F9_D4_readsCatalogue/F9_D4_allSites_noDups_H3K4ME1_intersect.bed
bedtools intersect -wa -wb -a F9_D4_TCP_readsCatalogue/F9_D4_TCP_allSites_noDups_final.bed -b H3K4ME1_Sub.bed  > F9_D4_TCP_readsCatalogue/F9_D4_TCP_allSites_noDups_H3K4ME1_intersect.bed
bedtools intersect -wa -wb -a F9_D4_PG_readsCatalogue/F9_D4_PG_allSites_noDups_final.bed -b H3K4ME1_Sub.bed > F9_D4_PG_readsCatalogue/F9_D4_PG_allSites_noDups_H3K4ME1_intersect.bed
