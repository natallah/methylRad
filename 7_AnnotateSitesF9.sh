# cp Results/F9_D4_D_vs_PG_D_vs_TCP_D/F9_D4_Min_UD_only_sites.bed FinalAnnotation
# cp Results/F9_UD_vs_F9_D4/F9_D4_diff_sites.bed FinalAnnotation
# cp Results/F9_UD_vs_F9_D4_PG/F9_D4_PG_diff_sites.bed FinalAnnotation
# cp Results/F9_UD_vs_F9_D4_TCP/F9_D4_TCP_diff_sites.bed FinalAnnotation

# ********************************* Annotate Sites for each samples ********************************* #
# # ************************ Convert all sites into bed files with tab delimited ************************ #
# awk 'NR>1 {print "chr"$1, $2, $2, $4, $6}' F9_UD_readsCatalogue/F9_UD_allSites_noDups_final.tsv | awk 'BEGIN{OFS="\t"}{print $1, $2, $3 + 1, $4, $5}' > F9_UD_readsCatalogue/F9_UD_allSites_noDups_final.bed
# awk 'NR>1 {print "chr"$1, $2, $2, $4, $6}' F9_D4_readsCatalogue/F9_D4_allSites_noDups_final.tsv | awk 'BEGIN{OFS="\t"}{print $1, $2, $3 + 1, $4, $5}' > F9_D4_readsCatalogue/F9_D4_allSites_noDups_final.bed
# awk 'NR>1 {print "chr"$1, $2, $2, $4, $6}' F9_D4_TCP_readsCatalogue/F9_D4_TCP_allSites_noDups_final.tsv | awk 'BEGIN{OFS="\t"}{print $1, $2, $3 + 1, $4, $5}' > F9_D4_TCP_readsCatalogue/F9_D4_TCP_allSites_noDups_final.bed
# awk 'NR>1 {print "chr"$1, $2, $2, $4, $6}' F9_D4_PG_readsCatalogue/F9_D4_PG_allSites_noDups_final.tsv | awk 'BEGIN{OFS="\t"}{print $1, $2, $3 + 1, $4, $5}' > F9_D4_PG_readsCatalogue/F9_D4_PG_allSites_noDups_final.bed

# Find overlapping features, print features from first file AND second file that overlap
bedtools intersect -wa -wb -a F9_UD_readsCatalogue/F9_UD_allSites_noDups_final.bed -b FinalAnnotation/LSD1occupied_enhancers.mm10.use.edit.bed  > F9_UD_readsCatalogue/F9_UD_allSites_noDups_LSD1_intersect.bed
bedtools intersect -wa -wb -a F9_D4_readsCatalogue/F9_D4_allSites_noDups_final.bed -b FinalAnnotation/LSD1occupied_enhancers.mm10.use.edit.bed  > F9_D4_readsCatalogue/F9_D4_allSites_noDups_LSD1_intersect.bed
bedtools intersect -wa -wb -a F9_D4_TCP_readsCatalogue/F9_D4_TCP_allSites_noDups_final.bed -b FinalAnnotation/LSD1occupied_enhancers.mm10.use.edit.bed  > F9_D4_TCP_readsCatalogue/F9_D4_TCP_allSites_noDups_LSD1_intersect.bed
bedtools intersect -wa -wb -a F9_D4_PG_readsCatalogue/F9_D4_PG_allSites_noDups_final.bed -b FinalAnnotation/LSD1occupied_enhancers.mm10.use.edit.bed  > F9_D4_PG_readsCatalogue/F9_D4_PG_allSites_noDups_LSD1_intersect.bed

# Find features that do NOT intersect with LSD1 occupied enhancers (these will be used in next steps):
bedtools intersect -a F9_UD_readsCatalogue/F9_UD_allSites_noDups_final.bed  -b FinalAnnotation/LSD1occupied_enhancers.mm10.use.edit.bed -v  > F9_UD_readsCatalogue/F9_UD_allSites_noDups_noLSD1_intersect.bed
bedtools intersect -a F9_D4_readsCatalogue/F9_D4_allSites_noDups_final.bed -b FinalAnnotation/LSD1occupied_enhancers.mm10.use.edit.bed -v  > F9_D4_readsCatalogue/F9_D4_allSites_noDups_noLSD1_intersect.bed
bedtools intersect -a F9_D4_TCP_readsCatalogue/F9_D4_TCP_allSites_noDups_final.bed -b FinalAnnotation/LSD1occupied_enhancers.mm10.use.edit.bed -v  > F9_D4_TCP_readsCatalogue/F9_D4_TCP_allSites_noDups_noLSD1_intersect.bed
bedtools intersect -a F9_D4_PG_readsCatalogue/F9_D4_PG_allSites_noDups_final.bed -b FinalAnnotation/LSD1occupied_enhancers.mm10.use.edit.bed -v  > F9_D4_PG_readsCatalogue/F9_D4_PG_allSites_noDups_noLSD1_intersect.bed

# Annotate with second enhancer dataset
bedtools intersect -wa -wb -a F9_UD_readsCatalogue/F9_UD_allSites_noDups_noLSD1_intersect.bed -b FinalAnnotation/ESC_J1.enhancers.use.edit.bed > F9_UD_readsCatalogue/F9_UD_allSites_noDups_ESC_J1_enhancers_intersect.bed
bedtools intersect -wa -wb -a F9_D4_readsCatalogue/F9_D4_allSites_noDups_noLSD1_intersect.bed -b FinalAnnotation/ESC_J1.enhancers.use.edit.bed > F9_D4_readsCatalogue/F9_D4_allSites_noDups_ESC_J1_enhancers_intersect.bed
bedtools intersect -wa -wb -a F9_D4_TCP_readsCatalogue/F9_D4_TCP_allSites_noDups_noLSD1_intersect.bed -b FinalAnnotation/ESC_J1.enhancers.use.edit.bed > F9_D4_TCP_readsCatalogue/F9_D4_TCP_allSites_noDups_ESC_J1_enhancers_intersect.bed
bedtools intersect -wa -wb -a F9_D4_PG_readsCatalogue/F9_D4_PG_allSites_noDups_noLSD1_intersect.bed -b FinalAnnotation/ESC_J1.enhancers.use.edit.bed > F9_D4_PG_readsCatalogue/F9_D4_PG_allSites_noDups_ESC_J1_enhancers_intersect.bed

# Find features that do NOT intersect with ESCJ1 occupied enhancers (these will be used in next steps):
bedtools intersect -a F9_UD_readsCatalogue/F9_UD_allSites_noDups_noLSD1_intersect.bed -b FinalAnnotation/ESC_J1.enhancers.use.edit.bed -v  > F9_UD_readsCatalogue/F9_UD_allSites_noDups_noenhancer_intersect.bed
bedtools intersect -a F9_D4_readsCatalogue/F9_D4_allSites_noDups_noLSD1_intersect.bed -b FinalAnnotation/ESC_J1.enhancers.use.edit.bed -v  > F9_D4_readsCatalogue/F9_D4_allSites_noDups_noenhancer_intersect.bed
bedtools intersect -a F9_D4_TCP_readsCatalogue/F9_D4_TCP_allSites_noDups_noLSD1_intersect.bed -b FinalAnnotation/ESC_J1.enhancers.use.edit.bed -v  > F9_D4_TCP_readsCatalogue/F9_D4_TCP_allSites_noDups_noenhancer_intersect.bed
bedtools intersect -a F9_D4_PG_readsCatalogue/F9_D4_PG_allSites_noDups_noLSD1_intersect.bed -b FinalAnnotation/ESC_J1.enhancers.use.edit.bed -v  > F9_D4_PG_readsCatalogue/F9_D4_PG_allSites_noDups_noenhancer_intersect.bed

# ********************************* Annotate Sites for each samples on Promoters ********************************* #
bedtools intersect -wa -wb -a F9_UD_readsCatalogue/F9_UD_allSites_noDups_noenhancer_intersect.bed -b promoters.bed > F9_UD_readsCatalogue/F9_UD_allSites_noDups_noenhancer_promoters_intersect.bed
bedtools intersect -wa -wb -a F9_D4_readsCatalogue/F9_D4_allSites_noDups_noenhancer_intersect.bed -b promoters.bed > F9_D4_readsCatalogue/F9_D4_allSites_noDups_noenhancer_promoters_intersect.bed
bedtools intersect -wa -wb -a F9_D4_TCP_readsCatalogue/F9_D4_TCP_allSites_noDups_noenhancer_intersect.bed -b promoters.bed > F9_D4_TCP_readsCatalogue/F9_D4_TCP_allSites_noDups_noenhancer_promoters_intersect.bed
bedtools intersect -wa -wb -a F9_D4_PG_readsCatalogue/F9_D4_PG_allSites_noDups_noenhancer_intersect.bed -b promoters.bed > F9_D4_PG_readsCatalogue/F9_D4_PG_allSites_noDups_noenhancer_promoters_intersect.bed

# Total overlap with LSD1
wc -l F9_UD_readsCatalogue/F9_UD_allSites_noDups_LSD1_intersect.bed
wc -l F9_D4_readsCatalogue/F9_D4_allSites_noDups_LSD1_intersect.bed
wc -l F9_D4_PG_readsCatalogue/F9_D4_PG_allSites_noDups_LSD1_intersect.bed
wc -l F9_D4_TCP_readsCatalogue/F9_D4_TCP_allSites_noDups_LSD1_intersect.bed


# Total overlap with ESC_J1
wc -l F9_UD_readsCatalogue/F9_UD_allSites_noDups_ESC_J1_enhancers_intersect.bed
wc -l F9_D4_readsCatalogue/F9_D4_allSites_noDups_ESC_J1_enhancers_intersect.bed
wc -l F9_D4_PG_readsCatalogue/F9_D4_PG_allSites_noDups_ESC_J1_enhancers_intersect.bed
wc -l F9_D4_TCP_readsCatalogue/F9_D4_TCP_allSites_noDups_ESC_J1_enhancers_intersect.bed

# Total all sites 
wc -l F9_UD_readsCatalogue/F9_UD_allSites_noDups_final.bed
wc -l F9_D4_readsCatalogue/F9_D4_allSites_noDups_final.bed
wc -l F9_D4_PG_readsCatalogue/F9_D4_PG_allSites_noDups_final.bed
wc -l F9_D4_TCP_readsCatalogue/F9_D4_TCP_allSites_noDups_final.bed

# ********************************* Annotate D4 - UD, D4_PG - UD, and D4_TCP - UD ********************************* #
# Intersect sites with count data
bedtools intersect -wa -a F9_D4_readsCatalogue/F9_D4_allSites_noDups_final.bed -b Results/F9_UD_vs_F9_D4/F9_D4_diff_sites.bed > Results/F9_UD_vs_F9_D4/F9_D4_diff_sites_counts.bed
bedtools intersect -wa -a F9_D4_PG_readsCatalogue/F9_D4_PG_allSites_noDups_final.bed -b Results/F9_UD_vs_F9_D4_PG/F9_D4_PG_diff_sites.bed > Results/F9_UD_vs_F9_D4_PG/F9_D4_PG_diff_sites_counts.bed
bedtools intersect -wa -a F9_D4_TCP_readsCatalogue/F9_D4_TCP_allSites_noDups_final.bed -b Results/F9_UD_vs_F9_D4_TCP/F9_D4_TCP_diff_sites.bed > Results/F9_UD_vs_F9_D4_TCP/F9_D4_TCP_diff_sites_counts.bed

# Annotate with LSD1 enhancers
bedtools intersect -wa -wb -a Results/F9_UD_vs_F9_D4/F9_D4_diff_sites_counts.bed -b FinalAnnotation/LSD1occupied_enhancers.mm10.use.edit.bed  > Results/F9_UD_vs_F9_D4/F9_D4_diff_sites_LSD1_intersect.bed
bedtools intersect -wa -wb -a Results/F9_UD_vs_F9_D4_PG/F9_D4_PG_diff_sites_counts.bed -b FinalAnnotation/LSD1occupied_enhancers.mm10.use.edit.bed  > Results/F9_UD_vs_F9_D4_PG/F9_D4_PG_diff_sites_LSD1_intersect.bed 
bedtools intersect -wa -wb -a Results/F9_UD_vs_F9_D4_TCP/F9_D4_TCP_diff_sites_counts.bed -b FinalAnnotation/LSD1occupied_enhancers.mm10.use.edit.bed  > Results/F9_UD_vs_F9_D4_TCP/F9_D4_TCP_diff_sites_LSD1_intersect.bed

# Minus annotation with LSD1 enhancers
bedtools intersect -a Results/F9_UD_vs_F9_D4/F9_D4_diff_sites_counts.bed -b FinalAnnotation/LSD1occupied_enhancers.mm10.use.edit.bed -v > Results/F9_UD_vs_F9_D4/F9_D4_diff_sites_noLSD1_intersect.bed
bedtools intersect -a Results/F9_UD_vs_F9_D4_PG/F9_D4_PG_diff_sites_counts.bed -b FinalAnnotation/LSD1occupied_enhancers.mm10.use.edit.bed  -v > Results/F9_UD_vs_F9_D4_PG/F9_D4_PG_diff_sites_noLSD1_intersect.bed 
bedtools intersect -a Results/F9_UD_vs_F9_D4_TCP/F9_D4_TCP_diff_sites_counts.bed -b FinalAnnotation/LSD1occupied_enhancers.mm10.use.edit.bed  -v > Results/F9_UD_vs_F9_D4_TCP/F9_D4_TCP_diff_sites_noLSD1_intersect.bed

# Annotate with ESC enhancer dataset
bedtools intersect -wa -wb -a Results/F9_UD_vs_F9_D4/F9_D4_diff_sites_noLSD1_intersect.bed -b FinalAnnotation/ESC_J1.enhancers.use.edit.bed > Results/F9_UD_vs_F9_D4/F9_D4_diff_sites_ESC_J1_intersect.bed
bedtools intersect -wa -wb -a Results/F9_UD_vs_F9_D4_PG/F9_D4_PG_diff_sites_noLSD1_intersect.bed  -b FinalAnnotation/ESC_J1.enhancers.use.edit.bed > Results/F9_UD_vs_F9_D4_PG/F9_D4_PG_diff_sites_ESC_J1_intersect.bed
bedtools intersect -wa -wb -a Results/F9_UD_vs_F9_D4_TCP/F9_D4_TCP_diff_sites_noLSD1_intersect.bed -b FinalAnnotation/ESC_J1.enhancers.use.edit.bed > Results/F9_UD_vs_F9_D4_TCP/F9_D4_TCP_diff_sites_ESC_J1_intersect.bed

# Minus annotate with ESC enhancer dataset
bedtools intersect -a Results/F9_UD_vs_F9_D4/F9_D4_diff_sites_noLSD1_intersect.bed -b FinalAnnotation/ESC_J1.enhancers.use.edit.bed -v > Results/F9_UD_vs_F9_D4/F9_D4_diff_sites_noenhancers_intersect.bed
bedtools intersect -a Results/F9_UD_vs_F9_D4_PG/F9_D4_PG_diff_sites_noLSD1_intersect.bed  -b FinalAnnotation/ESC_J1.enhancers.use.edit.bed -v > Results/F9_UD_vs_F9_D4_PG/F9_D4_PG_diff_sites_noenhancers_intersect.bed
bedtools intersect -a Results/F9_UD_vs_F9_D4_TCP/F9_D4_TCP_diff_sites_noLSD1_intersect.bed -b FinalAnnotation/ESC_J1.enhancers.use.edit.bed -v > Results/F9_UD_vs_F9_D4_TCP/F9_D4_TCP_diff_sites_noenhancers_intersect.bed

# Total overlap with LSD1
wc -l Results/F9_UD_vs_F9_D4/F9_D4_diff_sites_LSD1_intersect.bed
wc -l Results/F9_UD_vs_F9_D4_PG/F9_D4_PG_diff_sites_LSD1_intersect.bed 
wc -l Results/F9_UD_vs_F9_D4_TCP/F9_D4_TCP_diff_sites_LSD1_intersect.bed

# Total overlap with ESC_J1
wc -l Results/F9_UD_vs_F9_D4/F9_D4_diff_sites_ESC_J1_intersect.bed
wc -l Results/F9_UD_vs_F9_D4_PG/F9_D4_PG_diff_sites_ESC_J1_intersect.bed 
wc -l Results/F9_UD_vs_F9_D4_TCP/F9_D4_TCP_diff_sites_ESC_J1_intersect.bed

# Total no overlap
wc -l Results/F9_UD_vs_F9_D4/F9_D4_diff_sites_noenhancers_intersect.bed
wc -l Results/F9_UD_vs_F9_D4_PG/F9_D4_PG_diff_sites_noenhancers_intersect.bed
wc -l Results/F9_UD_vs_F9_D4_TCP/F9_D4_TCP_diff_sites_noenhancers_intersect.bed

