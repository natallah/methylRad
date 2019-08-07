# # cp ../Results/RA4_D_vs_PG_D_vs_TCP_D/RA4_UD_only_sites.bed ./
# # cp ../Results/UD_vs_RA4/RA4_diff_sites.bed ./
# # cp ../Results/UD_vs_RA4_PG/RA4_PG_diff_sites.bed ./
# # cp ../Results/UD_vs_RA4_TCP/RA4_TCP_diff_sites.bed ./

# Add 1KB to each LSD1 sites 
awk 'BEGIN{OFS="\t"}{print $1, $2 - 1000, $3 + 1000, $4, $5, $6, $7, $8, $9, $10, $11}' LSD1occupied_enhancers.mm10.use.bed > LSD1occupied_enhancers.mm10.use.edit.bed

# Add 1KB to each ESCJ1 sites 
awk 'BEGIN{OFS="\t"}{print $1, $2 - 1000, $3 + 1000, $4, $5}' ESC_J1.enhancers.use.bed > ESC_J1.enhancers.use.edit.bed

# # Find overlapping features, print features from first file AND second file that overlap
# bedtools intersect -wa -wb -a FinalAnnotation/RA4_UD_only_sites.bed -b LSD1occupied_enhancers.mm10.use.bed  > FinalAnnotation/RA4_UD_only_sites_LSD1_intersect.bed
# bedtools intersect -wa -wb -a FinalAnnotation/RA4_diff_sites.bed -b LSD1occupied_enhancers.mm10.use.bed  > FinalAnnotation/RA4_diff_sites_LSD1_intersect.bed
# bedtools intersect -wa -wb -a FinalAnnotation/RA4_PG_diff_sites.bed -b LSD1occupied_enhancers.mm10.use.bed  > FinalAnnotation/RA4_PG_diff_sites_LSD1_intersect.bed
# bedtools intersect -wa -wb -a FinalAnnotation/RA4_TCP_diff_sites.bed -b LSD1occupied_enhancers.mm10.use.bed  > FinalAnnotation/RA4_TCP_diff_sites_LSD1_intersect.bed

# Find overlapping features, print features from first file AND second file that overlap
bedtools intersect -wa -wb -a FinalAnnotation/RA4_UD_only_sites.bed -b FinalAnnotation/LSD1occupied_enhancers.mm10.use.edit.bed  > FinalAnnotation/RA4_UD_only_sites_LSD1_intersect.bed
bedtools intersect -wa -wb -a FinalAnnotation/RA4_diff_sites.bed -b FinalAnnotation/LSD1occupied_enhancers.mm10.use.edit.bed   > FinalAnnotation/RA4_diff_sites_LSD1_intersect.bed
bedtools intersect -wa -wb -a FinalAnnotation/RA4_PG_diff_sites.bed -b FinalAnnotation/LSD1occupied_enhancers.mm10.use.edit.bed  > FinalAnnotation/RA4_PG_diff_sites_LSD1_intersect.bed
bedtools intersect -wa -wb -a FinalAnnotation/RA4_TCP_diff_sites.bed -b FinalAnnotation/LSD1occupied_enhancers.mm10.use.edit.bed  > FinalAnnotation/RA4_TCP_diff_sites_LSD1_intersect.bed

# # Find features that do NOT intersect with LSD1 occupied enhancers (these will be used in next steps):
# bedtools intersect -a FinalAnnotation/RA4_UD_only_sites.bed -b LSD1occupied_enhancers.mm10.use.bed -v  > FinalAnnotation/RA4_UD_only_sites_noLSD1_intersect.bed
# bedtools intersect -a FinalAnnotation/RA4_diff_sites.bed -b LSD1occupied_enhancers.mm10.use.bed -v  > FinalAnnotation/RA4_diff_sites_noLSD1_intersect.bed
# bedtools intersect -a FinalAnnotation/RA4_PG_diff_sites.bed -b LSD1occupied_enhancers.mm10.use.bed -v  > FinalAnnotation/RA4_PG_diff_sites_noLSD1_intersect.bed
# bedtools intersect -a FinalAnnotation/RA4_TCP_diff_sites.bed -b LSD1occupied_enhancers.mm10.use.bed -v  > FinalAnnotation/RA4_TCP_diff_sites_noLSD1_intersect.bed

# Find features that do NOT intersect with LSD1 occupied enhancers (these will be used in next steps):
bedtools intersect -a FinalAnnotation/RA4_UD_only_sites.bed -b FinalAnnotation/LSD1occupied_enhancers.mm10.use.edit.bed -v  > FinalAnnotation/RA4_UD_only_sites_noLSD1_intersect.bed
bedtools intersect -a FinalAnnotation/RA4_diff_sites.bed -b FinalAnnotation/LSD1occupied_enhancers.mm10.use.edit.bed -v  > FinalAnnotation/RA4_diff_sites_noLSD1_intersect.bed
bedtools intersect -a FinalAnnotation/RA4_PG_diff_sites.bed -b FinalAnnotation/LSD1occupied_enhancers.mm10.use.edit.bed -v  > FinalAnnotation/RA4_PG_diff_sites_noLSD1_intersect.bed
bedtools intersect -a FinalAnnotation/RA4_TCP_diff_sites.bed -b FinalAnnotation/LSD1occupied_enhancers.mm10.use.edit.bed -v  > FinalAnnotation/RA4_TCP_diff_sites_noLSD1_intersect.bed

# # Annotate with second enhancer dataset
# bedtools intersect -wa -wb -a FinalAnnotation/RA4_UD_only_sites_noLSD1_intersect.bed -b ESC_J1.enhancers.use.bed > FinalAnnotation/RA4_UD_only_sites_ESC_J1_enhancers_intersect.bed
# bedtools intersect -wa -wb -a FinalAnnotation/RA4_diff_sites_noLSD1_intersect.bed -b ESC_J1.enhancers.use.bed > FinalAnnotation/RA4_diff_sites_ESC_J1_enhancers_intersect.bed
# bedtools intersect -wa -wb -a FinalAnnotation/RA4_PG_diff_sites_noLSD1_intersect.bed -b ESC_J1.enhancers.use.bed > FinalAnnotation/RA4_PG_diff_sites_ESC_J1_enhancers_intersect.bed
# bedtools intersect -wa -wb -a FinalAnnotation/RA4_TCP_diff_sites_noLSD1_intersect.bed -b ESC_J1.enhancers.use.bed > FinalAnnotation/RA4_TCP_diff_sites_ESC_J1_enhancers_intersect.bed

# Annotate with second enhancer dataset
bedtools intersect -wa -wb -a FinalAnnotation/RA4_UD_only_sites_noLSD1_intersect.bed -b FinalAnnotation/ESC_J1.enhancers.use.edit.bed > FinalAnnotation/RA4_UD_only_sites_ESC_J1_enhancers_intersect.bed
bedtools intersect -wa -wb -a FinalAnnotation/RA4_diff_sites_noLSD1_intersect.bed -b FinalAnnotation/ESC_J1.enhancers.use.edit.bed > FinalAnnotation/RA4_diff_sites_ESC_J1_enhancers_intersect.bed
bedtools intersect -wa -wb -a FinalAnnotation/RA4_PG_diff_sites_noLSD1_intersect.bed -b FinalAnnotation/ESC_J1.enhancers.use.edit.bed > FinalAnnotation/RA4_PG_diff_sites_ESC_J1_enhancers_intersect.bed
bedtools intersect -wa -wb -a FinalAnnotation/RA4_TCP_diff_sites_noLSD1_intersect.bed -b FinalAnnotation/ESC_J1.enhancers.use.edit.bed > FinalAnnotation/RA4_TCP_diff_sites_ESC_J1_enhancers_intersect.bed

# # Find features that do NOT intersect with ESCJ1 occupied enhancers (these will be used in next steps):
# bedtools intersect -a FinalAnnotation/RA4_UD_only_sites_noLSD1_intersect.bed -b ESC_J1.enhancers.use.bed -v  > FinalAnnotation/RA4_UD_only_sites_noenhancer_intersect.bed
# bedtools intersect -a FinalAnnotation/RA4_diff_sites_noLSD1_intersect.bed -b ESC_J1.enhancers.use.bed -v  > FinalAnnotation/RA4_diff_sites_noenhancer_intersect.bed
# bedtools intersect -a FinalAnnotation/RA4_PG_diff_sites_noLSD1_intersect.bed -b ESC_J1.enhancers.use.bed -v  > FinalAnnotation/RA4_PG_diff_sites_noenhancer_intersect.bed
# bedtools intersect -a FinalAnnotation/RA4_TCP_diff_sites_noLSD1_intersect.bed -b ESC_J1.enhancers.use.bed -v  > FinalAnnotation/RA4_TCP_diff_sites_noenhancer_intersect.bed

# Find features that do NOT intersect with ESCJ1 occupied enhancers (these will be used in next steps):
bedtools intersect -a FinalAnnotation/RA4_UD_only_sites_noLSD1_intersect.bed -b FinalAnnotation/ESC_J1.enhancers.use.edit.bed -v  > FinalAnnotation/RA4_UD_only_sites_noenhancer_intersect.bed
bedtools intersect -a FinalAnnotation/RA4_diff_sites_noLSD1_intersect.bed -b FinalAnnotation/ESC_J1.enhancers.use.edit.bed -v  > FinalAnnotation/RA4_diff_sites_noenhancer_intersect.bed
bedtools intersect -a FinalAnnotation/RA4_PG_diff_sites_noLSD1_intersect.bed -b FinalAnnotation/ESC_J1.enhancers.use.edit.bed -v  > FinalAnnotation/RA4_PG_diff_sites_noenhancer_intersect.bed
bedtools intersect -a FinalAnnotation/RA4_TCP_diff_sites_noLSD1_intersect.bed -b FinalAnnotation/ESC_J1.enhancers.use.edit.bed -v  > FinalAnnotation/RA4_TCP_diff_sites_noenhancer_intersect.bed

# ********************************* Annotate Sites in FinalAnnotation to Promoters ********************************* #
bedtools intersect -wa -wb -a FinalAnnotation/RA4_diff_sites_noenhancer_intersect.bed -b promoters.bed > FinalAnnotation/RA4_diff_sites_noenhancer_to_promotors.bed
bedtools intersect -wa -wb -a FinalAnnotation/RA4_PG_diff_sites_noenhancer_intersect.bed -b promoters.bed > FinalAnnotation/RA4_PG_diff_sites_noenhancer_to_promotors.bed
bedtools intersect -wa -wb -a FinalAnnotation/RA4_TCP_diff_sites_noenhancer_intersect.bed  -b promoters.bed > FinalAnnotation/RA4_TCP_diff_sites_noenhancer_to_promotors.bed
bedtools intersect -wa -wb -a FinalAnnotation/RA4_UD_only_sites_noenhancer_intersect.bed -b promoters.bed > FinalAnnotation/RA4_UD_sites_noenhancer_to_promotors.bed


# Total overlap with LSD1
wc -l FinalAnnotation/RA4_UD_only_sites_LSD1_intersect.bed
wc -l FinalAnnotation/RA4_diff_sites_LSD1_intersect.bed
wc -l FinalAnnotation/RA4_PG_diff_sites_LSD1_intersect.bed
wc -l FinalAnnotation/RA4_TCP_diff_sites_LSD1_intersect.bed

# Total overlap with ESC_J1
wc -l FinalAnnotation/RA4_UD_only_sites_ESC_J1_enhancers_intersect.bed
wc -l FinalAnnotation/RA4_diff_sites_ESC_J1_enhancers_intersect.bed
wc -l FinalAnnotation/RA4_PG_diff_sites_ESC_J1_enhancers_intersect.bed
wc -l FinalAnnotation/RA4_TCP_diff_sites_ESC_J1_enhancers_intersect.bed

# Total all sites 
wc -l FinalAnnotation/RA4_UD_only_sites.bed
wc -l FinalAnnotation/RA4_diff_sites.bed
wc -l FinalAnnotation/RA4_PG_diff_sites.bed
wc -l FinalAnnotation/RA4_TCP_diff_sites.bed

# ********************************* Annotate Sites for each samples ********************************* #
# ************************ Convert all sites into bed files with tab delimited ************************ #
awk 'NR>1 {print "chr"$1, $2, $2, $4, $6}' UD_readsCatalogue/UD_allSites_noDups_final.tsv | awk 'BEGIN{OFS="\t"}{print $1, $2, $3 + 1, $4, $5}' > UD_readsCatalogue/UD_allSites_noDups_final.bed
awk 'NR>1 {print "chr"$1, $2, $2, $4, $6}' RA4_readsCatalogue/RA4_allSites_noDups_final.tsv | awk 'BEGIN{OFS="\t"}{print $1, $2, $3 + 1, $4, $5}' > RA4_readsCatalogue/RA4_allSites_noDups_final.bed
awk 'NR>1 {print "chr"$1, $2, $2, $4, $6}' RA4_TCP_readsCatalogue/RA4_TCP_allSites_noDups_final.tsv | awk 'BEGIN{OFS="\t"}{print $1, $2, $3 + 1, $4, $5}' > RA4_TCP_readsCatalogue/RA4_TCP_allSites_noDups_final.bed
awk 'NR>1 {print "chr"$1, $2, $2, $4, $6}' RA4_PG_readsCatalogue/RA4_PG_allSites_noDups_final.tsv | awk 'BEGIN{OFS="\t"}{print $1, $2, $3 + 1, $4, $5}' > RA4_PG_readsCatalogue/RA4_PG_allSites_noDups_final.bed

# Find overlapping features, print features from first file AND second file that overlap
bedtools intersect -wa -wb -a UD_readsCatalogue/UD_allSites_noDups_final.bed -b FinalAnnotation/LSD1occupied_enhancers.mm10.use.edit.bed  > UD_readsCatalogue/UD_allSites_noDups_LSD1_intersect.bed
bedtools intersect -wa -wb -a RA4_readsCatalogue/RA4_allSites_noDups_final.bed -b FinalAnnotation/LSD1occupied_enhancers.mm10.use.edit.bed  > RA4_readsCatalogue/RA4_allSites_noDups_LSD1_intersect.bed
bedtools intersect -wa -wb -a RA4_TCP_readsCatalogue/RA4_TCP_allSites_noDups_final.bed -b FinalAnnotation/LSD1occupied_enhancers.mm10.use.edit.bed  > RA4_TCP_readsCatalogue/RA4_TCP_allSites_noDups_LSD1_intersect.bed
bedtools intersect -wa -wb -a RA4_PG_readsCatalogue/RA4_PG_allSites_noDups_final.bed -b FinalAnnotation/LSD1occupied_enhancers.mm10.use.edit.bed  > RA4_PG_readsCatalogue/RA4_PG_allSites_noDups_LSD1_intersect.bed

# Find features that do NOT intersect with LSD1 occupied enhancers (these will be used in next steps):
bedtools intersect -a UD_readsCatalogue/UD_allSites_noDups_final.bed  -b FinalAnnotation/LSD1occupied_enhancers.mm10.use.edit.bed -v  > UD_readsCatalogue/UD_allSites_noDups_noLSD1_intersect.bed
bedtools intersect -a RA4_readsCatalogue/RA4_allSites_noDups_final.bed -b FinalAnnotation/LSD1occupied_enhancers.mm10.use.edit.bed -v  > RA4_readsCatalogue/RA4_allSites_noDups_noLSD1_intersect.bed
bedtools intersect -a RA4_TCP_readsCatalogue/RA4_TCP_allSites_noDups_final.bed -b FinalAnnotation/LSD1occupied_enhancers.mm10.use.edit.bed -v  > RA4_TCP_readsCatalogue/RA4_TCP_allSites_noDups_noLSD1_intersect.bed
bedtools intersect -a RA4_PG_readsCatalogue/RA4_PG_allSites_noDups_final.bed -b FinalAnnotation/LSD1occupied_enhancers.mm10.use.edit.bed -v  > RA4_PG_readsCatalogue/RA4_PG_allSites_noDups_noLSD1_intersect.bed

# Annotate with second enhancer dataset
bedtools intersect -wa -wb -a UD_readsCatalogue/UD_allSites_noDups_noLSD1_intersect.bed -b FinalAnnotation/ESC_J1.enhancers.use.edit.bed > UD_readsCatalogue/UD_allSites_noDups_ESC_J1_enhancers_intersect.bed
bedtools intersect -wa -wb -a RA4_readsCatalogue/RA4_allSites_noDups_noLSD1_intersect.bed -b FinalAnnotation/ESC_J1.enhancers.use.edit.bed > RA4_readsCatalogue/RA4_allSites_noDups_ESC_J1_enhancers_intersect.bed
bedtools intersect -wa -wb -a RA4_TCP_readsCatalogue/RA4_TCP_allSites_noDups_noLSD1_intersect.bed -b FinalAnnotation/ESC_J1.enhancers.use.edit.bed > RA4_TCP_readsCatalogue/RA4_TCP_allSites_noDups_ESC_J1_enhancers_intersect.bed
bedtools intersect -wa -wb -a RA4_PG_readsCatalogue/RA4_PG_allSites_noDups_noLSD1_intersect.bed -b FinalAnnotation/ESC_J1.enhancers.use.edit.bed > RA4_PG_readsCatalogue/RA4_PG_allSites_noDups_ESC_J1_enhancers_intersect.bed

# Find features that do NOT intersect with ESCJ1 occupied enhancers (these will be used in next steps):
bedtools intersect -a UD_readsCatalogue/UD_allSites_noDups_noLSD1_intersect.bed -b FinalAnnotation/ESC_J1.enhancers.use.edit.bed -v  > UD_readsCatalogue/UD_allSites_noDups_noenhancer_intersect.bed
bedtools intersect -a RA4_readsCatalogue/RA4_allSites_noDups_noLSD1_intersect.bed -b FinalAnnotation/ESC_J1.enhancers.use.edit.bed -v  > RA4_readsCatalogue/RA4_allSites_noDups_noenhancer_intersect.bed
bedtools intersect -a RA4_TCP_readsCatalogue/RA4_TCP_allSites_noDups_noLSD1_intersect.bed -b FinalAnnotation/ESC_J1.enhancers.use.edit.bed -v  > RA4_TCP_readsCatalogue/RA4_TCP_allSites_noDups_noenhancer_intersect.bed
bedtools intersect -a RA4_PG_readsCatalogue/RA4_PG_allSites_noDups_noLSD1_intersect.bed -b FinalAnnotation/ESC_J1.enhancers.use.edit.bed -v  > RA4_PG_readsCatalogue/RA4_PG_allSites_noDups_noenhancer_intersect.bed

# ********************************* Annotate Sites for each samples on Promoters ********************************* #
bedtools intersect -wa -wb -a UD_readsCatalogue/UD_allSites_noDups_noenhancer_intersect.bed -b promoters.bed > UD_readsCatalogue/UD_allSites_noDups_noenhancer_promoters_intersect.bed
bedtools intersect -wa -wb -a RA4_readsCatalogue/RA4_allSites_noDups_noenhancer_intersect.bed -b promoters.bed > RA4_readsCatalogue/RA4_allSites_noDups_noenhancer_promoters_intersect.bed
bedtools intersect -wa -wb -a RA4_TCP_readsCatalogue/RA4_TCP_allSites_noDups_noenhancer_intersect.bed -b promoters.bed > RA4_TCP_readsCatalogue/RA4_TCP_allSites_noDups_noenhancer_promoters_intersect.bed
bedtools intersect -wa -wb -a RA4_PG_readsCatalogue/RA4_PG_allSites_noDups_noenhancer_intersect.bed -b promoters.bed > RA4_PG_readsCatalogue/RA4_PG_allSites_noDups_noenhancer_promoters_intersect.bed

# Total overlap with LSD1
wc -l UD_readsCatalogue/UD_allSites_noDups_LSD1_intersect.bed
wc -l RA4_readsCatalogue/RA4_allSites_noDups_LSD1_intersect.bed
wc -l RA4_TCP_readsCatalogue/RA4_TCP_allSites_noDups_LSD1_intersect.bed
wc -l RA4_PG_readsCatalogue/RA4_PG_allSites_noDups_LSD1_intersect.bed

# Total overlap with ESC_J1
wc -l UD_readsCatalogue/UD_allSites_noDups_ESC_J1_enhancers_intersect.bed
wc -l RA4_readsCatalogue/RA4_allSites_noDups_ESC_J1_enhancers_intersect.bed
wc -l RA4_TCP_readsCatalogue/RA4_TCP_allSites_noDups_ESC_J1_enhancers_intersect.bed
wc -l RA4_PG_readsCatalogue/RA4_PG_allSites_noDups_ESC_J1_enhancers_intersect.bed

# Total all sites 
wc -l UD_readsCatalogue/UD_allSites_noDups_final.bed
wc -l RA4_readsCatalogue/RA4_allSites_noDups_final.bed
wc -l RA4_TCP_readsCatalogue/RA4_TCP_allSites_noDups_final.bed
wc -l RA4_PG_readsCatalogue/RA4_PG_allSites_noDups_final.bed

# ********************************* Annotate RA4 - UD, RA4_PG - UD, and RA4_TCP - UD ********************************* #
# Intersect sites with count data
bedtools intersect -wa -a RA4_readsCatalogue/RA4_allSites_noDups_final.bed -b Results/UD_vs_RA4/RA4_diff_sites.bed > Results/UD_vs_RA4/RA4_diff_sites_counts.bed
bedtools intersect -wa -a RA4_PG_readsCatalogue/RA4_PG_allSites_noDups_final.bed -b Results/UD_vs_RA4_PG/RA4_PG_diff_sites.bed > Results/UD_vs_RA4_PG/RA4_PG_diff_sites_counts.bed
bedtools intersect -wa -a RA4_TCP_readsCatalogue/RA4_TCP_allSites_noDups_final.bed -b Results/UD_vs_RA4_TCP/RA4_TCP_diff_sites.bed > Results/UD_vs_RA4_TCP/RA4_TCP_diff_sites_counts.bed

# Annotate with LSD1 enhancers
bedtools intersect -wa -wb -a Results/UD_vs_RA4/RA4_diff_sites_counts.bed -b FinalAnnotation/LSD1occupied_enhancers.mm10.use.edit.bed  > Results/UD_vs_RA4/RA4_diff_sites_LSD1_intersect.bed
bedtools intersect -wa -wb -a Results/UD_vs_RA4_PG/RA4_PG_diff_sites_counts.bed -b FinalAnnotation/LSD1occupied_enhancers.mm10.use.edit.bed  > Results/UD_vs_RA4_PG/RA4_PG_diff_sites_LSD1_intersect.bed 
bedtools intersect -wa -wb -a Results/UD_vs_RA4_TCP/RA4_TCP_diff_sites_counts.bed -b FinalAnnotation/LSD1occupied_enhancers.mm10.use.edit.bed  > Results/UD_vs_RA4_TCP/RA4_TCP_diff_sites_LSD1_intersect.bed

# Minus annotation with LSD1 enhancers
bedtools intersect -a Results/UD_vs_RA4/RA4_diff_sites_counts.bed -b FinalAnnotation/LSD1occupied_enhancers.mm10.use.edit.bed -v > Results/UD_vs_RA4/RA4_diff_sites_noLSD1_intersect.bed
bedtools intersect -a Results/UD_vs_RA4_PG/RA4_PG_diff_sites_counts.bed -b FinalAnnotation/LSD1occupied_enhancers.mm10.use.edit.bed  -v > Results/UD_vs_RA4_PG/RA4_PG_diff_sites_noLSD1_intersect.bed 
bedtools intersect -a Results/UD_vs_RA4_TCP/RA4_TCP_diff_sites_counts.bed -b FinalAnnotation/LSD1occupied_enhancers.mm10.use.edit.bed  -v > Results/UD_vs_RA4_TCP/RA4_TCP_diff_sites_noLSD1_intersect.bed

# Annotate with ESC enhancer dataset
bedtools intersect -wa -wb -a Results/UD_vs_RA4/RA4_diff_sites_noLSD1_intersect.bed -b FinalAnnotation/ESC_J1.enhancers.use.edit.bed > Results/UD_vs_RA4/RA4_diff_sites_ESC_J1_intersect.bed
bedtools intersect -wa -wb -a Results/UD_vs_RA4_PG/RA4_PG_diff_sites_noLSD1_intersect.bed  -b FinalAnnotation/ESC_J1.enhancers.use.edit.bed > Results/UD_vs_RA4_PG/RA4_PG_diff_sites_ESC_J1_intersect.bed
bedtools intersect -wa -wb -a Results/UD_vs_RA4_TCP/RA4_TCP_diff_sites_noLSD1_intersect.bed -b FinalAnnotation/ESC_J1.enhancers.use.edit.bed > Results/UD_vs_RA4_TCP/RA4_TCP_diff_sites_ESC_J1_intersect.bed

# Minus annotate with ESC enhancer dataset
bedtools intersect -a Results/UD_vs_RA4/RA4_diff_sites_noLSD1_intersect.bed -b FinalAnnotation/ESC_J1.enhancers.use.edit.bed -v > Results/UD_vs_RA4/RA4_diff_sites_noenhancers_intersect.bed
bedtools intersect -a Results/UD_vs_RA4_PG/RA4_PG_diff_sites_noLSD1_intersect.bed  -b FinalAnnotation/ESC_J1.enhancers.use.edit.bed -v > Results/UD_vs_RA4_PG/RA4_PG_diff_sites_noenhancers_intersect.bed
bedtools intersect -a Results/UD_vs_RA4_TCP/RA4_TCP_diff_sites_noLSD1_intersect.bed -b FinalAnnotation/ESC_J1.enhancers.use.edit.bed -v > Results/UD_vs_RA4_TCP/RA4_TCP_diff_sites_noenhancers_intersect.bed

# Total overlap with LSD1
wc -l Results/UD_vs_RA4/RA4_diff_sites_LSD1_intersect.bed
wc -l Results/UD_vs_RA4_PG/RA4_PG_diff_sites_LSD1_intersect.bed 
wc -l Results/UD_vs_RA4_TCP/RA4_TCP_diff_sites_LSD1_intersect.bed

# Total overlap with ESC_J1
wc -l Results/UD_vs_RA4/RA4_diff_sites_ESC_J1_intersect.bed
wc -l Results/UD_vs_RA4_PG/RA4_PG_diff_sites_ESC_J1_intersect.bed 
wc -l Results/UD_vs_RA4_TCP/RA4_TCP_diff_sites_ESC_J1_intersect.bed

# # ********************************* Find Sites in RA4 minus UD, TCP, and PG that are in promoters and enhancers ********************************* #
# # First intersect with enhancers in RA4
# bedtools intersect -wa -wb -a Results/RA4_D_vs_PG_D_vs_TCP_D/RA4_UD_only_sites.bed -b RA4_readsCatalogue/RA4_allSites_noDups_noenhancer_promoters_intersect.bed > Results/RA4_D_vs_PG_D_vs_TCP_D/RA4_UD_only_sites_annoPromotersAll.bed
# cut -f1-3 --complement Results/RA4_D_vs_PG_D_vs_TCP_D/RA4_UD_only_sites_annoPromotersAll.bed > Results/RA4_D_vs_PG_D_vs_TCP_D/RA4_UD_only_sites_annoPromoters.bed
# bedtools intersect -wa -wb -a Results/RA4_D_vs_PG_D_vs_TCP_D/RA4_UD_only_sites.bed -b RA4_readsCatalogue/RA4_allSites_noDups_ESC_J1_enhancers_intersect.bed > Results/RA4_D_vs_PG_D_vs_TCP_D/RA4_UD_only_sites_annoESCAll.bed
# cut -f1-3 --complement Results/RA4_D_vs_PG_D_vs_TCP_D/RA4_UD_only_sites_annoESCAll.bed > Results/RA4_D_vs_PG_D_vs_TCP_D/RA4_UD_only_sites_annoESC.bed

# # ********************************* Find Sites in RA4_PG - UD - RA4 and Annotate with LSD1 and ESC ********************************* #
# bedtools intersect -a Results/UD_vs_RA4_PG/RA4_PG_diff_sites.bed -b RA4_readsCatalogue/RA4_allSites_noDups_final.bed  -v > Results/RA4_PG_Minus_UD_Minus_RA4/RA4_PG_Minus_UD_Minus_RA4_sites.bed
# bedtools intersect -wa -wb -a Results/RA4_PG_Minus_UD_Minus_RA4/RA4_PG_Minus_UD_Minus_RA4_sites.bed -b FinalAnnotation/LSD1occupied_enhancers.mm10.use.bed > Results/RA4_PG_Minus_UD_Minus_RA4/RA4_PG_Minus_UD_Minus_RA4_sites_annoLSD1.bed
# bedtools intersect -a Results/RA4_PG_Minus_UD_Minus_RA4/RA4_PG_Minus_UD_Minus_RA4_sites.bed -b FinalAnnotation/LSD1occupied_enhancers.mm10.use.bed -v > Results/RA4_PG_Minus_UD_Minus_RA4/RA4_PG_Minus_UD_Minus_RA4_sites_noLSD1.bed
# bedtools intersect -wa -wb -a Results/RA4_PG_Minus_UD_Minus_RA4/RA4_PG_Minus_UD_Minus_RA4_sites_noLSD1.bed -b FinalAnnotation/ESC_J1.enhancers.use.bed > Results/RA4_PG_Minus_UD_Minus_RA4/RA4_PG_Minus_UD_Minus_RA4_sites_annoESC_J1.bed
# bedtools intersect -a Results/RA4_PG_Minus_UD_Minus_RA4/RA4_PG_Minus_UD_Minus_RA4_sites_noLSD1.bed -b FinalAnnotation/ESC_J1.enhancers.use.bed -v > Results/RA4_PG_Minus_UD_Minus_RA4/RA4_PG_Minus_UD_Minus_RA4_sites_noenhancers.bed

# # ********************************* Find Sites in RA4_TCP - UD - RA4 and Annotate with LSD1 and ESC ********************************* #
# bedtools intersect -a Results/UD_vs_RA4_TCP/RA4_TCP_diff_sites.bed -b RA4_readsCatalogue/RA4_allSites_noDups_final.bed  -v > Results/RA4_TCP_Minus_UD_Minus_RA4/RA4_TCP_Minus_UD_Minus_RA4_sites.bed
# bedtools intersect -wa -wb -a Results/RA4_TCP_Minus_UD_Minus_RA4/RA4_TCP_Minus_UD_Minus_RA4_sites.bed -b FinalAnnotation/LSD1occupied_enhancers.mm10.use.bed > Results/RA4_TCP_Minus_UD_Minus_RA4/RA4_TCP_Minus_UD_Minus_RA4_sites_annoLSD1.bed
# bedtools intersect -a Results/RA4_TCP_Minus_UD_Minus_RA4/RA4_TCP_Minus_UD_Minus_RA4_sites.bed -b FinalAnnotation/LSD1occupied_enhancers.mm10.use.bed -v > Results/RA4_TCP_Minus_UD_Minus_RA4/RA4_TCP_Minus_UD_Minus_RA4_sites_noLSD1.bed
# bedtools intersect -wa -wb -a Results/RA4_TCP_Minus_UD_Minus_RA4/RA4_TCP_Minus_UD_Minus_RA4_sites_noLSD1.bed -b FinalAnnotation/ESC_J1.enhancers.use.bed > Results/RA4_TCP_Minus_UD_Minus_RA4/RA4_TCP_Minus_UD_Minus_RA4_sites_annoESC_J1.bed
# bedtools intersect -a Results/RA4_TCP_Minus_UD_Minus_RA4/RA4_TCP_Minus_UD_Minus_RA4_sites_noLSD1.bed -b FinalAnnotation/ESC_J1.enhancers.use.bed -v > Results/RA4_TCP_Minus_UD_Minus_RA4/RA4_TCP_Minus_UD_Minus_RA4_sites_noenhancers.bed

# # **************** Annotate with LSD1 and ESC sites for comparison between RA4_TCP minus UD and RA4_PG minus UD ************************* #
# # Sites only in RA4_PG minus UD
# bedtools intersect -wa -wb -a Results/RA4_D_vs_PG_D_vs_TCP_D/RA4_PG_Min_UD_vs_RA4_TCP_Min_UD/RA4_PG_Min_UD_only_sites.bed -b FinalAnnotation/LSD1occupied_enhancers.mm10.use.bed > Results/RA4_D_vs_PG_D_vs_TCP_D/RA4_PG_Min_UD_vs_RA4_TCP_Min_UD/RA4_PG_Min_UD_only_sites_annoLSD1.bed
# bedtools intersect -a Results/RA4_D_vs_PG_D_vs_TCP_D/RA4_PG_Min_UD_vs_RA4_TCP_Min_UD/RA4_PG_Min_UD_only_sites.bed -b FinalAnnotation/LSD1occupied_enhancers.mm10.use.bed -v > Results/RA4_D_vs_PG_D_vs_TCP_D/RA4_PG_Min_UD_vs_RA4_TCP_Min_UD/RA4_PG_Min_UD_only_sites_noLSD1.bed
# bedtools intersect -wa -wb -a Results/RA4_D_vs_PG_D_vs_TCP_D/RA4_PG_Min_UD_vs_RA4_TCP_Min_UD/RA4_PG_Min_UD_only_sites_noLSD1.bed -b FinalAnnotation/ESC_J1.enhancers.use.bed > Results/RA4_D_vs_PG_D_vs_TCP_D/RA4_PG_Min_UD_vs_RA4_TCP_Min_UD/RA4_PG_Min_UD_only_sites_annoESC1_J1.bed
# bedtools intersect -a Results/RA4_D_vs_PG_D_vs_TCP_D/RA4_PG_Min_UD_vs_RA4_TCP_Min_UD/RA4_PG_Min_UD_only_sites_noLSD1.bed -b FinalAnnotation/ESC_J1.enhancers.use.bed -v > Results/RA4_D_vs_PG_D_vs_TCP_D/RA4_PG_Min_UD_vs_RA4_TCP_Min_UD/RA4_PG_Min_UD_only_sites_noenhancers.bed
# bedtools intersect -wa -wb -a Results/RA4_D_vs_PG_D_vs_TCP_D/RA4_PG_Min_UD_vs_RA4_TCP_Min_UD/RA4_PG_Min_UD_only_sites_noenhancers.bed -b promoters.bed -v > Results/RA4_D_vs_PG_D_vs_TCP_D/RA4_PG_Min_UD_vs_RA4_TCP_Min_UD/RA4_PG_Min_UD_only_sites_noenhancers_promoters.bed

# # Sites only in RA4_TCP minus UD
# bedtools intersect -wa -wb -a Results/RA4_D_vs_PG_D_vs_TCP_D/RA4_PG_Min_UD_vs_RA4_TCP_Min_UD/RA4_TCP_Min_UD_only_sites.bed -b FinalAnnotation/LSD1occupied_enhancers.mm10.use.bed > Results/RA4_D_vs_PG_D_vs_TCP_D/RA4_PG_Min_UD_vs_RA4_TCP_Min_UD/RA4_TCP_Min_UD_only_sites_annoLSD1.bed
# bedtools intersect -a Results/RA4_D_vs_PG_D_vs_TCP_D/RA4_PG_Min_UD_vs_RA4_TCP_Min_UD/RA4_TCP_Min_UD_only_sites.bed -b FinalAnnotation/LSD1occupied_enhancers.mm10.use.bed -v > Results/RA4_D_vs_PG_D_vs_TCP_D/RA4_PG_Min_UD_vs_RA4_TCP_Min_UD/RA4_TCP_Min_UD_only_sites_noLSD1.bed
# bedtools intersect -wa -wb -a Results/RA4_D_vs_PG_D_vs_TCP_D/RA4_PG_Min_UD_vs_RA4_TCP_Min_UD/RA4_TCP_Min_UD_only_sites_noLSD1.bed -b FinalAnnotation/ESC_J1.enhancers.use.bed > Results/RA4_D_vs_PG_D_vs_TCP_D/RA4_PG_Min_UD_vs_RA4_TCP_Min_UD/RA4_TCP_Min_UD_only_sites_annoESC1_J1.bed
# bedtools intersect -a Results/RA4_D_vs_PG_D_vs_TCP_D/RA4_PG_Min_UD_vs_RA4_TCP_Min_UD/RA4_TCP_Min_UD_only_sites_noLSD1.bed -b FinalAnnotation/ESC_J1.enhancers.use.bed -v > Results/RA4_D_vs_PG_D_vs_TCP_D/RA4_PG_Min_UD_vs_RA4_TCP_Min_UD/RA4_TCP_Min_UD_only_sites_noenhancers.bed
# bedtools intersect -wa -wb -a Results/RA4_D_vs_PG_D_vs_TCP_D/RA4_PG_Min_UD_vs_RA4_TCP_Min_UD/RA4_TCP_Min_UD_only_sites_noenhancers.bed -b promoters.bed -v > Results/RA4_D_vs_PG_D_vs_TCP_D/RA4_PG_Min_UD_vs_RA4_TCP_Min_UD/RA4_TCP_Min_UD_only_sites_noenhancers_promoters.bed

# # Sites in intersection of RA4_TCP minus UD and RA4_PG minus UD
# bedtools intersect -wa -wb -a Results/RA4_D_vs_PG_D_vs_TCP_D/RA4_PG_Min_UD_vs_RA4_TCP_Min_UD/RA4_TCP_minus_UD_intersect_RA4_PG_minus_UD_sites.bed -b FinalAnnotation/LSD1occupied_enhancers.mm10.use.bed > Results/RA4_D_vs_PG_D_vs_TCP_D/RA4_PG_Min_UD_vs_RA4_TCP_Min_UD/RA4_TCP_minus_UD_intersect_RA4_PG_minus_UD_sites_annoLSD1.bed
# bedtools intersect -a Results/RA4_D_vs_PG_D_vs_TCP_D/RA4_PG_Min_UD_vs_RA4_TCP_Min_UD/RA4_TCP_minus_UD_intersect_RA4_PG_minus_UD_sites.bed -b FinalAnnotation/LSD1occupied_enhancers.mm10.use.bed -v > Results/RA4_D_vs_PG_D_vs_TCP_D/RA4_PG_Min_UD_vs_RA4_TCP_Min_UD/RA4_TCP_minus_UD_intersect_RA4_PG_minus_UD_sites_noLSD1.bed
# bedtools intersect -wa -wb -a Results/RA4_D_vs_PG_D_vs_TCP_D/RA4_PG_Min_UD_vs_RA4_TCP_Min_UD/RA4_TCP_minus_UD_intersect_RA4_PG_minus_UD_sites_noLSD1.bed -b FinalAnnotation/ESC_J1.enhancers.use.bed > Results/RA4_D_vs_PG_D_vs_TCP_D/RA4_PG_Min_UD_vs_RA4_TCP_Min_UD/RA4_TCP_minus_UD_intersect_RA4_PG_minus_UD_sites_annoESC1_J1.bed
# bedtools intersect -a Results/RA4_D_vs_PG_D_vs_TCP_D/RA4_PG_Min_UD_vs_RA4_TCP_Min_UD/RA4_TCP_minus_UD_intersect_RA4_PG_minus_UD_sites_noLSD1.bed -b FinalAnnotation/ESC_J1.enhancers.use.bed -v > Results/RA4_D_vs_PG_D_vs_TCP_D/RA4_PG_Min_UD_vs_RA4_TCP_Min_UD/RA4_TCP_minus_UD_intersect_RA4_PG_minus_UD_sites_noenhancers.bed
# bedtools intersect -wa -wb -a Results/RA4_D_vs_PG_D_vs_TCP_D/RA4_PG_Min_UD_vs_RA4_TCP_Min_UD/RA4_TCP_minus_UD_intersect_RA4_PG_minus_UD_sites_noenhancers.bed -b promoters.bed -v > Results/RA4_D_vs_PG_D_vs_TCP_D/RA4_PG_Min_UD_vs_RA4_TCP_Min_UD/RA4_TCP_minus_UD_intersect_RA4_PG_minus_UD_sites_noenhancers_promoters.bed

# # **************** Annotate with LSD1 and ESC sites in the overlap of RA4_UD, RA4_PG_UD, and RA4_TCP_UD ************************* #
# bedtools intersect -wa -wb -a Results/RA4_D_vs_PG_D_vs_TCP_D/RA4_vs_PG_vs_TCP_intersection_sites.bed -b FinalAnnotation/LSD1occupied_enhancers.mm10.use.bed > Results/RA4_D_vs_PG_D_vs_TCP_D/RA4_vs_PG_vs_TCP_intersection_sites_annoLSD1.bed
# bedtools intersect -a Results/RA4_D_vs_PG_D_vs_TCP_D/RA4_vs_PG_vs_TCP_intersection_sites.bed -b FinalAnnotation/LSD1occupied_enhancers.mm10.use.bed -v > Results/RA4_D_vs_PG_D_vs_TCP_D/RA4_vs_PG_vs_TCP_intersection_sites_noLSD1.bed
# bedtools intersect -wa -wb -a Results/RA4_D_vs_PG_D_vs_TCP_D/RA4_vs_PG_vs_TCP_intersection_sites_noLSD1.bed -b FinalAnnotation/ESC_J1.enhancers.use.bed > Results/RA4_D_vs_PG_D_vs_TCP_D/RA4_vs_PG_vs_TCP_intersection_sites_annoESC1_J1.bed
# bedtools intersect -a Results/RA4_D_vs_PG_D_vs_TCP_D/RA4_vs_PG_vs_TCP_intersection_sites_noLSD1.bed -b FinalAnnotation/ESC_J1.enhancers.use.bed -v > Results/RA4_D_vs_PG_D_vs_TCP_D/RA4_vs_PG_vs_TCP_intersection_sites_noenhancers.bed







