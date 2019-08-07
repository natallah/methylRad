# # ************************ Convert gtf files to bed files ***********************
# awk 'OFS="\t" {if($3 == "gene"){print "chr"$1,$4-1,$5,$10,$14,$6,$7}}' Mus_musculus.GRCm38.93.chr.gtf | tr -d '";' > Mus_musculus.GRCm38.93.chr.bed

# # ************************ Convert all sites into bed files with tab delimited ************************
# awk 'NR>1 {print "chr"$1, $2, $2}' UD_readsCatalogue/UD_allSites.tsv | awk 'BEGIN{OFS="\t"}{print $1, $2, $3 + 1}' > UD_readsCatalogue/UD_allSites.bed
# awk 'NR>1 {print "chr"$1, $2, $2}' RA4_readsCatalogue/RA4_allSites.tsv | awk 'BEGIN{OFS="\t"}{print $1, $2, $3 + 1}' > RA4_readsCatalogue/RA4_allSites.bed
# awk 'NR>1 {print "chr"$1, $2, $2}' RA4_TCP_readsCatalogue/RA4_TCP_allSites.tsv | awk 'BEGIN{OFS="\t"}{print $1, $2, $3 + 1}' > RA4_TCP_readsCatalogue/RA4_TCP_allSites.bed
# awk 'NR>1 {print "chr"$1, $2, $2}' RA4_PG_readsCatalogue/RA4_PG_allSites.tsv | awk 'BEGIN{OFS="\t"}{print $1, $2, $3 + 1}' > RA4_PG_readsCatalogue/RA4_PG_allSites.bed

# ************************ Use bedtools intersect to annotate all sites ************************
bedtools intersect -a UD_readsCatalogue/UD_allSites.bed -b Mus_musculus.GRCm38.93.chr.bed -wb | cut -f1,2,3,7,8,10 > UD_readsCatalogue/UD_allSitesAnno.bed 
bedtools intersect -a RA4_readsCatalogue/RA4_allSites.bed -b Mus_musculus.GRCm38.93.chr.bed -wb | cut -f1,2,3,7,8,10 > RA4_readsCatalogue/RA4_allSitesAnno.bed 
bedtools intersect -a RA4_TCP_readsCatalogue/RA4_TCP_allSites.bed -b Mus_musculus.GRCm38.93.chr.bed -wb | cut -f1,2,3,7,8,10 > RA4_TCP_readsCatalogue/RA4_TCP_allSitesAnno.bed 
bedtools intersect -a RA4_PG_readsCatalogue/RA4_PG_allSites.bed -b Mus_musculus.GRCm38.93.chr.bed -wb | cut -f1,2,3,7,8,10 > RA4_PG_readsCatalogue/RA4_PG_allSitesAnno.bed 

# ************************ Use bedtools intersect to annotate all sites comparisons overlap ************************
bedtools intersect -a Results/UD_vs_RA4/UD_diff_sites.bed -b Mus_musculus.GRCm38.93.chr.bed -wb | cut -f1,2,3,7,8,10 > Results/UD_vs_RA4/UD_diff_sitesAnno.bed
bedtools intersect -a Results/UD_vs_RA4/RA4_diff_sites.bed -b Mus_musculus.GRCm38.93.chr.bed -wb | cut -f1,2,3,7,8,10 > Results/UD_vs_RA4/RA4_diff_sitesAnno.bed
bedtools intersect -a Results/UD_vs_RA4/Intersection_sites.bed -b Mus_musculus.GRCm38.93.chr.bed -wb | cut -f1,2,3,7,8,10 > Results/UD_vs_RA4/Intersection_sitesAnno.bed

bedtools intersect -a Results/UD_vs_RA4_PG/UD_diff_sites.bed -b Mus_musculus.GRCm38.93.chr.bed -wb | cut -f1,2,3,7,8,10 > Results/UD_vs_RA4_PG/UD_diff_sitesAnno.bed
bedtools intersect -a Results/UD_vs_RA4_PG/RA4_PG_diff_sites.bed -b Mus_musculus.GRCm38.93.chr.bed -wb | cut -f1,2,3,7,8,10 > Results/UD_vs_RA4_PG/RA4_PG_diff_sitesAnno.bed
bedtools intersect -a Results/UD_vs_RA4_PG/Intersection_sites.bed -b Mus_musculus.GRCm38.93.chr.bed -wb | cut -f1,2,3,7,8,10 > Results/UD_vs_RA4_PG/Intersection_sitesAnno.bed

bedtools intersect -a Results/UD_vs_RA4_TCP/UD_diff_sites.bed -b Mus_musculus.GRCm38.93.chr.bed -wb | cut -f1,2,3,7,8,10 > Results/UD_vs_RA4_TCP/UD_diff_sitesAnno.bed
bedtools intersect -a Results/UD_vs_RA4_TCP/RA4_TCP_diff_sites.bed -b Mus_musculus.GRCm38.93.chr.bed -wb | cut -f1,2,3,7,8,10 > Results/UD_vs_RA4_TCP/RA4_TCP_diff_sitesAnno.bed
bedtools intersect -a Results/UD_vs_RA4_TCP/Intersection_sites.bed -b Mus_musculus.GRCm38.93.chr.bed -wb | cut -f1,2,3,7,8,10 > Results/UD_vs_RA4_TCP/Intersection_sitesAnno.bed

bedtools intersect -a Results/RA4_vs_RA4_PG/RA4_diff_sites.bed -b Mus_musculus.GRCm38.93.chr.bed -wb | cut -f1,2,3,7,8,10 > Results/RA4_vs_RA4_PG/RA4_diff_sitesAnno.bed
bedtools intersect -a Results/RA4_vs_RA4_PG/RA4_PG_diff_sites.bed -b Mus_musculus.GRCm38.93.chr.bed -wb | cut -f1,2,3,7,8,10 > Results/RA4_vs_RA4_PG/RA4_PG_diff_sitesAnno.bed
bedtools intersect -a Results/RA4_vs_RA4_PG/Intersection_sites.bed -b Mus_musculus.GRCm38.93.chr.bed -wb | cut -f1,2,3,7,8,10 > Results/RA4_vs_RA4_PG/Intersection_sitesAnno.bed

bedtools intersect -a Results/RA4_vs_RA4_TCP/RA4_diff_sites.bed -b Mus_musculus.GRCm38.93.chr.bed -wb | cut -f1,2,3,7,8,10 > Results/RA4_vs_RA4_TCP/RA4_diff_sitesAnno.bed
bedtools intersect -a Results/RA4_vs_RA4_TCP/RA4_TCP_diff_sites.bed -b Mus_musculus.GRCm38.93.chr.bed -wb | cut -f1,2,3,7,8,10 > Results/RA4_vs_RA4_TCP/RA4_TCP_diff_sitesAnno.bed
bedtools intersect -a Results/RA4_vs_RA4_TCP/Intersection_sites.bed -b Mus_musculus.GRCm38.93.chr.bed -wb | cut -f1,2,3,7,8,10 > Results/RA4_vs_RA4_TCP/Intersection_sitesAnno.bed

bedtools intersect -a Results/RA4_PG_vs_RA4_TCP/RA4_PG_diff_sites.bed -b Mus_musculus.GRCm38.93.chr.bed -wb | cut -f1,2,3,7,8,10 > Results/RA4_PG_vs_RA4_TCP/RA4_PG_diff_sitesAnno.bed
bedtools intersect -a Results/RA4_PG_vs_RA4_TCP/RA4_TCP_diff_sites.bed -b Mus_musculus.GRCm38.93.chr.bed -wb | cut -f1,2,3,7,8,10 > Results/RA4_PG_vs_RA4_TCP/RA4_TCP_diff_sitesAnno.bed
bedtools intersect -a Results/RA4_PG_vs_RA4_TCP/Intersection_sites.bed -b Mus_musculus.GRCm38.93.chr.bed -wb | cut -f1,2,3,7,8,10 > Results/RA4_PG_vs_RA4_TCP/Intersection_sitesAnno.bed

bedtools intersect -a Results/RA_PG_TCP_vs_UD/RA_diff_sites.bed -b Mus_musculus.GRCm38.93.chr.bed -wb | cut -f1,2,3,7,8,10 > Results/RA_PG_TCP_vs_UD/RA_diff_sitesAnno.bed
bedtools intersect -a Results/RA_PG_TCP_vs_UD/UD_diff_sites.bed -b Mus_musculus.GRCm38.93.chr.bed -wb | cut -f1,2,3,7,8,10 > Results/RA_PG_TCP_vs_UD/UD_sitesAnno.bed
bedtools intersect -a Results/RA_PG_TCP_vs_UD/Intersection_sites.bed -b Mus_musculus.GRCm38.93.chr.bed -wb | cut -f1,2,3,7,8,10 > Results/RA_PG_TCP_vs_UD/Intersection_sitesAnno.bed

# ************************ Use bedtools nearest to annotate all sites comparisons ************************
bedtools closest -D ref -a Results/UD_vs_RA4/UD_diff_sites.bed -b Mus_musculus.GRCm38.93.chr.bed -wb | cut -f1,2,3,7,8,10,11 > Results/UD_vs_RA4/UD_diff_sitesAnnoNearest.bed
bedtools closest -D ref -a Results/UD_vs_RA4/RA4_diff_sites.bed -b Mus_musculus.GRCm38.93.chr.bed -wb | cut -f1,2,3,7,8,10,11 > Results/UD_vs_RA4/RA4_diff_sitesAnnoNearest.bed
bedtools closest -D ref -a Results/UD_vs_RA4/Intersection_sites.bed -b Mus_musculus.GRCm38.93.chr.bed -wb | cut -f1,2,3,7,8,10,11 > Results/UD_vs_RA4/Intersection_sitesAnnoNearest.bed

bedtools closest -D ref -a Results/UD_vs_RA4_PG/UD_diff_sites.bed -b Mus_musculus.GRCm38.93.chr.bed -wb | cut -f1,2,3,7,8,10,11 > Results/UD_vs_RA4_PG/UD_diff_sitesAnnoNearest.bed
bedtools closest -D ref -a Results/UD_vs_RA4_PG/RA4_PG_diff_sites.bed -b Mus_musculus.GRCm38.93.chr.bed -wb | cut -f1,2,3,7,8,10,11 > Results/UD_vs_RA4_PG/RA4_PG_diff_sitesAnnoNearest.bed
bedtools closest -D ref -a Results/UD_vs_RA4_PG/Intersection_sites.bed -b Mus_musculus.GRCm38.93.chr.bed -wb | cut -f1,2,3,7,8,10,11 > Results/UD_vs_RA4_PG/Intersection_sitesAnnoNearest.bed

bedtools closest -D ref -a Results/UD_vs_RA4_TCP/UD_diff_sites.bed -b Mus_musculus.GRCm38.93.chr.bed -wb | cut -f1,2,3,7,8,10,11 > Results/UD_vs_RA4_TCP/UD_diff_sitesAnnoNearest.bed
bedtools closest -D ref -a Results/UD_vs_RA4_TCP/RA4_TCP_diff_sites.bed -b Mus_musculus.GRCm38.93.chr.bed -wb | cut -f1,2,3,7,8,10,11 > Results/UD_vs_RA4_TCP/RA4_TCP_diff_sitesAnnoNearest.bed
bedtools closest -D ref -a Results/UD_vs_RA4_TCP/Intersection_sites.bed -b Mus_musculus.GRCm38.93.chr.bed -wb | cut -f1,2,3,7,8,10,11 > Results/UD_vs_RA4_TCP/Intersection_sitesAnnoNearest.bed

bedtools closest -D ref -a Results/RA4_vs_RA4_PG/RA4_diff_sites.bed -b Mus_musculus.GRCm38.93.chr.bed -wb | cut -f1,2,3,7,8,10,11 > Results/RA4_vs_RA4_PG/RA4_diff_sitesAnnoNearest.bed
bedtools closest -D ref -a Results/RA4_vs_RA4_PG/RA4_PG_diff_sites.bed -b Mus_musculus.GRCm38.93.chr.bed -wb | cut -f1,2,3,7,8,10,11 > Results/RA4_vs_RA4_PG/RA4_PG_diff_sitesAnnoNearest.bed
bedtools closest -D ref -a Results/RA4_vs_RA4_PG/Intersection_sites.bed -b Mus_musculus.GRCm38.93.chr.bed -wb | cut -f1,2,3,7,8,10,11 > Results/RA4_vs_RA4_PG/Intersection_sitesAnnoNearest.bed

bedtools closest -D ref -a Results/RA4_vs_RA4_TCP/RA4_diff_sites.bed -b Mus_musculus.GRCm38.93.chr.bed -wb | cut -f1,2,3,7,8,10,11 > Results/RA4_vs_RA4_TCP/RA4_diff_sitesAnnoNearest.bed
bedtools closest -D ref -a Results/RA4_vs_RA4_TCP/RA4_TCP_diff_sites.bed -b Mus_musculus.GRCm38.93.chr.bed -wb | cut -f1,2,3,7,8,10,11 > Results/RA4_vs_RA4_TCP/RA4_TCP_diff_sitesAnnoNearest.bed
bedtools closest -D ref -a Results/RA4_vs_RA4_TCP/Intersection_sites.bed -b Mus_musculus.GRCm38.93.chr.bed -wb | cut -f1,2,3,7,8,10,11 > Results/RA4_vs_RA4_TCP/Intersection_sitesAnnoNearest.bed

bedtools closest -D ref -a Results/RA4_PG_vs_RA4_TCP/RA4_PG_diff_sites.bed -b Mus_musculus.GRCm38.93.chr.bed -wb | cut -f1,2,3,7,8,10,11 > Results/RA4_PG_vs_RA4_TCP/RA4_PG_diff_sitesAnnoNearest.bed
bedtools closest -D ref -a Results/RA4_PG_vs_RA4_TCP/RA4_TCP_diff_sites.bed -b Mus_musculus.GRCm38.93.chr.bed -wb | cut -f1,2,3,7,8,10,11 > Results/RA4_PG_vs_RA4_TCP/RA4_TCP_diff_sitesAnnoNearest.bed
bedtools closest -D ref -a Results/RA4_PG_vs_RA4_TCP/Intersection_sites.bed -b Mus_musculus.GRCm38.93.chr.bed -wb | cut -f1,2,3,7,8,10,11 > Results/RA4_PG_vs_RA4_TCP/Intersection_sitesAnnoNearest.bed

bedtools closest -D ref -a Results/RA_PG_TCP_vs_UD/RA_diff_sites.bed -b Mus_musculus.GRCm38.93.chr.bed -wb | cut -f1,2,3,7,8,10,11 > Results/RA_PG_TCP_vs_UD/RA_diff_sitesAnnoNearest.bed
bedtools closest -D ref -a Results/RA_PG_TCP_vs_UD/UD_diff_sites.bed -b Mus_musculus.GRCm38.93.chr.bed -wb | cut -f1,2,3,7,8,10,11 > Results/RA_PG_TCP_vs_UD/UD_sitesAnnoNearest.bed
bedtools closest -D ref -a Results/RA_PG_TCP_vs_UD/Intersection_sites.bed -b Mus_musculus.GRCm38.93.chr.bed -wb | cut -f1,2,3,7,8,10,11 > Results/RA_PG_TCP_vs_UD/Intersection_sitesAnnoNearest.bed


# ************************ Use bedtools intersect to annotate all sites comparisons overlap ************************
bedtools intersect -a Results/RA4_D_vs_PG_D_vs_TCP_D/RA4_PG_UD_only_sites.bed -b Mus_musculus.GRCm38.93.chr.bed -wb | cut -f1,2,3,7,8,10 > Results/RA4_D_vs_PG_D_vs_TCP_D/RA4_PG_UD_only_sitesAnno.bed
bedtools intersect -a Results/RA4_D_vs_PG_D_vs_TCP_D/RA4_TCP_UD_only_sites.bed -b Mus_musculus.GRCm38.93.chr.bed -wb | cut -f1,2,3,7,8,10 > Results/RA4_D_vs_PG_D_vs_TCP_D/RA4_TCP_UD_only_sitesAnno.bed
bedtools intersect -a Results/RA4_D_vs_PG_D_vs_TCP_D/RA4_UD_only_sites.bed -b Mus_musculus.GRCm38.93.chr.bed -wb | cut -f1,2,3,7,8,10 > Results/RA4_D_vs_PG_D_vs_TCP_D/RA4_UD_only_sitesAnno.bed

bedtools intersect -a Results/RA4_D_vs_PG_D_vs_TCP_D/RA4_PG_vs_TCP_intersection_sites.bed -b Mus_musculus.GRCm38.93.chr.bed -wb | cut -f1,2,3,7,8,10 > Results/RA4_D_vs_PG_D_vs_TCP_D/RA4_PG_vs_TCP_intersection_sitesAnno.bed
bedtools intersect -a Results/RA4_D_vs_PG_D_vs_TCP_D/RA4_UD_vs_PG_intersection_sites.bed -b Mus_musculus.GRCm38.93.chr.bed -wb | cut -f1,2,3,7,8,10 > Results/RA4_D_vs_PG_D_vs_TCP_D/RA4_UD_vs_PG_intersection_sitesAnno.bed
bedtools intersect -a Results/RA4_D_vs_PG_D_vs_TCP_D/RA4_UD_vs_TCP_intersection_sites.bed -b Mus_musculus.GRCm38.93.chr.bed -wb | cut -f1,2,3,7,8,10 > Results/RA4_D_vs_PG_D_vs_TCP_D/RA4_UD_vs_TCP_intersection_sitesAnno.bed
bedtools intersect -a Results/RA4_D_vs_PG_D_vs_TCP_D/RA4_vs_PG_vs_TCP_intersection_sites.bed -b Mus_musculus.GRCm38.93.chr.bed -wb | cut -f1,2,3,7,8,10 > Results/RA4_D_vs_PG_D_vs_TCP_D/RA4_vs_PG_vs_TCP_intersection_sitesAnno.bed

# # ************************ Convert promoters to bed files ************************
awk 'NR>1 {print $2, $3, $4, $5, $6, $7, $8}' promoters.txt | awk 'BEGIN{OFS="\t"}{print $1, $2, $3, $4, $5, $6, $7}' > promoters.bed















