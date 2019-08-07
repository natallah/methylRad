# methylRad
- Usage:
        python3 1_genomeCountSites.py -i genome_ref.fasta

1b_CountTotSites.sh

- Counts the number of sites found in Reference Genome for each pattern

2_extractFields.sh

- Extract necessary fields from bamfiles and save in a .tsv file

2b_indexBamFile.sh

- Index sorted bamfiles

3_conv2CountMat.py

- Convert .tsv file to count matrix. Results are stored in the provided output folder name
- Output:
        (pattern_name)_sitesLoc.tsv - contains sites and counts, regardless of CGAR strings
        (pattern_name)_sitesLocRaw.tsv - contains all reads with insertion and deletion in the CGAR strings.
                                                                         Use to edit the sites accordingly.
- Usage:
        python3 3_conv2CountMat.py -i (name_of_file.tsv) -o (name_of_folder_for_results)
        ex: python3 3_conv2CountMat.py -i UD_MAPQ10_coord_sorted.tsv -o UD_readsCatalogue

3_conv2CountMat.sh

- Script to run 3_conv2CountMat.py on all samples

4_CompareSites.py

- Compare sites found in 3_conv2CountMat.py to the sites found in 1_genomeCountSites.py to find mismatches.
- Prints out number of sites in and out of the reference genome.
- Edit sites according to the CGAR strings insertion, deletions, substitution
- Output:
        (pattern_name)_finalSites.tsv - contains edited sites and counts
- Usage:
        python3 4_compareSites.py -i (folder_name_of_samples_with_tsv_files)
        ex: python3 4_CompareSites.py -i UD_readsCatalogue

4_CompareSites.sh

- Script to run 4_CompareSites.py on all samples

4a_RemoveDupSites.py

- Collect all sites from every samples into one dictionary. This collection still has duplicates between patterns.
        Results saved to SaveData/SampleSites.pkl and (sample_name)_readsCatalogue/(sample_name)_allSites.tsv

- Find duplicates between patterns, combine counts, and save results to (sample_folder)/(sample_name)_allSites_noDups.tsv

- Include RPM and save to (sample_folder)/(sample_name)_allSites_noDups_final.tsv
- Usage:
        python3 4a_RemoveDupSites.py

4b_VizCounts.py

- Count number of overlapping sites in each sample vs the sites found in the reference genome
- Note: still uses unedited sites. Need to redo this to use no duplicated sites and edited sites
- Usage:
        python3 4b_VizCounts.py

4b_runVizCounts.py
- Script to run 4b_VizCounts.py

5_CompareSamples.py

- Find overlapping sites and differences for every pair of samples as well as RA4_TCP + RA4_PG - UD
        Save results to SaveData/SampleSitesCompare.pkl and SaveData/SampleSitesCompareList.pkl
        Bed files are saved in their respective comparisons under Results.
        Output example: Results/UD_vs_RA4
                                        RA4_diff_sites.bed - contains sites that are in RA4 and not in UD
                                        UD_diff_sites.bed  - contains sites that are in UD and not in RA4
                                        Intersection_sites.bed - contains sites that are in both UD and RA4
                                        * naming convention are similar for all comparisons

- Find overlapping sites and differences amongst all RA4_TCP - UD, RA4_PG - UD, and RA4 - UD.
        Creates Venn diagram of all sets which is found in Figures/RA4_UD_diffVenn.pdf
        Save results as a bed file under Results/
        Output example: Results/RA4_D_vs_PG_D_vs_TCP_D/RA4_UD_only - are sites found only in RA4 - UD
                                        Results/RA4_D_vs_PG_D_vs_TCP_D/RA4_UD_vs_TCP_intersection - are sites found in both RA4 - UD and RA4_TCP - UD
                                        Results/RA4_D_vs_PG_D_vs_TCP_D/RA4_PG_UD_only - are sites found in RA4_PG - UD only
                                        Results/RA4_D_vs_PG_D_vs_TCP_D/RA4_PG_vs_TCP_intersection - are sites found in both RA4_PG - UD and RA4_TCP - UD
                                        Results/RA4_D_vs_PG_D_vs_TCP_D/RA4_TCP_UD_only - are sites found in RA4_TCP - UD only
                                        Results/RA4_D_vs_PG_D_vs_TCP_D/RA4_UD_vs_PG_intersection - are sites found in both RA4 - UD and RA4_PG - UD
                                        Results/RA4_D_vs_PG_D_vs_TCP_D/RA4_vs_PG_vs_TCP_intersection - are sites found in all RA4 - UD, RA4_TCP - UD and RA4_PG - UD

Usage:
        Line-by-line from python console.
        * Needs a little cleaning up

5a_VizCompareSamples.py

- Create Venn diagrams for site differences between samples

- Usage:
        python3 5a_VizCompareSamples.py

6_AnnotateWithGenome.sh

- Convert gtf files to bed files
- Convert (sample_name)_allSites.tsv into bed files and change chromosome type from "1" to "chr1"
- Use BEDTools intersect to annotate (sample_comparisons)_allSites.bed to gene using Mus_musculus.GRCm38.93.chr.bed
- Use BEDTools closest to annotate (sample_comparisons)_allSites.bed to nearest gene using Mus_musculus.GRCm38.93.chr.bed
- Use BEDTools intersect to annotate Results/RA4_D_vs_PG_D_vs_TCP_D for sets using Mus_musculus.GRCm38.93.chr.bed
- Convert promoters.txt to bed file

7_AnnotateSites.py

- In progress, annotate sites using pybedtools instead of BEDTools. Maybe useful for creating a package later.

7_AnnotateSites.sh

- Adds 1KB to LSD1occupied_enhancers.mm10.use.bed and save as LSD1occupied_enhancers.mm10.use.edit.bed
- Adds 1KB to ESC_J1.enhancers.use.bed and save as ESC_J1.enhancers.use.edit.bed
- Order of annotations: LSD1, ESC, promotors.
- Annotate to LSD1 sites in RA4 - UD, RA4_PG - UD, RA4_TCP - UD, and sites found only in RA4 - UD when compare to the rest
        (RA4_UD_only_sites.bed). Store reminder of sites in FinalAnnotation/(sample_name_diff)_sites_noLSD1*
- Annotate FinalAnnotation/(sample_name_diff)_sites_noLSD1* to ESC1
- Annotate to promoters
- Annotate (sample_name)_allSites_noDups_final.tsv to LSD1, ESC, promoters.
- Annotate with count data to get foldChanges

8_makePieCharts.R
- Produces pie chart from annotation using ChIPSeeker

9_FoldChange.py
- Compute Log2 fold changes and save results under Results/FoldChanges
- Usage:
        Line-by-line from python console

9a_FoldChange.sub
- Script to submit 9_FoldChange.py

10_AdditionalAnalysis.py
- Find common genes between sites only in RA4 annotated to promoters vs enhancers
- Sites analysis between RA4_PG - UD and RA4_TCP - UD
- Count the number of patterns in each samples
