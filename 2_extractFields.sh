module load bioinfo
module load samtools

samtools view ../bamfiles/F9_UD_MAPQ10_coord_sorted.bam | cut -f 3,4,6,10 > tabFiles/F9_UD_MAPQ10_coord_sorted.tsv
samtools view ../bamfiles/F9_D4_MAPQ10_coord_sorted.bam | cut -f 3,4,6,10 > tabFiles/F9_D4_MAPQ10_coord_sorted.tsv
samtools view ../bamfiles/F9_D4_PG_MAPQ10_coord_sorted.bam | cut -f 3,4,6,10 > tabFiles/F9_D4_PG_MAPQ10_coord_sorted.tsv
samtools view ../bamfiles/F9_D4_TCP_MAPQ10_coord_sorted.bam | cut -f 3,4,6,10 > tabFiles/F9_D4_TCP_MAPQ10_coord_sorted.tsv