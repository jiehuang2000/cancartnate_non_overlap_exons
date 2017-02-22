# cancartnate_non_overlap_exons

## This is an instruction of cancartnating the RNA-seq read depth data for exons from the same gene, and exclude genes with overlaped exons.

### control group: Wild type cell; 
### treatment group: healthy cell with SETD2 gene knockedout. 


from hg19 bed files of 5’UTRExon, CodingExon, 3’UTRExon, extract all genes with non-overlapped exons, and concatenate all exons

###1. make the bed file(in zinba format) : /proj/dllab/Jie/Catherine/RNA_seq/GENES_WITH_NON_OVERLAP_EXON.bed

includes exons that are shared by only one gene.

first few lines:
NM_000014_utr5_35_0_chr12_9268446_r	chr12	9268445	9268558	- 0
NM_000014_cds_0_0_chr12_9220419_r	chr12	9220418	9220435	- 0
NM_000014_cds_1_0_chr12_9220779_r	chr12	9220778	9220820	- 0
NM_000014_cds_2_0_chr12_9221336_r	chr12	9221335	9221438	- 0
NM_000014_cds_3_0_chr12_9222341_r	chr12	9222340	9222409	- 0


###2. from zinba to coord: 

library(zinba)
coord.sbpc(coordfile="BEDFILE_OF_REGIONS.bed", inputfile="SAMPLE_WIG_FILE.wig", outputfile="OUTPUT.coord", twobitfile="/proj/dllab/Catherine/hg19/hg19.2bit")
 
Make a Rscript file, and run on killdevil using /proj/.test/roach/FAIRE/bin/Rscript yourscript.R
 
mine is: bsub /proj/.test/roach/FAIRE/bin/Rscript /proj/dllab/Jie/Catherine/RNA_seq/zinba_to_coord.R



###3. each file is too large, so break down into chromosomes (repeat for each chr):


grep 'chr21' /proj/dllab/Jie/Catherine/RNA_seq/WT1.coord        > /proj/dllab/Jie/Catherine/RNA_seq/chromosome/ch21_WT1.coord
grep 'chr21' /proj/dllab/Jie/Catherine/RNA_seq/WT2.coord        > /proj/dllab/Jie/Catherine/RNA_seq/chromosome/ch21_WT2.coord
grep 'chr21' /proj/dllab/Jie/Catherine/RNA_seq/Setd2_del1.coord > /proj/dllab/Jie/Catherine/RNA_seq/chromosome/ch21_Setd2_del1.coord
grep 'chr21' /proj/dllab/Jie/Catherine/RNA_seq/Setd2_del2.coord > /proj/dllab/Jie/Catherine/RNA_seq/chromosome/ch21_Setd2_del2.coord


###4. fill in NA

Once done, this needs to be filled in with NA values to make a dataframe. Austin has a script to do this. (can use this to do multiple chromosomes: /proj/dllab/Jie/Catherine/RNA_seq/fill_in_na.R)


bsub perl /proj/dllab/Austin/Scripts/fill_out_uneven_coord.pl /proj/dllab/Jie/Catherine/RNA_seq/chromosome/ch21_WT1.coord        /proj/dllab/Jie/Catherine/RNA_seq/chromosome/ch21_WT1_filledNA.coord
bsub perl /proj/dllab/Austin/Scripts/fill_out_uneven_coord.pl /proj/dllab/Jie/Catherine/RNA_seq/chromosome/ch21_WT2.coord        /proj/dllab/Jie/Catherine/RNA_seq/chromosome/ch21_WT2_filledNA.coord
bsub perl /proj/dllab/Austin/Scripts/fill_out_uneven_coord.pl /proj/dllab/Jie/Catherine/RNA_seq/chromosome/ch21_Setd2_del1.coord /proj/dllab/Jie/Catherine/RNA_seq/chromosome/ch21_Setd2_del1_filledNA.coord
bsub perl /proj/dllab/Austin/Scripts/fill_out_uneven_coord.pl /proj/dllab/Jie/Catherine/RNA_seq/chromosome/ch21_Setd2_del2.coord /proj/dllab/Jie/Catherine/RNA_seq/chromosome/ch21_Setd2_del2_filledNA.coord


###5. cancatenate all exons

/proj/dllab/Jie/Catherine/RNA_seq/merge_exon_coord_LSF.R

raw read depth plots were generated in the folder: proj/dllab/Jie/Catherine/RNA_seq/chromosome/plots
the read depth matrix is in the folder: proj/dllab/Jie/Catherine/RNA_seq/chromosome/gene_matrix
