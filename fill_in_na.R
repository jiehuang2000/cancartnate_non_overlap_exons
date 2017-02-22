# bsub Rscript /proj/dllab/Jie/Catherine/RNA_seq/fill_in_na.R 
 
for (i in 1:9){
system(paste0('bsub perl /proj/dllab/Austin/Scripts/fill_out_uneven_coord.pl /proj/dllab/Jie/Catherine/RNA_seq/chromosome/ch',i,'_WT1.coord        /proj/dllab/Jie/Catherine/RNA_seq/chromosome/ch',i,'_WT1_filledNA.coord'))
system(paste0('bsub perl /proj/dllab/Austin/Scripts/fill_out_uneven_coord.pl /proj/dllab/Jie/Catherine/RNA_seq/chromosome/ch',i,'_WT2.coord        /proj/dllab/Jie/Catherine/RNA_seq/chromosome/ch',i,'_WT2_filledNA.coord'))
system(paste0('bsub perl /proj/dllab/Austin/Scripts/fill_out_uneven_coord.pl /proj/dllab/Jie/Catherine/RNA_seq/chromosome/ch',i,'_Setd2_del1.coord /proj/dllab/Jie/Catherine/RNA_seq/chromosome/ch',i,'_Setd2_del1_filledNA.coord'))
system(paste0('bsub perl /proj/dllab/Austin/Scripts/fill_out_uneven_coord.pl /proj/dllab/Jie/Catherine/RNA_seq/chromosome/ch',i,'_Setd2_del2.coord /proj/dllab/Jie/Catherine/RNA_seq/chromosome/ch',i,'_Setd2_del2_filledNA.coord'))
}