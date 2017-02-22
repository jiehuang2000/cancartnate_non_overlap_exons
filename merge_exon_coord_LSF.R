# bsub Rscript /proj/dllab/Jie/Catherine/RNA_seq/merge_exon_coord_LSF.R 20 21

args = commandArgs(trailingOnly = TRUE)

chr_start = as.numeric(as.character(args[1]))
chr_end = as.numeric(as.character(args[2]))

print(chr_start)
print(chr_end)


# for (i in chr_start:chr_end){
# system(paste0('bsub perl /proj/dllab/Austin/Scripts/fill_out_uneven_coord.pl /proj/dllab/Jie/Catherine/RNA_seq/chromosome/ch',i,'_WT1.coord        /proj/dllab/Jie/Catherine/RNA_seq/chromosome/ch',i,'_WT1_filledNA.coord'))
# system(paste0('bsub perl /proj/dllab/Austin/Scripts/fill_out_uneven_coord.pl /proj/dllab/Jie/Catherine/RNA_seq/chromosome/ch',i,'_WT2.coord        /proj/dllab/Jie/Catherine/RNA_seq/chromosome/ch',i,'_WT2_filledNA.coord'))
# system(paste0('bsub perl /proj/dllab/Austin/Scripts/fill_out_uneven_coord.pl /proj/dllab/Jie/Catherine/RNA_seq/chromosome/ch',i,'_Setd2_del1.coord /proj/dllab/Jie/Catherine/RNA_seq/chromosome/ch',i,'_Setd2_del1_filledNA.coord'))
# system(paste0('bsub perl /proj/dllab/Austin/Scripts/fill_out_uneven_coord.pl /proj/dllab/Jie/Catherine/RNA_seq/chromosome/ch',i,'_Setd2_del2.coord /proj/dllab/Jie/Catherine/RNA_seq/chromosome/ch',i,'_Setd2_del2_filledNA.coord'))
# }







bed = read.table('/proj/dllab/Jie/Catherine/RNA_seq/GENES_WITH_NON_OVERLAP_EXON.bed', header=F, sep='\t')
library(dplyr)

cat_exon = function(data){
        # this function will take the coord data for each exon (sorted by position in gene), and cancatenate all exons in each gene
		# it return two lists:
		# the "complete_gene" list is the gene id and the read depth for each bp of all exons in the same gene
		# the "complete_gene_exon_info" list gives the length of each exon for each gene.
		colnames(bed)=c('NAME', colnames(data)[2:6])
		bed= mutate(bed, coord=paste0(CHROM,':', START,'-',STOP))
		#head(bed)

		data2= mutate(data, coord=paste0(CHROM,':', START,'-',STOP))


		combined = merge(bed, data2, by = "coord", all.x=FALSE, all.y=TRUE, sort=FALSE )
		remove(data);remove(data2);remove(bed)
		#combined[1:10,1:9]
		stopifnot(sum(is.na(combined$START.x))==0)  # check everything in data2 file matches with the bed file

		gene_id=c()
		for (i in 1:nrow(combined)) {
		gene_id = c(gene_id,  unlist(strsplit(as.character(combined[i,2]), "_"))[2])
		 }

		combined = cbind(gene_id, combined)
		combined[1:10,1:20]


		gene_start = which(colnames(combined)=="Position1")
		L = ncol(combined)



		#cancatenate all exons

		complete_gene = vector(mode="list", length = length(unique(gene_id)))
		names(complete_gene) <- unique(gene_id)
		#head(complete_gene)

		complete_gene_exon_info = vector(mode="list", length = length(unique(gene_id)))
		names(complete_gene_exon_info) <- unique(gene_id)

		i=1
		id = combined$gene_id[i]
		seq = combined[i,gene_start:L]; seq = seq[!is.na(seq)]
		len = length(seq)

		for (i in 2:nrow(combined)){
		  #print(i)
		  seq_temp = combined[i,gene_start:L]; seq_temp = seq_temp[!is.na(seq_temp)]
		  len_temp = length(seq_temp)
		  if (combined$gene_id[i] == combined$gene_id[i-1]) {
			seq = c(seq, seq_temp)
			len = c(len, len_temp)
		  } else { 
			pointer = which( names(complete_gene) == combined$gene_id[i-1])
			complete_gene[[pointer]] = seq
			complete_gene_exon_info[[pointer]] = len
			#print(names(complete_gene)[pointer])
			seq = seq_temp
			len = len_temp
		  }
		  if (i==nrow(combined)) {
			pointer = which( names(complete_gene) == combined$gene_id[i])
			complete_gene[[pointer]] = seq
			#print(names(complete_gene)[pointer])
		  }
		}

		result = list(complete_gene=complete_gene, complete_gene_exon_info = complete_gene_exon_info)
		return(result)
} #end of cat_exon function




for (chr_num in chr_start:chr_end){

		print('now processing chromosome:')
		print(chr_num)

		wt1 = read.table(paste0('/proj/dllab/Jie/Catherine/RNA_seq/chromosome/ch',chr_num,'_WT1_filledNA.coord'), sep='\t', header=T)
		dim(wt1)
		result_temp = cat_exon(wt1) # time consuming
		complete_gene.wt1 = result_temp$complete_gene
		complete_gene_exon_info.wt1 = result_temp$complete_gene_exon_info

		wt2 = read.table(paste0('/proj/dllab/Jie/Catherine/RNA_seq/chromosome/ch',chr_num,'_WT2_filledNA.coord'), sep='\t', header=T)
		dim(wt2)
		result_temp = cat_exon(wt2) # time consuming
		complete_gene.wt2 = result_temp$complete_gene
		complete_gene_exon_info.wt2 = result_temp$complete_gene_exon_info

		mut1 = read.table(paste0('/proj/dllab/Jie/Catherine/RNA_seq/chromosome/ch',chr_num,'_Setd2_del1_filledNA.coord'), sep='\t', header=T)
		dim(mut1)
		result_temp = cat_exon(mut1) # time consuming
		complete_gene.mut1 = result_temp$complete_gene
		complete_gene_exon_info.mut1 = result_temp$complete_gene_exon_info

		mut2 = read.table(paste0('/proj/dllab/Jie/Catherine/RNA_seq/chromosome/ch',chr_num,'_Setd2_del2_filledNA.coord'), sep='\t', header=T)
		dim(mut2)
		result_temp = cat_exon(mut2) # time consuming
		complete_gene.mut2 = result_temp$complete_gene
		complete_gene_exon_info.mut2 = result_temp$complete_gene_exon_info




		for (i in 1: length(complete_gene)){
			  max_val = max(complete_gene.mut1[[i]], complete_gene.mut2[[i]], complete_gene.wt1[[i]], complete_gene.wt2[[i]])
			  if (max_val <5) next
			  png(paste0('/proj/dllab/Jie/Catherine/RNA_seq/chromosome/plots/Chr',chr_num,'_',i,'_NM_',names(complete_gene)[i],'.png'), width=8*200, height=4*200, res=200, pointsize = 8)
			  plot(complete_gene.mut1[[i]], type='l', col='red', ylim=c(0, max_val+10),
				   xlab='Gene length (bp)', ylab='Read Depth',
				   main=paste0('NM_',names(complete_gene)[i]))
			  lines(complete_gene.mut2[[i]], type='l', col='red')
			  lines(complete_gene.wt1[[i]], type='l')
			  lines(complete_gene.wt2[[i]], type='l')
			  abline(v=cumsum(complete_gene_exon_info.wt1[[i]]), col='green4')
			  abline(v=cumsum(complete_gene_exon_info.mut1[[i]]), col='yellow')
			  legend('topright', c('WT', 'Setd2 del', 'Exon joint'), col=c('black', 'red','yellow'), lty=1)
			  dev.off()
		}


		#max(sapply(complete_gene_exon_info, sum)) # longest gene length


		# put the result into a matrix 
		list2matrix = function(complete_gene, complete_gene_exon_info, note){
				x=do.call(rbind.data.frame, complete_gene)
				y=sapply(complete_gene_exon_info, sum)
				x = cbind(y,x)
				colnames(x)=NULL
				colnames(x)[1] = "Gene_Length"
				#x[1:10,1:20]
				write.table(x, paste0('/proj/dllab/Jie/Catherine/RNA_seq/chromosome/gene_matrix/Chr',chr_num,'_complete_gene_matrix_',note,'.txt'), sep='\t', col.names = F, row.names = T, quote = F)
		}

		list2matrix(complete_gene.wt1, complete_gene_exon_info.wt1, 'wt1')
		list2matrix(complete_gene.wt2, complete_gene_exon_info.wt2, 'wt2')
		list2matrix(complete_gene.mut1, complete_gene_exon_info.mut1, 'mut1')
		list2matrix(complete_gene.mut2, complete_gene_exon_info.mut2, 'mut2')


}