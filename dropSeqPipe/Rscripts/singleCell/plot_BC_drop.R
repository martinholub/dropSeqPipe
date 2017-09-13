library(ggplot2)
library(reshape2)
library(grid)
library(gridExtra)
samples = snakemake@config$Samples
data = data.frame(matrix(nrow=length(samples), ncol=4))
data[,4] = samples
colnames(data) = c('BC_drop','UMI_drop','Total_reads','Sample')
fastqc_summary = read.table('summary/fastqc.txt', header = TRUE, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
BC_length = snakemake@params$BC_length
UMI_length = snakemake@params$UMI_length
for(i in 1:length(samples)){
  sample = samples[i]
  data[i,3] = as.integer(unique(as.character(fastqc_summary['Total Sequences',grep(sample, colnames(fastqc_summary))])))
  
  BC_drop_data = read.table(file.path('logs',paste0(sample,'_CELL_barcode.txt')), header = TRUE)
  BC_drop_data = BC_drop_data[-1,] # delete 0 error line
  ids = which(snakemake@params$min_num_below_BC < BC_drop_data$num_failed_bases)
  num_reads_BC_dropped = sum(BC_drop_data$num_barcodes[ids])
  data[i,1] = num_reads_BC_dropped
  
  UMI_drop_data = read.table(file.path('logs',paste0(sample,'_UMI_barcode.txt')), header = TRUE)
  UMI_drop_data = UMI_drop_data[-1,]# delete 0 error line
  ids = which(snakemake@params$min_num_below_UMI < UMI_drop_data$num_failed_bases)
  num_reads_UMI_dropped = sum(UMI_drop_data$num_barcodes[ids])
  data[i,2] = num_reads_UMI_dropped
}

data_long = melt(data, 'Sample')
# Keep the order of the barcodes using factor and levels.
p1 = ggplot(data_long, aes(x=Sample, y = value, fill = variable))
p1 = p1 + geom_bar(stat = 'identity')
p1 = p1 + theme(axis.text.x=element_text(angle = 90, hjust = 0))
p1 = p1 + labs(title=paste('BC and UMI drop comparison.\nMin Cell BC quality =',snakemake@params$min_BC_quality,'& Max number bellow =',snakemake@params$min_num_below_BC,'\nMin UMI BC quality =',snakemake@params$min_UMI_quality,'& Max number bellow =',snakemake@params$min_num_below_UMI), x='Samples', y='Number of reads')
p1 = p1 + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) 
p1 = p1 + scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))
data$BC_drop_pct = data$BC_drop/data$Total_reads
data$UMI_drop_pct = data$UMI_drop/data$Total_reads
data$not_dropped = (data$Total_reads-data$BC_drop-data$UMI_drop)/data$Total_reads
data_long_pct = melt(data,'Sample', measure.vars = c('BC_drop_pct', 'UMI_drop_pct'))
p2 = ggplot(data_long_pct, aes(x=Sample, y = value, fill = variable))
p2 = p2 + geom_bar(stat = 'identity')
p2 = p2 + theme(axis.text.x=element_text(angle = 90, hjust = 0))
p2 = p2 + labs(x='Samples', y='Percentage of reads')
p2 = p2 + scale_y_continuous(labels = scales::percent)
p2 = p2 + scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))
# This allows to align the main plots so that we can relate both directly with the label from the bottom one.
gp1 <- ggplotGrob(p1)
gp2 <- ggplotGrob(p2)
pdf(file = snakemake@output[[1]], width = 8, height = 7)
grid::grid.newpage()
grid::grid.draw(rbind(gp1, gp2))
dev.off()