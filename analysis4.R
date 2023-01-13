size = read.table('errc_transcriptome_Bo.E.5H_all_rnaseq_combined_pass.trimmed_stranded.fastq.gmap_with_alignments.txt', header = T, sep = '\t')

hist(size$align_identity, breaks = 50, xlab = 'Percentage alignment identity')
summary(size$align_identity)
plot(size$align_identity,size$read_identity)
plot(size$align_identity,size$read_length)
