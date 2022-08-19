###############################
#A) Download relevant packages

library(BiocManager)
library(GenomicRanges)
library(Rsamtools)
library(BSgenome)
library(Biostrings)
library(Bios2cor)
library(rtracklayer)
library(Gviz)
library(stats)

###############################
#B) Define & import BED # 

single_array_BED <- file.path("DIRECTORY", "FOLDER", "BED FILE.bed")
peaks <- import.bed(single_array_BED)

###############################
#C) Import bigwig & define histone parameters #

options(ucscChromosomeNames=FALSE)
#Always include if you're not using genome/chromosomes from UCSC database

#1  Dmel-MSL2 to Dmel HL
#2  Dvir-MSL2 ChIP-seq to Dmel HL 
#3  Dvir-MSL2 DIP-seq to Dmel HL 
#4  Non-OE MSL2 from S2 
#5  Non-OE MSL2 from salivary glands (Figueiredo 2014) 

#1
dmel_msl2_to_dmelHL_rep1_bw <- import.bw("DIRECTORY/FOLDER/BIGWIG FILE 1 REPLICATE 1.bigwig", as = "GRanges")
dmel_msl2_to_dmelHL_rep2_bw <- import.bw("DIRECTORY/FOLDER/BIGWIG FILE 1 REPLICATE 2.bigwig", as = "GRanges")

#2
dvir_msl2_chip_to_dmelHL_rep1_bw <- import.bw("DIRECTORY/FOLDER/BIGWIG FILE 2 REPLICATE 1.bigwig", as = "GRanges")
dvir_msl2_chip_to_dmelHL_rep2_bw <- import.bw("DIRECTORY/FOLDER/BIGWIG FILE 2 REPLICATE 2.bigwig", as = "GRanges")

#3
dvir_msl2_dip_to_dmelHL_rep1_bw <- import.bw("DIRECTORY/FOLDER/BIGWIG FILE 3 REPLICATE 1.bigwig", as = "GRanges")
dvir_msl2_dip_to_dmelHL_rep2_bw <- import.bw("DIRECTORY/FOLDER/BIGWIG FILE 3 REPLICATE 2.bigwig", as = "GRanges")

#4
nonOE_msl2_s2_rep1_bw <- import.bw("DIRECTORY/FOLDER/BIGWIG FILE 4 REPLICATE 1.bigwig", as = "GRanges")
nonOE_msl2_s2_rep2_bw <- import.bw("DIRECTORY/FOLDER/BIGWIG FILE 4 REPLICATE 2.bigwig", as = "GRanges")
nonOE_msl2_s2_rep3_bw <- import.bw("DIRECTORY/FOLDER/BIGWIG FILE 4 REPLICATE 3.bigwig", as = "GRanges")
nonOE_msl2_s2_rep4_bw <- import.bw("DIRECTORY/FOLDER/BIGWIG FILE 4 REPLICATE 4.bigwig", as = "GRanges")

#5
nonOE_msl2_salivary_glands_bw <- import.bw("DIRECTORY/FOLDER/BIGWIG FILE 5 REPLICATE 1.bigwig", as = "GRanges")

track_bed <- AnnotationTrack(peaks, id = c("H2B","H4","H2A","H3","H1"), cex = 0.5 ,fill = "darkblue", col = "darkblue", name = '')
#Histone track parameters

############################################################# 
#D) Create single GRanges BIGWIG #

#In order to produce ranges & avg, all replicates (BIGWIGs) from a single experiment must be generated into one GRanges object.
seqi1 <- seqinfo(dmel_msl2_to_dmelHL_rep1_bw)
seqi2 <- seqinfo(dmel_msl2_to_dmelHL_rep2_bw)

seqi3 <- seqinfo(dvir_msl2_chip_to_dmelHL_rep1_bw)
seqi4 <- seqinfo(dvir_msl2_chip_to_dmelHL_rep2_bw)

seqi5 <- seqinfo(dvir_msl2_dip_to_dmelHL_rep1_bw)
seqi6 <- seqinfo(dvir_msl2_dip_to_dmelHL_rep2_bw)

seqi7 <- seqinfo(nonOE_msl2_s2_rep1_bw)
seqi8 <- seqinfo(nonOE_msl2_s2_rep2_bw)
seqi9 <- seqinfo(nonOE_msl2_s2_rep3_bw)
seqi10 <- seqinfo(nonOE_msl2_s2_rep4_bw)

seqi11 <- seqinfo(nonOE_msl2_salivary_glands_bw)

#Alter smoothness of plot by tweaking tilewidth.
combined_bigwig_1 <- tileGenome(seqi1, tilewidth = 30, cut.last.tile.in.chrom = TRUE)
overlap1 <- findOverlaps(combined_bigwig_1,dmel_msl2_to_dmelHL_rep1_bw)
overlap2 <- findOverlaps(combined_bigwig_1,dmel_msl2_to_dmelHL_rep2_bw)


combined_bigwig_2 <- tileGenome(seqi3, tilewidth = 30, cut.last.tile.in.chrom = TRUE)
overlap3 <- findOverlaps(combined_bigwig_2,dvir_msl2_chip_to_dmelHL_rep1_bw)
overlap4 <- findOverlaps(combined_bigwig_2,dvir_msl2_chip_to_dmelHL_rep2_bw)

combined_bigwig_3 <- tileGenome(seqi5, tilewidth = 30, cut.last.tile.in.chrom = TRUE)
overlap5 <- findOverlaps(combined_bigwig_3,dvir_msl2_dip_to_dmelHL_rep1_bw)
overlap6 <- findOverlaps(combined_bigwig_3,dvir_msl2_dip_to_dmelHL_rep2_bw)

combined_bigwig_4 <- tileGenome(seqi7, tilewidth = 30, cut.last.tile.in.chrom = TRUE)
overlap7 <- findOverlaps(combined_bigwig_4,nonOE_msl2_s2_rep1_bw)
overlap8 <- findOverlaps(combined_bigwig_4,nonOE_msl2_s2_rep2_bw)
overlap9 <- findOverlaps(combined_bigwig_4,nonOE_msl2_s2_rep3_bw)
overlap10 <- findOverlaps(combined_bigwig_4,nonOE_msl2_s2_rep4_bw)

combined_bigwig_5 <- tileGenome(seqi11, tilewidth = 30, cut.last.tile.in.chrom = TRUE)
overlap11 <- findOverlaps(combined_bigwig_5,nonOE_msl2_salivary_glands_bw)


combined_bigwig_1$score_1 = NA_real_
combined_bigwig_1$score_2 = NA_real_
combined_bigwig_1[queryHits(overlap1)]$score_1 = dmel_msl2_to_dmelHL_rep1_bw[subjectHits(overlap1)]$score
combined_bigwig_1[queryHits(overlap2)]$score_2 = dmel_msl2_to_dmelHL_rep2_bw[subjectHits(overlap2)]$score

combined_bigwig_2$score_1 = NA_real_
combined_bigwig_2$score_2 = NA_real_
combined_bigwig_2[queryHits(overlap3)]$score_1 = dvir_msl2_chip_to_dmelHL_rep1_bw[subjectHits(overlap3)]$score
combined_bigwig_2[queryHits(overlap4)]$score_2 = dvir_msl2_chip_to_dmelHL_rep2_bw[subjectHits(overlap4)]$score

combined_bigwig_3$score_1 = NA_real_
combined_bigwig_3$score_2 = NA_real_
combined_bigwig_3[queryHits(overlap5)]$score_1 = dvir_msl2_dip_to_dmelHL_rep1_bw[subjectHits(overlap5)]$score
combined_bigwig_3[queryHits(overlap6)]$score_2 = dvir_msl2_dip_to_dmelHL_rep2_bw[subjectHits(overlap6)]$score

combined_bigwig_4$score_1 = NA_real_
combined_bigwig_4$score_2 = NA_real_
combined_bigwig_4$score_3 = NA_real_
combined_bigwig_4$score_4 = NA_real_
combined_bigwig_4[queryHits(overlap7)]$score_1 = nonOE_msl2_s2_rep1_bw[subjectHits(overlap7)]$score
combined_bigwig_4[queryHits(overlap8)]$score_2 = nonOE_msl2_s2_rep2_bw[subjectHits(overlap8)]$score
combined_bigwig_4[queryHits(overlap9)]$score_3 = nonOE_msl2_s2_rep3_bw[subjectHits(overlap9)]$score
combined_bigwig_4[queryHits(overlap10)]$score_4 = nonOE_msl2_s2_rep4_bw[subjectHits(overlap10)]$score

combined_bigwig_5$score_1 = NA_real_
combined_bigwig_5[queryHits(overlap11)]$score_1 = nonOE_msl2_salivary_glands_bw[subjectHits(overlap11)]$score

#############################################################
#E) Plot #

track_1 <- DataTrack(range = combined_bigwig_1,  name = "Dmel-MSL2 to Dmel HL", col = "purple",col.name = "black", type = c("a","confint"), baseline = 0, col.baseline= "black", lty.baseline = "dashed", lwd.baseline = 1, fontsize = 15)
track_2 <- DataTrack(range = combined_bigwig_2, name = "Dvir-MSL2 ChIP-seq to Dmel HL", col = "purple",  col.name = "black", type = c("a","confint"), baseline = 0, col.baseline= "black", lty.baseline = "dashed", lwd.baseline = 1)
plotTracks(c(GenomeAxisTrack(), track_1, track_2, track_bed), yTicksAt = c(0,2, 4,6,8,10,12), background.title="transparent",cex = 0.7,cex.axis=0.5,cex.main=2, sizes = c(1,4,4,1),cex.title = 0.5, featureAnnotation = "id", ylim = c(0, 13),col.title = "black",  type = c("a","confint"), from = 1, to = 5048, fontcolor = "black", col.axis = "black", fontsize = 15, margin = 40, innerMargin = 0, baseline = 0, col.baseline= "black", lty.baseline = "dashed", lwd.baseline = 1)

track_3 <- DataTrack(range = combined_bigwig_3, name = "Dvir-MSL2 DIP-seq to Dmel HL ",col = "purple", col.name = "black", type = c("a","confint"), baseline = 0, col.baseline= "black", lty.baseline = "dashed", lwd.baseline = 1)
plotTracks(c(GenomeAxisTrack(), track_3, track_bed), sizes = c(1,8,1), yTicksAt = c(0,2, 4,6,8,10,12), background.title="transparent", ylim = c(-1, 13), featureAnnotation = "id",cex = 0.7,cex.axis=0.5,cex.main=2, cex.title = 0.6 ,col.title = "black",  type = c("a","confint"), from = 1, to = 5048, fontcolor = "black", col.axis = "black", fontsize = 15, margin = 40, innerMargin = 0, baseline = 0, col.baseline= "black", lty.baseline = "dashed", lwd.baseline = 1)

track_4 <- DataTrack(range = combined_bigwig_4,  name = "Non-OE MSL2 from S2", col = "purple",  col.name = "black", type = c("a","confint"), baseline = 0, col.baseline= "black", lty.baseline = "dashed", lwd.baseline = 1, fontsize = 15)
track_5 <- DataTrack(range = combined_bigwig_5,  name = "Non-OE MSL2 from salivary glands (Figueiredo 2014)", col = "purple",  col.name = "black", type = c("a","confint"), baseline = 0, col.baseline= "black", lty.baseline = "dashed", lwd.baseline = 1, fontsize = 15)
plotTracks(c(GenomeAxisTrack(), track_4, track_5, track_bed), ylim = c(-1.5,1),sizes = c(1,4,4,1), background.title="transparent", featureAnnotation = "id",cex = 0.7,cex.axis=0.5,cex.main=2, cex.title = 0.6 ,col.title = "black",  type = c("a","confint"), from = 1, to = 5048, fontcolor = "black", col.axis = "black", fontsize = 15, margin = 40, innerMargin = 0, baseline = 0, col.baseline= "black", lty.baseline = "dashed", lwd.baseline = 1)

###############################################################
#sizes => ratio of heights between GenomeAxisTrack()(maps bases), track_n (data track), and track_bed (histone track) 
#type => presents the average and confidence intervals of the data
#baseline => refers to dashed horizontal line 
