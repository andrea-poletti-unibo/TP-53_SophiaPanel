
# DEFINE THE PROJECT in the "analisi_in_corso" folder

project <- "TP53"

#################################################### BATCH 1 ######################################################## 
library(ASCAT)
directories <- list.dirs(paste0("D:/analisi_in_corso/", project, "/ASCAT"), recursive = F, full.names = F)

i = directories[1]

#_________PART 1: set up dir tree and files_________
#step 0: choose working directory
# setwd(choose.dir()) 
setwd(paste0("D:/analisi_in_corso/", project, "/ASCAT/", i))

#step 1: identify multi-sample logR and BAF files (pre-generated and already present in the working directory) 
file_logR<-list.files(pattern="logR")
file_BAF <-list.files(pattern="BAF")

#step 2: creation of sub-directories for the different GISTIC phases outputs
dir.create("1.rawdata_plots")
dir.create("2.Sep_plots")
dir.create("3.segment_data")
dir.create("4.segmented_plots")
dir.create("5.ASCAT_output_plots")

#_________PART 2: running ASCAT_________
#step 0: loading data
ascat.bc = ascat.loadData(file_logR, file_BAF, 
                          chrs = c(1:22))

#step 1: plotting raw logR and BAF data 
ascat.plotRawData(ascat.bc,
                  img.dir = "1.rawdata_plots")

#step 2: predict genotypes areas from BAF data
gg<-ascat.predictGermlineGenotypes(ascat.bc, 
                                   platform = "AffySNP6",
                                   img.dir = "2.Sep_plots")

#step 3: running ASPCF segmentation on both logR and BAF
ascat.bc = ascat.aspcf(ascat.bc, 
                       ascat.gg=gg, 
                       penalty=70, #default = 25 ( bigger => less segments)
                       out.dir= "3.segment_data") 

#step 4: plotting segmented data
ascat.plotSegmentedData(ascat.bc,
                        img.dir="4.segmented_plots")

#step 5: running ASCAT analysis and plotting outputs
ascat.output = ascat.runAscat(ascat.bc,
                              gamma=0.55, # default= 0.55
                              y_limit=5, # default= 5
                              circos = NA,
                              rho_manual = NA, # rho= aberrant cell fraction
                              psi_manual = NA, # psi= ploidy
                              img.dir = "5.ASCAT_output_plots")

#_________ part 3: exporting data results_________
# #visualize results
# ascat.output$aberrantcellfraction
# ascat.output$ploidy
# ascat.output$psi
# ascat.output$goodnessOfFit
# ascat.output$failedarrays
# ascat.output$nonaberrantarrays

# step 1: export a csv table of results
exportRes<- data.frame(ascat.output$goodnessOfFit,
                       ascat.output$aberrantcellfraction,
                       ascat.output$ploidy,
                       ascat.output$psi)

write.csv(exportRes, paste0( i, "_ASCAT.results_table.csv"))

#step 2: export the lists of failed arrays and non-aberrant arrays
write.csv(ascat.output$failedarrays, paste0( i,"_failed_arrays_list.csv"))
write.csv(ascat.output$nonaberrantarrays, paste0( i,"_nonAberrant_arrays_list.csv"))

#step 3: export a tsv file with the segments for every sample (adjusted-fitted and adjusted-raw)
library(readr)
segmentsAdj<-ascat.output$segments
segmentsAdj$CNvalue<- segmentsAdj$nMajor+segmentsAdj$nMinor
segmentsAdj$length<- segmentsAdj$endpos - segmentsAdj$startpos
write_tsv(segmentsAdj, paste0( i,"_adj_fitted_segments.tsv"))

segmentsAdjRaw<-ascat.output$segments_raw
segmentsAdjRaw$CNvalue<- segmentsAdjRaw$nMajor + segmentsAdjRaw$nMinor
segmentsAdjRaw$CNvalueRaw<- segmentsAdjRaw$nAraw + segmentsAdjRaw$nBraw
segmentsAdjRaw$length<- segmentsAdjRaw$endpos - segmentsAdjRaw$startpos
write_tsv(segmentsAdjRaw, paste0( i,"adj_fitted_and_raw_segments.tsv"))

# export ascat.output as a .Rdata single object - use readRDS() function to restore/load it in a new R session
saveRDS(ascat.output, paste0( i,"ascat_output.RData"))

# my_output<-readRDS("ascat_output.RData")


#################################################### BATCH 2 ######################################################## 
library(ASCAT)
directories <- list.dirs(paste0("D:/analisi_in_corso/", project, "/ASCAT"), recursive = F, full.names = F)

i = directories[2]

#_________PART 1: set up dir tree and files_________
#step 0: choose working directory
# setwd(choose.dir()) 
setwd(paste0("D:/analisi_in_corso/", project, "/ASCAT/", i))

#step 1: identify multi-sample logR and BAF files (pre-generated and already present in the working directory) 
file_logR<-list.files(pattern="logR")
file_BAF <-list.files(pattern="BAF")

#step 2: creation of sub-directories for the different GISTIC phases outputs
dir.create("1.rawdata_plots")
dir.create("2.Sep_plots")
dir.create("3.segment_data")
dir.create("4.segmented_plots")
dir.create("5.ASCAT_output_plots")

#_________PART 2: running ASCAT_________
#step 0: loading data
ascat.bc = ascat.loadData(file_logR, file_BAF, 
                          chrs = c(1:22))

#step 1: plotting raw logR and BAF data 
ascat.plotRawData(ascat.bc,
                  img.dir = "1.rawdata_plots")

#step 2: predict genotypes areas from BAF data
gg<-ascat.predictGermlineGenotypes(ascat.bc, 
                                   platform = "AffyCytoScanHD",
                                   img.dir = "2.Sep_plots")

#step 3: running ASPCF segmentation on both logR and BAF
ascat.bc = ascat.aspcf(ascat.bc, 
                       ascat.gg=gg, 
                       penalty=70, #default = 25 ( bigger => less segments)
                       out.dir= "3.segment_data") 

#step 4: plotting segmented data
ascat.plotSegmentedData(ascat.bc,
                        img.dir="4.segmented_plots")

#step 5: running ASCAT analysis and plotting outputs
ascat.output = ascat.runAscat(ascat.bc,
                              gamma=0.55, # default= 0.55
                              y_limit=5, # default= 5
                              circos = NA,
                              rho_manual = NA, # rho= aberrant cell fraction
                              psi_manual = NA, # psi= ploidy
                              img.dir = "5.ASCAT_output_plots")

#_________ part 3: exporting data results_________
# #visualize results
# ascat.output$aberrantcellfraction
# ascat.output$ploidy
# ascat.output$psi
# ascat.output$goodnessOfFit
# ascat.output$failedarrays
# ascat.output$nonaberrantarrays

# step 1: export a csv table of results
exportRes<- data.frame(ascat.output$goodnessOfFit,
                       ascat.output$aberrantcellfraction,
                       ascat.output$ploidy,
                       ascat.output$psi)

write.csv(exportRes, paste0( i, "_ASCAT.results_table.csv"))

#step 2: export the lists of failed arrays and non-aberrant arrays
write.csv(ascat.output$failedarrays, paste0( i,"_failed_arrays_list.csv"))
write.csv(ascat.output$nonaberrantarrays, paste0( i,"_nonAberrant_arrays_list.csv"))

#step 3: export a tsv file with the segments for every sample (adjusted-fitted and adjusted-raw)
library(readr)
segmentsAdj<-ascat.output$segments
segmentsAdj$CNvalue<- segmentsAdj$nMajor+segmentsAdj$nMinor
segmentsAdj$length<- segmentsAdj$endpos - segmentsAdj$startpos
write_tsv(segmentsAdj, paste0( i,"_adj_fitted_segments.tsv"))

segmentsAdjRaw<-ascat.output$segments_raw
segmentsAdjRaw$CNvalue<- segmentsAdjRaw$nMajor + segmentsAdjRaw$nMinor
segmentsAdjRaw$CNvalueRaw<- segmentsAdjRaw$nAraw + segmentsAdjRaw$nBraw
segmentsAdjRaw$length<- segmentsAdjRaw$endpos - segmentsAdjRaw$startpos
write_tsv(segmentsAdjRaw, paste0( i,"adj_fitted_and_raw_segments.tsv"))

# export ascat.output as a .Rdata single object - use readRDS() function to restore/load it in a new R session
saveRDS(ascat.output, paste0( i,"ascat_output.RData"))

# my_output<-readRDS("ascat_output.RData")


#################################################### BATCH 3 ######################################################## 
library(ASCAT)
directories <- list.dirs(paste0("D:/analisi_in_corso/", project, "/ASCAT"), recursive = F, full.names = F)

i = directories[3]

#_________PART 1: set up dir tree and files_________
#step 0: choose working directory
# setwd(choose.dir()) 
setwd(paste0("D:/analisi_in_corso/", project, "/ASCAT/", i))

#step 1: identify multi-sample logR and BAF files (pre-generated and already present in the working directory) 
file_logR<-list.files(pattern="logR")
file_BAF <-list.files(pattern="BAF")

#step 2: creation of sub-directories for the different GISTIC phases outputs
dir.create("1.rawdata_plots")
dir.create("2.Sep_plots")
dir.create("3.segment_data")
dir.create("4.segmented_plots")
dir.create("5.ASCAT_output_plots")

#_________PART 2: running ASCAT_________
#step 0: loading data
ascat.bc = ascat.loadData(file_logR, file_BAF, 
                          chrs = c(1:22))

#step 1: plotting raw logR and BAF data 
ascat.plotRawData(ascat.bc,
                  img.dir = "1.rawdata_plots")

#step 2: predict genotypes areas from BAF data
gg<-ascat.predictGermlineGenotypes(ascat.bc, 
                                   platform = "AffyCytoScanHD",
                                   img.dir = "2.Sep_plots")

#step 3: running ASPCF segmentation on both logR and BAF
ascat.bc = ascat.aspcf(ascat.bc, 
                       ascat.gg=gg, 
                       penalty=70, #default = 25 ( bigger => less segments)
                       out.dir= "3.segment_data") 

#step 4: plotting segmented data
ascat.plotSegmentedData(ascat.bc,
                        img.dir="4.segmented_plots")

#step 5: running ASCAT analysis and plotting outputs
ascat.output = ascat.runAscat(ascat.bc,
                              gamma=0.55, # default= 0.55
                              y_limit=5, # default= 5
                              circos = NA,
                              rho_manual = NA, # rho= aberrant cell fraction
                              psi_manual = NA, # psi= ploidy
                              img.dir = "5.ASCAT_output_plots")

#_________ part 3: exporting data results_________
# #visualize results
# ascat.output$aberrantcellfraction
# ascat.output$ploidy
# ascat.output$psi
# ascat.output$goodnessOfFit
# ascat.output$failedarrays
# ascat.output$nonaberrantarrays

# step 1: export a csv table of results
exportRes<- data.frame(ascat.output$goodnessOfFit,
                       ascat.output$aberrantcellfraction,
                       ascat.output$ploidy,
                       ascat.output$psi)

write.csv(exportRes, paste0( i, "_ASCAT.results_table.csv"))

#step 2: export the lists of failed arrays and non-aberrant arrays
write.csv(ascat.output$failedarrays, paste0( i,"_failed_arrays_list.csv"))
write.csv(ascat.output$nonaberrantarrays, paste0( i,"_nonAberrant_arrays_list.csv"))

#step 3: export a tsv file with the segments for every sample (adjusted-fitted and adjusted-raw)
library(readr)
segmentsAdj<-ascat.output$segments
segmentsAdj$CNvalue<- segmentsAdj$nMajor+segmentsAdj$nMinor
segmentsAdj$length<- segmentsAdj$endpos - segmentsAdj$startpos
write_tsv(segmentsAdj, paste0( i,"_adj_fitted_segments.tsv"))

segmentsAdjRaw<-ascat.output$segments_raw
segmentsAdjRaw$CNvalue<- segmentsAdjRaw$nMajor + segmentsAdjRaw$nMinor
segmentsAdjRaw$CNvalueRaw<- segmentsAdjRaw$nAraw + segmentsAdjRaw$nBraw
segmentsAdjRaw$length<- segmentsAdjRaw$endpos - segmentsAdjRaw$startpos
write_tsv(segmentsAdjRaw, paste0( i,"adj_fitted_and_raw_segments.tsv"))

# export ascat.output as a .Rdata single object - use readRDS() function to restore/load it in a new R session
saveRDS(ascat.output, paste0( i,"ascat_output.RData"))

# my_output<-readRDS("ascat_output.RData")





