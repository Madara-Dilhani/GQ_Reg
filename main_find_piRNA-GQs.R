########################################################################################################
## title: "Find G-Quadraplex structures in given set of sequences"
## author: "Madara Hetti-Arachchilage; mhettiar 'at' illinois.edu"
## date: "August 3, 2017"
########################################################################################################

##### Importing libraries and custom scripts 
# Loading R source code which has defined custom functions needed for later steps
# Loading R/Biconductor packages (checks if the requested package is already installed in 
# the environment. If not, the package is installed from the appropriate source and then 
# loaded into the workspace
# load packages
require(seqinr)
library(Biostrings)
########################################################################################################

seq_piRNA <- read.fasta(file = "Input/3/C.elegans/piR_cel_v1.0.fa", seqtype = "DNA", forceDNAtolower = TRUE)
pattern <- "G{2,3}?[A|C|G|T|U|N]{1,7}?G{2,3}?[A|C|G|T|U|N]{1,7}?G{2,3}?[A|C|G|T|U|N]{1,7}?G{2,3}?"


gq_pos <- data.frame("Sequence_ID" = character(0), 
                     "GQ_StartPos" = integer(0), 
                     "Matched_GQ_sequence" = character(0), 
                     "sequence_length" = character(0))

gq_count <- data.frame(matrix(NA, nrow=1, ncol=length(seq_piRNA)))
colnames(gq_count) <- names(seq_piRNA)

# Get GQ Positions per each sequence
for (x in 1: nrow(seq_piRNA)){
  
  seq = paste(seq_piRNA[x][[1]], collapse="")
  gq_table <- gquad(seq, names(seq_piRNA)[x], pattern)
  gq_pos <- rbind(gq_pos, gq_table)
}

output = "./Output/3/piR_cel_v1.0_GQ_Pos.csv"
write.csv(x = GQ_Pos, quote = F, row.names = T)
##################################################################
