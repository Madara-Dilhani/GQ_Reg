########################################################################################################
# This function predicts G quadruplex motif(s) in 'x' (nucleotide sequence(s)). 
# Nucleotide sequence should be provided in fasta format.
########################################################################################################
gquad <- function(seq, seq_name, pattern){
  GQ = NULL
  GQ_temp <- gregexpr(pattern, seq, ignore.case = TRUE)
  sequence_position <- GQ_temp[[1]][1:length(GQ_temp[[1]])]
  
  if (sequence_position > 0){
  GQ_temp_pos <- regmatches(seq, GQ_temp)
  sequence <- GQ_temp_pos[[1]][1:length(GQ_temp_pos[[1]])]
  sequence_length <- nchar(sequence)
  
  df_temp <- c(seq_name, sequence_position, sequence, sequence_length)
  GQ <- rbind(GQ,df_temp)
  colnames(GQ) <- c("Sequence_ID", "GQ_StartPos", "Matched_GQ_sequence", "sequence_length")
  return(GQ)
  } else {return(GQ)}
}
########################################################################################################