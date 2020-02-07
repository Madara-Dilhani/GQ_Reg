########################################################################################################
## title: "Find G-Quadraplex structures in given set of sequences"
## author: "Madara Hetti-Arachchilage; mhettiar 'at' illinois.edu"
## date: "August 3, 2017"
########################################################################################################

##### 1. Select given list of sequences from fasta files and write into relevant file

require(seqinr)
unique_seq <- read.csv("Input/1/Most_diverged_unique_seq_pairs.csv", header = TRUE, sep = ",") 
filenames <- list.files("Input/1", pattern="*.fas", full.names=TRUE)
ldf <- lapply(filenames, read.fasta)

for (i in 1: length(ldf)) {
  for(p in length(ldf[[i]]): 1){
    if (is.element(strsplit(attributes(ldf[[i]][p])$name, "[.]")[[1]][1], unique_seq[, 2]) ==FALSE){
      ldf[[i]][p] <- NULL
    }
  }
}

# specifing output file name
filename_out <- gsub("/Input", "/Output", filenames)
filename_out <- gsub(".fas", "-unique.fas", filenames)

# write new fasta file 
for(i in 1: length(filename_out)){
  write.fasta(sequences = ldf[[i]], names = names(ldf[[i]]), 
              file.out = filename_out[i], 
              open = "w")
}
##########################################################################################################

##### 2. Convert accession IDs from one format to another; 
##### eg.NCBI refseq_mrna ids to ensembl transcript ids/ensembl gene ids

BiocManager::install("biomaRt")
library("biomaRt")
ensembl<-  useMart("ensembl", dataset="hsapiens_gene_ensembl")

val <- read.csv("Input/2/NCBI_ref_seq.csv", header = TRUE)
values<- val[,1]
seq <- getBM(attributes=c("refseq_mrna", "ensembl_transcript_id", "ensembl_gene_id", "hgnc_symbol"), 
             filters = "refseq_mrna", values = values, mart= ensembl)
write.csv(x = seq, file = "Output/2/NCBI_ref_seq_converted_to_ensemble.csv")
############################################################################################################

##### This example input uniprot_swissprot accession ids and 
##### output protein sequences & coding sequences with ensemble gene ids, ensembl transcript ids

library(biomaRt)
df = data.frame(ensembl_transcript_id = numeric(), ensembl_gene_id = numeric(),  
                uniprot_swissprot = numeric(), peptide = character(), coding = character(),	
                Peptide_Pattern = character(), Seq_Num_by_Pattern= numeric())
df2 = data.frame(element= numeric(), start= numeric(), End= numeric())

# listing folders 
filenames <- list.files("/Input/3", full.names=TRUE, include.dirs = TRUE)

# loop over each folder containg data for each species
for (i in 1: length(filenames)){
  
txtfiles <- list.files(filenames[i], full.names=TRUE)
for (c in 1: length(txtfiles)){
  data_PIR_raw = read.table(file = txtfiles[c], sep = "\t", quote = "", skip = 3)
  if (i == 1){
    ensembl = useMart(dataset= "dmelanogaster_gene_ensembl", biomart = "ensembl")
  } else if (i == 2){
    ensembl = useMart(dataset= "hsapiens_gene_ensembl", biomart = "ensembl")
  }
  
  # Get sequence attributes for given accessions as in first column for this example
  seq_AC = getBM(attributes = c("ensembl_gene_id","ensembl_transcript_id", "uniprotswissprot"),
                 filters = "uniprotswissprot", values = as.character(data_PIR_raw[,1]), mart = ensembl)
  
  #get protein sequence, coding sequence for matching transcripts from ensemble
  seq_peptide = getSequence(id=seq_AC[,2],  type="ensembl_transcript_id", 
                            seqType="peptide", mart=ensembl)
  
  seq_coding = getSequence(id=seq_AC[,2],  type="ensembl_transcript_id",
                           seqType="coding", mart=ensembl)
  
  
  #combine protein and dna sequence before writting
  merge_Ac_Pep = merge(seq_AC,seq_peptide, by="ensembl_transcript_id")
  merge_Pep_CD = merge(merge_Ac_Pep, seq_coding, by="ensembl_transcript_id") 
  
  #read peptide match tool file to extract seq string
  con=file(filename,open="r")
  line=readLines(con) 
  pattern = unlist(strsplit(line[1], split = "\\s"))
  close(con)
  
  string <- data.frame(matrix(nrow=dim(merge_Pep_CD)[1], ncol=1))
  string_rep <- rep( pattern[3], each=dim(merge_Pep_CD)[1])
  string[,1] <- string_rep
  colnames(string) <- "Peptide_Pattern"
  
  #Writting final data to a excel file
  number_rep <- 1: dim(merge_Pep_CD)[1]
  final_data <- cbind(merge_Pep_CD, string,number_rep)
  
  df <- rbind(df,final_data)
  sy <- easyGregexpr(pattern[3], final_data[,4])
  df2 <- rbind(df2,sy)
  
}
write.csv(x = df2, file= paste(gsub(filenames[i], "/Input", "/Output"),"/ensemble_sequences_downloads.csv"), sep="")
}

