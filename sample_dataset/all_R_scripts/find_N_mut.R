## ---------------------------------------------------------------------------------------

# Script: Mutational Analysis Pipeline for SARS-CoV-2 Nucleocapsid Protein
# Purpose: To detect amino acid mutations in the Nucleocapsid region
# Authors: Tamara Alrawahneh (CEPHR, UCD), Gabo Gonzalez (NVRL, UCD)

# Last edits: 2022/10/27
# Comments: nuc_to_aa_matrix takes into account deletions by translating codon by codon


## ---------------------------------------------------------------------------------------
# Set working directory ------------------------------------------------------------------


# Load packages to be used ---------------------------------------------------------------
if(!require(seqinr)){
      install.packages("seqinr")
      library(seqinr)
}

if(!require(ape)){
      install.packages("ape")
      library(ape)
}


# Load functions -------------------------------------------------------------------------
# Function to read sequences in batches
# read_n_seqs reads n fasta format sequences, pattern indicating line is ">"
# use read_n_Seqs to read sequence fasta file
read_n_seqs <- function(reading_control, n = 100, output_file, previous_line = NULL){ 
      writing_control <- file(output_file, "w")
      seqs_read <-  0 
      my_line <-  ""
      
      if(!is.null(previous_line)){ #if previous line is not missing, write it to the control
            write(previous_line, file = writing_control)
            seqs_read <-  1 #count this as how many sequences read
      }
      
      while(seqs_read <= n && !is.null(my_line)){ #while loop to read lines one by one
            my_line <-  readLines(con = reading_control, n = 1)
            
            if(!is.null(my_line) && length(grep(pattern = ">", my_line)) > 0){ 
                  seqs_read <- seqs_read + 1 #update the sequences read by 1, while will keep reading until n = 100
            }
            
            if(!is.null(my_line) && seqs_read <= n){ #write the lines to the writing control
                  write(my_line, file = writing_control)
            }
            #THIS WORKS
            if(length(my_line) == 0){
                  print(paste0("reached end of sequence FASTA file"))
                  my_line = NULL
                  break
            }
      }
      close(writing_control)
      return(my_line)
}

# Write chosen sequences to a file with the refseq
write_with_ref <- function(filename, seqs.list, refseq.file) {
      refseq_string <- read.fasta(file =refseq.file, as.string = TRUE)
      con <- file(filename, "w") 
      write(paste0(">", names(refseq_string)[1]), file= con)
      write(refseq_string[[1]], file= con)
      for(i.seq in 1:length(seqs.list)) {
            write(paste0(">", names(seqs.list)[i.seq]), file= con)
            write(paste(seqs.list[[i.seq]], collapse = ""), file= con)
      }
      close(con)
}


# From unaligned to aligned positions to account for dynamic alignment
# to account for gaps in reference sequence
unaligned_to_aligned <- function(unaligned_pos, aligned_seq){
      i.unaligned <- 0
      i <- 0
      while(i.unaligned < unaligned_pos){
            i = i + 1
            if(aligned_seq[i] != '-'){ # if there is a gap in the aligned_seq (refseq)
                  i.unaligned <- i.unaligned + 1
            }
      }
      return(i) # returns aligned position according to reference (without gaps)
}



## GET AA FROM NUC 
nuc_to_aa_matrix <- function(nuc_mat){
      
      # get length of prot
      prot_length <- ncol(nuc_mat)/3
      #create empty matrix
      aa_mat <- matrix(nrow=0, ncol=prot_length)
      
      # loop over nucleotide matrix to translate by codon
      for(i.row in 1:nrow(nuc_mat)){
            
            sequence_ID <- rownames(nuc_mat)[i.row]
            
            # initiate codon start and end positions
            c_start <- 0
            c_end <- 0
            # initiate empty vec to store all amino acids
            all_aa <- c()
            # initiate empty vector to combine nucleotides in case of gaps
            to_keep <- c()
            
            # begin extracting codons
            while(c_end < ncol(nuc_mat)){
                  
                  # for first nucleotide:
                  if(c_start==0){
                        c_start <-c_start+1
                        # continue
                  } else {
                        c_start <- c_start+3
                  }
                  c_end <- c_end+3
                  codon <- nuc_mat[i.row,c_start:c_end]
                  
                  
                  # check for nucleotide deletions indicated by "-" (from alignment)
                  # in CODON
                  gap_no <- length(which(codon=="-"))
                  #print(which(n_nuc[i.row,]!="-"))
                  
                  # if no gaps, translate normally
                  if(gap_no==0){
                        aa <- translate(codon)
                  }
                  
                  # if there are gaps, investigate:
                  if(gap_no>0){
                        # if 3 gaps all(codon=="-")
                        if(gap_no==3){
                              # assign the amino acid as del
                              aa <- "-"
                        }
                        # if 3 gaps are dispersed between two codons: do more work:
                        if(gap_no<3){
                              # assign THIS aa as del
                              aa <- "-"
                              nuc_ind <- which(codon != "-")
                              keep_nuc <- codon[nuc_ind]

                              # store the nucleotides to keep, to follow up with 
                              to_keep <- c(to_keep, keep_nuc)
                        }

                        # take this codon and translate it
                        if(length(to_keep)==3){
                              aa <- translate(to_keep)
                              #empty when 3
                              to_keep  <- c()
                        }
                  }
                  # Add all the amino acids together for each row
                  all_aa <- c(all_aa,aa)
                  
            }
            #rbind to form aa matrix
            aa_mat <- rbind(aa_mat, all_aa)
      }
      # keep IDs by keeping rownames
      rownames(aa_mat) <- rownames(nuc_mat)
      return(aa_mat)
}

# WRITE SEQUENCES FROM LIST TO FASTA
list_to_fasta <- function(filename, seqs.list, append = TRUE) {
      
      if(append == TRUE){
            con <- file(filename, "a")
      } else {
            con <- file(filename, "w") 
      }
      for(i.seq in 1:length(seqs.list)) {
            write(paste0(">", names(seqs.list)[i.seq]), file= con)
            write(paste(seqs.list[[i.seq]], collapse = ""), file= con)
      }
      close(con)
}

# WRITE SEQUENCES FROM LIST TO FASTA
matrix_to_fasta <- function(seq_mat, fasta_name, add_ref=T, append=T) {
      if(append == T){
            output_fasta <- file(fasta_name, "a")
      } else {
            output_fasta <- file(fasta_name, "w")
            
      }
      #OPTION TO ADD REF
      if(add_ref == T){
            
            for(i.row in 1:nrow(seq_mat)) { #begins from 2nd row
                  write(paste0(">", rownames(seq_mat)[i.row]), file = output_fasta)
                  write(paste(seq_mat[i.row,], collapse = ""), file = output_fasta)
            }
            
      } else {
            
            for(i.row in 2:nrow(seq_mat)) { #begins from 2nd row
                  write(paste0(">", rownames(seq_mat)[i.row]), file = output_fasta)
                  write(paste(seq_mat[i.row,], collapse = ""), file = output_fasta)
            }
            
      }
      
      
      close(output_fasta)
}


# Mutation detection function
find_AA_mutations <- function(aa_mat, output_file, append=T, header=F){
      header_names <- paste("Sequence_ID", "AA_Position", "Reference_AA", "AA_Mutation", "Type", sep=",")
      
      if(append==T){
            mutation_con <- file(output_file, "a")

      } else {
            mutation_con <- file(output_file, "w")
            #write(header_names, mutation_con)
      }
      
      if(header==T){
            write(header_names, mutation_con)
      }
      
      # create an empty matrix with number of cols of output
      # cols will be: ID, position, reference aa, mutation, type of mutation
      # mutations_mat <- matrix(nrow=0, ncol = 5)
      for(i.row in 2:nrow(aa_mat)) {
            sequence_ID <- rownames(aa_mat)[i.row]
            aa_position <- which(aa_mat[1,] != aa_mat[i.row,])
            if(length(aa_position) > 0){
                  for(j.mut in aa_position){
                        ref_aa <- aa_mat[1,][j.mut]
                        mutation_aa <- aa_mat[i.row,][j.mut]
                        
                        if(mutation_aa == "-"){
                              type <- "deletion"
                        }
                        if(mutation_aa == "X"){
                              type <- "missing"
                        }
                        aa_letters <- as.character(unique(SEQINR.UTIL$CODON.AA$L))
                        if(mutation_aa %in% aa_letters){
                              type <- "substitution"
                        }
                        #mutations_mat <- rbind(mutations_mat, c(sequence_ID, j.mut, ref_aa, mutation_aa, type))
                        output <- paste(sequence_ID, j.mut, ref_aa, mutation_aa, type, sep=",")
                        #print(output)
                        write(output, file=mutation_con)
                        
                  }
            }
      }
      #colnames(mutations_mat) <- c("ID", "aa_position", "ref_aa", "mutation_aa", "type")
      #return(mutations_mat)
      close(mutation_con)
}


##########################################################################################
args <-  commandArgs(trailingOnly=TRUE)

# input wd path
wd <- args[1]
#setwd(wd)

# results dir
results_wd <- paste0(wd, "/mutations_output")
dir.create(results_wd)
setwd(results_wd)

# input fasta name
read_big_fasta <- file(paste0(wd, "/", args[2]), "r")

# input refseq path
refseq_path <- args[3]


# get ORF positions (file names here are for NP, but can be changed)
start_pos <- as.numeric(args[4]) #28274
end_pos <- as.numeric(args[5]) #29533

# set simultaneous seq number
simultaneous_seqs <- 100


# open log file to track steps
log_file <- file("find_mut_log.txt", "w")


# begin record 
rec <- paste("Mutation analysis begins in batches of:", simultaneous_seqs, "//", Sys.time())
print(rec)
write(rec, log_file)

# set counter
counter <- 1
batch_no <- paste0("batch ", counter,":")
rec <- paste("Processing batch", counter, "begins")
print(rec)
write(rec, log_file)

start_time <- Sys.time()
prev_seq <- read_n_seqs(read_big_fasta , n = simultaneous_seqs, "batch.fasta", previous_line = NULL)

Last_part <- FALSE
while(!is.null(prev_seq) || Last_part){
      
      start_batch <- Sys.time()
      
      # READ AND ALIGNMENT STEPS #########################################################
      # read batch file as fasta (seqinr)
      seqs <- read.fasta(file = "batch.fasta", whole.header = TRUE)
      
      # write to reference
      write_with_ref("input_seqs.fasta", seqs, refseq_path)
      
      # align with MAFFT
      time_now <- Sys.time()
      cmd <- "mafft --thread -1 --quiet input_seqs.fasta > aligned_seqs.fasta"
      system(cmd)
      
      ## RECORD TIME: ALIGNMENT TO REF
      time_taken <- round(as.numeric(difftime(Sys.time(), time_now, units="secs")), digits = 3)
      rec <- paste(batch_no, "Alignment with MAFFT took", time_taken, "secs")
      print(rec)
      write(rec, log_file)
      
      # read alignment as matrix
      time_now <- Sys.time()
      genomes_aligned <- read.dna(file = "aligned_seqs.fasta", as.character = TRUE, as.matrix = TRUE, format="fasta")
      ref_genome <- genomes_aligned[1,]
      
      
      #EXTRACT NUCLEOCAPSID REGION #######################################################
      # get aligned positions of Nucleocapsid regions
      # taking into consideration dynamic alignment
      start_pos <- unaligned_to_aligned(28274, ref_genome)
      end_pos <- unaligned_to_aligned(29533, ref_genome)
      
      # extract N nucleotide matrix
      n_nuc <- genomes_aligned[,c(start_pos:end_pos)]

      ## RECORD TIME: EXTRACTING N REGION
      time_taken <- round(as.numeric(difftime(Sys.time(), time_now, units="secs")), digits = 3)
      rec <- paste(batch_no, "Extracting N nucloetide matrix took", time_taken, "secs")
      print(rec)
      write(rec, log_file)
      
      ## ANALYSE FOR MISSING DATA: INSERTIONS +EXCESS GAPS ###############################
      # check for insertions: gaps in reference: write in log
      time_now <- Sys.time()
      ref_gaps <- which(n_nuc[1,]=="-")
      
      #IF INSERTIONS FOUND
      if(length(ref_gaps)>0){
            rec <- paste(batch_no, "Reference contains gaps: check for insertions")
            print(rec)
            write(rec,log_file)
            
            
            # maintain file: if first batch: write a new file
            if(counter==1){
                  ins_con <- file("insertions_detected.csv", "w")
            }
            
            # if in following batches,allow to append
            #file.exists("insertions_detected.csv")
            if(counter>1){
                  ins_con <- file("insertions_detected.csv", "a")
            }
            
            # take the N nucleotide sequence to investigate later..
            for(i.row in 1:nrow(n_nuc)){
                  sequence_ID <- rownames(n_nuc)[i.row]
                  
                  t <- which(n_nuc[i.row,c(ref_gaps)]!="-")
                  
                  ref_ID <- rownames(n_nuc)[1]
                  ref_nuc <- paste(n_nuc[1,], collapse = "")
                  #print(t)
                  if(length(t)>0){
                        # write reference
                        write(paste0(">",ref_ID),ins_con)
                        write(ref_nuc,ins_con)
                        
                        # write sequence with insertion
                        write(paste0(">",sequence_ID),ins_con)
                        nuc_string <- paste(n_nuc[i.row,], collapse = "")
                        write(nuc_string,ins_con)
                        
                  }
            }
            
            close(ins_con)

            
            # remove the columns to maintain n_nuc length for analysis
            n_nuc <- n_nuc[,-c(ref_gaps)]
            
            
      }
      
      # if deletions, write to missing data con
      for(i.row in 1:nrow(n_nuc)){
            sequence_ID <- rownames(n_nuc)[i.row]
            nuc_length <- length(which(n_nuc[i.row,]!="-"))
            gap_ind <- which(n_nuc[i.row,]=="-")
            md_out <- paste(sequence_ID, nuc_length, length(gap_ind), sep=",")
            
            # if gaps greater than the reported for omicron NP: 
            # common del at: 89 90 91 92 93 94 95 96 97
            
            # if consecutive gaps..at terminals?
            
            # IF EXCESSIVE DELETIONS FOUND
            if(length(gap_ind)>9){
                  
                  if(!file.exists("excessive_dels.csv")){
                        md_con <- file("excessive_dels.csv", "w")
                        md_header <- paste("ID", "nuc_length", "del_length",sep=",")
                        write(md_header, md_con)
                        
                  } else {
                        md_con <- file("excessive_dels.csv", "a")
                  }
                  
                  write(md_out, md_con)
                  close(md_con)
                  
            }
            
      }
      
      ## RECORD TIME: CHECKING FOR GAPS
      time_taken <- round(as.numeric(difftime(Sys.time(), time_now, units="secs")), digits = 3)
      rec <- paste(batch_no, "Checking for gaps took", time_taken, "secs")
      print(rec)
      write(rec, log_file)
      
      # TRANSLATE AND FIND MUTATIONS #####################################################

      # translate nucleotide to amino acid sequences
      # convert nucleotide matrix to amino acid matrix : takes into account deletions
      time_now <- Sys.time()
      n_aa_mat <- nuc_to_aa_matrix(n_nuc)
      
      ## RECORD TIME: TRANSLATING MATRIX 
      time_taken <- round(as.numeric(difftime(Sys.time(), time_now, units="secs")), digits = 3)
      rec <- paste(batch_no, "Translating to amino-acid matrix took", time_taken, "secs")
      print(rec)
      write(rec, log_file)

      # FIND MUTATIONS
      # For first batch, add the reference, overwrite existing files.
      if(counter==1){
            time_now <- Sys.time()
            # store all N nucleotide sequences to a file
            matrix_to_fasta(n_nuc, "all_N_nuc.fasta", add_ref = T, append=F)
            
            # store all N AA sequences to a file
            matrix_to_fasta(n_aa_mat, "all_N_peptides.fasta", add_ref = T, append=F)
            
            # find mutations in N peptide matrix
            find_AA_mutations(n_aa_mat, "N_aa_mutations.csv", append=F, header=T)
            
            ## RECORD TIME: FINDING N MUTATIONS
            time_taken <- round(as.numeric(difftime(Sys.time(), time_now, units="secs")), digits = 3)
            rec <- paste(batch_no, "N mutation analysis took", time_taken, "secs")
            print(rec)
            write(rec, log_file)
      } 
      
      # After first batch, append to existing
      if(counter>1){
            time_now <- Sys.time()
            # store all N nucleotide sequences to a file
            matrix_to_fasta(n_nuc, "all_N_nuc.fasta", add_ref = F, append=T)
            
            # store all N AA sequences to a file
            matrix_to_fasta(n_aa_mat, "all_N_peptides.fasta", add_ref = F, append=T)
            
            # find mutations in N peptide matrix
            find_AA_mutations(n_aa_mat, "N_aa_mutations.csv", append=T, header=F)
            
            ## RECORD TIME: FINDING N MUTATIONS
            time_taken <- round(as.numeric(difftime(Sys.time(), time_now, units="secs")), digits = 3)
            rec <- paste(batch_no, "N mutation analysis took", time_taken, "secs")
            print(rec)
            write(rec, log_file)
      }
      
      # MOVE TO NEXT BATCH ###############################################################
      end_batch <- Sys.time()
      
      time_taken <- round(as.numeric(difftime(end_batch, start_batch, units="secs")), digits = 3)
      rec <- paste(batch_no, "Complete batch processing took", time_taken, "secs")
      print(rec)
      write(rec, log_file)

      if(Last_part){
            break
      }

      counter <- counter + 1
      batch_no <- paste0("batch ", counter,":")
      rec <- paste("Processing batch", counter, "begins")
      print(rec)
      write(rec, log_file)
      
      # analyse next batch of seqs
      prev_seq <- read_n_seqs(read_big_fasta, n = simultaneous_seqs, "batch.fasta", prev_seq)
      
      # IF END OF FILE: END ##############################################################
      if(!Last_part && is.null(prev_seq)){
            Last_part <- TRUE
      }
      
}
# close the files
close(read_big_fasta)
#close(md_con)

# RECORD END #############################################################################
rec <- paste("Mutation analysis ends in batches of:", simultaneous_seqs, "//", Sys.time())
print(rec)
write(rec, log_file)
end_time <- Sys.time()
time_taken <- round(as.numeric(difftime(end_time, start_time, units="mins")), digits = 3)
rec <- paste("Mutation analysis took", time_taken, "mins")
print(rec)
write(rec, log_file)

#zip temporary files: with r.utils #######################################################
time_now <- Sys.time()
suppressPackageStartupMessages(library(R.utils))
gzip("batch.fasta")
gzip("input_seqs.fasta")
gzip("aligned_seqs.fasta")
gzip("all_N_nuc.fasta")
gzip("all_N_peptides.fasta")
time_taken <- round(as.numeric(difftime(Sys.time(), time_now, units="secs")), digits = 3)
rec <- paste("Compressing temp files took", time_taken, "secs")
print(rec)
write(rec, log_file)
close(log_file)

