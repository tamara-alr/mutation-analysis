## ---------------------------------------------------------------------------------------

# Script: File chunking
# Purpose: To chunk large file (fasta or metadata) into smaller files
# Authors: Tamara Alrawahneh (CEPHR, UCD), Gabriel Gonzalez (NVRL, UCD)
# Date: 2022/08/01
# Last edits: 2022/08/25
# Comments: each chunk (fasta, metadata) is stored in a chunk_index folder. 
# Chunk pairs (fasta, metadata) are taken separately without extracting using vector IDs.
# chunk fasta IDs and chunk data IDs are checked, if not identical, reported in error file

## ---------------------------------------------------------------------------------------
#ARGUMENTS
# args[1] is working directory
# args[2] is .fasta file
# args[3] is .tsv file
# args[4] is number of sequences in chunk

## ---------------------------------------------------------------------------------------
# libraries
# if want to zip output chunks for storage, will require R.utils

#if(!require(R.utils)){
      #install.packages("R.utils")
      #library(R.utils)
#}

## FUNCTIONS -----------------------------------------------------------------------------
# BATCH FROM TSV FILE
batch_from_tsv <- function(reading_control, n = 1000, output_file, previous_line = NULL, header = NULL, back_to_wd){
      
      # create folder to store chunk in
      folder_name <- sub("^([^.]*).*", "\\1", output_file) # name folder based on file name
      dir.create(folder_name, recursive = TRUE) # create directory
      setwd(folder_name) # setwd with same name
      
      # assign writing control 
      writing_control <- file(output_file, "w")
      seqs_read <-  0  #initiate number of sequences read
      my_line <-  ""   #initiate my_line
      
      # add metdata header
      if(!is.null(header)){
            write(header, file = writing_control)
      }
      
      if(!is.null(previous_line)){ #if previous line is not missing, write it to the control
            write(previous_line, file = writing_control)
            seqs_read <-  1 #count this as how many sequences read
      }
      
      while(seqs_read <= n && !is.null(my_line)){ #while loop to read lines one by one
            my_line <-  readLines(con = reading_control, n = 1)
            
            if(!is.null(my_line) && length(grep(pattern = "hCoV", my_line, ignore.case = T)) > 0){ 
                  seqs_read <- seqs_read + 1 #update the sequences read by 1, while will keep reading until n = 100
            }
            
            if(!is.null(my_line) && seqs_read <= n){ #write the lines to the writing control
                  write(my_line, file = writing_control)
            }
            #THIS WORKS
            if(length(my_line) == 0){
                  print(paste0("reached end of sequence metadata file"))
                  my_line = NULL
                  break
            }
      }
      # end
      setwd(back_to_wd) #back to wd
      close(writing_control) 
      return(my_line) #returns last line
}

# BATCH FROM FASTA FILE
batch_from_fasta <- function(reading_control, n = 1000, output_file, previous_line = NULL, back_to_wd){ 
      
      folder_name <- sub("^([^.]*).*", "\\1", output_file) # name folder based on file name
      setwd(folder_name) # setwd with same name
      
      # open writing control
      writing_control <- file(output_file, "w")
      seqs_read <-  0  #initiate number of sequences read
      my_line <-  ""   #initiate my_line
      
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
      # end
      setwd(back_to_wd) #back to wd
      close(writing_control) 
      return(my_line) #returns last line
}


# GET DATA IDS FROM TSV
# returns vector called data.id
get_data_ids <- function(data_file){
      data_lines <- readLines(data_file)
      split_lines <- strsplit(data_lines, split="\t", fixed=TRUE)
      data.id <- c()
      for(line in split_lines){
            id_ind <- grep("^hCoV", line, ignore.case = T)
            if(length(id_ind) > 0){
                  data.id <- c(data.id, line[1])
            }
      }
      assign("data.id" ,data.id, envir = globalenv())
}

# GET FASTA IDS FROM FASTA
# returns vector called fasta.id
get_fasta_ids <- function(fasta_file){
      fasta_lines <- readLines(fasta_file)
      fasta.id <- c()
      for(line in fasta_lines){
            fasta_ind <- grep("^>", line)
            if(length(fasta_ind) > 0){
                  new_format  <- gsub(">", "",line)
                  new_format <- gsub("|\\|(.*?)\\d{4}-\\d{2}-\\d{2}$", "", new_format)
                  fasta.id <- c(fasta.id,new_format)
            }
            
      }
      assign("fasta.id", fasta.id, envir = globalenv())
      
}

# DO VECTORS MATCH?
# do IDs match returns a character answer which can be checked
do_IDs_match <- function(vec1, vec2){
      if(all(vec1==vec2)){
            ans <- "IDs matched in order"
      }
      if(!all(vec1==vec2) && all(sort(vec1)==sort(vec2))){
            index <- which(vec1 != vec2)
            mismatch <- paste(index, collapse = ",")
            ans <- paste("IDs matched but not in order, order mixed at position(s):" , mismatch)
            
      }
      if(!all(vec1==vec2) && !all(sort(vec1)==sort(vec2))){
            index <- which(vec1 != vec2)
            mismatch <- paste(index, collapse = ",")
            ans <- paste("Mismatches found, mismatches at position(s): ", mismatch)
            #pos1 <- setdiff(vec1,vec2)
            #pos2 <- setdiff(vec2, vec1)
            #if(length(pos2)>0 && length(pos1)>0){
                  #ans <- paste("Mismatches found, mismatches at position(s): ", mismatch, "/", 
                               #pos1, "ID only present in 1st set",  pos2, "ID only present in 2nd set")
                  
            #}
      }
      return(ans)
}


# ----------------------------------------------------------------------------------------
# set args to run in cmd
args <-  commandArgs(trailingOnly=TRUE)
wd <- args[1]
setwd(wd)


# open log files 
time_log <- file("chunks_log.txt", "w")
error_log <- file("error_log.txt", "w")


# set simultaneous seqs
simultaneous_seqs <- as.numeric(args[4])

# BEGIN FILE CHUNKING
start_time <- Sys.time()
rec <- paste("Files chunking begins", paste0("(sequences per chunk: ", simultaneous_seqs,")"),start_time)
print(rec)
write(rec, time_log)
write("------------------", file= time_log)


# open reading controls
read_big_fasta <- file(args[2], "r")
read_big_data <- file(args[3], "r")
header_names <- readLines(read_big_data, n=1)

# set index. folder name, file names
file_index <- 1
folder_name <-  paste0("chunk", file_index)
fasta_name <- paste0("chunk", file_index, ".fasta")
tsv_name <- paste0("chunk", file_index, ".tsv")

# note: prev_line refers to tsv line. prev_seq refers to fasta format sequence

# BEGIN
rec <- paste("Extracting chunk", file_index)
print(rec)
write(rec, time_log)

# begin with first batch METADATA
time_now <- Sys.time()
prev_line <- batch_from_tsv(read_big_data, n = simultaneous_seqs, tsv_name, previous_line = NULL, header = header_names, wd)
# print updates and time taken
mins <- round(as.numeric(difftime(Sys.time(), time_now, units="mins")), digits = 2)
rec <- paste("Metadata chunk", file_index, "complete", "took:", mins, "mins")
print(rec)
write(rec, file=time_log)

# begin with first batch FASTA
time_now <-  Sys.time()
prev_seq <- batch_from_fasta(read_big_fasta, n = simultaneous_seqs, fasta_name, previous_line = NULL, wd)
# print updates and time taken
mins <- round(as.numeric(difftime(Sys.time(), time_now, units="mins")), digits = 2)
rec <- paste("FASTA chunk", file_index, "complete", "took:", mins, "mins")
print(rec)
write(rec, file=time_log)

# go to chunk_index folder
changewd <- paste0(wd, "/chunk", file_index)
setwd(changewd)
# get fasta ids and data ids
time_now <-  Sys.time()
get_fasta_ids(fasta_name)
get_data_ids(tsv_name)

# go back to wd
setwd(wd)

# check if IDs match
match <- "IDs matched in order"
check <- do_IDs_match(fasta.id, data.id)
if(check == match){
      print(check)
} else {
      # write output to error file
      write(paste("Mismatches in chunk", file_index, check), error_log)
}
# go back to main folder
mins <- round(as.numeric(difftime(Sys.time(), time_now, units="mins")), digits = 2)
rec <- paste("Checking IDs in chunk", file_index, "took", mins, "mins", check)
print(rec)
write(rec, file=time_log)
write("------------------", file= time_log)


last_batch <- FALSE
while(!is.null(prev_line) && !is.null(prev_seq) || last_batch){

      # set index. folder name, file names
      file_index <- file_index + 1
      folder_name <-  paste0("chunk", file_index)
      fasta_name <- paste0("chunk", file_index, ".fasta")
      tsv_name <- paste0("chunk", file_index, ".tsv")
      
      if(last_batch){
            break
      }
      
      rec <- paste("Extracting chunk", file_index)
      print(rec)
      write(rec, time_log)
      
      # begin with first batch METADATA
      time_now <- Sys.time()
      prev_line <- batch_from_tsv(read_big_data, n = simultaneous_seqs, tsv_name, previous_line = prev_line, header = header_names, wd)
      # print updates and time taken
      mins <- round(as.numeric(difftime(Sys.time(), time_now, units="mins")), digits = 2)
      rec <- paste("Metadata chunk", file_index, "complete", "took:", mins, "mins")
      print(rec)
      write(rec, file=time_log)
      
      # begin with first batch FASTA
      time_now <-  Sys.time()
      prev_seq <- batch_from_fasta(read_big_fasta, n = simultaneous_seqs, fasta_name, previous_line = prev_seq, wd)
      # print updates and time taken
      mins <- round(as.numeric(difftime(Sys.time(), time_now, units="mins")), digits = 2)
      rec <- paste("FASTA chunk", file_index, "complete", "took:", mins, "mins")
      print(rec)
      write(rec, file=time_log)
      
      # go to chunk_index folder
      changewd <- paste0(wd, "/chunk", file_index)
      setwd(changewd)
      # get fasta ids and data ids
      time_now <-  Sys.time()
      get_fasta_ids(fasta_name)
      get_data_ids(tsv_name)
      
      # go back to wd
      setwd(wd)
      
      # check if IDs match
      match <- "IDs matched in order"
      check <- do_IDs_match(fasta.id, data.id)
      if(check == match){
            print(check)
      } else {
            # write output to error file
            write(paste("Mismatches in chunk", file_index, check), error_log)
      }
      # go back to main folder
      mins <- round(as.numeric(difftime(Sys.time(), time_now, units="mins")), digits = 2)
      rec <- paste("Checking IDs in chunk", file_index, "took", mins, "mins", check)
      print(rec)
      write(rec, file=time_log)
      write("------------------", file= time_log)
      
      if(!last_batch && is.null(prev_line) && is.null(prev_seq) ){
            last_batch <- TRUE
      }
      
}
end_time <- Sys.time()
rec <- paste("Files chunking ends", end_time)
print(rec)
write(rec, time_log)

mins <- round(as.numeric(difftime(end_time, start_time, units="mins")), digits = 2)
print(paste("complete FASTA and Metadata file chunking took: ", mins, " mins"))

#close connections
close(read_big_data)
close(read_big_fasta)
close(time_log)
close(error_log)


# remove BigChunk FASTA
#unlink(args[2])
# zip BigChunk TSV
#R.utils::gzip(args[3])



