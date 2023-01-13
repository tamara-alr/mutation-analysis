## ---------------------------------------------------------------------------------------

# Script: Pipeline for sequence custom filtering (metadata and fasta parsing)
# Purpose: To filter and extract complete high coverage sequences and good quality sequences
# by custom function by batch processing two parallel files (TSV & FASTA)

# Authors: Tamara Alrawahneh (CEPHR, UCD), Gabriel Gonzalez (NVRL, UCD)
# Last edits: 2022/09/14

## ---------------------------------------------------------------------------------------
#ARGUMENTS
# args[1] is working directory
# args[2] is .fasta file
# args[3] is .tsv file
# args[4] is number of sequences to process

## LIBRARIES -----------------------------------------------------------------------------
# load seqinr to read.fasta

if(!require(seqinr)){
      install.packages("seqinr")
      library(seqinr)
}

# load R.utils to compress temp files at end

if(!require(R.utils)){
      install.packages("R.utils")
      library(R.utils)
}

## ---------------------------------------------------------------------------------------


## FUNCTIONS -----------------------------------------------------------------------------
# six functions
# READ_N_DATA, READ_N_SEQS : used to read data lines or sequences in batches
# NB: OPEN READING CONTROL OUTSIDE LOOP, CLOSE AT END

# read_n_data reads n lines, returns last line, pattern indicating line is "hCov" 
# use read_n_data to read metadata TSV file
read_n_data <- function(reading_control, n = 100, output_file, previous_line = NULL){ 
      
      writing_control <- file(output_file, "w")
      seqs_read <-  0 
      my_line <-  ""
      
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
      close(writing_control)
      return(my_line)
}

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

# EXTRACT SEQUENCES USING A PATTERN VECTOR - FASTA OUTPUT
# NB: two versions to use in loop: APPEND = TRUE to grow file, APPEND = FALSE to overwrite batch
extract_fasta_from_vec <- function(big_file, output_file, pattern.vec, append = TRUE){ 
      
      reading_control <- file(big_file, "r")
      if(append == TRUE){
            writing_control <- file(output_file, "a")
      } else {
            writing_control <- file(output_file, "w") 
      }
      
      my_line <-  "" 
      seq_count <-  0 #initiate sequence counter
      
      write_to_file <-  FALSE #initiate writing to file as FALSE
      my_line <-  readLines(con = reading_control, n = 1) 
      
      
      while(!is.null(my_line)){ 
            
            if(length(grep(pattern = ">", my_line)) > 0){
                  
                  index <-  grep(pattern = gsub("^>|\\|(.*?)\\d{4}-\\d{2}-\\d{2}$", "", my_line), pattern.vec, fixed = F)
                  
                  if(length(index) > 0){
                        #print(paste(index, my_line, sep="  "))
                        # remove ID once found
                        pattern.vec <-  pattern.vec[-index]
                        write_to_file <- TRUE
                        seq_count <-  seq_count + 1
                  } else {
                        write_to_file <- FALSE
                  }
            }
            if(write_to_file){ 
                  write(my_line, file = writing_control)
            }
            my_line <-  readLines(con = reading_control, n = 1)
            if(length(my_line) == 0){
                  break
            }
      }
      
      close(reading_control)
      close(writing_control)
      
}


# EXTRACT DATA USING A PATTERN VECTOR - TSV OUTPUT
# modified to readLines matching beginning of ID vec
# NB: two versions to use in loop: APPEND = TRUE to grow file, APPEND = FALSE to overwrite batch

extract_data_from_vec <- function(big_file, output_file, pattern.vec, append = TRUE){ 
      reading_control <- file(big_file, "r")
      if(append == TRUE){
            writing_control <- file(output_file, "a")
      } else {
            writing_control <- file(output_file, "w") 
      }
      my_line <-  ""
      seq_count <-  0 #initiate sequence counter
      write_to_file <-  FALSE #initiate writing to file as FALSE
      my_line <-  readLines(con = reading_control, n = 1) #read line by line
      #store <- c() #  debugging, vec to store IDs
      while(!is.null(my_line)){
            if(length(grep(pattern = "hCoV", my_line, ignore.case = T)) > 0){
      
                  # Split lines 
                  split_lines <- strsplit(my_line, split="\t", fixed=TRUE)
                  # take first split, header names
                  for(line in split_lines){
                        index <-  grep(pattern = line[1], pattern.vec, fixed = TRUE)
                        # get which headers match the pattern vector
                        if(length(index) > 0){
                              #print(paste(index, my_line, sep="  "))
                              # remove ID once found
                              pattern.vec = pattern.vec[-index]
                              write_to_file <- TRUE
                              seq_count <-  seq_count + 1
                              #store <- c(store, index)
                        } else {
                              write_to_file <- FALSE
                        }
                        
                        if(write_to_file){
                              #store <- c(store,my_line)
                              write(my_line, file = writing_control)
                        }
                  }
            }
            my_line <-  readLines(con = reading_control, n = 1)
            if(length(my_line) == 0){
                  break
            }
      }
      #return(length(store))
      close(reading_control)
      close(writing_control)
}


# CUSTOM FILTER FUNCTION 
# to remove sequences with missing data (Ns and gaps)
filter_seqs <- function(seq.vec){
      count_nuc <- length(which(seq.vec == "a" | seq.vec == "c" | 
                                      seq.vec == "g" | seq.vec == "t"))
      return((count_nuc / length(seq.vec)) >= 0.9999)
}


# WRITE SEQUENCES FROM LIST TO FASTA
write_to_fasta <- function(filename, seqs.list, append = TRUE) {
      
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

# GET DATA IDS FROM TSV
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

## ---------------------------------------------------------------------------------------
#Rscript check.R C:/Users/alraw/Documents/test_files/BigChunk86/chunk15 chunk15.fasta chunk15.tsv

## FILTER SEQUENCES ANALYSIS -------------------------------------------------------------

# set wd
args <-  commandArgs(trailingOnly=TRUE)
wd <- args[1]
setwd(wd)

# set simultaneous sequences
#simultaneous_seqs <- 1000
simultaneous_seqs <- as.numeric(args[4])

# writing controls
log_file <- file("time_log.txt", "w") #open log file to record times
error_log <- file("filter_errors.txt", "w")

# Open reading controls
#read_big_fasta <- file("10k_sequences.fasta", "r")
#read_big_data <- file("10k_metadata.tsv", "r")
read_big_fasta <- file(args[2], "r")
read_big_data <- file(args[3], "r")


# write header to output metadata files
#header <- readLines("10k_metadata.tsv", n=1)
header <- readLines(args[3], n=1)
write(header, file = "chc_data.tsv")
write(header, file = "filtered_data.tsv")


# begin record
rec <- paste("Filtering of metadata and FASTA files begins in batches of:", simultaneous_seqs, "//", Sys.time())
print(rec)
write(rec, log_file)

# initiate counter and batch number
counter <- 1
batch_no <- paste0("batch ", counter,":")
start_time <- Sys.time()
rec <- paste("processing batch", counter, "begins")
print(rec)
write(rec, log_file)

# FIRST METADATA BATCH
time_now <- Sys.time()
prev_line <- read_n_data(read_big_data, n = simultaneous_seqs, "batch.tsv", previous_line = NULL)
minutes <-  round(as.numeric(difftime(Sys.time(), time_now, units="secs")), digits = 3)
rec <- paste(batch_no,"extracting batch metadata took:", minutes, "secs")
print(rec)
write(rec, log_file)

# FIRST FASTA BATCH
time_now <- Sys.time()
prev_seq <- read_n_seqs(read_big_fasta, n = simultaneous_seqs, "batch.fasta", previous_line = NULL)
minutes <-  round(as.numeric(difftime(Sys.time(), time_now, units="secs")), digits = 3)
rec <- paste(batch_no,"extracting batch fasta took:", minutes, "secs")
print(rec)
write(rec, log_file)
#write("------------------", file= time_log)

last_batch <- FALSE
while(!is.null(prev_line) &&!is.null(prev_seq) || last_batch){
      
      ### GET CHC IDS ####################################################################
      # extract complete high coverage (chc) sequeneces
      rec <- paste(batch_no, "extracting CHC IDs") 
      #print(rec)
      time_now <- Sys.time()
      # readlines and store chc IDs as vector
      # read data lines and split lines to filter by header
      data_lines <- readLines("batch.tsv")
      split_lines <- strsplit(data_lines, split="\t", fixed=TRUE)
      chc <- c() # initiate empty vector for chc
      with.date <- c() 
      for(line in split_lines){
            # find which lines have full collection date recorded
            have_date <- grep(pattern = "\\d{4}-\\d{2}-\\d{2}", line[4])
            # no_date <- grep(pattern = "\\d{4}$", line[4])
            # find chc only in lines with collection date:
            if(length(have_date) > 0){
                  #do i need to store the with.date part?
                  with.date <- c(with.date, line[have_date])
                  ind <- grep(pattern = "True", line[c(18,19)])
                  if(length(ind) == 2){
                        chc_ids <- line[ind[1]]
                        chc <- c(chc, chc_ids)
                  }
                  
            }
            
      }
      
      # remove duplicates
      unique_chc <- unique(chc)
      dup <- chc[duplicated(chc)]
      
      if(length(dup)>0){
            chc <- chc[!(duplicated(chc) | duplicated(chc, fromLast = TRUE))]
            rec <- paste(batch_no, length(dup), "duplicate(s) removed")
            print(rec)
            write(rec, file = log_file)
            #print(length(chc))
      }
      
      # print to screen number of sequences removed
      seqs_removed <- simultaneous_seqs - length(with.date)
      if(seqs_removed > 0){
            rec <- paste(batch_no, seqs_removed, "sequence(s) missing collection date removed")
            print(rec)
            write(rec, file = log_file)
      }
      
      # if there are no CHC sequences, move to next batch
      if(length(chc) == 0){
            rec <- paste("no chc sequences found: move to next batch")
            print(rec)
            write(rec, log_file)
      }
      
      # if there are chc sequences, move ahead with filtering
      if(length(chc)>0){
            # record time taken to get chc vector
            minutes <-  round(as.numeric(difftime(Sys.time(), time_now, units="secs")), digits = 3)
            rec <- paste(batch_no,"extracting", length(chc), "CHC IDs took:", minutes, "secs")
            print(rec)
            write(rec, file = log_file)
            #print(paste(batch_no,  "number of CHC sequnces found:",length(chc), "/", simultaneous_seqs))
            
            ### GET CHC DATA ###################################################################
            
            # use chc pattern.vec to extract metadata, begins now
            rec <- paste(batch_no, "extracting CHC metadata")
            #print(rec)
            time_now <- Sys.time()
            # write chc_data which will keep growing
            extract_data_from_vec("batch.tsv", "chc_data.tsv", chc, append = TRUE)
            # write chc_batch which will be overwritten - purpose is to use it to temporarily 
            extract_data_from_vec("batch.tsv", "chc_batch.tsv", chc, append = FALSE)
            # record time taken to get chc metadata - writing TSV
            minutes <-  round(as.numeric(difftime(Sys.time(), time_now, units="secs")), digits = 3)
            rec <- paste(batch_no,"extracting CHC metadata took:", minutes, "secs")
            print(rec)
            write(rec, file = log_file)
            
            ### GET CHC FASTA ##################################################################
            # use chc pattern.vec to extract FASTA sequences, begins now
            rec <- paste(batch_no,"extracting CHC FASTA sequences")
            #print(rec)
            time_now <- Sys.time()
            # write chc_data which will keep growing
            extract_fasta_from_vec("batch.fasta", "chc_seqs.fasta", chc, append = TRUE)
            # write chc_batch which will be overwritten - purpose is to use it to temporarily 
            extract_fasta_from_vec("batch.fasta", "chc_batch.fasta", chc, append = FALSE)
            # record time taken to get chc sequences - writing FASTA
            minutes <-  round(as.numeric(difftime(Sys.time(), time_now, units="secs")), digits = 3)
            rec <- paste(batch_no,"extracting CHC FASTA sequences took:", minutes, "secs")
            print(rec)
            write(rec, file = log_file)
            
            
            ### GET FILTERED FASTA  ############################################################
            # custom filtering step begins now
            rec <- paste(batch_no, "custom-filtering CHC FASTA sequences")
            #print(rec)
            time_now <- Sys.time()
            
            # using the chc_batch file, read the sequences as list using seqinr
            #NB! ADD whole.header = TRUE, some IDs have spaces. TRUE avoids truncating names after space
            seqs <- read.fasta(file = "chc_batch.fasta", whole.header = TRUE)
            
            good_quality_index <- which(sapply(seqs, filter_seqs)) # take index of good quality seqs
            good_seqs <- seqs[good_quality_index] # subset good quality seqs
            good_seqs_num <- length(good_seqs) # record numbers
            seqs_num <- length(seqs)
            
            if(length(good_seqs) == 0){
                  rec <- paste("no good quality sequences found: move to next batch")
                  print(rec)
                  write(rec, log_file)
            }
            
            # if good seqs > 0 : write to fasta and tsv
            if(length(good_seqs)> 0){
                  # record time taken to get chc sequences - writing FASTA
                  minutes <-  round(as.numeric(difftime(Sys.time(), time_now, units="secs")), digits = 3)
                  rec <- paste(batch_no,"custom-filtering CHC sequences took:", minutes, "secs")
                  print(rec)
                  write(rec, file = log_file)
                  
                  # record number of filtered sequences
                  rec <- paste(batch_no,"number of sequences filtered:", good_seqs_num, "/", seqs_num)
                  print(rec)
                  write(rec, file = log_file)
                  
                  # write filtered seqs to fasta
                  rec <- paste(batch_no, "writing filtered FASTA sequences")
                  #print(rec)
                  time_now <- Sys.time()
                  # filtered_seqs.fasta will keep growing, append=T
                  write_to_fasta("filtered_seqs.fasta", good_seqs, append = TRUE)
                  # filtered_batch.fasta will be recycled, append=F
                  write_to_fasta("filtered_batch.fasta", good_seqs, append = FALSE) # will be used to get fastaIDs
                  
                  # record time taken to write filtered FASTA
                  minutes <-  round(as.numeric(difftime(Sys.time(), time_now, units="secs")), digits = 3)
                  rec <- paste(batch_no,"writing filtered sequences to FASTA took:", minutes, "secs")
                  print(rec)
                  write(rec, file = log_file)
                  
                  ### GET FILTERED DATA  #############################################################
                  # extract fiiltered IDs to extract filtered sequence metadata
                  time_now <- Sys.time()
                  rec <- paste(batch_no, "extracting filtered sequence metadata")
                  #print(rec)
                  # get fasta IDs function returns the IDs in metadata header format
                  fasta_IDs <- get_fasta_ids("filtered_batch.fasta")
                  #rec <- paste("FASTA IDs", batch_no, length(fasta_IDs))
                  #print(rec)
                  #write(rec, file = log_file)
                  
                  # extract filtered tsv 
                  extract_data_from_vec("chc_batch.tsv", "filtered_data.tsv", fasta_IDs, append = TRUE)
                  minutes <-  round(as.numeric(difftime(Sys.time(), time_now, units="secs")), digits = 3)
                  rec <- paste(batch_no, "extracting filtered sequence metadata took:", minutes, "secs")
                  print(rec)
                  write(rec, file = log_file)
                  
                  
            }
            
            
      }
      
      # ENDS #
      rec <- paste("processing batch", counter, "complete")
      print(rec)
      write(rec, file = log_file)
      write("----------------------------------------------------------------", file= log_file)
      
      
      ### MOVE TO NEXT BATCH #############################################################
      if(last_batch){
            break
      }
      
      counter <- counter + 1
      batch_no <- paste0("batch ", counter,":")
      rec <- paste("processing batch", counter, "begins")
      print(rec)
      write(rec, file = log_file)
      
      # NEXT METADATA BATCH
      time_now <- Sys.time()
      prev_line <- read_n_data(read_big_data, n = simultaneous_seqs, "batch.tsv", previous_line = prev_line)
      minutes <-  round(as.numeric(difftime(Sys.time(), time_now, units="secs")), digits = 3)
      rec <- paste(batch_no,"extracting batch metadata took:", minutes, "secs")
      print(rec)
      write(rec, log_file)
      
      # NEXT FASTA BATCH
      time_now <- Sys.time()
      prev_seq <- read_n_seqs(read_big_fasta, n = simultaneous_seqs, "batch.fasta", previous_line = prev_seq)
      minutes <-  round(as.numeric(difftime(Sys.time(), time_now, units="secs")), digits = 3)
      rec <- paste(batch_no,"extracting batch fasta took:", minutes, "secs")
      print(rec)
      write(rec, log_file)
      
      if(!last_batch && is.null(prev_line) && is.null(prev_seq)){
            last_batch <- TRUE
            
      }
      
}

end_time <- Sys.time()
rec <- paste("Filtering of metadata and FASTA files ends", "//", end_time)
print(rec)
write(rec, file = log_file)
rec <- paste("Filtering took:", round(as.numeric(difftime(end_time, start_time, units="mins")), digits = 2), "mins")
print(rec)
write(rec, file = log_file)


#find alternative to counting sequences? when file gets bigger might be an issue

total_filtered <- length(get_data_ids("filtered_data.tsv"))
total_chc <- length(get_data_ids("chc_data.tsv"))
rec <- paste("total number of selected sequences:", total_filtered)
print(rec)
write(rec, file = log_file)

rec <- paste("total number of CHC sequences:", total_chc)
print(rec)
write(rec, file = log_file)




#check if filtered_IDs match
check_filtered_IDs <- function(path){
      setwd(path)
      #print(path)
      get_fasta_ids("filtered_seqs.fasta")
      get_data_ids("filtered_data.tsv")
      do_IDs_match(fasta.id, data.id)
      
}

check <- check_filtered_IDs(args[1])
match <- "IDs matched in order"
if(check == match){
      rec <- paste("Filtered", check, "in FASTA and TSV files")
      print(rec)
      write(rec, log_file)
      
} else {
      #write output to error file
      write(paste("Mismatches in filtered output", file_index, check), error_log)
      write(paste("data IDs:",length(data.id),",fasta IDs:", length(fasta.id)), error_log)
}

# close connections
close(read_big_fasta)
close(read_big_data)
close(log_file)


# zip temp files
library(R.utils)
#gzip(args[2], remove=T);
gzip("batch.tsv")
gzip("batch.fasta")
#unlink("chc_data.tsv")
#unlink("chc_seqs.fasta")
#unlink("selected_IDs.txt")
#gzip("filtered_data.tsv")
gzip("filtered_batch.fasta")
#gzip("filtered_seqs.fasta")
gzip("chc_batch.fasta")
gzip("chc_batch.tsv")

