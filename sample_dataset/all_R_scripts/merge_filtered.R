## ---------------------------------------------------------------------------------------

# Script: File merging of filtered sequences
# Purpose: Get chunk directories, check file sizes and errors, merge and store in new dir
# Authors: Tamara Alrawahneh (CEPHR, UCD), Gabriel Gonzalez (NVRL, UCD)
# Last edits: 2022/09/07

## ---------------------------------------------------------------------------------------
## LIBRARIES -----------------------------------------------------------------------------

## ---------------------------------------------------------------------------------------
## FUNCTIONS -----------------------------------------------------------------------------
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

# to check if filtered data and filtered sequences have exact IDs
check_filtered_IDs <- function(path){
      setwd(path)
      #print(path)
      get_fasta_ids("filtered_seqs.fasta")
      get_data_ids("filtered_data.tsv")
      do_IDs_match(fasta.id, data.id)
      
}

# to check if chc sequences and chc data have exact IDs
check_chc_IDs <- function(path){
      setwd(path)
      #print(path)
      get_fasta_ids("chc_seqs.fasta")
      get_data_ids("chc_data.tsv")
      return(do_IDs_match(fasta.id, data.id))
      
}


# MERGE FASTA: reads lines in all paths and writes to another file
# returns the total number of sequences found
merge_fasta <- function(paths_vec, file_name, append=FALSE){
      # open file to write all
      if(append==FALSE){
            writing_control <- file(file_name, "w")
            
      } else {
            writing_control <- file(file_name, "a")
      }
      
      fasta_log <- file("fasta_log.txt", "w")
      
      # initiate to store number of sequences in fasta
      count_seqs <- c()
      
      # loop over paths and readlines
      for(f.path in paths_vec){
            #print(f.path)
            f.lines <- readLines(f.path, skipNul = T)
            #count sequences ">"
            index <- grep(pattern = ">", f.lines)
            count_seqs <- c(count_seqs, index)
            print(paste(f.path, length(index)))
            write(paste(f.path, length(index)), fasta_log)
            
            #write to merged file
            write(f.lines, writing_control)
      }
      
      write(length(count_seqs), fasta_log)
      close(writing_control)
      close(fasta_log)
      return(length(count_seqs))
      
}


# MERGE DATA: reads lines in all paths and writes to another file
# returns the total number of sequences found
merge_data <- function(paths_vec, file_name, append=FALSE){
      # open file to write all
      if(append==FALSE){
            writing_control <- file(file_name, "w")
            
      } else {
            writing_control <- file(file_name, "a")
      }
      
      data_log <- file("metadata_log.txt", "w")
      
      header <- readLines(paths_vec[1], n=1)
      write(header, writing_control)
      
      # initiate to store number of sequences in fasta
      count_lines <- c()
      
      for(d.path in paths_vec){
            
            all_lines <- readLines(d.path)
            index <- grep(pattern = "hCoV", all_lines, ignore.case = T)
            print(paste(d.path, length(index)))
            write(paste(d.path, length(index)), data_log)
            d.line <- all_lines[index]
            
            # write tp writing con
            write(d.line, writing_control)
            count_lines <- c(count_lines, index)
            
            
      }
      write(length(count_lines), data_log)
      close(writing_control)
      close(data_log)
      return(length(count_lines))
}


## ---------------------------------------------------------------------------------------


## set arguments and WD
args <-  commandArgs(trailingOnly=TRUE)
start_path <- args[1]
setwd(start_path)


# open connection to record info
merge_log <- file("merge_log.txt", "w")
rec <- paste("Checking and storing file paths begins", Sys.time())
print(rec)
write(rec, merge_log)



# get small chunk paths
error_path <-  paste0(start_path, "/error_log.txt")
if(file.size(error_path)>0){
      rec <- paste("errors found in path", error_path)
      print(rec)
      #write(rec, merge_log)
} else {
      print("no errors from file chunking")
}



chunk_paths <- list.files(start_path, pattern="^chunk[0-9]*$", full.names = T)
# record how many chunks found
rec <- paste(length(chunk_paths), "chunk paths found")
print(rec)
write(rec, merge_log)




# get FASTA & TSV paths for filtered sequences
fasta_paths <- c()
data_paths <- c()
no_fasta_paths <- c()
no_data_paths <- c()
filter_error_exists <- c()
fasta_size <- c()
data_size <- c()



# check if filter run had errors, only take files that were filtered

for(path in chunk_paths){
      # check if error
      filter_errors <- paste0(path, "/", "filter_errors.txt")

      if(file.size(filter_errors) > 0){
            print(filter_errors)
            filter_error_exists <- c(filter_error_exists, filter_errors)
      }
      
      # FASTA PATHS
      fasta <- paste0(path, "/", "filtered_seqs.fasta")
      if(file.size(fasta) > 0){
            #print(paste(fasta, file.size(fasta)))
            fasta_size <- c(fasta_size, file.size(fasta))
            fasta_paths <- c(fasta_paths, fasta)
      } else {
            #store missing info
            no_fasta_paths <- c(no_fasta_paths, fasta)
            
      }
      
      # TSV PATHS
      data <- paste0(path, "/", "filtered_data.tsv")
      if(file.size(data) > 0){
            #print(paste(data, file.size(data)))
            data_size <- c(data_size, file.size(data))
            data_paths <- c(data_paths, data)
      } else {
            #store missing info
            no_data_paths <- c(no_data_paths, fasta)

      }
}

# record filtered fasta numbers
rec <- paste(length(fasta_paths), "filtered fasta paths found")
print(rec)
write(rec, merge_log)

# record filtwred tsv numbers
rec <- paste(length(data_paths), "filtered metadata paths found")
print(rec)
write(rec, merge_log)

# record if any missing
if(length(no_fasta_paths)>0){
      rec <- paste(length(no_fasta_paths), "paths missing filtered fasta", no_fasta_paths)
      print(rec)
      write(rec, merge_log)
}

if(length(no_data_paths)){
      rec <- paste(length(no_data_paths), "paths missing filtered metadata", no_data_paths)
      print(rec)
      write(rec, merge_log)
}

# get combined size 
total_fasta_size <- sum(fasta_size)
rec <- paste("combined size of all FASTA files:", total_fasta_size, "bytes")
print(rec)
write(rec, merge_log)

total_data_size <- sum(data_size)
rec <- paste("combined size of all TSV files:", total_data_size, "bytes")
print(rec)
write(rec, merge_log)


## MERGE DATA ############################################################################
results_wd <- paste0(start_path, "/merged_filtered")

if(dir.exists(results_wd)){
      setwd(results_wd)
} else{
      dir.create(results_wd)
      setwd(results_wd)
}

## MERGE DATA ############################################################################
time_now <- Sys.time()
rec <- paste("Merging", length(data_paths),"filtered TSV files begins", Sys.time())
print(rec)
write(rec, merge_log)

# merge all filtered_data in data_paths into one merged_data
merge_data(data_paths, "merged_filtered_data.tsv", append=FALSE)
minutes <-  round(as.numeric(difftime(Sys.time(), time_now, units="mins")), digits = 3)
rec <- paste("Merging",length(fasta_paths), "filtered TSV files ends", Sys.time(), ", took", minutes, "mins" )
print(rec)
write(rec, merge_log)


## MERGE FASTA ##########################################################################
time_now <- Sys.time()
rec <- paste("Merging", length(fasta_paths),"filtered FASTA files begins", Sys.time())
print(rec)
write(rec, merge_log)

# merge all filtered_seqs in fasta_paths into one merged_seqs 
merge_fasta(fasta_paths, "merged_filtered_seqs.fasta", append=FALSE)
minutes <-  round(as.numeric(difftime(Sys.time(), time_now, units="mins")), digits = 3)
rec <- paste("Merging",length(fasta_paths), "filtered FASTA files ends", Sys.time(), ", took", minutes, "mins" )
print(rec)
write(rec, merge_log)
close(merge_log)
