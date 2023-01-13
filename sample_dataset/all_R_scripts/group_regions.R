## ---------------------------------------------------------------------------------------

# Script: Group merged sequences by region, processes merged in batches of 1000
# Authors: Tamara Alrawahneh (CEPHR, UCD)
# Last edits: 2022/09/14

## ---------------------------------------------------------------------------------------
#ARGUMENTS
# args[1] is working directory
# args[2] is .fasta file
# args[3] is .tsv file


## ---------------------------------------------------------------------------------------

args <-  commandArgs(trailingOnly=TRUE)



## FUNCTIONS -----------------------------------------------------------------------------

# read_n_data reads n lines, returns last line, pattern indicating line is "hCov" 
# use read_n_data to read metadata TSV file
read_n_data <- function(reading_control, n = 1000, output_file, previous_line = NULL){ 
      
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
read_n_seqs <- function(reading_control, n = 1000, output_file, previous_line = NULL){ 
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


# extract metadata subset with sepecic pattern
extract_specific_data <- function(input_data_file, n = 1000, output_file, pattern_name, append=TRUE){ 
      
      reading_control <- file(input_data_file, "r")
      
      if(append == TRUE){
            writing_control <- file(output_file, "a")
      } else {
            writing_control <- file(output_file, "w") 
      }
      
      my_line <-  "" 
      seq_count <-  0 #initiate sequence counter
      
      write_to_file <-  FALSE #initiate writing to file as FALSE
      
      store_ids <- c()
      
      while(!is.null(my_line) && seq_count < n){ 
            #print(seq_count) # for observation purposes
            
            if(length(grep(pattern = "hCoV", my_line, ignore.case = T)) > 0){
                  split_lines <- strsplit(my_line, split="\t", fixed=TRUE)
                  for(line in split_lines){
                        index <-  grep(pattern = pattern_name, line[5])
                        if(length(index) > 0){
                              #print(paste(index, my_line, sep="  "))
                              write_to_file <- TRUE
                              seq_count <-  seq_count + 1
                              #store <- c(store, index)
                              store_ids <- c(store_ids, line[1])
                        } else {
                              
                              write_to_file <- FALSE
                        }
                        if(write_to_file){
                              write(my_line, file= writing_control)
                        }
                        
                  }
                  
            }
            my_line <-  readLines(con = reading_control, n = 1)
            
            if(length(my_line) == 0){
                  break
            }
            
      }
      close(reading_control)
      close(writing_control)
      return(store_ids)
      
}

# extract metadata subset with sepecic pattern
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
                        #pattern.vec = pattern.vec[-index]
                        write_to_file <- TRUE
                        seq_count <-  seq_count + 1
                        
                        
                  #} else if(seq_count == length(pattern.vec)){
                        #break
                        
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



##########################################################################################

start_path <- args[1]
setwd(start_path)
start_time <- Sys.time()

log_file <- file("region_log.txt", "w")

# open connections
read_metadata <- file(args[2], "r")
header <- readLines(args[2], n=1)


regions <- c("Africa", "Asia", "Europe", "North_America", "Oceania", "South_America")
for(region in regions){
      
      region_dir <- paste0(start_path, "/", region)
      
      if(dir.exists(region_dir)){
            unlink(region_dir, recursive=T)
            dir.create(region)
      } else {
            dir.create(region_dir)
      }
      
      region_tsv <- paste0(region_dir, "/", region, "_data.tsv")
      write(header,region_tsv) 
      
}

read_fasta <- file(args[3], "r")


counter <- 1
batch_no <- paste0("batch ", counter,":")
rec <- paste("processing batch", counter, "begins", Sys.time())
print(rec)
write(rec, log_file)


prev_line <- read_n_data(read_metadata, n = 1000, "batch.tsv", previous_line = NULL)
prev_seq <- read_n_seqs(read_fasta, n = 1000, "batch.fasta", previous_line = NULL)


last_batch <- FALSE
while(!is.null(prev_line) && !is.null(prev_seq)|| last_batch){
      
      # loop through regions
      for(region in regions){
            # get region dir, tsv name, fasta na,e
            region_dir <- paste0(start_path, "/", region)
            region_tsv <- paste0(region_dir, "/", region, "_data.tsv")
            region_fasta <- paste0(region_dir, "/", region, "_seqs.fasta")
            
            print(paste("Extracting sequences from", region))
            
            region_patt <- gsub("_", " ", region)
            # extract metadata by region name, see function, selects row 5 where region is specified
            region_vec <- extract_specific_data("batch.tsv", n=1000, region_tsv, pattern_name = region_patt, append=TRUE)
            # extract sequences region name, see function, selects row 5 where region is specified
            extract_fasta_from_vec("batch.fasta", region_fasta, region_vec, append = TRUE)
      }

      
      if(last_batch){
            break
      }
      
      counter <- counter + 1
      batch_no <- paste0("batch ", counter,":")
      rec <- paste("processing batch", counter, "begins", Sys.time())
      print(rec)
      write(rec, log_file)
      
      
      prev_line <- read_n_data(read_metadata, n = 1000, "batch.tsv", previous_line = prev_line)
      prev_seq <- read_n_seqs(read_fasta, n = 1000, "batch.fasta", previous_line = prev_seq)
      
      if(!last_batch && is.null(prev_line) && is.null(prev_seq)){
            last_batch <- TRUE
      }
}



end_time <- Sys.time()

rec <- paste("Grouping by region complete", Sys.time())
print(rec)
write(rec, log_file)


minutes <-  round(as.numeric(difftime(end_time, start_time, units="mins")), digits = 3)
rec <- paste("Grouping by region took:", minutes, "mins")
print(rec)
write(rec, log_file)



close(read_metadata)
close(read_fasta)
close(log_file)

unlink("batch.tsv")
unlink("batch.fasta")
