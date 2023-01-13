## ---------------------------------------------------------------------------------------

# Script: Assign epiweeks and extract lineages
# Authors: Tamara Alrawahneh (CEPHR, UCD)
# Last edits: 2022/09/14

## ---------------------------------------------------------------------------------------
#ARGUMENTS
# args[1] is working directory
# args[2] is .fasta file
# args[3] is .tsv file
  
## ---------------------------------------------------------------------------------------

## load functions ########################################################################

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

create_epi_list <- function(start, end, list_name, wk_name){
      start_date <- as.Date(start)
      end_date <- as.Date(end)
      #get_name <- deparse(substitute(list_name))
      
      list_name <- list()
      
      wk_ind <- 0
      name <- paste0(wk_name, "_", as.character(wk_ind))
      date_count <- 0
      
      while(!is.null(date_count)){
            add_date <- start_date + date_count
            list_name[[""]] <- add_date
            date_count <- date_count + 1
            for(date in 1:length(list_name)){
                  #print(names(date_count))
                  names(list_name)[date] <- name
            }
            
            if(add_date == end_date){
                  break
            }
            
      }
      
      start_wk <- 1
      end_wk <- start_wk + 6
      wk_ind <- 1
      
      while(wk_ind <= length(list_name)/7){
            
            names(list_name)[start_wk:end_wk] <- paste0(wk_name, "_", as.character(wk_ind))
            
            start_wk <- start_wk + 7
            end_wk <- start_wk + 6
            wk_ind <- wk_ind + 1
            
      }
      #return(list_name)
      
      #assign(get_name, list_name, envir = globalenv())
      return(list_name)
      
      
}

get_week <- function(input_list, input_date){
      
      if(input_date < input_list[[1]] | input_date > input_list[[length(input_list)]] ){
            get_name <- "no_epi"
      }
      
      for(i in 1:length(input_list)){
            date_here <- input_list[[i]] == input_date
            if(date_here){
                  get_name <- names(input_list)[i]
            }
      }
      return(get_name)
      
}


##########################################################################################
args <-  commandArgs(trailingOnly=TRUE)
wd <- args[1]
setwd(wd)


epi_20 <- create_epi_list("2019-12-29", "2021-01-02", epi_20, "20_EPI")
epi_21 <- create_epi_list("2021-01-03", "2022-01-01", epi_21, "21_EPI")
epi_22 <- create_epi_list("2022-01-02", "2022-04-09", epi_22, "22_EPI")
all_epi <- c(epi_20, epi_21, epi_22)


read_tsv <- file(args[2], "r")
n_seqs <- args[3]


epi_log <- file("Epi_data.csv", "w")
header <- paste("ID", "Collection_Date", "Epiweek", "Lineage", sep =",")
write(header, epi_log)
time_log <- file("epi_log.txt", "w")

start_time <- Sys.time()
rec <- paste("processing collection dates beings, in batches of", n_seqs, "//", start_time)
print(rec)
write(rec, time_log)

#
counter <- 1
batch_no <- paste0("batch ", counter,":")
rec <- paste("processing batch", counter, "begins", Sys.time())
print(rec)
#write(rec, time_log)

#
prev_line <- read_n_data(read_tsv, n = n_seqs, "batch.tsv", previous_line = NULL)
last_batch <- FALSE

while(!is.null(prev_line)|| last_batch){
      time_now <- Sys.time()
      
      data_lines <- readLines("batch.tsv")
      split_lines <- strsplit(data_lines, split="\t", fixed=TRUE)
      for(line in split_lines){
            date_ind <- grep(pattern = "\\d{4}-\\d{2}-\\d{2}", line[4])
            if(length(date_ind)>0){
                  get_date <- line[4]
                  get_id <- line[1]
                  epi_wk <- get_week(all_epi, get_date)
                  lineage <- line[12]
                  output <- paste(get_id, get_date, epi_wk, lineage, sep =",")
                  write(output, epi_log)
            }
      }
      
      minutes <-  round(as.numeric(difftime(Sys.time(), time_now, units="secs")), digits = 3)
      rec <- paste("processing batch", counter, "took", minutes, "secs")
      print(rec)
      write(rec, time_log)
      
      if(last_batch){
            break
      }
      
      counter <- counter + 1
      batch_no <- paste0("batch ", counter,":")
      rec <- paste("processing batch", counter, "begins", Sys.time())
      print(rec)
      
      prev_line <- read_n_data(read_tsv, n = n_seqs, "batch.tsv", previous_line = prev_line)
      
      if(!last_batch && is.null(prev_line)){
            last_batch <- TRUE
      }
}




minutes <-  round(as.numeric(difftime(Sys.time(), start_time, units="mins")), digits = 3)
rec <- paste("processing collection dates took",  minutes, "mins")
print(rec)
write(rec, time_log)


close(read_tsv)
close(epi_log)
close(time_log)
unlink("batch.tsv")


