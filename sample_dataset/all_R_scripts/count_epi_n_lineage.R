## ---------------------------------------------------------------------------------------

# Script: Count sequences per epiweek and lineage
# Authors: Tamara Alrawahneh (CEPHR, UCD)
# Last edits: 2022/09/14

## ---------------------------------------------------------------------------------------
# generate epiweeks
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
epi_20 <- create_epi_list("2019-12-29", "2021-01-02", epi_20, "20_EPI")
epi_21 <- create_epi_list("2021-01-03", "2022-01-01", epi_21, "21_EPI")
epi_22 <- create_epi_list("2022-01-02", "2022-04-09", epi_22, "22_EPI")
all_epi <- c(epi_20, epi_21, epi_22)
epi_names <- unique(names(all_epi))

# set wd
args <-  commandArgs(trailingOnly=TRUE)
wd <- args[1]
setwd(wd)

# create output dir
results_wd <- paste0(wd, "/epi_lineage_counts")
dir.create(results_wd)


count_epi <- function(data_epi, epi_week_names, output_file){
      # read epiweeks
      count_con <- file(output_file, "w")
      input_epi <- read.csv(data_epi)$Epiweek
      tot_seqs <- length(input_epi)
      write(paste("Epiweek", "Count", "Frequency", sep=","), count_con)
      # count per epiweek and get frequency
      print("counting sequences per epiweek")
      for(epi in epi_week_names){
            count <- length(which(input_epi==epi))
            freq <- count/tot_seqs
            out <- paste(epi,count, freq, sep=",")
            # write to file
            write(out,count_con)
      }
      close(count_con)

}
count_lineage <- function(input_data,output_file){
      count_con <- file(output_file, "w")
      # create a lineage lvec from unique_lineages
      input_lineages <- read.csv(input_data)$Lineage
      unique_lineages <- unique(input_lineages)
      print(paste(length(unique_lineages), "unique lineages in dataset"))
      unique_lineages <- unique_lineages[order(unique_lineages)]
      tot_seqs <- length(input_lineages)
      # count from unique lineages and get frequency
      write(paste("Lineage", "Count", "Frequency", sep=","), count_con)
      print("counting sequences per unique lineage")
      for(lin in unique_lineages){
            count <- length(which(input_lineages==lin))
            freq <- count/tot_seqs
            # write to file
            out <- paste(lin,count, freq, sep=",")
            write(out,count_con)
      }
      close(count_con)

}

input_file <- paste0(wd, "/Epi_data.csv")
count_epi(input_file, epi_names, paste0(results_wd, "/epiweek_count.csv"))
count_lineage(input_file, paste0(results_wd, "/lineage_count.csv"))
