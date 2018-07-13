#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
term <- as.character(args[1])
variant_type <- as.character(args[2])
significance <- as.character(args[3])
risk_factor <- as.character(args[4])
options(scipen = 999)

Clinvar_Retrieval_Pipe <- function(term, variant_type = NULL, significance = 1, risk_factor = NULL){
  # term: 
  # Key words you are using to search Clinvar. Please use logical statements "AND", "OR", and "NOT" 
  # to seperate your terms. Paretheses can be used for multiple logical statements such as:
  # '((HIV AND Black Man) OR (Gonorrhea White Woman)) NOT Georgia
  # variant_type: 
  # 1 = "single nucleotide variant", 2 = "Duplication" (Will add more as they arise),
  # 3 = "Deletion"
  # significance:
  # 1 = Not subsetted on significance (default)   
  # 2 = data %in% c("Benign", "Likely benign", "Uncertain significance", "Likely pathogenic", "Pathogenic")
  # 3 = data %in% c("Uncertain significance", "Likely pathogenic", "Pathogenic")
  # 4 = data %in% c("Likely pathogenic", "Pathogenic")
  # risk_factor:
  # c("Yes", "yes", "YES") = you are only interested in variants that have been labeled as a risk factor
  # to some degree.
  # c("No", "no", "NO") = you are not interested in variants that have been labeled as a risk factor
  # to any degree.
  # NULL = Bring back risk factors and non risk factors (default)
  
  library(httr)
  library(jsonlite)
  library(dplyr)
  library(rentrez)
  
  # function to determine character postion of short string inside a longer string
  string_finder <- function(full_phrase_vector, short_phrase, group = 1){
    # group: if the str of intrest is repeated more than once, put the repeat number in the string. 
    placement <- matrix(".", length(full_phrase_vector), 2)
    for(j in 1:length(full_phrase_vector)){
      #   if(length(grep(short_phrase, full_phrase_vector[j])) == 0){
      #     next
      #   }
      # splitting string into single characters
      full_phrase <- substring(full_phrase_vector[j], seq(1, nchar(full_phrase_vector[j]), 1), seq(1, nchar(full_phrase_vector[j]), 1))
      # finding location last element in short phrase
      for(i in 1:nchar(short_phrase)){
        if(i == 1){
          n <- which(full_phrase == substr(short_phrase, i, i))
        }else{
          x <- which(full_phrase == substr(short_phrase, i, i))
          n <- x[x %in% I(n + 1)]
        }
      }
      
      if(length(n) == 0){
        next
      }else{
        if(length(n) > 1){
          n2 <- n[group]
          # finding first location of short phrase
          n1 <- n[group] - nchar(short_phrase) + 1
        }else{
          n2 <- n
          # finding first location of short phrase
          n1 <- n - nchar(short_phrase) + 1
        }
        
      }
      
      placement[j, 1] <- n1
      placement[j, 2] <- n2
      
    }
    # returns the placement of the first/last integer in the phrase
    placement <- data.frame(placement)
    colnames(placement) <- c("Start", "Finish")
    placement[, 1] <- as.numeric(as.character(placement[, 1]))
    placement[, 2] <- as.numeric(as.character(placement[, 2]))
    return(placement)
  }
  
  
  
  # Doing converstions for relevant inputs
  variant_label <- function(variant_type){
    if(is.null(variant_type)){
      return(NULL)
    }else{
      if(variant_type == 1){
        variant_type <- "single nucleotide variant"
      }else{
        if(variant_type == 2){
          variant_type <- "Duplication"
        }else{
          if(variant_type == 3){
            variant_type <- "Deletion"
          }else{
            variant_type <- NULL
            print("Invalid number for variant type. Will not filter for it.")
          }
        }
      }
    }
    
    return(variant_type)
  }
  variant_type <- variant_label(variant_type)
  
  
  # Searching the database
  r_search <- entrez_search(db = "clinvar", term = term, retmax = 50)
  
  # Returning, arguably the most important thing, the ids 
  #r_search$ids
  
  # Getting summary from ids
  taxize_summ <- entrez_summary(db = "clinvar", id = r_search$ids)
  
  # Function to create full dataset for all terms related to search words
  for(i in 1:length(taxize_summ)){
    # Function to extract inforation from messy_data
    extracter_clinvar <- function(identifier, group, skip = "no"){
      # Identifer: the minimum number of character up until the numeric information you are looking for
      # inlcuidng " " 
      # group: the reccorrence postiion of the phrase you are looking for
      # skip: if yes, dont start search until after the first seen ",  " (2 spacse are intentional)
      
      n <- as.numeric(string_finder(full_phrase_vector = messy_data, short_phrase = identifier, group = group)[2]) + 1
      
      # Skips to next comma if needed
      if(skip %in% c("No", "no", "NO")){
        # Do nothing
      }else{
        if(skip %in% c("Yes", "yes", "YES")){
          skip_list <- substring(messy_data, seq(1, nchar(messy_data), 1), seq(1, nchar(messy_data), 1))
          loop <- 1
          while(loop == 1){
            if(skip_list[n] != ","){
              n = n + 1
            }else{
              n = n + 2
              loop <- 0
            }
          }
        }else{
          print(paste("Check spelling on skip = '", skip, "'", sep = ""))
        }
      }
      
      result <- NULL
      while(!is.na(suppressWarnings(as.numeric(messy_data_long[n]))) | 
            messy_data_long[n] == "X"){
        
        result <- paste(result, messy_data_long[n], sep = "")
        n <- n + 1
      }
      
      return(result)
    }
    
    # Organizing data so it's easier to access
    messy_data <- as.character(taxize_summ[[i]]$variation_set["variation_loc"])
    messy_data_long <- substring(messy_data, seq(1, nchar(messy_data), 1), seq(1, nchar(messy_data), 1))
    messy_data_long <- messy_data_long[-which(messy_data_long == "\"")]
    messy_data <- paste(messy_data_long, collapse = "")
    if(nchar(messy_data) == 0){
      next
    }
    
    # Defining build
    build <- "GRCh37"
    
    # Data pulled that include "current" have more characters than those that do not
    # since we are using the 37 build, "current" refers to the 38 build and we are looking for
    # the most recent "previous"
    check_1 <- grep("current", messy_data)
    
    if(length(check_1) == 0){
      skip = "No"
    }else{
      c_spot <- string_finder(full_phrase_vector = messy_data, short_phrase = "current", group = 1)[2]
      p_spot1 <- string_finder(full_phrase_vector = messy_data, short_phrase = "previous", group = 1)[1]
      
      if(p_spot1 > c_spot){ # Checks if the previous build (the one we are looking for) comes after the current build
        skip = "Yes"
      }else{
        print("Have Octavious check if 'postion' comes before 'chromosome' in the data you read in.")
      }
      
    }
    
    # Finding the Chromosome
    chr <- extracter_clinvar("chr = ", 1, skip = skip)
    
    # Finding the starts/stops
    position_start <- extracter_clinvar("start = ", 1, skip = skip)
    position_stop <- extracter_clinvar("stop = ", 1, skip = skip)
    
    position_inner_start <- extracter_clinvar("inner_start = ", 1, skip = skip)
    position_inner_stop <- extracter_clinvar("inner_stop = ", 1, skip = skip)
    
    position_outer_start <- extracter_clinvar("outer_start = ", 1, skip = skip)
    position_outer_stop <- extracter_clinvar("outer_stop = ", 1, skip = skip)
    
    # Finding Gene and strand
    if(length(taxize_summ[[i]]$genes) == 0){
      gene_strand <- NA
    }else{
      gene_strand <- paste(taxize_summ[[i]]$genes[, "symbol"]," (", taxize_summ[[i]]$genes[, "strand"], ")" , collapse = ", ", sep = "")
    }
    
    # Grabbing conditions 
    conditions <- paste(taxize_summ[[i]]$trait_set["trait_name"], collapse = ", ")
    remove <- c("\"", "\\)", "c\\(")
    for(i1 in 1:length(remove)){
      conditions <- gsub(remove[i1], "", conditions)
    }
    
    # Grabbing significance
    significance_ind <- as.character(taxize_summ[[i]]$clinical_significance["description"])
    
    # Grabbing variant type
    vt <- as.character(taxize_summ[[i]]$variation_set["variant_type"])
    for(i1 in 1:length(remove)){
      vt <- gsub(remove[i1], "", vt)
    }
    
    # building dataframe
    entries <- c(build, paste(chr, position_start, sep = ":"), chr, position_start, position_stop,
                 gene_strand, conditions, significance_ind, vt)
    if(i == 1){
      df <- matrix(entries, 1, length(entries))
      colnames(df) <- c("Build", "CHR:Position", "CHR", "Position_Start", "Position_Stop", "Gene_Strand",
                        "Conditions", "Significance", "Variant_Type")
    }else{
      df1 <- matrix(entries, 1, length(entries))
      colnames(df1) <- c("Build", "CHR:Position", "CHR", "Position_Start", "Position_Stop", "Gene_Strand",
                         "Conditions", "Significance", "Variant_Type")
      
      df <- rbind(df, df1)
    }
  }
  
  # Subsetting to variant type
  if(!is.null(variant_type)){
    df <- df[grep(variant_type, df[, "Variant_Type"]), ]
  }
  
  # Subsetting to significance
  if(significance == 2){
    s_n <- unique(sort(unlist(lapply(c("Benign", "Likely benign", "Uncertain significance", "uncertain significance", "Likely pathogenic", "Pathogenic"), grep, 
                                     df[, "Significance"]))))
    df <- df[s_n, ]
  }else{
    if(significance == 3){
      s_n <- unique(sort(unlist(lapply(c("uncertain significance", "Uncertain significance", "Likely pathogenic", "Pathogenic"), grep, 
                                       df[, "Significance"]))))
      df <- df[s_n, ]
    }else{
      if(significance == 4){
        s_n <- unique(sort(unlist(lapply(c("Likely pathogenic", "Pathogenic"), grep, 
                                         df[, "Significance"]))))
        df <- df[s_n, ]
      }else{
        if(significance != 1){
          significance = 1
          print(paste("significance = ", significance, "is not a valid entry. Will default to 1."))
        }
      }
    }
  }
  
  # Subsetting to risk_factor
  if(!is.null(risk_factor)){
    if(risk_factor %in% c("Yes", "yes", "YES")){
      n_r <- unique(c(grep("risk factor", df[, "Significance"]), grep("Risk factor", df[, "Significance"])))
      df <- df[n_r, ]
    }else{
      if(risk_factor %in% c("No", "no", "NO")){
        n_r <- unique(c(grep("risk factor", df[, "Significance"]), grep("Risk factor", df[, "Significance"])))
        df <- df[-n_r, ]
      }else{
        risk_factor <- NULL
        print(paste("risk_factor = ", risk_factor, "is not a valid entry (check spelling). Will default to NULL."))
      }
    }
  }
  
  return(df)
}

# Outputs Results
results <- Clinvar_Retrieval_Pipe(term, variant_type, significance, risk_factor)
write.csv(results, "Clinvar_variables.tsv", row.names = FALSE, quote = FALSE)
