# Default arguments
default_args <- list(
  obj_nullmodel.file="",
  category = "disruptive_missense",
  chr = 1,
  transcript = "standard",
  gds.path = "",
  gene_name = "",
  rare_maf_cutoff = 0.01,
  QC_label = "annotation/info/QC_label3",
  variant_type = "variant",
  geno_missing_imputation = "mean",  
  Annotation_dir = "annotation/info/FunctionalAnnotation",
  Annotation_name_catalog.file = "/data/Annotation_name_catalog_V2.csv",
  Set_Based_Analysis.file="/data/Set_Based_Analysis_Variant_List_20250224.R",
  seq.path = "/data/sequences.csv",
  Use_annotation_weights = FALSE,
  protein_sequence_selection_strategy = "one-step",
  outfile = "")


num_args = c("chr")

args <- commandArgs(TRUE)
print(args)

args_list <- default_args

for (arg in args) {
  
  key_value <- strsplit(arg, "=")[[1]]
  
  if (length(key_value) == 2) {
    arg_name = sub("--", "", key_value[1])
    
    if (arg_name %in% names(default_args)) {
      if (arg_name %in% num_args) {
        value <- suppressWarnings(as.numeric(key_value[2]))
        if (!is.na(value)) {
          args_list[[arg_name]] <- value
        } else {
          warning(paste("Invalid value: ", arg_name, ". Numeric value required."))
        }
      }  else {
        args_list[[arg_name]] <- key_value[2]
      }
    } else {
      warning(paste("Unknown argument: ", arg_name))
    }
  }
}

print(sessionInfo())
print(args_list)


############## load source code
library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
library(STAAR)
library(STAARpipeline)
library(STAARpipelineSummary)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(data.table)
library(stringr)

source(args_list$Set_Based_Analysis.file)
obj_nullmodel <- get(load(args_list$obj_nullmodel.file))
category <- args_list$category
chr <- args_list$chr
transcript <- args_list$transcript
genofile <- seqOpen(args_list$gds.path)
gene_name <- args_list$gene_name
rare_maf_cutoff <- args_list$rare_maf_cutoff
QC_label <- args_list$QC_label
variant_type <- args_list$variant_type
geno_missing_imputation <- args_list$geno_missing_imputation
Annotation_dir <- args_list$Annotation_dir
Annotation_name_catalog <- read.csv(args_list$Annotation_name_catalog.file)
Use_annotation_weights <- args_list$Use_annotation_weights
seq_file <- read.csv(args_list$seq.path, stringsAsFactors = FALSE)
protein_sequence_selection_strategy <- args_list$protein_sequence_selection_strategy
outfile <- args_list$outfile

seq <- ""

if (transcript == "standard") {
  matching_row <- seq_file[seq_file$gene_symbol == gene_name & seq_file$Canonical == "MANE select", ]
  if (nrow(matching_row) == 0) {
    stop(paste("No MANE select transcript found for gene:", gene_name))
  } else if (nrow(matching_row) > 0) {
    matching_row <- matching_row[1, , drop = FALSE]
  }
  transcript <- matching_row$transcript
  seq <- matching_row$sequence
} else {
  matching_row <- seq_file[seq_file$transcript == transcript & seq_file$gene_symbol == gene_name, ]
  if (nrow(matching_row) == 0) {
    stop(paste("No matching transcript found for transcript:", transcript, "and gene:", gene_name))
  } else if (nrow(matching_row) > 0) {
    warning(paste("Multiple matches found for transcript:", transcript, "and gene:", gene_name, ". Selecting the first one."))
    matching_row <- matching_row[1, , drop = FALSE]
  }
  seq <- matching_row$sequence
}
output_file <- paste0(gene_name, "_", transcript, "_seq.txt")
write.table(seq, file = output_file)

Annotation_name <- c("rs_num","GENCODE.Category","GENCODE.Info","GENCODE.EXONIC.Category","GENCODE.EXONIC.Info","am.pathogenicity","MetaSVM",
                     "GeneHancer","CAGE","DHS",
                     "CADD","LINSIGHT","FATHMM.XF","aPC.EpigeneticActive","aPC.EpigeneticRepressed","aPC.EpigeneticTranscription",
                     "aPC.Conservation","aPC.LocalDiversity","aPC.Mappability","aPC.TF","aPC.Protein")


df <- Gene_Centric_Coding_Info(category=category,chr=chr,genofile=genofile,obj_nullmodel=obj_nullmodel,gene_name=gene_name,known_loci=NULL,
                                          rare_maf_cutoff=rare_maf_cutoff,QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                          Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,Annotation_name=Annotation_name)
write.csv(df,paste0(outfile, gene_name,"_",category,".csv"))
df <- as.data.frame(df)

matches <- grepl(transcript, df$GENCODE.EXONIC.Info)
print(matches)
print(sum(matches))
if (sum(matches) == 0) {
  stop(paste("No variants found for transcript:", transcript))
}
df_subset <- df[matches, ]
print(df_subset)
write.csv(df_subset,paste0(outfile, gene_name,"_",category,"_",transcript,".csv"))

exonic_info <- as.character(df_subset$GENCODE.EXONIC.Info)
result_dict <- list()

for (i in seq_along(exonic_info)) {
  line <- exonic_info[i]
  split_line <- strsplit(line, ",")[[1]]    
  for (acid in split_line) {
    if (grepl(transcript, acid)) {
      matches <- regmatches(acid, regexpr("p\\.[A-Za-z]\\d+[A-Za-z]+", acid))
      if (length(matches) > 0) {
        matches <- sub("^p\\.", "", matches)
        pos <- sub(".*?(\\d+).*", "\\1", matches)
        if (is.null(result_dict[[pos]])) {
          result_dict[[pos]] <- i
        } else {
          result_dict[[pos]] <- c(result_dict[[pos]], i)
        }
      }
    }
  }
}

combination_count <- 1
for (position in names(result_dict)) {
  unique_count <- length(result_dict[[position]])
  combination_count <- combination_count * unique_count
}

cat("Transcript", transcript, "combination count:", combination_count, "\n")


STAAR_Alphafold_one_step_analysis <- function(result_dict, df_subset, genofile, obj_nullmodel,
                                              Annotation_name_catalog, Annotation_name,
                                              gene_name, transcript, seq) {
  

  single_value_dict <- list()
  multi_value_dict <- list()

  for (key in names(result_dict)) {
    if (length(result_dict[[key]]) == 1) {
      single_value_dict[[key]] <- result_dict[[key]]
    } else if (length(result_dict[[key]]) >= 2) {
      multi_value_dict[[key]] <- result_dict[[key]]
    }
  }

  if (length(multi_value_dict) > 0) {
    
    p_values <- list()
    
    for (key in names(multi_value_dict)) {
      values <- multi_value_dict[[key]]
      
      for (val_to_remove in values) {
        all_values <- unlist(multi_value_dict)
        single_values <- unlist(single_value_dict)
        combined_set <- all_values[all_values != val_to_remove]
        combined_set <- c(single_values, combined_set)
        
        filtered_rows <- df_subset[combined_set, ]
        variant_list <- filtered_rows[, 3:6]
        
        results_info <- Set_Based_Analysis_Variant_List(genofile=genofile,variant_list=variant_list,obj_nullmodel=obj_nullmodel,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=2,
                                                QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                                Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                                Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
                                                mask_name="self_defined_variant_list",silent=FALSE)
        
        results_info <- as.data.frame(results_info)
        if (obj_nullmodel$use_SPA == FALSE) {
          p_value <- results_info$`STAAR-O`[[1]]
        } else if (obj_nullmodel$use_SPA == TRUE) {
          p_value <- results_info$`STAAR-B`[[1]]
        }
        
        p_values[[paste(key, val_to_remove, sep = "_")]] <- p_value
      }
    }
    
    for (key in names(multi_value_dict)) {
      key_pattern <- paste0("^", key, "_")
      key_p_values <- p_values[grep(key_pattern, names(p_values))]
      
      if (length(key_p_values) > 0) {
        max_p_value <- max(unlist(key_p_values))
        max_key_value <- names(key_p_values)[which.max(unlist(key_p_values))]
        
        max_value <- as.numeric(sub(key_pattern, "", max_key_value))
        
        single_value_dict[[key]] <- max_value
      }
    }

  }

  all_values <- unlist(single_value_dict)

  filtered_rows <- df_subset[all_values, ]

  output_file <- paste0(gene_name, "_", transcript, "_onestep_variant.txt")
  write.table(filtered_rows, file = output_file, sep = "\t", row.names = FALSE, quote = TRUE)

  variant_list <- filtered_rows[, 3:6]

  results_info <- Set_Based_Analysis_Variant_List(genofile=genofile,variant_list=variant_list,obj_nullmodel=obj_nullmodel,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=2,
                                            QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                            Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                            Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
                                            mask_name="self_defined_variant_list",silent=FALSE)
  results_info <- as.data.frame(results_info)

  if (obj_nullmodel$use_SPA == FALSE) {
    p_value <- results_info$`STAAR-O`[[1]]
  } else if (obj_nullmodel$use_SPA == TRUE) {
    p_value <- results_info$`STAAR-B`[[1]]
  }
  
  combined_df <- do.call(cbind, results_info)

  output_file <- paste0(gene_name, "_", transcript, "_onestep.txt")
  write.table(combined_df, file = output_file, sep = "\t", row.names = FALSE, quote = TRUE)

  exonic_info <- filtered_rows$GENCODE.EXONIC.Info

  aa_change <- c()

  for (info in exonic_info) {
    variants <- unlist(strsplit(info, ","))
    for (variant in variants) {
      if (grepl(transcript, variant)) {
        match <- str_match(variant, "p\\.([A-Z][0-9]+[A-Z]+)")[,2]
        if (!is.na(match)) {
          aa_change <- c(aa_change, match)
        }
      }
    }
  }

  for (change in aa_change) {
    match <- str_match(change, "([A-Z])([0-9]+)([A-Z]+)")
    if (!is.na(match[1])) {
      original_aa <- match[2]
      position <- as.integer(match[3])
      mutated_aa <- match[4]
      
      if (substr(seq, position, position) == original_aa) {
        substr(seq, position, position) <- mutated_aa
      } else {
        cat(sprintf("Position %d: original amino acid %s does not match %s. No mutation applied.\n",
                    position, substr(seq, position, position), original_aa))
      }
    }
  }

  output_file <- paste0(gene_name, "_", transcript, "_onestep_seq.txt")
  write.table(seq, file = output_file)

}


STAAR_Alphafold_stepwise_analysis <- function(result_dict, df_subset, genofile, obj_nullmodel,
                                              Annotation_name_catalog, Annotation_name,
                                              gene_name, transcript, seq) {

  single_value_dict <- list()
  multi_value_dict <- list()

  for (key in names(result_dict)) {
    if (length(result_dict[[key]]) == 1) {
      single_value_dict[[key]] <- result_dict[[key]]
    } else if (length(result_dict[[key]]) >= 2) {
      multi_value_dict[[key]] <- result_dict[[key]]
    }
  }

  while (length(multi_value_dict) > 0) {
    
    p_values <- list()
    
    for (key in names(multi_value_dict)) {
      values <- multi_value_dict[[key]]
      
      for (val_to_remove in values) {
        all_values <- unlist(multi_value_dict)
        single_values <- unlist(single_value_dict)
        combined_set <- all_values[all_values != val_to_remove]
        combined_set <- c(single_values, combined_set)
        filtered_rows <- df_subset[combined_set, ]
        variant_list <- filtered_rows[, 3:6]
        
        results_info <- Set_Based_Analysis_Variant_List(genofile=genofile,variant_list=variant_list,obj_nullmodel=obj_nullmodel,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=2,
                                                QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                                Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                                Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
                                                mask_name="self_defined_variant_list",silent=FALSE)
        
        results_info <- as.data.frame(results_info)
        if (obj_nullmodel$use_SPA == FALSE) {
          p_value <- results_info$`STAAR-O`[[1]]
        } else if (obj_nullmodel$use_SPA == TRUE) {
          p_value <- results_info$`STAAR-B`[[1]]
        }
        
        p_values[[paste(key, val_to_remove, sep = "_")]] <- p_value
      }
    }
    
    max_p_value <- max(unlist(p_values))
    max_key_value <- names(p_values)[which.max(unlist(p_values))]
    
    max_key <- strsplit(max_key_value, "_")[[1]][1]
    max_value <- as.numeric(strsplit(max_key_value, "_")[[1]][2])
    
    single_value_dict[[max_key]] <- max_value
    
    multi_value_dict[[max_key]] <- NULL
 }

  all_values <- unlist(single_value_dict)

  filtered_rows <- df_subset[all_values, ]

  output_file <- paste0(gene_name, "_", transcript, "_stepwise_variant.txt")
  write.table(filtered_rows, file = output_file, sep = "\t", row.names = FALSE, quote = TRUE)

  variant_list <- filtered_rows[, 3:6]

  results_info <- Set_Based_Analysis_Variant_List(genofile=genofile,variant_list=variant_list,obj_nullmodel=obj_nullmodel,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=2,
                                            QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                            Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                            Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
                                            mask_name="self_defined_variant_list",silent=FALSE)
  results_info <- as.data.frame(results_info)
  if (obj_nullmodel$use_SPA == FALSE) {
    p_value <- results_info$`STAAR-O`[[1]]
  } else if (obj_nullmodel$use_SPA == TRUE) {
    p_value <- results_info$`STAAR-B`[[1]]
  }

  combined_df <- do.call(cbind, results_info)
  output_file <- paste0(gene_name, "_", transcript, "_stepwise.txt")
  write.table(combined_df, file = output_file, sep = "\t", row.names = FALSE, quote = TRUE)

  exonic_info <- filtered_rows$GENCODE.EXONIC.Info

  aa_change <- c()

  for (info in exonic_info) {
    variants <- unlist(strsplit(info, ","))
    for (variant in variants) {
      if (grepl(transcript, variant)) {
        match <- str_match(variant, "p\\.([A-Z][0-9]+[A-Z]+)")[,2]
        if (!is.na(match)) {
          aa_change <- c(aa_change, match)
        }
      }
    }
  }

  for (change in aa_change) {
    match <- str_match(change, "([A-Z])([0-9]+)([A-Z]+)")
    if (!is.na(match[1])) {
      original_aa <- match[2]
      position <- as.integer(match[3])
      mutated_aa <- match[4]
      
      if (substr(seq, position, position) == original_aa) {
        substr(seq, position, position) <- mutated_aa
      } else {
        cat(sprintf("Position %d: original amino acid %s does not match %s. No mutation applied.\n",
                    position, substr(seq, position, position), original_aa))
      }
    }
  }

  output_file <- paste0(gene_name, "_", transcript, "_stepwise_seq.txt")
  write.table(seq, file = output_file)
}


STAAR_Alphafold_all_combinations_analysis <- function(result_dict, df_subset, genofile, obj_nullmodel,
                                                     Annotation_name_catalog, Annotation_name,
                                                     gene_name, transcript, seq) {
  combination_list <- expand.grid(result_dict)
  p_values <- list()
  
  for (i in 1:nrow(combination_list)) {
    current_combination <- as.numeric(combination_list[i, ])
    filtered_rows <- df_subset[current_combination, ]
    variant_list <- filtered_rows[, 3:6]
    
    results_info <- Set_Based_Analysis_Variant_List(genofile=genofile,variant_list=variant_list,obj_nullmodel=obj_nullmodel,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=2,
                                            QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                            Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                            Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
                                            mask_name="self_defined_variant_list",silent=FALSE)
    results_info <- as.data.frame(results_info)
    if (obj_nullmodel$use_SPA == FALSE) {
      p_value <- results_info$`STAAR-O`[[1]]
    } else if (obj_nullmodel$use_SPA == TRUE) {
      p_value <- results_info$`STAAR-B`[[1]]
    }
    p_values[[paste(current_combination, collapse = "_")]] <- p_value
    }
  
  min_p_value <- min(unlist(p_values))
  min_key_value <- names(p_values)[which.min(unlist(p_values))]
  min_combination <- as.numeric(strsplit(min_key_value, "_")[[1]])
  
  filtered_rows <- df_subset[min_combination, ]
  output_file <- paste0(gene_name, "_", transcript, "_all_combinations_variant.txt")
  write.table(filtered_rows, file = output_file, sep = "\t", row.names = FALSE, quote = TRUE)
  variant_list <- filtered_rows[, 3:6]
  results_info <- Set_Based_Analysis_Variant_List(genofile=genofile,variant_list=variant_list,obj_nullmodel=obj_nullmodel,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=2,
                                                 QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                                 Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                                 Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
                                                 mask_name="self_defined_variant_list",silent=FALSE)
  results_info <- as.data.frame(results_info)
  if (obj_nullmodel$use_SPA == FALSE) {
    p_value <- results_info$`STAAR-O`[[1]]
  } else if (obj_nullmodel$use_SPA == TRUE) {
    p_value <- results_info$`STAAR-B`[[1]]
  }
  
  combined_df <- do.call(cbind, results_info)
  output_file <- paste0(gene_name, "_", transcript, "_allcombinations.txt")
  write.table(combined_df, file = output_file, sep = "\t", row.names = FALSE, quote = TRUE)
  
  exonic_info <- filtered_rows$GENCODE.EXONIC.Info

  aa_change <- c()

  for (info in exonic_info) {
    variants <- unlist(strsplit(info, ","))
    for (variant in variants) {
      if (grepl(transcript, variant)) {
        match <- str_match(variant, "p\\.([A-Z][0-9]+[A-Z]+)")[,2]
        if (!is.na(match)) {
          aa_change <- c(aa_change, match)
        }
      }
    }
  }

  for (change in aa_change) {
    match <- str_match(change, "([A-Z])([0-9]+)([A-Z]+)")
    if (!is.na(match[1])) {
      original_aa <- match[2]
      position <- as.integer(match[3])
      mutated_aa <- match[4]
      
      if (substr(seq, position, position) == original_aa) {
        substr(seq, position, position) <- mutated_aa
      } else {
        cat(sprintf("Position %d: original amino acid %s does not match %s. No mutation applied.\n",
                    position, substr(seq, position, position), original_aa))
      }
    }
  }

  output_file <- paste0(gene_name, "_", transcript, "_all_combinations_seq.txt")
  write.table(seq, file = output_file)
}

if (protein_sequence_selection_strategy == "one-step") {
  STAAR_Alphafold_one_step_analysis(
    result_dict = result_dict,
    df_subset = df_subset,
    genofile = genofile,
    obj_nullmodel = obj_nullmodel,
    Annotation_name_catalog = Annotation_name_catalog,
    Annotation_name = Annotation_name,
    gene_name = gene_name,
    transcript = transcript,
    seq = seq
  )
} else if (protein_sequence_selection_strategy == "stepwise") {
  STAAR_Alphafold_stepwise_analysis(
    result_dict = result_dict,
    df_subset = df_subset,
    genofile = genofile,
    obj_nullmodel = obj_nullmodel,
    Annotation_name_catalog = Annotation_name_catalog,
    Annotation_name = Annotation_name,
    gene_name = gene_name,
    transcript = transcript,
    seq = seq
  )
} else if (protein_sequence_selection_strategy == "best-subset-selection") {
  STAAR_Alphafold_all_combinations_analysis(
    result_dict = result_dict,
    df_subset = df_subset,
    genofile = genofile,
    obj_nullmodel = obj_nullmodel,
    Annotation_name_catalog = Annotation_name_catalog,
    Annotation_name = Annotation_name,
    gene_name = gene_name,
    transcript = transcript,
    seq = seq
  )
}
