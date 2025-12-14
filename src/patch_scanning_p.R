
# Default arguments
default_args <- list(
  Set_Based_Analysis.file="/data/Set_Based_Analysis_Variant_List_20250224.R",
  obj_nullmodel.file="",
  variant_file="",
  rare_maf_cutoff = 0.01,
  Annotation_dir = "annotation/info/FunctionalAnnotation",
  Annotation_name_catalog.file = "/data/Annotation_name_catalog_V2.csv",
  QC_label = "annotation/info/QC_label3",
  geno_missing_imputation = "mean",
  Use_annotation_weights = FALSE,
  chr = 1,
  gene_name = "PCSK9",
  category = "disruptive_missense",
  variant_type = "variant",
  gds.path = "",
  outfile = "")

num_args = c("chr")
args <- commandArgs(TRUE)

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

############## load source code
suppressPackageStartupMessages({
library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
library(STAAR)
library(STAARpipeline)
library(STAARpipelineSummary)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
})

source(args_list$Set_Based_Analysis.file)
combination_data <- read.table(file("stdin"), sep="\t", header=TRUE, stringsAsFactors=FALSE)
obj_nullmodel <- get(load(args_list$obj_nullmodel.file))
variant_file <- read.table(args_list$variant_file, sep = "\t", stringsAsFactors = FALSE, header = TRUE)
Annotation_dir <- args_list$Annotation_dir
Annotation_name_catalog <- read.csv(args_list$Annotation_name_catalog.file)
QC_label <- args_list$QC_label
geno_missing_imputation <- args_list$geno_missing_imputation
Use_annotation_weights <- args_list$Use_annotation_weights
rare_maf_cutoff <- args_list$rare_maf_cutoff
gene_name <- args_list$gene_name
chr <- args_list$chr
genofile <- seqOpen(args_list$gds.path)
category <- args_list$category
variant_type <- args_list$variant_type

Annotation_name <- c("rs_num","GENCODE.Category","GENCODE.Info","GENCODE.EXONIC.Category","GENCODE.EXONIC.Info","MetaSVM",
                     "GeneHancer","CAGE","DHS",
                     "CADD","LINSIGHT","FATHMM.XF","aPC.EpigeneticActive","aPC.EpigeneticRepressed","aPC.EpigeneticTranscription",
                     "aPC.Conservation","aPC.LocalDiversity","aPC.Mappability","aPC.TF","aPC.Protein")


all_results <- list()
for (i in 1:nrow(combination_data)) {

  values <- as.numeric(unlist(strsplit(combination_data[i,3], ",")))
  
  filtered_rows <- variant_file[values, ]

  if (nrow(filtered_rows) <= 1) {
    next
  }

  
  variant_list <- filtered_rows[, 3:6]
  print(variant_list)


  results_info <- tryCatch({Set_Based_Analysis_Variant_List(genofile=genofile,variant_list=variant_list,obj_nullmodel=obj_nullmodel,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=2,
                                              QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                              Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                              Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
                                              mask_name="self_defined_variant_list",silent=FALSE)

                            }, error=function(e) {
                              if (grepl("Number of rare variant in the set is less than 2", e$message)) {
                                message("Skipping combination ", i, " due to insufficient rare variants.")
                                return(NULL)
                              } else {
                                stop(e)
                              }
                            })

  if (is.null(results_info) || nrow(as.data.frame(results_info)) == 0) {next}

  combination_row <- combination_data[i, 1:3]
  results_info <- as.data.frame(results_info)
  combined_row <- cbind(combination_row, results_info)

  all_results[[i]] <- combined_row

}

if (obj_nullmodel$use_SPA == FALSE) {
  final_results <- do.call(rbind, lapply(all_results, function(df) {
    data.frame(
      `central_id` = as.integer(df$central_id), 
      `neighbors` = as.character(df$neighbors),
      `variant_indices` = as.character(df$variant_indices), 
      Chr = as.character(df$Chr[[1]]),
      Category = as.character(df$Category[[1]]),
      `#SNV` = as.integer(df$`#SNV`[[1]]), 
      cMAC = as.numeric(df$cMAC[[1]]),  
      `STAAR-O` = as.numeric(df$`STAAR-O`[[1]])  
    )
  }))
} else if (obj_nullmodel$use_SPA == TRUE) {
  final_results <- do.call(rbind, lapply(all_results, function(df) {
    data.frame(
      `central_id` = as.integer(df$central_id), 
      `neighbors` = as.character(df$neighbors),
      `variant_indices` = as.character(df$variant_indices), 
      Chr = as.character(df$Chr[[1]]),
      Category = as.character(df$Category[[1]]),
      `#SNV` = as.integer(df$`#SNV`[[1]]), 
      cMAC = as.numeric(df$cMAC[[1]]),  
      `STAAR-B` = as.numeric(df$`STAAR-B`[[1]])  
    )
  }))
}

write.table(final_results, file = stdout(), sep = "\t", row.names = FALSE, quote = FALSE)

