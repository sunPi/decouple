#---- Requirements ----
# library("DESeq2")
# library("S4Vectors")
pkgs <- c("docopt", "decoupleR", "pheat0map", "here", "readxl", 
          "ggplot2", "org.Hs.eg.db", "AnnotationDbi", "mesocore", "SummarizedExperiment")

# "AUCell", "doMC", "ranger", "viper", "progeny", "dorothea", "OmnipathR"
suppressMessages(mesocore::handleRequirements(pkgs))

#---- Source Functions ----
# Get the path to the script being executed
args <- commandArgs(trailingOnly = FALSE)
script.path <- normalizePath(sub("--file=", "", args[grep("--file=", args)]))
script.dir <- dirname(script.path)

# Source the functions.R script located in ./R/ relative to the script
source(file.path(script.dir, "functions.R"))

#---- Header ----
"Mesothelioma AI Pipeline - Decouple Transcription Factors

Usage: decouple_TF.R [options]

Options:
  -h --help                     Show this screen.
  -d --data_path=<PATH>         Path to the counts data set (RNA-seq data).
  -e --experiment_design=<PATH> Path to the table with experiment design.
  -o --outfolder=<PATH>         Folder where the results are saved.
  -r --reference=<STRING>       In contrast experiments, which is the baseline (reference) mapping.
  -m --methods=<STRING>         If specified, runs decouple only on the selected method. Otherwise runs all methods.
  -s --deseq_exp=<PATH>         Specify a pre-ranked list of genes for GSEA.
  -v --verbose=<BOOLEAN>        If set to 1, prints messages verbously.
  -V --version
"-> doc

#---- Arguments ----
arguments <- docopt(doc, quoted_args = TRUE, help = TRUE)
print(arguments)

cts.path  <- normalizePath(arguments$data_path)
exp.path  <- normalizePath(arguments$experiment_design)
outfolder <- normalizePath(arguments$outfolder, mustWork = F)
reference <- arguments$reference
methods   <- arguments$methods
deseq.res <- readFile(arguments$deseq_exp)

if(as.integer(arguments$verbose) == 1){
  verbose <- T
} else{
  verbose <- F
}
if(arguments$methods == "NULL"){
  methods <- NULL
} 

if(is.null(methods)){
  message("No specific method selected, running decouple on all available methods...")
  
  methods  <- c("fgsea", "gsva", "mdt", "mlm", "ora", "udt", "viper", "wmean", "wsum")
  do.split <- T
  
} else{
  methods <- arguments$methods
  do.split <- F
}

cts       <- readFile(cts.path)
exp       <- readFile(exp.path)
dir.create(outfolder, T, T) # Creates a directory based on the outfoolder path provided

# Load Arguments
# cts       <- mesocore::readFile("/home/jr453/Documents/Cancer_Studies_PhD/Studies/Study_MEDUSA_Tumour_Volumes/data/filtered/M54_cts.csv")
# exp       <- mesocore::readFile("/home/jr453/Documents/Cancer_Studies_PhD/Studies/Study_MEDUSA_Tumour_Volumes/data/filtered/M54_exp.csv")
# deseq.res <- mesocore::readFile("/home/jr453/Documents/Cancer_Studies_PhD/Studies/Study_MEDUSA_Tumour_Volumes/results/GSEA/h.all.v2024.1.Hs.symbols/Q4_DEG.xlsx")
# outfolder <- "/home/jr453/Documents/Cancer_Studies_PhD/Studies/Study_MEDUSA_Tumour_Volumes/decoupleR/TEST"
# reference <- "Q1"
# top.n <- 500
# methods <- NULL

#----- Main -----
### Preprocess
processed <- preprocess(cts, exp, reference)

### Run Decouple for TF Inference
decouple_TF(processed$normalized.counts, processed$deg, methods, do.split, outfolder)

