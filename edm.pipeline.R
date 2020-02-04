require("ExomeDepth")
require("dplyr")
require("tools")
require("yaml")
require("EDM")
require("clustermq")
require("utils")
require("data.table")

library(data.table)
library(utils)
library(ExomeDepth)
library(optparse,quietly = T)
library(tools,quietly = T)
library(yaml,quietly = T)
library(purrr,quietly = T)
library(dplyr,quietly = T)
library(clustermq)
library(EDM)

rm(list = ls())
capture <- commandArgs(trailingOnly = TRUE)

opt1 = list(make_option(c("--input"), type = "character", default = "input.yaml", dest = "input"))


user_input <- function(name, argv) {
  return(name %in% names(argv))
}

argv <- parse_args(OptionParser(option_list = opt1))

if (user_input("input", argv)) {
  input = argv$input
  
  if(file.exists(input)){
    input.yaml <- read_yaml(input)
  }else{
    message('Input YAML not found.\n')
    break;
  } 
} else {
  message('Cannot proceed without input yaml file. Please use "--input" flag .\n')
}


## read  manifest
if(is.null(input.yaml[["manifest"]])){
  message('Cannot proceed without a  manifest. Please look at the sample manifest provided in the documentation.\n')
  break;
}else if(!is.null(input.yaml[["manifest"]]) & file.exists(input.yaml$manifest)){
  message('Reading the manifest...')
  manifest <- read.table(input.yaml$manifest, header = T, sep ="\t", stringsAsFactors = F)
  names(manifest) <- c("bam","sampleID","sex")
  
  
  nSamples <- NROW(manifest)
  nFemales <- NROW(manifest[manifest$sex == "M",])
  nMales <- NROW(manifest[manifest$sex == "F",])
  
  message(paste0('Found:\n * ',
                 nSamples,' individuals (',
                 nMales,' Males and ',nFemales,' females.)\n'))
  
}else if(sum(is.na(manifest)) > 0){
  # check for empty entries in the manifest
  message('Entries in manifest cannot be empty. Please check your manifest before running.\n');
  break;
}


# Check for duplicate individual IDs
if(sum(duplicated(sort(manifest$sampleID))) > 0){
  message('Found these duplicates in individual IDs\n', manifest$sampleID[duplicated(manifest$sampleID)] );
  break;
}


## check bam files
# message('Checking the bam files...\n');
#  if(sum(file.exists(as.character(manifest$bam))) < length(manifest$bam)){
#        message('List of missing bam files:\n\n');
#        print(paste(manifest$bam[!file.exists(as.character(manifest$bam))],"\n"))
#        stop('Looks like some of your bam files are missing. Check your bam paths.\n\n\n')
#  }

## check for index files
#  message('Checking the index files for bams...\n');
#  manifest$bai <- NA;
#  for (i in 1:length(manifest$bam)){
#    if(file.exists(gsub(pattern = "bam","bai",manifest$bam[i]))){
#      manifest$bai[i] <- gsub(pattern = "bam","bai",manifest$bam[i])
#    } else if(file.exists(paste0(as.character(manifest$bam[i]),".bai"))){
#      manifest$bai[i] <- paste0(as.character(manifest$bam[i]),".bai")
#    } else{
#      message(paste0('Index file not found for ',manifest$bam[i],"\n"));
#      break;
#    }
#  }

## clustermq template
if(is.null(input.yaml[["scheduler"]])){
  message('Cannot proceed without a Scheduler type\n')
  break;
}else{ 
  options(
    clustermq.scheduler = paste0(input.yaml$scheduler),
    clustermq.template = paste0("template/",input.yaml$scheduler,"_template"),
    clustermq.defaults = list(conda="edm_env")
  )
}


## bed-file
if(is.null(input.yaml[["bed.file"]])){
  data(exons.hg19,package = "ExomeDepth")
  data(exons.hg19.X,package = "ExomeDepth")
  exons <- rbind(exons.hg19,exons.hg19.X)
  exons$chromosome <- paste0("chr",exons$chromosome)
  i = 1
  while(length(which(duplicated(exons$name))) >0) { 
    exons$name[duplicated(exons$name)] <- paste0( exons$name[duplicated(exons$name)],"_dup",i)
    i = i+1 
  }
  write.table(exons,"bed_file.bed",row.names =F,sep="\t",quote=F, col.names=F) 
  system("bigWigAverageOverBed wgEncodeDukeMapabilityUniqueness35bp.bigWig bed_file.bed map_coverage.tab")
  message('Mappability results written to map_coverage.tab  \n')
  
  foo <- read.table("map_coverage.tab",stringsAsFactors=F)
  exons %>% filter(name %in% (foo %>% filter(V5 >= 0.7) %>% pull(1))) -> exons_1
  exons %>% filter(name %in% (foo %>% filter(V5 < 0.7) %>% pull(1))) -> exons_excluded
  
  exons_1$chromosome <- gsub("chr","",exons_1$chromosome) 
  write.table(exons_1[,1:3],"bed_file.bed",row.names =F,sep="\t",quote=F,col.names=F)
  message('Exons with low mappability written to excluded_exons.bed  \n')
  write.table(exons_excluded,"excluded_exons.bed",row.names =F,sep="\t",quote=F,col.names=F)
  exon_path <- "bed_file.bed"
} else if(is.na(input.yaml[["bed.file"]])){
  data(exons.hg19,package = "ExomeDepth")
  data(exons.hg19.X,package = "ExomeDepth")
  exons <- rbind(exons.hg19,exons.hg19.X)
  exons$chromosome <- paste0("chr",exons$chromosome)
  i = 1
  while(length(which(duplicated(exons$name))) >0) {
    exons$name[duplicated(exons$name)] <- paste0( exons$name[duplicated(exons$name)],"_dup",i)
    i = i+1
  }
  write.table(exons,"bed_file.bed",row.names =F,sep="\t",quote=F, col.names=F)
  
  system("bigWigAverageOverBed wgEncodeDukeMapabilityUniqueness35bp.bigWig bed_file.bed map_coverage.tab")
  message('Mappability results written to map_coverage.tab  \n')
  
  foo <- read.table("map_coverage.tab",stringsAsFactors=F)
  exons %>% filter(name %in% (foo %>% filter(V5 >= 0.7) %>% pull(1))) -> exons_1
  exons %>% filter(name %in% (foo %>% filter(V5 < 0.7) %>% pull(1))) -> exons_excluded
  
  exons_1$chromosome <- gsub("chr","",exons_1$chromosome)
  write.table(exons_1[,1:3],"bed_file.bed",row.names =F,sep="\t",quote=F,col.names=F)
  message('Exons with low mappability written to excluded_exons.bed  \n')
  write.table(exons_excluded,"excluded_exons.bed",row.names =F,sep="\t",quote=F,col.names=F)
  exon_path <- "bed_file.bed"
  
} else {
  if(tools::file_ext(input.yaml$bed.file)=="rds"){
    exons = readRDS(input.yaml$bed.file)
    if (dim(exons)[2] == 3){
      message('Found only 3 columns in ',input.yaml$bed.file,' file. Need 4th Column.\n')
    } else {
      names(exons) <- c("chromosome","start","end","name")
      exons$chromosome <- gsub("chr","",exons$chromosome)
      exons$chromosome <- paste0("chr",exons$chromosome)
      i = 1
      while(length(which(duplicated(exons$name))) >0) {
        exons$name[duplicated(exons$name)] <- paste0( exons$name[duplicated(exons$name)],"_dup",i)
        i = i+1
      }
      write.table(exons,"bed_file.bed",row.names =F,sep="\t",quote=F, col.names=F)
      
      system("bigWigAverageOverBed wgEncodeDukeMapabilityUniqueness35bp.bigWig bed_file.bed map_coverage.tab")
      message('Mappability results written to map_coverage.tab  \n')
      foo <- read.table("map_coverage.tab",stringsAsFactors=F)
      exons %>% filter(name %in% (foo %>% filter(V5 >= 0.7) %>% pull(1))) -> exons_1
      exons %>% filter(name %in% (foo %>% filter(V5 < 0.7) %>% pull(1))) -> exons_excluded
      
      exons_1$chromosome <- gsub("chr","",exons_1$chromosome)
      write.table(exons_1[,1:3],"bed_file.bed",row.names =F,sep="\t",quote=F,col.names=F)
      message('Exons with low mappability written to excluded_exons.bed  \n')
      write.table(exons_excluded,"excluded_exons.bed",row.names =F,sep="\t",quote=F,col.names=F)
      exon_path <- "bed_file.bed"
    }
    
  } else{
    exons <- read.table(input.yaml$bed.file,header = F)
    if (dim(exons)[2] == 3){
      message('Found only 3 columns in ',input.yaml$bed.file,' file. Need 4th Column.\n')
    } else {
      names(exons) <- c("chromosome","start","end","name")
      exons$chromosome <- gsub("chr","",exons$chromosome)
      exons$chromosome <- paste0("chr",exons$chromosome)
      i = 1
      while(length(which(duplicated(exons$name))) >0) {
        exons$name[duplicated(exons$name)] <- paste0( exons$name[duplicated(exons$name)],"_dup",i)
        i = i+1
      }
      write.table(exons,"bed_file.bed",row.names =F,sep="\t",quote=F, col.names=F)
      
      system("bigWigAverageOverBed wgEncodeDukeMapabilityUniqueness35bp.bigWig bed_file.bed map_coverage.tab")
      message('Mappability results written to map_coverage.tab  \n')
      foo <- read.table("map_coverage.tab",stringsAsFactors=F)
      exons %>% filter(name %in% (foo %>% filter(V5 >= 0.7) %>% pull(1))) -> exons_1
      exons %>% filter(name %in% (foo %>% filter(V5 < 0.7) %>% pull(1))) -> exons_excluded
      
      exons_1$chromosome <- gsub("chr","",exons_1$chromosome)    
      write.table(exons_1[,1:3],"bed_file.bed",row.names =F,sep="\t",quote=F,col.names=F)
      message('Exons with low mappability written to excluded_exons.bed  \n')
      write.table(exons_excluded,"excluded_exons.bed",row.names =F,sep="\t",quote=F,col.names=F)
      exon_path <- "bed_file.bed"
    }
  }
}


system(paste0("mkdir -p ",input.yaml$outputDir,"/coverage"))
system(paste0("mkdir -p ",input.yaml$outputDir,"/results"))
system(paste0("mkdir -p ",input.yaml$outputDir,"/logs"))

## get counts

source("get.counts.R")

## gather counts
source("gather.counts.R")

## call variants
source("call.variants.R")

## move log files
system(paste0("mv compute_coverage* ",input.yaml$outputDir,"/logs/."))
system(paste0("mv gather_coverage* ",input.yaml$outputDir,"/logs/."))
system(paste0("mv call_variants* ",input.yaml$outputDir,"/logs/."))

