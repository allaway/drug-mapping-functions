options(java.parameters = "-Xmx8g" ) 
library(shiny)
library(DT)
library(rJava)
library(rcdk)
library(fingerprint)
library(enrichR)
library(webchem)
library(plyr)
library(tidyverse)
library(synapser)
synLogin()

db <- readRDS(synGet("syn11712148")$path) %>% 
  filter(!is.na(hugo_gene)) %>% 
  select(internal_id, common_name, hugo_gene, mean_pchembl, n_quantitative, n_qualitative)

fp.db <- readRDS(synGet("syn11693143")$path)

is.smiles <- function(x, verbose = TRUE) { ##corrected version from webchem
  if (!requireNamespace("rcdk", quietly = TRUE)) {
    stop("rcdk needed for this function to work. Please install it.",
         call. = FALSE)
  }
  # x <- 'Clc(c(Cl)c(Cl)c1C(=O)O)c(Cl)c1Cl'
  if (length(x) > 1) {
    stop('Cannot handle multiple input strings.')
  }
  out <- try(rcdk::parse.smiles(x), silent = TRUE)
  if (inherits(out[[1]], "try-error") | is.null(out[[1]])) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}

parseInputFingerprint <- function(input) {
  test_smiles <- is.smiles(input)
  if(is.smiles(input==TRUE)){
    input.mol <- parse.smiles(as.character(input))
    lapply(input.mol, do.typing)
    lapply(input.mol, do.aromaticity)
    lapply(input.mol, do.isotopes)
    fp.inp <- lapply(input.mol, get.fingerprint, type = "extended")
  }else{
    print('Please input a valid SMILES string.')
  }
}


convertDrugToSmiles <- function(input) {
  filt <- filter(db.names, common_name == input) %>% dplyr::select(smiles)
  filt
}

getTargetList <- function(selectdrugs) {
  targets <- filter(db, common_name %in% selectdrugs) %>% 
    arrange(-n_quantitative) %>%
    select(common_name, hugo_gene, mean_pchembl, n_quantitative, n_qualitative)
  
  if (nrow(targets) > 1) {
    targets
  } else {
    print("none found")
  }
}

similarityFunction <- function(input) {
  input <- input
  fp.inp <- parseInputFingerprint(input)

    sims <- lapply(fp.inp, function(i) {
      sim <- lapply(fp.db, function(j) {
        distance(i, j)
      })
      bar <- ldply(sim)
      colnames(bar) <- c("match", "similarity")
      bar
  })
  sims <- ldply(sims)
}

getSimMols <- function(sims, sim.thres) {
  sims2 <- sims %>% dplyr::filter(similarity >= sim.thres) %>% arrange(-similarity)
  sims2$internal_id <- as.character(sims2$match)
  sims2$`Tanimoto Similarity` <- signif(sims2$similarity, 3)
  targets <- left_join(sims2, db) %>% 
    dplyr::select(common_name, `Tanimoto Similarity`) %>% 
    distinct()
}

getMolsFromGenes <- function(inp.gene) {
  genes <- trimws(unlist(strsplit(inp.gene,",")))
  mols <- db %>% 
    mutate(hugo_gene = as.character(hugo_gene)) %>% 
    mutate(keep = hugo_gene %in% genes) %>% 
    group_by(internal_id, keep) %>% 
    mutate(count = n()) %>% 
    filter(keep == TRUE, count >= length(genes)) %>% 
    ungroup() %>% 
    distinct() %>% 
    select(-keep, -count)
}
