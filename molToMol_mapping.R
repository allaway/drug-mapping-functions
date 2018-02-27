options(java.parameters = "-Xmx8g" ) 
library(rJava)
library(rcdk)
library(plyr)
library(tidyverse)
library(fingerprint)

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
    input.mol <- parse.smiles(as.character(input))
    lapply(input.mol, do.typing)
    lapply(input.mol, do.aromaticity)
    lapply(input.mol, do.isotopes)
    fp.inp <- lapply(input.mol, get.fingerprint, type = "circular")
}


mapDrugSets <- function(names.a, smiles.a, names.b, smiles.b, similarity.cutoff){ 

  temp <- lapply(smiles.a, function(x){
    suppressMessages(foo <- isTRUE(is.smiles(x)))
    if(foo == FALSE){
      stop(paste0("SMILES ",x, " in set A is invalid."))
    }
  })
  
  temp <- lapply(smiles.b, function(x){
    if(isTRUE(is.smiles(x)) == FALSE){
      stop(paste0("SMILES ",x, " in set B is invalid."))
    }
  })
  
  fing.a <- parseInputFingerprint(smiles.a)
  names(fing.a) <- names.a
  fing.b <- parseInputFingerprint(smiles.b)
  names(fing.b) <- names.b
  
  sims <- lapply(fing.a, function(i) {
    sim <- lapply(fing.b, function(j) {
      distance(i, j)
    })
    bar <- ldply(sim) 
    colnames(bar) <- c("match", "similarity")
    bar %>% filter(similarity >= similarity.cutoff)
  })
  sims <- ldply(sims)
  colnames(sims) <- c("set A", "set B", "similarity")
  sims
}

# #sample data
# names.a <- c("a","b","c")
# smiles.a <- c("CCCC", "CCNC", "CNNC")
# names.b <- c("1", "2", "3")
# smiles.b <- c("CNNC", "CCCC", "CNOC")

## names.a, names.b = character vectors of names of compounds
## smiles.a, smiles.b = character vectors of SMILES of compounds, matched to names
## similarity.cutoff = tanimoto distance cutoff, this uses circular fingerprint which is fairly stringent (i.e. will have low tanimoto similarity unless identical)
## can confuse opposite chiral centers as identical
res <- mapDrugSets(names.a, smiles.a, names.b, smiles.b, 0.99)



