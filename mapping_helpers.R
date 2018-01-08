library(rJava)
library(rcdk)
library(fingerprint)
library(plyr)
library(tidyverse)
library(synapseClient)
synapseLogin()

evo <- readRDS(synGet("syn11658938")@filePath)
evo$Structure_ID <- as.character(evo$Structure_ID)
evo$Common_Name <- as.character(evo$Common_Name)
evo <- evo %>% filter(N_quantitative >= N_inactive | N_qualitative >= N_inactive | N_DGIDB > 0)

db.genes <- unique(evo$Hugo_Gene)

fp.evo <- readRDS(synGet("syn11658941")@filePath)[unique(evo$Original_molecule_SMILES)]

##converts SMILES string to fingerprint
parseInputFingerprint <- function(input) {
  input.mol <- parse.smiles(input)
  lapply(input.mol, do.typing)
  lapply(input.mol, do.aromaticity)
  lapply(input.mol, do.isotopes)
  fp.inp <- lapply(input.mol, get.fingerprint, type = "extended")
}

getTargetList <- function(selectdrugs) {
  targets <- filter(evo, Common_Name %in% selectdrugs) %>% dplyr::select(Common_Name, Hugo_Gene, MedianActivity_nM, N_quantitative, N_qualitative, 
                                                                         N_inactive, N_DGIDB, Confidence_Score) %>% arrange(-N_quantitative)
  if (nrow(targets) > 1) {
    targets
  } else {
    print("none found")
  }
}

getSimMols <- function(input, sim.thres) {
  input <- input
  fp.inp <- parseInputFingerprint(input)
  
  sims <- lapply(fp.inp, function(i) {
    sim <- sapply(fp.evo, function(j) {
      distance(i, j)
    })
    bar <- as.data.frame(sim)
    bar$match <- rownames(bar)
    bar
  })
  sims <- ldply(sims)
  sims2 <- sims %>% filter(sim >= sim.thres) %>% arrange(-sim)
  sims2$Original_molecule_SMILES <- as.character(sims2$match)
  sims2$`Tanimoto Similarity` <- signif(sims2$sim, 3)
  targets <- left_join(sims2, evo, by = "Original_molecule_SMILES") %>% dplyr::select(Common_Name, `Tanimoto Similarity`) %>% distinct()
}

dbs <- c("GO_Molecular_Function_2017", "GO_Cellular_Component_2017", "GO_Biological_Process_2017", 
         "KEGG_2016")

getGeneOntologyfromTargets <- function(selectdrugs) {
  selectdrugs <- selectdrugs
  targets <- getTargetList(selectdrugs)
  
  if (nrow(targets) > 1) {
    enriched <- enrichr(as.vector(targets$Hugo_Gene), dbs)
  } else {
    print("no targets")
  }
  
}

getMolsFromGenes <- function(inp.gene) {
  mols <- filter(evo, Hugo_Gene == inp.gene) %>% 
    select(Structure_ID,Supplier_Molname, Common_Name, MedianActivity_nM, N_quantitative, N_qualitative, N_inactive, Original_molecule_SMILES) %>% 
    distinct()
  
}

getSmiles <- function(input.name) {
  input.name <- input.name
  input.name <- URLencode(input.name)
  query <- as.vector(cir_query(input.name, representation = "smiles", first = TRUE))
  query
}