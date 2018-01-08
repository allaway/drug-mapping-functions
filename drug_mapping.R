library(parallel)
source("mapping_helpers.R") ##some helper functions lifted from shiny app

hugo_list <- read.table("all-symbols.tsv", sep = "\t") %>% 
  sample_n(20)

###### TARGETS TO DRUGS ###########################################################

###hugo_list = character vector of HGNC Symbols/HUGO IDs 
###parallelized = boolean indicating whether to use mclapply or not. Killing mclapply prematurely can cause zombie processes

getDrugsfromGenes <- function(hugo_list, parallelized = T){
  evo_overlap <- intersect(hugo_list, unique(evo$Hugo_Gene))
  print(paste0(length(evo_overlap), " gene products are druggable:"))
  print(c(evo_overlap))
  foo <- evo %>% filter(Hugo_Gene %in% hugo_list)
  print(paste0(length(unique(foo$Common_Name)), 
               " drugs found from ", 
               length(evo_overlap), 
               " genes."))
  
  n <- length(unique(evo$Hugo_Gene))
  
  if(parallelized == T){
    fisher.list <- mclapply(unique(foo$Common_Name), function(x){
      evo_overlap <- evo_overlap
      bar <- filter(evo, Common_Name == x)
      a <- n - length(union(unique(bar$Hugo_Gene), unique(evo_overlap)))
      b <- length(setdiff(unique(bar$Hugo_Gene), unique(evo_overlap)))
      c <- length(setdiff(unique(evo_overlap), unique(bar$Hugo_Gene)))
      d <- length(intersect(unique(bar$Hugo_Gene), unique(evo_overlap)))
      res <- fisher.test(matrix(c(a,b,c,d), nrow=2), alternative="greater")
      c("p.value" = res$p.value, res$estimate, "alt" = res$alternative)
    }, mc.cores = detectCores())
  }
  
  if(parallelized == F){
    fisher.list <- lapply(unique(foo$Common_Name), function(x){
      evo_overlap <- evo_overlap
      bar <- filter(evo, Common_Name == x)
      a <- n - length(union(unique(bar$Hugo_Gene), unique(evo_overlap)))
      b <- length(setdiff(unique(bar$Hugo_Gene), unique(evo_overlap)))
      c <- length(setdiff(unique(evo_overlap), unique(bar$Hugo_Gene)))
      d <- length(intersect(unique(bar$Hugo_Gene), unique(evo_overlap)))
      res <- fisher.test(matrix(c(a,b,c,d), nrow=2), alternative="greater")
      c("p.value" = res$p.value, res$estimate, "alt" = res$alternative)
    })
  }
  
  names(fisher.list) <- unique(foo$Common_Name)
  res <- ldply(fisher.list)
  res$bh2 <- p.adjust(res$p.value, method = "BH")
  res <- res %>% set_names(c("drug", "p.value", "odds ratio", "alt", "bh")) %>% select(drug, p.value, bh, `odds ratio`, alt)
  list("drug_fisher_test" = res, "drug_target_matrix" = foo)
}


start_time <- Sys.time()
list <- getDrugsfromGenes(hugo_list$V1, parallelized = F)
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
list <- getDrugsfromGenes(hugo_list$V1, parallelized = T)
end_time <- Sys.time()
end_time - start_time



###### DRUGS TO TARGETS ###########################################################
structs <- read.table("drug_structures.txt", header = T, sep = "\t", comment.char = "") #overlapping FIMM/OHSU drugs annotated with structures
input_struct <- structs$smiles

###drugnames = vector containing drug names
###input_struct = equal length, same order vector as drugnames containing 1 SMILES structure for each drug name
###tanimoto_threshold = 0 to 1 range for similarity, allows flexiiblity to map to structural analogs or similar molecules, where 1 is identical
   ###note that in testing this I have found that a similarity of 1 can still sometimes occur for nearly identical "large" small molecules (statins are an example)
   ###recommend 0.95 as a starting point
###whether to use mclapply or not. Killing mclapply prematurely can cause zombie processes

getGenesfromDrugs <- function(drugnames, input_struct, tanimoto_threshold, parallelized = T){
  if(parallelized == T){
    foo <- mclapply(input_struct, function(x){
      structure <- as.character(x)
      getSimMols(structure, tanimoto_threshold)
      }, mc.cores = detectCores())
    names(foo) <- drugnames
    bar <- ldply(foo)
    colnames(bar) <- (c("input_drug", "Common_Name", "Tanimoto Similarity"))
  }
  
  if(parallelized == F){
    foo <- lapply(input_struct, function(x){
      structure <- as.character(x)
      getSimMols(structure, tanimoto_threshold)
    })
    names(foo) <- drugnames
    bar <- ldply(foo)
    colnames(bar) <- (c("input_drug", "Common_Name", "Tanimoto Similarity"))
  }
  
  bar <- inner_join(evo, bar, by = "Common_Name")
  
  
}

start_time <- Sys.time()
getGenesfromDrugs(structs$drug[1:5], structs$smiles[1:5], 0.95, parallelized = F)
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
getGenesfromDrugs(structs$drug[1:5], structs$smiles[1:5], 0.95, parallelized = T)
end_time <- Sys.time()
end_time - start_time

genes<-getGenesfromDrugs(structs$drug[1:5], structs$smiles[1:5], 0.95, parallelized = F)

