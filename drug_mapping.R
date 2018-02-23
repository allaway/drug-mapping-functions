library(parallel)
source("mapping_helpers.R") ##some helper functions lifted from shiny app

hugo_list <- read.table("all-symbols.tsv", sep = "\t", header = T) %>% sample_n(100)

###### TARGETS TO DRUGS ###########################################################

###hugo_list = character vector of HGNC Symbols/HUGO IDs 
###parallelized = boolean indicating whether to use mclapply or not. Killing mclapply prematurely can cause zombie processes

getDrugsfromGenes <- function(hugo_list, parallelized = T){
  db_overlap <- intersect(hugo_list, unique(db$hugo_gene))
  print(paste0(length(db_overlap), " gene products are druggable:"))
  print(c(db_overlap))
  foo <- db %>% filter(hugo_gene %in% hugo_list)
  print(paste0(length(unique(foo$common_name)), 
               " drugs found from ", 
               length(db_overlap), 
               " genes."))
  
  n <- length(unique(db$hugo_gene))
  
  if(parallelized == T){
    fisher.list <- mclapply(unique(foo$common_name), function(x){
      bar <- filter(db, common_name == x)
      a <- n - length(union(unique(bar$hugo_gene), unique(db_overlap)))
      b <- length(setdiff(unique(bar$hugo_gene), unique(db_overlap)))
      c <- length(setdiff(unique(db_overlap), unique(bar$hugo_gene)))
      d <- length(intersect(unique(bar$hugo_gene), unique(db_overlap)))
      res <- fisher.test(matrix(c(a,b,c,d), nrow=2), alternative="greater")
      c("p.value" = res$p.value, res$estimate, "alt" = res$alternative)
    }, mc.cores = detectCores())
  }
  
  if(parallelized == F){
    fisher.list <- lapply(unique(foo$common_name), function(x){
      bar <- filter(db, common_name == x)
      a <- n - length(union(unique(bar$hugo_gene), unique(db_overlap)))
      b <- length(setdiff(unique(bar$hugo_gene), unique(db_overlap)))
      c <- length(setdiff(unique(db_overlap), unique(bar$hugo_gene)))
      d <- length(intersect(unique(bar$hugo_gene), unique(db_overlap)))
      res <- fisher.test(matrix(c(a,b,c,d), nrow=2), alternative="greater")
      c("p.value" = res$p.value, res$estimate, "alt" = res$alternative)
    })
  }
  
  names(fisher.list) <- unique(foo$common_name)
  res <- ldply(fisher.list)
  res$bh2 <- p.adjust(res$p.value, method = "BH")
  res <- res %>% set_names(c("drug", "p.value", "odds ratio", "alt", "bh")) %>% select(drug, p.value, bh, `odds ratio`, alt)
  list("drug_fisher_test" = res, "drug_target_matrix" = foo)
}
# 
# 
# start_time <- Sys.time()
# list <- getDrugsfromGenes(hugo_list$symbol, parallelized = F)
# end_time <- Sys.time()
# end_time - start_time
# 
# start_time <- Sys.time()
# list <- getDrugsfromGenes(hugo_list$V1, parallelized = T)
# end_time <- Sys.time()
# end_time - start_time
# 


###### DRUGS TO TARGETS ###########################################################
# structs <- read.table("drug_structures.txt", header = T, sep = "\t", comment.char = "") #overlapping FIMM/OHSU drugs annotated with structures
# input_struct <- structs$smiles

###drugnames = vector containing drug names
###input_struct = equal length, same order vector as drugnames containing 1 SMILES structure for each drug name
###tanimoto_threshold = 0 to 1 range for similarity, allows flexiiblity to map to structural analogs or similar molecules, where 1 is identical
   ###note that in testing this I have found that a similarity of 1 can still sometimes occur for nearly identical "large" small molecules (statins are an example)
   ###recommend 0.95 as a starting point
###parallelized = boolean, whether to use mclapply or not. Killing mclapply prematurely can cause zombie processes

getGenesfromDrugs <- function(drugnames, input_struct, tanimoto_threshold, parallelized = T){
  if(parallelized == T){
    foo <- mclapply(input_struct, function(x){
      structure <- as.character(x)
      sims <- similarityFunction(structure)
      getSimMols(sims, tanimoto_threshold)
      }, mc.cores = detectCores())
    names(foo) <- drugnames
    bar <- ldply(foo)
    colnames(bar) <- (c("input_drug", "common_name", "Tanimoto Similarity"))
  }
  
  if(parallelized == F){
    foo <- lapply(input_struct, function(x){
      structure <- as.character(x)
      sims <- similarityFunction(structure)
      getSimMols(sims, tanimoto_threshold)
    })
    names(foo) <- drugnames
    bar <- ldply(foo)
    colnames(bar) <- (c("input_drug", "common_name", "Tanimoto Similarity"))
  }
  
  bar <- inner_join(db, bar, by = "common_name")
  
  
}

# start_time <- Sys.time()
# getGenesfromDrugs(structs$drug[1:5], structs$smiles[1:5], 0.95, parallelized = F)
# end_time <- Sys.time()
# end_time - start_time
# 
# start_time <- Sys.time()
# getGenesfromDrugs(structs$drug[1:5], structs$smiles[1:5], 0.95, parallelized = T)
# end_time <- Sys.time()
# end_time - start_time
# 
# genes<-getGenesfromDrugs(structs$drug[1:5], structs$smiles[1:5], 0.95, parallelized = T)

