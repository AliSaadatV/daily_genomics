### This code shows how to use random walk with restart to rank candidate genes based on a set of known genes (called seed genes)
### Input: is a data frame which contains CHROM, POS, REF, ALT, SYMBOL, Gene (ensembl-id)
### Output is the ranked genes

library(tidyverse)
library(RandomWalkRestartMH) 
library(STRINGdb)

data <- read_delim("toy_df.tsv")

#get STRINGdb
string_db <- STRINGdb$new( version="11.5", species=9606, score_threshold=500, input_directory="")

#map ensembl-id to STRING-id
string_id_mapped <- string_db$map(as.data.frame(data),
                                  "Gene", removeUnmappedRows = F ) %>% 
  dplyr::select(SYMBOL, Gene, STRING_id) %>% 
  distinct(Gene, .keep_all = T)
undefined <- string_id_mapped %>% filter(is.na(STRING_id))
string_id_mapped <- string_id_mapped %>% filter(!is.na(STRING_id))

#read graph and clean them.
string_db_graph <- string_db$get_graph()
string_db_graph <- igraph::simplify(string_db_graph, remove.multiple = TRUE, remove.loops = TRUE) 
STRING_Multiplex <- create.multiplex(list(PPI=string_db_graph))
AdjMatrix_STRING_PATH <- compute.adjacency.matrix(STRING_Multiplex)
AdjMatrixNorm_STRING_PATH <- normalize.multiplex.adjacency(AdjMatrix_STRING_PATH)

#get all nodes names (used later)
all_nodes_STRING_id <- igraph::V(string_db_graph)$name

#prepare seeds: either use OpenTargets/DisGeNET or use your own set of gene seeds.
#let's imagine there are 3 genes related to our phenotype: IFIH1, DDX58, IRF7
seeds <- c("IFIH1", "DDX58", "IRF7")
seeds_id <- string_db$map(as.data.frame(seeds), "seeds", removeUnmappedRows = T) %>%
  distinct(STRING_id, .keep_all = T)

#random walk with restart
RWR_STRING_Results <- Random.Walk.Restart.Multiplex(AdjMatrixNorm_STRING_PATH,
                                                    STRING_Multiplex,seeds_id$STRING_id)
RWR_STRING_Results_Scores <- left_join(string_id_mapped,
                                      RWR_STRING_Results$RWRM_Results, by=c("STRING_id"="NodeNames"))

#now we order the data according to the Score from random walk
data_ordered <- data %>% 
  left_join(RWR_STRING_Results_Scores) %>% 
  arrange(desc(Score))
