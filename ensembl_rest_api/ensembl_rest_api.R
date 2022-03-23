### This code uses ensembl API to extract number of het/hom samples for query variants from 1000G project
### Documentation can be found at: http://rest.ensembl.org/documentation/info/variation_id
### Input: is a data frame which contains CHROM, POS, REF, ALT, and dbSNP rs_id of variants (AF of 1000G variants is usually available, if annotated with Annovar)
### Output is number of heterozygous and homozygous samples for each variant

library(tidyverse)
library(jsonlite)
library(httr)
library(xml2)

data <- read_delim("toy_df.tsv")
server <- "http://rest.ensembl.org"
data$kg_eur_hom_samples <- 0
data$kg_eur_het_samples <- 0
data$kg_eur_samples <- 0
for(i in 1:nrow(data)){
  #check if rs_ids exists and AF of variant in 1000G is not 0 or NA
  if(is.na(data$avsnp147[i]) | data$ALL.sites.2015_08[i]==0 | is.na(data$ALL.sites.2015_08[i])) {next}
  
  #send request
  ext <- str_c("/variation/human/", data$avsnp147[i], "?population_genotypes=1")
  r <- GET(paste(server, ext, sep = ""), content_type("application/json"))
  #if(status_code(r)==400) {next}
  stop_for_status(r)
  temp_df <- head(fromJSON(toJSON(content(r))))$population_genotypes
  
  #if there is no information, go next
  if(length(temp_df)==0) {next}   
  
  #extract 1000G phae3 european ancestry
  temp_df <- temp_df %>%
    filter(population=="1000GENOMES:phase_3:EUR")
  
  #count homozygous
  count_hom <- temp_df %>%
    filter(genotype == str_c(data$ALT[i], "|", data$ALT[i])) %>%
    select(count) %>%
    unlist(use.names = F)
  
  #some variants do NOT have homozygous in population, so add count_hom with care
  if(length(count_hom)>=1) {data$kg_eur_hom_samples[i] <- count_hom}
  
  #count heterozygous
  count_het <- temp_df %>%
    filter(genotype %in% c(str_c(data$REF[i], "|", data$ALT[i]), str_c(data$ALT[i], "|", data$REF[i]))) %>%
    select(count) %>%
    unlist(use.names = F)
  
  #add heterozygous information
  if(length(count_het)>=1) {data$kg_eur_het_samples[i] <- count_het}
  
  #get total number of samples for the variant
  data$kg_eur_samples[i] <- sum(unlist(temp_df$count))
  
  #be kind to ensembl server, sleep for 0.02 seconds!
  Sys.sleep(0.02)                                                               
}
