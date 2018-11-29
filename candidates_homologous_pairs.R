#Libraries---------------------------------------------
if (!require(tidyverse)) install.packages('tidyverse')
if (!require(igraph)) install.packages('igraph')
if (!require(SDMTools)) install.packages('SDMTools')
library(tidyverse)
library(igraph)
library(SDMTools)

if(!require(IRanges)){
    source("https://bioconductor.org/biocLite.R")
    biocLite("IRanges")
}
library(IRanges)

#run if needed-----------------------------------------
  sessionInfo()
#------------------------------------------------------
  
#Parameters--------------------------------------------
#may be edited to your whishes
filename_IN = file.choose() #is a tabular output of blast
filter_pcID = 70
filter_alignLength = 100
#------------------------------------------------------

#IMPORTING DATA----------------------------------------
blastp_out = read.table(
    file = filename_IN,
    header = F, stringsAsFactors = F, 
    comment.char = "#")
colnames(blastp_out) = c(
    "query",
    "subject",
    "pc_id",
    "alignment_length",
    "mistmatches",
    "gap_openings",
    "q.start",
    "q.end",
    "s.start",
    "s.end",
    "e-value",
    "bit_score"
  )
#conversion to tible for better performances
#and use of tydiverse functions
blastp_out = as.tbl(blastp_out)

# ======================================

# First filtering based on length and id % and get rid of "autoblast"
blastp_out = blastp_out %>% filter(query != subject)
blastp_out = blastp_out %>% filter(pc_id >= filter_pcID,
                                   alignment_length >= filter_alignLength )


# Second filtering based on couples of query - subject
blastp_out$couple_name = apply(blastp_out, 1, function(x){
                            ids = sort(c(x[1], x[2]))
                            paste(ids, collapse = "-")
                        })
blastp_out = blastp_out %>% mutate(pc_times_length = pc_id * alignment_length)
blastp_out = blastp_out %>% group_by(couple_name) %>% slice(which.max(pc_times_length)) %>% ungroup()

# similarity graph construction, using igraph
sim_graph = graph.data.frame(blastp_out[,1:3], 
                             directed = F)
#pc_id is added as an attribute of the graph and latter
#used as weights

# family contructions by edge betweeness clutering
families = cluster_edge_betweenness(sim_graph,
                         weights = E(sim_graph)$pc_id,
                         directed = F,
                         edge.betweenness = F,
                         merges = TRUE,
                         bridges = TRUE,
                         modularity = TRUE,
                         membership = TRUE)
V(sim_graph)$fam = families$membership 
#family data can then be exported
