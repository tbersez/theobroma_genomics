#Libraries
if (!require(tidyverse)) install.packages('tidyverse')
if (!require(igraph)) install.packages('igraph')
if (!require(SDMTools)) install.packages('SDMTools')
library(tidyverse)
library(igraph)
library(SDMTools)

#Parameters
filename_IN = file.choose()
filter_pcID = 70
filter_alignLength = 100

#IMPORTING DATA
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

blastp_out = as.tbl(blastp_out)

# First filtering based on length and id % and get rid of autoblast
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

# similarity graph construction

sim_graph = graph.data.frame(blastp_out[,1:3], 
                             directed = F)

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

#----Ploting the biggest sub graph (may take some time)----

c = clusters(sim_graph)
biggest_sub = induced.subgraph(sim_graph,
                               c$membership == order(-c$csize)[1])
cols = heat.colors(length(unique(round(E(biggest_sub)$pc_id) -
                     min(round(E(biggest_sub)$pc_id)))))
par(cex = 2)
plot(biggest_sub,
     vertex.label = NA,
     vertex.size = 5,
     vertex.color = V(sim_graph)$fam,
     layout = layout_as_tree(biggest_sub),
     edge.color = cols[round(E(biggest_sub)$pc_id) -
       min(round(E(biggest_sub)$pc_id))],
     main = 'Similarity subgraph'
     )
legend.gradient(cbind(x =c(-1.4,-1.5,-1.5,-1.4), y =c(1.0,1.0,0.8,0.8)),
                cols = cols,
                limits = c(100, filter_pcID),
                title = "Percentage of identity")

     