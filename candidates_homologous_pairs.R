<<<<<<< HEAD
#Libraries
if (!require(tidyverse)) install.packages('tidyverse')
if (!require(igraph)) install.packages('igraph')
library(tidyverse)
library(igraph)

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

=======
library(tidyverse)
library(igraph)

find_id = function(blast_tab,
                   query,
                   subject)
  {
  res = blast_tab$pc_id[which(blast_tab[,1] == query 
                         && blast_tab[,2] == subject)]
  if(length(res) != 0){
    return(res)
  }
  return(0)
}

#Parameters
filename_IN = file.choose()
filter_pcID = 90
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

#hierarchical clustering
dis = dist(blastp_out$pc_id,
     method = "euclidian")
tree = hclust(dis,
              method = "centroid")
plot(tree)

#similarity graph construction

vertex = unique(c(blastp_out$query,blastp_out$subject))

sim_graph = graph.empty (length(vertex), directed = FALSE) %>%
  set_vertex_attr("name", value = vertex) %>%
  add.edges(c(rbind(blastp_out$query, blastp_out$subject))) %>%
  set_edge_attr("weight", value = blastp_out$pc_id) #pas dans le bon ordre?


>>>>>>> 93363317fa646675c7e9140735d156a2668e5fff
