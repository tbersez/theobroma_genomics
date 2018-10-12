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
similarity = matrix(data = NA,
                    nrow = length(vertex),
                    ncol = length(vertex)
                    )
colnames(similarity) = vertex
rownames(similarity) = vertex
rm(vertex)#for memory management

#NON COMPUTIONNALY EFFECTIVE !!
#filling the similarity matrix with percent of identity
#for(i in 1:length(rownames(similarity))){
#  for(j in 1:length(colnames(similarity))){
#    qu = colnames(similarity)[i]
#    subj = rownames(similarity)[j]
#    similarity[i,j] = find_id(blastp_out,
#                              qu,
#                              subj)
#  }
#}

#for test only
for(i in 1:10){
  for(j in 1:10){
    qu = colnames(similarity)[i]
    subj = rownames(similarity)[j]
    similarity[i,j] = find_id(blastp_out,
                              qu,
                              subj)
  }
}

