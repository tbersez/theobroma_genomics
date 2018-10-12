library(tidyverse)

#Parameters
filename_IN = file.choose()
filter_pcID = 60
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

blastp_out %>% group_by(couple_name) %>% slice(which.max(pc_times_length)) %>% ungroup()



