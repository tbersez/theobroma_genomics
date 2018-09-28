#IMPORTING DATA
setwd("m2_saclay/comparative_genomics/theobroma_genomics/")
  blastp_out = read.table(
    file = "./results_BLASTp.tsv",
    header = F)
  colnames(blastp_out) = c(
    "query",
    "subject",
    "%id",
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
#FILTERING
    blastp_out = subset(x = blastp_out, query != subject)
    blastp_out = subset(x = blastp_out, alignment_length >= 100)
    blastp_out = subset(x = blastp_out, `%id` >= 75)
#EXPORTING
    write.table(blastp_out,
                file = "blastp_filtered",
                sep = '\t')
    