
if(!require("plotrix")){
    install.packages("plotrix")
}
library("plotrix")



data_families  = read.table(file = "results_from_galaxy_FTAG_finder/families_filter30-70.tabular",
                            stringsAsFactors = FALSE)

# to get the size of families from the raw file, one need to call the table() function twice.
# First to get the size of each family, then to get the number of families with a given size.

family_sizes = table(table(data_families$V2))


gap.barplot(family_sizes[-1], gap = c(500,1480), ytics = c(0,50,100,150,250,480,1500,1568),
            xaxlab = names(family_sizes[-1]), las =2, col = rep("chocolate4", 62), 
            xlab = "Family size", ylab = "Number of families")

# ==================================================================================================

library(ggplot2)

data_KaKs = read.csv(file = "ka_ks_values", sep = ";", header = T)

ggplot(data_KaKs) + 
    geom_density(mapping = aes(x = ratio), bg = "blue",alpha = 0.5) + 
    xlim(0, 2) +
    geom_vline(xintercept = 1, linetype = "dotted") +
    theme_bw()

# add  other species in the plot later, when data will have been uploaded







