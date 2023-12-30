# # Set the CRAN mirror directly in R
# options(repos = 'http://cran.us.r-project.org')

# # Install the venn package if you haven't already
# install.packages("venn")

# Load the venn library
library(venn)

# Read the TSV file
data <- read.table("/Users/u1983390/Projects/NewToolsProject/exp/Effectors/cdhit/AllCandidateEffectorSets.clstr_heatmapdata.tsv", header = TRUE, sep = "\t", row.names = 1)

# Create sets for the Venn diagram
sets <- lapply(data, function(col) rownames(data)[col == 1])

# Display the Venn diagram
venn(sets, category.names = c("SY-2", "S16", "S32", "S6"))
