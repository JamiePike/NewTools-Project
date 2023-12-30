# # Set the CRAN mirror directly in R
# options(repos = 'http://cran.us.r-project.org')

# # Install the venn package if you haven't already
# install.packages("venn")

# Load the venn library
library(ggplot2)
library(ggvenn)


# Read the TSV file
data <- read.table("/Users/u1983390/Projects/NewToolsProject/exp/Effectors/cdhit/AllCandidateEffectorSets.clstr_heatmapdata.tsv", header = TRUE, sep = "\t", row.names = 1)

# Create sets for the Venn diagram
sets <- lapply(data, function(col) rownames(data)[col == 1])

# Convert sets to a format suitable for ggvenn
ggvenn_data <- list(
  S16 = sets[[1]],
  S32 = sets[[2]],
  S6 = sets[[3]],
  SY_2 = sets[[4]]
)

# Create the Venn diagram using ggplot2 and ggvenn
ggvenn_plot <- ggvenn(ggvenn_data, fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4)

# Save the Venn diagram
ggsave("/Volumes/Jamie_EXT/Projects/NewToolsProject/exp/Effectors/cdhit/sharedCandEffsVenn.png", plot = ggvenn_plot, width = 8, height = 6, units = "in")
