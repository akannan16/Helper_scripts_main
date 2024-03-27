# Install and load the "gplots" package
install.packages("gplots")
library(gplots)

# Read the data from the file Coverage_Profile that cointains coverage info from 25X to 1000X in the form of %)

# Select the columns you want to include in the heatmap
subset_columns <- c("25x", "50x", "100x", "200x", "500x", "1000x")

# Create a subset of the data frame with the selected columns
heatmap_data_subset <- Coverage_Profile[, c("ID", subset_columns)]

# Set the DDM_Patient_ID column as row names
rownames(heatmap_data_subset) <- heatmap_data_subset$ID
heatmap_data_subset <- heatmap_data_subset[, -1]

# Convert the data frame subset to a numeric matrix
heatmap_matrix <- as.matrix(heatmap_data_subset[, subset_columns])

# Create the heatmap with color key only
heatmap_plot <- heatmap.2(heatmap_matrix, Rowv = NA, Colv = NA, col = colorRampPalette(brewer.pal(11, "RdYlGn"))(100),
                          scale = "none", trace = "none", density.info = "none",
                          xlab = "Coverage", labRow = Coverage_Profile$ID)

# Add a color legend
library("RColorBrewer")
legend("right", legend = c("Low", "High"), fill = colorRampPalette(brewer.pal(11, "RdYlGn"))(100), cex = 0.8)
