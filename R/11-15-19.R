###################################################################################################################################################
# Logan Nguyen
# Fehl Lab
# 10-26-19
# mgfs
######################################################################################################################################################

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("proteomics")


demo(graphics)

### Reading files into R ############################################################################################################################

rm(list = ls()) # clearing the enviroment
setwd("/home/logan/R") # setting working directory

# Loading libraries
library(readxl) # for loading xlsx files
library(ggplot2) # for plotting
library(RColorBrewer) ## Color palettes

# Loading in the file path, specifically the first sheet
filename <- "/home/logan/Documents/Research/Fehl Lab/Data/9-25-19/ja8b07328_si_002.xlsx"
# data <- read.csv(filename, sep = "\t", header = TRUE) # for .txt files and .csv files
data <- read_excel(filename, sheet = "24h", na = "---") # for .xls and .xlsx files

# Renaming the first columns
names(data)[4:13] <- c("Da", "Db", "Dc", "Dd", "De", "osmi2a", "osmi2b", "osmi2c", "osmi2d", "osmi2e")

# Deleting the first row as we don't need it anymore
data <- data[-c(1), ]

# Converting the 4 character columns to the numeric data type
data <- transform(data, Da = as.numeric(Da))
data <- transform(data, Db = as.numeric(Db))
data <- transform(data, Dc = as.numeric(Dc))
data <- transform(data, Dd = as.numeric(Dd))
data <- transform(data, De = as.numeric(De))
data <- transform(data, osmi2a = as.numeric(osmi2a))
data <- transform(data, osmi2b = as.numeric(osmi2b))
data <- transform(data, osmi2c = as.numeric(osmi2c))
data <- transform(data, osmi2d = as.numeric(osmi2d))
data <- transform(data, osmi2e = as.numeric(osmi2e))


## Basic Data Info ############################################################################################################################
head(data) # looking at the first 6 rows of our data
summary(data) # basic stats on our data
sapply(data, class) # checking the class of each of our columns
dim(data) # dimensions
typeof(data)


# Some basic plots

ggplot(diamonds)  # if only the dataset is known.
ggplot(diamonds, aes(x=carat))  # if only X-axis is known. The Y-axis can be specified in respective geoms.
ggplot(diamonds, aes(x=carat, y=price)) + geom_point() # if both Xand Y axes are fixed for all layers.


ggplot(diamonds, aes(x=carat, color=cut))  # Each category of the 'cut' variable will now have a distinct  color, once a geom is added.
ggplot(diamonds, aes(x=carat), color="steelblue")

ggplot(diamonds, aes(x=carat, y=price, color=cut)) + geom_point() + geom_smooth() # Adding scatterplot geom (layer1) and smoothing geom (layer2).
ggplot(diamonds) + geom_point(aes(x=carat, y=price, color=cut)) + geom_smooth(aes(x=carat, y=price, color=cut)) 
ggplot(diamonds) + geom_point(aes(x=carat, y=price, color=cut)) + geom_smooth(aes(x=carat, y=price)) # Remove color from geom_smooth
ggplot(diamonds, aes(x=carat, y=price)) + geom_point(aes(color=cut)) + geom_smooth()  # same but simpler
ggplot(diamonds, aes(x=carat, y=price, color=cut, shape=color)) + geom_point()

### Volcano Plot ############################################################################################################################ 
theme_set(theme_bw())
theme(legend.position = "bottom")

# We need to plot -log(p value) vs fold change, so load in these vector columns
data$'-logpvalue' <- -log(data$p.value) # creating a new column and modify pvalue as we want to plot -log(pvalue)
pvalue <- data$'-logpvalue'
foldchange <- data$fold_change

# Plotting
volcano <- ggplot(data, aes(x = foldchange, y = pvalue, color = data$Quantified.peptides)) + geom_point()
volcano
volcano <- volcano + xlim(0.25, 1.75) + ylim(0, 16) + scale_color_gradient(limits = c(0, 30))
volcano <- volcano + theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))
volcano <- volcano + labs(x = "Fold Change", y = "-log(p value)") + ggtitle("24 hr")
volcano <- volcano + labs(color = "# of peptides")

volcano <- volcano + geom_hline(yintercept = 1.3, linetype = "dashed", color = "red")
volcano <- volcano + geom_vline(xintercept = 0.8, linetype = "dashed", color = "red")
volcano <- volcano + geom_vline(xintercept = 1.2, linetype = "dashed", color = "red")

volcano <- volcano + geom_text(aes(label = ifelse(foldchange >= 1.2 & pvalue >= 1.3, as.character(data$Gene.Symbol), '')), hjust = 0, vjust = 0)
volcano <- volcano + geom_text(aes(label = ifelse(foldchange <= 0.8 & pvalue >= 1.3, as.character(data$Gene.Symbol), '')), hjust = 0, vjust = 0)


volcano
dev.off()

### Heatmap ##########################################################################################################################################
library(ComplexHeatmap)
library(plyr)

# Prep the matrix
data.matrix <- as.matrix(data[ , c(4:13)]) # convert the TMT relative abundance into a matrix
genes <- as.vector(data$Gene.Symbol) # grab the gene names
rownames(data.matrix) <- genes
data.matrix <- t(data.matrix) # transpose/flip the matrix

# Making the heatmap
heat <- Heatmap(
    data.matrix,
    name = "Relative abundance (24 hr)",
    
    column_title = "Gene Symbol",
    column_title_side = "top",
    column_title_gp = gpar(cex = 2),
    
    column_names_side = "top",
    column_dend_side = "top",
    column_names_gp = gpar(cex = 0.2), # 2 font size
    
    column_dend_height = unit(5, "cm"), # set the height
    column_km = 10,
    clustering_distance_columns = "maximum",
    clustering_method_columns = "ward.D",
    
    cluster_rows = FALSE,
    show_row_dend = FALSE,
    row_title = "Tandem mass tag",
    row_names_side = "left"
)
dev.off()

### Saving Plots ###################################################################################################################################################

# .png
png("Projects/9-25-19/volcano.png")
volcano
dev.off()
png("Projects/9-25-19/volcano_hr.png", 1000, 1000) # High resolution
volcano
dev.off()
png("Projects/9-25-19/heat.png")
heat
dev.off()
png("Projects/9-25-19/heat_hr.png", 1000, 1000) # High resolution
heat
dev.off()

# .tiff
tiff("Projects/9-25-19/volcano.tiff")
volcano
dev.off()
tiff("Projects/9-25-19/heat.tiff")
heat
dev.off()

# .jpg
jpeg("Projects/9-25-19/volcano.jpg")
volcano
dev.off()
jpeg("Projects/9-25-19/heat,jpg")
heat
dev.off()

# .pdf
pdf("Projects/9-25-19/volcano.pdf")
volcano
dev.off()
pdf("Projects/9-25-19/heat.pdf")
heat
dev.off()







###################################################################################################################################################
# Logan Nguyen
# Fehl Lab
# 9-25-19
# Practicing basic R skills
######################################################################################################################################################


### Reading files into R ############################################################################################################################

rm(list = ls()) # clearing the enviroment
setwd("/home/logan/R") # setting working directory

# Loading libraries
library(readxl) # for loading xlsx files
library(ggplot2) # for plotting

# Loading in the file path, specifically the first sheet
filename <- "/home/logan/Documents/Research/Fehl Lab/Data/9-25-19/ja8b07328_si_002.xlsx"
# data <- read.csv(filename, sep = "\t", header = TRUE) # for .txt files and .csv files
data <- read_excel(filename, sheet = "24h", na = "---") # for .xls and .xlsx files

# Renaming the first columns
names(data)[4:13] <- c("Da", "Db", "Dc", "Dd", "De", "osmi2a", "osmi2b", "osmi2c", "osmi2d", "osmi2e")

# Deleting the first row as we don't need it anymore
data <- data[-c(1), ]

# Converting the 4 character columns to the numeric data type
data <- transform(data, Da = as.numeric(Da))
data <- transform(data, Db = as.numeric(Db))
data <- transform(data, Dc = as.numeric(Dc))
data <- transform(data, Dd = as.numeric(Dd))
data <- transform(data, De = as.numeric(De))
data <- transform(data, osmi2a = as.numeric(osmi2a))
data <- transform(data, osmi2b = as.numeric(osmi2b))
data <- transform(data, osmi2c = as.numeric(osmi2c))
data <- transform(data, osmi2d = as.numeric(osmi2d))
data <- transform(data, osmi2e = as.numeric(osmi2e))


## Basic Data Info ############################################################################################################################
head(data) # looking at the first 6 rows of our data
summary(data) # basic stats on our data
sapply(data, class) # checking the class of each of our columns
dim(data) # dimensions
typeof(data)

### Volcano Plot ############################################################################################################################ 
theme_set(theme_bw())
theme(legend.position = "bottom")

# We need to plot -log(p value) vs fold change, so load in these vector columns
data$'-logpvalue' <- -log(data$p.value) # creating a new column and modify pvalue as we want to plot -log(pvalue)
pvalue <- data$'-logpvalue'
foldchange <- data$fold_change

# Plotting
volcano <- ggplot(data, aes(x = foldchange, y = pvalue, color = data$Quantified.peptides)) + geom_point()
volcano <- volcano + xlim(0.25, 1.75) + ylim(0, 16) + scale_color_gradient(limits = c(0, 30))
volcano <- volcano + theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))
volcano <- volcano + labs(x = "Fold Change", y = "-log(p value)") + ggtitle("24 hr")
volcano <- volcano + labs(color = "# of peptides")

volcano <- volcano + geom_hline(yintercept = 1.3, linetype = "dashed", color = "red")
volcano <- volcano + geom_vline(xintercept = 0.8, linetype = "dashed", color = "red")
volcano <- volcano + geom_vline(xintercept = 1.2, linetype = "dashed", color = "red")

volcano <- volcano + geom_text(aes(label = ifelse(foldchange >= 1.2 & pvalue >= 1.3, as.character(data$Gene.Symbol), '')), hjust = 0, vjust = 0)
volcano <- volcano + geom_text(aes(label = ifelse(foldchange <= 0.8 & pvalue >= 1.3, as.character(data$Gene.Symbol), '')), hjust = 0, vjust = 0)


volcano
dev.off()

### Heatmap ##########################################################################################################################################
library(ComplexHeatmap)
library(plyr)

# Prep the matrix
data.matrix <- as.matrix(data[ , c(4:13)]) # convert the TMT relative abundance into a matrix
genes <- as.vector(data$Gene.Symbol) # grab the gene names
rownames(data.matrix) <- genes
data.matrix <- t(data.matrix) # transpose/flip the matrix

# Making the heatmap
heat <- Heatmap(
    data.matrix,
    name = "Relative abundance (24 hr)",
    
    column_title = "Gene Symbol",
    column_title_side = "top",
    column_title_gp = gpar(cex = 2),
    
    column_names_side = "top",
    column_dend_side = "top",
    column_names_gp = gpar(cex = 0.2), # 2 font size
    
    column_dend_height = unit(5, "cm"), # set the height
    column_km = 10,
    clustering_distance_columns = "maximum",
    clustering_method_columns = "ward.D",
    
    cluster_rows = FALSE,
    show_row_dend = FALSE,
    row_title = "Tandem mass tag",
    row_names_side = "left"
)
dev.off()

### Saving Plots ###################################################################################################################################################

# .png
png("Projects/9-25-19/volcano.png")
volcano
dev.off()
png("Projects/9-25-19/volcano_hr.png", 1000, 1000) # High resolution
volcano
dev.off()
png("Projects/9-25-19/heat.png")
heat
dev.off()
png("Projects/9-25-19/heat_hr.png", 1000, 1000) # High resolution
heat
dev.off()

# .tiff
tiff("Projects/9-25-19/volcano.tiff")
volcano
dev.off()
tiff("Projects/9-25-19/heat.tiff")
heat
dev.off()

# .jpg
jpeg("Projects/9-25-19/volcano.jpg")
volcano
dev.off()
jpeg("Projects/9-25-19/heat,jpg")
heat
dev.off()

# .pdf
pdf("Projects/9-25-19/volcano.pdf")
volcano
dev.off()
pdf("Projects/9-25-19/heat.pdf")
heat
dev.off()

demo(graphics)