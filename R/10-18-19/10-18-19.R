### Script Information ##################################################################################################################################################
# Logan Nguyen
# Fehl Lab
# 10-18-19
# R Bloggers Tutorial

# 6 samples
# Parental = intensity data from breast cancer line SKBR3
# Resistant = intensity data from drug-resistant cell line made from culturing Parental cell line
# ...in the prescence of an EGFR inhibitor
# Publication can be found here: https://europepmc.org/abstract/MED/26883193

### Reading files into R ############################################################################################################################
rm(list = ls()) # clearing the enviroment
setwd("/home/logan/R/10-18-19") # setting working directory

# Loading libraries
library(readxl) # for loading xlsx files
library(ggplot2) # for plotting

# Read in the 1787 by 79 dataframe. Proteins in rows, information (eg. abundance!) by columns
raw = read.delim("proteinGroups.txt", stringsAsFactors = FALSE, colClasses = "character")
grep("^LFQ.intensity", names(raw), value = TRUE)

## [1] "LFQ.intensity.Parental_bR1"  "LFQ.intensity.Parental_bR2" 
## [3] "LFQ.intensity.Parental_bR3"  "LFQ.intensity.Resistant_bR1"
## [5] "LFQ.intensity.Resistant_bR2" "LFQ.intensity.Resistant_bR3"

### Data Cleaning ############################################################################################################################
# Removing False Hits (potential contaminants, reverse proteins, proteins only identified by site)
library(dplyr)
df = dplyr::filter(raw, Potential.contaminant != "+") %>%
    dplyr::filter(Reverse != "+") %>%
    dplyr::filter(Only.identified.by.site != "+")

# note: Q-value column is the probability that the protein is a false hit (typical cutoff is 0.1)

summary(as.numeric(df$Q.value))
##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
## 0.0000000 0.0000000 0.0000000 0.0007193 0.0012759 0.0091429


### Define Protein and Gene ID's ############################################################################################################################
head(df$Protein.IDs)
## [1] "sp|A0AV96|RBM47_HUMAN" "sp|A0FGR8|ESYT2_HUMAN" "sp|A1L0T0|ILVBL_HUMAN"
## [4] "sp|A4D1S0|KLRG2_HUMAN" "sp|A5YKK6|CNOT1_HUMAN" "sp|A6NDG6|PGP_HUMAN"

head(df$Fasta.headers)
## [1] ">sp|A0AV96|RBM47_H  UMAN RNA-binding protein 47 OS=Homo sapiens GN=RBM47 PE=1 SV=2"                               
## [2] ">sp|A0FGR8|ESYT2_HUMAN Extended synaptotagmin-2 OS=Homo sapiens GN=ESYT2 PE=1 SV=1"                             
## [3] ">sp|A1L0T0|ILVBL_HUMAN Acetolactate synthase-like protein OS=Homo sapiens GN=ILVBL PE=1 SV=2"                   
## [4] ">sp|A4D1S0|KLRG2_HUMAN Killer cell lectin-like receptor subfamily G member 2 OS=Homo sapiens GN=KLRG2 PE=1 SV=3"
## [5] ">sp|A5YKK6|CNOT1_HUMAN CCR4-NOT transcription complex subunit 1 OS=Homo sapiens GN=CNOT1 PE=1 SV=2"             
## [6] ">sp|A6NDG6|PGP_HUMAN Glycerol-3-phosphate phosphatase OS=Homo sapiens GN=PGP PE=1 SV=1"


# Loading in the file path, specifically the first sheet
filename <- "/home/logan/R/"
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