if (!require("devtools")) install.packages("devtools", dependencies = TRUE)
if (!require("pacman")) install.packages("pacman", dependencies = TRUE)
if (!require("plyr")) install.packages("plyr", dependencies = TRUE)
if (!require("dplyr")) install.packages("dplyr", dependencies = TRUE)
if (!require("psych")) install.packages("psych", dependencies = TRUE)
if (!require("car")) install.packages("car", dependencies = TRUE)
if (!require("MASS")) install.packages("MASS", dependencies = TRUE)
if (!require("kableExtra")) install.packages("kableExtra", dependencies = TRUE)
if (!require("rmarkdown")) install.packages("rmarkdown", dependencies = TRUE)
if (!require("knitr")) install.packages("knitr", dependencies = TRUE)
if (!require("pastecs")) install.packages("pastecs", dependencies = TRUE)
if (!require("paran")) install.packages("paran", dependencies = TRUE)
if (!require("ggplot2")) install.packages("ggplot2", dependencies = TRUE)
if (!require("gridExtra")) install.packages("gridExtra", dependencies = TRUE)
if (!require("ggrepel")) install.packages("ggrepel", dependencies = TRUE)
if (!require("ggpubr")) install.packages("ggpubr", dependencies = TRUE)
if (!require("ggbiplot")) install_github("vqv/ggbiplot")

library(pacman)

p_load(devtools, plyr, dplyr, psych, car, MASS, kableExtra, rmarkdown, knitr,
       pastecs, paran, pcaMethods, ggplot2, gridExtra,ggrepel, ggpubr, ggbiplot)

# Load the full data set
load("morpho_complete.Rdata")

# Eliminate factor variables & untransformed weights from full dataset
morpho <- morpho_complete[,-c(1,2,6)]

morpho.stats <- round(stat.desc(morpho, basic = FALSE), digits = 3)
morpho.stats

# Visualize the pairwise correlations of the variables
pairs.panels(
  morpho,
  main = "Wolf Spider Morphometrics - Correlation Summary",
  gap = 0, # set to zero for no gap between plot panels
  lm = TRUE, # draw linear regression lines for pairs
  stars = TRUE, # display significance
  bg = c("red", "blue")[morpho_complete$sex], # color based on sex 
  pch = 21) # data point shape

# Define a standardization function
standardize <- function(x) {(x - mean(x))/sd(x)}

# Eliminate factor variables & untransformed weights from the scaled data
my.scaled.data <- as.data.frame(apply(morpho, 2, standardize))

# Example plot showing the linear relationship between the first two variables
ggplot(my.scaled.data, aes(interoc, cwidth)) +
  geom_point(size = 2) +
  geom_smooth(method = 'lm') +
  ggtitle("Plot of Interoccular Distance and Carapace Width") +
  theme(plot.title = element_text(hjust = 0.5)) # center plot title

# Calculate correlation matrix
my.cor <- cor(my.scaled.data)
my.cor

# Save the eigenvalues of the correlation matrix
my.eigen <- eigen(my.cor)

# Rename matrix rows and columns for easier interpretation
rownames(my.eigen$vectors) <- c("interoc", "cwidth", "clength", "T.weight")
colnames(my.eigen$vectors) <- c("PC1", "PC2", "PC3", "PC4")

eigen.table <- as.data.frame(my.eigen$vectors)
eigen.table

eigen.values <- data.frame(PC = c("PC1", "PC2", "PC3", "PC4"),
                           eigenvalues = my.eigen$values)
eigen.values

#the sum of the eigenvalues equals the total variance of the scaled data
sum(my.eigen$values)

sum(var(my.scaled.data[,1]),
    var(my.scaled.data[,2]),
    var(my.scaled.data[,3]),
    var(my.scaled.data[,4]))

# Amount of variation captured by the PCs in the dataset
pc1.var <- 100*round(my.eigen$values[1]/sum(my.eigen$values), digits = 3)
pc2.var <- 100*round(my.eigen$values[2]/sum(my.eigen$values), digits = 3)
pc3.var <- 100*round(my.eigen$values[3]/sum(my.eigen$values), digits = 3)
pc4.var <- 100*round(my.eigen$values[4]/sum(my.eigen$values), digits = 3)

pc <- data.frame(PC = c("PC1", "PC2", "PC3", "PC4"),
                 Percentage = c(pc1.var, pc2.var, pc3.var, pc4.var))
pc

# The total variation should sum to ~100% depending on rounding error
sum(pc1.var, pc2.var, pc3.var, pc4.var)

# Compute PCA scores
loadings <- my.eigen$vectors
my.scaled.matrix <- as.matrix(my.scaled.data)
# %*% is matrix multiplication
scores <- my.scaled.matrix %*% loadings
sd <- sqrt(my.eigen$values)
rownames(loadings) <- colnames(my.scaled.data)
head(scores)

# Accomplishing the PCA with native R functions
pca_morpho <- prcomp(morpho, center = TRUE, scale. = TRUE)

# Show the variables in the class "prcomp"
ls(pca_morpho)

# Rotations are often referred to as the "loadings" of a PCA.
pca_loadings <- pca_morpho$rotation
pca_loadings

# Summary output of the PCA
pca_summary <- summary(pca_morpho)$importance %>%
  as.data.frame()
pca_summary

# Orthogonality of PCs
pairs.panels(
  pca_morpho$x,
  main = "PCA Correlation Summary",
  gap = 0, # set to zero for no gap between plot panels
  lm = TRUE, # draw linear regression lines for pairs
  stars = TRUE, # display significance
  bg = c("red", "blue")[morpho_complete$sex], # color based on sex 
  pch = 21) # data point shape

# Introduction to the biplot
# Plot without ID labels
g <- ggbiplot(pca_morpho,
              obs.scale = 1,
              var.scale = 1,
              groups = morpho_complete$sex,
              ellipse = TRUE,
              circle = FALSE,
              # draw ellipse around points (+/-) 1 standard deviation
              ellipse.prob = 0.68) +
  scale_color_discrete(name = '') +
  theme(legend.direction = "horizontal", legend.position = "top",
        aspect.ratio = 1)

# Plot with ID labels
h <- ggbiplot(pca_morpho,
              obs.scale = 1,
              var.scale = 1,
              groups = morpho_complete$sex,
              ellipse = TRUE,
              circle = FALSE,
              # draw ellipse around points (+/-) 1 standard deviation
              ellipse.prob = 0.68) +
  scale_color_discrete(name = '') +
  theme(legend.direction = "horizontal", legend.position = "top",
        aspect.ratio = 1) +
  geom_text_repel(aes(label = morpho_complete$id))

figure <- ggarrange(g, h, ncol = 2)

annotate_figure(figure,
               top = text_grob("PCA Biplots - Morphometric Data",
                               color = "black", face = "bold", size = 14),
               fig.lab.face = "bold")

# Choosing points from the biplot for comparison
filter(morpho_complete, id == "336" | id == "339")
paste("minimum carapace length =",
      min(morpho_complete$clength),
      "maximum carapace length =",
      max(morpho_complete$clength),
      sep = ' ')

filter(morpho_complete, id == "360" | id == "55")
paste("minimum weight =",
      min(morpho_complete$weight.mg),
      "maximum weight =",
      max(morpho_complete$weight.mg),
      sep = ' ')

# How many PCs can I retain to explain "enough" variation?
# Compute eigenvalues from PCA object: Kaiser Criterion
pca_morpho$sdev ^ 2

# Scree plot
screeplot(pca_morpho, type = "line", main = "Scree Plot - Morphometric Data")

# Horn’s parallel analysis of principal components/factors
paran(morpho, iterations = 5000, centile = 0, quietly = TRUE, 
    status = FALSE, all = TRUE, cfa = TRUE, graph = TRUE, color = TRUE, 
    col = c("black", "red", "blue"), lty = c(1, 2, 3), lwd = 1, legend = TRUE, 
    file = "", width = 640, height = 640, grdevice = "png", seed = 0)

# Oak Woods dataset
load("OakWoods.Rdata")
boxplot(oak2)

# Save data frame as a matrix and eliminate the first column (character)
oak.matrix <- as.matrix(oak2[,-1])
pca_oak_robust <- pca(oak.matrix, method = "robustPca", nPcs = 27,
                      center = TRUE, scaled = TRUE, scores = TRUE)

# Parallel anaysis of Oak Woods dataset
paran(oak.matrix, iterations = 5000, centile = 0, quietly = TRUE, 
    status = FALSE, all = TRUE, cfa = TRUE, graph = TRUE, color = TRUE, 
    col = c("black", "red", "blue"), lty = c(1, 2, 3), lwd = 1, legend = TRUE, 
    file = "", width = 640, height = 640, grdevice = "png", seed = 0)