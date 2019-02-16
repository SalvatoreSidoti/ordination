library(RCurl)

morpho_doc <- getURL("https://docs.google.com/spreadsheets/d/1NkklI07JRY142-0iBB0r1j8dHW3ECd5FeipxMfWWsWI/pub?gid=0&single=true&output=csv")

morpho <- read.csv(textConnection(morpho_doc),
                   header = TRUE, na.strings=c("", "NA"))

# na.strings = c("", "NA")) will insert "NA" in blank cells

# Covert weights g >>> mg
morpho$weight.mg <- 1000*(morpho$weight)

# Eliminate columns: light, weight, motionless
morpho <- morpho[,-c(3,4,8)]

# Remove rows with missing data
# Equalize the number of males and females
morpho <- morpho[-c(14,15,20:25,27,28),]

males <- subset(morpho, sex == "m")
females <- subset(morpho, sex == "f")

### Data Imputation ###
browseURL("https://bioconductor.org/packages/release/bioc/vignettes/pcaMethods/inst/doc/missingValues.pdf")

library(pcaMethods)

# MALES 
# Run pca
morpho_pca <- pca(males[,-c(1,2)], nPcs = 3, method = "ppca")

# Impute missing data
imputed_males <- as.data.frame(completeObs(morpho_pca))

# Round values to 3 decimal places
imputed_males <- round(imputed_males, digits = 3)

# Glue original columns (1,2) back
pieces_males <- list(males[,1:2], imputed_males)

males_complete <- do.call(cbind.data.frame, pieces_males)

# FEMALES 
# Run pca
morpho_pca <- pca(females[,-c(1,2)], nPcs = 3, method = "ppca")

# Impute missing data
imputed_females <- as.data.frame(completeObs(morpho_pca))

# Round values to 3 decimal places
imputed_females <- round(imputed_females, digits = 3)

# Glue original columns (1,2) back
pieces_females <- list(females[,1:2], imputed_females)

females_complete <- do.call(cbind.data.frame, pieces_females)

final_stitch <- list(females_complete, males_complete)

morpho_complete <- do.call(rbind.data.frame, final_stitch)

png(filename = "correlation-summary.png",
    width = 1920, height = 1080, units = "px", pointsize = 36,
    bg = "white", family = "Liberation Sans", restoreConsole = TRUE,
    type = "cairo-png")
pairs.panels(morpho_complete[,-c(1,2,7)], # eliminate factor variables from plots
             main = "Wolf Spider Morphometrics - Correlation Summary",
             gap = 0, # set to zero for no gap between plot panels
             lm = TRUE, # draw linear regression lines for pairs
             stars = TRUE, # display significance
             bg = c("red", "blue")[morpho_complete$sex], # color based on sex 
             pch = 21)
dev.off()

shapiro.test(morpho_complete$weight.mg)

morpho_lm <- lm(weight.mg ~ sex, morpho_complete)
shapiro.test(resid(morpho_lm))

morpho_boxcox <- boxcox(morpho_lm)

best_lambda <- morpho_boxcox$x[which(
  morpho_boxcox$y == max(morpho_boxcox$y))]

best_lambda

# Create a new vector of transformed weight measures:
morpho_complete$T.weight <- (morpho_complete$weight ^ best_lambda - 1)/best_lambda

morpho_complete

morpho_lm_trans <- lm(T.weight ~ sex, data = morpho_complete)

shapiro.test(resid(morpho_lm_trans))

shapiro.test(morpho_complete$T.weight)

png(filename = "histogram-1.png",
    width = 1080, height = 1080, units = "px", pointsize = 36,
    bg = "white", family = "Liberation Sans", restoreConsole = TRUE,
    type = "cairo-png")
hist(morpho_complete$weight.mg, breaks = "FD", col = "light blue")
dev.off()

png(filename = "histogram-2.png",
    width = 1080, height = 1080, units = "px", pointsize = 36,
    bg = "white", family = "Liberation Sans", restoreConsole = TRUE,
    type = "cairo-png")
hist(morpho_complete$T.weight, breaks = "FD", col = "light blue")
dev.off()

shapiro.test(morpho$interoc)
shapiro.test(morpho$cwidth)
shapiro.test(morpho$clength)

save(morpho_complete,
     file = "d:/biology/OSU/spiderpapers/loco/summer 2016/morpho_complete.Rdata")

# Save a local copy to working directory
save(morpho_complete,
     file = "morpho_complete.Rdata")

png(filename = "correlation-summary-2.png",
    width = 1920, height = 1080, units = "px", pointsize = 36,
    bg = "white", family = "Liberation Sans", restoreConsole = TRUE,
    type = "cairo-png")
pairs.panels(morpho_complete[,-c(1,2,6)], # eliminate factor variables from plots
             main = "Wolf Spider Morphometrics - Correlation Summary",
             gap = 0, # set to zero for no gap between plot panels
             lm = TRUE, # draw linear regression lines for pairs
             stars = TRUE, # display significance
             bg = c("red", "blue")[morpho_complete$sex], # color based on sex 
             pch = 21)
dev.off()