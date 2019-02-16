install.packages('RCurl')
install.packages('plyr')
install.packages('boot')
install.packages('xtable')
install.packages('lubridate')
install.packages('nortest')
install.packages('pwr')
install.packages('lsr')
install.packages('ggplot2')
install.packages('psych')
install.packages('FactoMineR')
install.packages('devtools')
install.packages('ggfortify')

library(RCurl)
library(plyr)
library(boot)
library(xtable)
library(lubridate)
library(nortest)
library(pwr)
library(lsr)
library(ggplot2)
library(psych)
library(FactoMineR)
library(devtools)
library(ggfortify)

morpho_doc <- getURL("https://docs.google.com/spreadsheets/d/1NkklI07JRY142-0iBB0r1j8dHW3ECd5FeipxMfWWsWI/pub?gid=0&single=true&output=csv")

morpho <- read.csv(textConnection(morpho_doc),
                   header = TRUE, na.strings=c("", "NA"))
# na.strings=c("", "NA")) will insert "NA" in blank cells

morpho <- morpho[,-c(3,8:10)]

save(morpho,
     file = "d:/biology/OSU/spiderpapers/loco/summer 2016/morpho_complete.Rdata")

load("d:/biology/OSU/spiderpapers/loco/summer 2016/morpho.Rdata")

# Subset 'morpho' by motionless spiders (TRUE) versus those that moved (FALSE)
morpho_move <- subset(morpho, motionless!="TRUE")
morpho_still <- subset(morpho, motionless!="FALSE")

save(morpho_move,
     file = "f:/biology/OSU/spider papers/loco/summer 2016/morpho_move.Rdata")

load("d:/biology/OSU/spiderpapers/loco/summer 2016/morpho_move.Rdata")

# Subset by metric variables
metrics <- morpho_move[,c(4:7)]

# Grams to milligrams
metrics$weight <- 1000*(metrics$weight)

### Numerical Summary ### ----
df <- summary(metrics)
df_table <- xtable(df, digits = 4)
print.xtable(df_table, file = "F:\\biology\\OSU\\spider papers\\loco\\summer 2016\\df.html", type = "html")

pdf("F:/biology/OSU/ddig/box_weight.pdf",
    width = 3, height = 6, useDingbats = FALSE)
boxplot(metrics$weight, ylim = c(40,140), main = "Weight (mg)")
dev.off()

# ln() transfom weight variable (non-normal to normal)
metrics$weight <- log(metrics$weight)

pdf("F:/biology/OSU/ddig/box_weight_trans.pdf", width = 3, height = 6, useDingbats = FALSE)
boxplot(metrics$weight, main = "log[Weight (mg)]")
dev.off()

pdf("F:/biology/OSU/ddig/box.pdf", width = 6, height = 6, useDingbats = FALSE)
boxplot(metrics[,c(2,3,4)], main = "Linear metrics")
dev.off()

cor_box <- metrics[,c(2,3,4)]
boxplot(cor_box, main = "Linear metrics (mm)")

# Requires summarySE function
browseURL("http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/#Helper functions")

summarySE(metrics, measurevar = "weight")
summarySE(metrics, measurevar = "interoc")
summarySE(metrics, measurevar = "cwidth")
summarySE(metrics, measurevar = "clength")

# metrics matrix with p-values
cor.test.p <- function(x)
{
  FUN <- function(x, y) cor.test(x, y)[["p.value"]]
  z <- outer(
    colnames(x), 
    colnames(x), 
    Vectorize(function(i,j) FUN(x[,i], x[,j]))
  )
  dimnames(z) <- list(colnames(x), colnames(x))
  z
}

cor_matrix <- as.data.frame(cor.test.p(metrics))
df <- xtable(cor_matrix, digits = -4)
print.xtable(df, file = "F:\\biology\\OSU\\spider papers\\loco\\summer 2016\\df.html",
             type = "html")

# Scatterplot matix with metrics and p-values
# https://www.r-bloggers.com/scatter-plot-matrices-in-r/
panel.cor <- function(x, y, digits = 2, cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  # metrics coefficient
  r <- cor(x, y)
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste("r= ", txt, sep = "")
  text(0.5, 0.6, txt)
  
  # p-value calculation
  p <- cor.test(x, y)$p.value
  txt2 <- format(c(p, 0.123456789), digits = digits)[1]
  txt2 <- paste("p= ", txt2, sep = "")
  if(p<0.01) txt2 <- paste("p= ", "<0.01", sep = "")
  text(0.5, 0.4, txt2)
}

windows(7,7)

pdf("F:/biology/OSU/ddig/corplot_raw.pdf", width = 6, height = 6, useDingbats = FALSE)
pairs(metrics, upper.panel = panel.cor)
dev.off()

### PCA ###
pca_morpho <- princomp(metrics, cor = TRUE)

pca_morpho_extra <- PCA(metrics, graph = FALSE)

pca_loadings <- xtable(unclass(loadings(pca_morpho)), digits = 4)
print.xtable(pca_loadings,
             file = "F:\\biology\\OSU\\spider papers\\loco\\summer 2016\\pca_loadings.html",
             type = "html")

loadings(pca_morpho)
# gives the coefficients of the linear combination
# subtitute the observed values into the linear combination to obtain the PCA "scores"
# aka eigenvalues. The eigenvalues can be though of as the weights, or the magnitude of
# the eigenvector

windows(7,7)

pdf(file = "F:/biology/OSU/ddig/biplot.pdf")
biplot(pca_morpho,
       main = "Biplot: Variables & Indivuduals",
       xlab = "PC1",
       ylab = "PC2")
abline(h = 0, col = "blue")
abline(v = 0, col = "green")
dev.off()

# Nice scree plot function using ggplot2
ggscreeplot <- function(pcobj, type = c('pev', 'cev')) 
{
  type <- match.arg(type)
  d <- pcobj$sdev^2
  yvar <- switch(type, 
                 pev = d / sum(d), 
                 cev = cumsum(d) / sum(d))
  
  yvar.lab <- switch(type,
                     pev = 'proportion of explained variance',
                     cev = 'cumulative proportion of explained variance')
  
  df <- data.frame(PC = 1:length(d), yvar = yvar)
  
  ggplot(data = df, aes(x = PC, y = yvar)) + 
    xlab('principal component number') + ylab(yvar.lab) +
    geom_point() + geom_path()
}

# Scree plot
pdf("F:/biology/OSU/ddig/screeplot.pdf", width = 6, height = 6)
ggscreeplot(pca_morpho) +
  ggtitle("Scree Plot: Principle Components")
theme_bw() +
  theme(panel.grid.major.y = element_line(color = "light grey")) +
  theme(panel.grid.major.x = element_line(color = "light grey")) +
  theme(plot.title = element_text(size=16,face="bold", hjust=0.5, vjust=1)) +
  theme(axis.title.x = element_text(size=14,face="bold", hjust=0.5, vjust=-0.2)) +
  theme(axis.title.y = element_text(size=14,face="bold", vjust=1)) +
  theme(axis.text.y = element_text(angle=0,size=12,face="bold", hjust=.5)) +
  theme(axis.text.x = element_text(angle=0,size=12,face="bold", hjust=.5))
dev.off()

pca_coord <- xtable(pca_morpho_extra$var$coord, digits = 7)
print.xtable(pca_coord,
             file = "F:\\biology\\OSU\\spider papers\\loco\\summer 2016\\pca_coord.html",
             type = "html")

pca_contrib <- xtable(pca_morpho_extra$var$contrib/100, digits = 7)
print.xtable(pca_contrib,
             file = "F:\\biology\\OSU\\spider papers\\loco\\summer 2016\\pca_contrib.html",
             type = "html")

pca_contrib <- data.frame(
  metric = c("weight", "interoc", "cwidth", "clength"),
  C1 = c(0.1406505, 0.2338513, 0.2838265, 0.3416717),
  C2 = c(0.808481362, 0.068396475, 0.121259561, 0.001862602))

pca_contrib_c1 <- c(pca_contrib$C1)

pdf("F:/biology/OSU/ddig/barplot1.pdf", width = 5, height = 6)
ggplot(pca_contrib, aes(metric, C1)) +
  geom_bar(fill = "dark grey", width = 0.8, stat = "identity") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.36)) +
  # above moves bars flush with x-axis with headroom on y-axis
  ggtitle("Proportion of Influence - PC1") +
  ylab("Proportion") +
  xlab("Metric") +
  geom_text(aes(label = c("0.137", "0.236", "0.288", "0.339")),
            vjust = -0.5, size = 5, color="black") +
  theme_bw() +
  theme(panel.grid.major.y = element_line(color = "light grey")) +
  theme(panel.grid.major.x = element_line(color = "light grey")) +
  theme(plot.title = element_text(size=18,face="bold", hjust=0.5, vjust=1)) +
  theme(axis.title.x = element_text(size=16,face="bold", hjust=0.5, vjust=-0.2)) +
  theme(axis.title.y = element_text(size=16,face="bold", vjust=1)) +
  theme(axis.text.y = element_text(angle=0,size=14,face="bold", hjust=.5)) +
  theme(axis.text.x = element_text(angle=0,size=14,face="bold", hjust=.5))
dev.off()

pdf("F:/biology/OSU/ddig/barplot2.pdf", width = 5, height = 6)
ggplot(pca_contrib, aes(metric, C2)) +
  geom_bar(fill = "dark grey", width = 0.8, stat = "identity") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.86)) +
  # above moves bars flush with x-axis with headroom on y-axis
  ggtitle("Proportion of Influence - PC2") +
  ylab("Proportion") +
  xlab("Metric") +
  geom_text(aes(label = c("0.822", "0.052", "0.122", "0.004")),
            vjust = -0.5, size = 5, color="black") +
  theme_bw() +
  theme(panel.grid.major.y = element_line(color = "light grey")) +
  theme(panel.grid.major.x = element_line(color = "light grey")) +
  theme(plot.title = element_text(size=18,face="bold", hjust=.5, vjust=1)) +
  theme(axis.title.x = element_text(size=16,face="bold", hjust=.5, vjust=-0.2)) +
  theme(axis.title.y = element_text(size=16,face="bold", vjust=1)) +
  theme(axis.text.y = element_text(angle=0,size=14,face="bold", hjust=.5)) +
  theme(axis.text.x = element_text(angle=0,size=14,face="bold", hjust=.5))
dev.off()

# weighted geometric mean - higest correlation with PC1
wgm <- data.frame(
  weight = exp(metrics$weight) ^ pca_contrib_c1[1], # back transform from log scale
  interoc = metrics$interoc ^ pca_contrib_c1[2],
  cwidth = metrics$cwidth ^ pca_contrib_c1[3],
  clength = metrics$clength ^ pca_contrib_c1[4])

wgm$size <- apply(wgm, 1, prod)
wgm$size_scaled <- scale(wgm$size)
wgm_prod <- apply(wgm, 1, prod)
wgm_calc <- scale(wgm_prod)

pc1 <- prcomp(metrics, scale.=TRUE)$x[,1] #scaled "scores" of the PCA

pdf("F:/biology/OSU/ddig/wgm.pdf", width = 7, height = 7, useDingbats = FALSE)
plot(pc1, wgm_calc, xlab = "PC1 Scores", ylab = "metrics - Weighted Geometric Mean", main = "metrics vs. PC1")
reg4 <- lm(wgm_calc ~ pc1)
abline(reg4, col = "red")
dev.off()

cor(cbind(pc1, wgm_calc))
cor.test(pc1,wgm_calc)

m <- mean(wgm$size)
std <- sd(wgm$size)

pdf("F:/biology/OSU/ddig/wgm_dist.pdf", width = 7, height = 7, useDingbats = FALSE)
h <- hist(wgm$size, main = "Distribution of Size Scores", xlab = "Size Scores",
          breaks = "FD", xaxs="i", yaxs="i", las = 1, ylim = c(0,1.6), freq = FALSE)
curve(dnorm(x, mean=m, sd=std), xlim = c(3,4.4),
      col="black", lwd = 2, add = TRUE)
dev.off()

shapiro.test(wgm$size)

pdf("F:/biology/OSU/ddig/wgm_dist_qq.pdf", width = 7, height = 7, useDingbats = FALSE)
qqnorm(wgm$size, main = "Normal Q-Q Plot - Size Scores")
qqline(wgm$size, col = "red")
dev.off()

# geomtric mean
gm <- apply(metrics, 1, prod)^(1/4)
pc1 <- prcomp(metrics, scale.=TRUE)$x[, 1] #scaled "scores" of the PCA
plot(pc1, gm)
reg1 <- lm(gm~pc1)
abline(reg1, col = "red")
cor(cbind(pc1, gm))
cor.test(pc1,gm)

# arithmetic mean
am <- apply(metrics, 1, mean)
pc1 <- prcomp(metrics, scale.=TRUE)$x[, 1] #scaled "scores" of the PCA
plot(pc1, am)
reg2 <- lm(am~pc1)
abline(reg2, col = "red")
cor(cbind(pc1, am))
cor.test(pc1,am)

# scaled arithmetic mean
sam <- scale(metrics)
sam <- apply(metrics, 1, mean)
pc1 <- prcomp(metrics, scale.=TRUE)$x[, 1] #scaled "scores" of the PCA
plot(pc1, sam)
reg3 <- lm(sam~pc1)
abline(reg3, col = "red")
cor(cbind(pc1, sam))
cor.test(pc1, sam)