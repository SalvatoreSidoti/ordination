

### PCA ###
load("morpho_complete.Rdata")

morpho <- morpho_complete[,-c(1,2,6)]

pca_morpho <- princomp(morpho, cor = TRUE)

pca_morpho_extra <- PCA(morpho, graph = FALSE)

pca_loadings <- unclass(loadings(pca_morpho)) %>%
  kable(pca_loadings, caption = "Table X: PCA Loadings") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"),
                full_width = F)

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