if (!require("devtools")) install.packages("devtools", dependencies = TRUE)
if (!require("pacman")) install.packages("pacman", dependencies = TRUE)
if (!require("plyr")) install.packages("plyr", dependencies = TRUE)
if (!require("psych")) install.packages("psych", dependencies = TRUE)
if (!require("car")) install.packages("car", dependencies = TRUE)
if (!require("MASS")) install.packages("MASS", dependencies = TRUE)
if (!require("kableExtra")) install.packages("kableExtra", dependencies = TRUE)
if (!require("rmarkdown")) install.packages("rmarkdown", dependencies = TRUE)
if (!require("knitr")) install.packages("knitr", dependencies = TRUE)
if (!require("pastecs")) install.packages("pastecs", dependencies = TRUE)
if (!require("ggplot2")) install.packages("ggplot2", dependencies = TRUE)
if (!require("ggbiplot")) install_github("vqv/ggbiplot")

p_load(devtools, pacman, plyr, psych, car, MASS, kableExtra,
       rmarkdown, knitr, pastecs, ggplot2, ggbiplot)

load("morpho_complete.Rdata")

