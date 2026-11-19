# --- Minimal reinstall of packages after R upgrade ---

#Run in previous R
installed <- as.data.frame(installed.packages())
#save a copy
write.csv(installed, 'previously_installed.csv')


#Run in new R
previously_installed <- read.csv('previously_installed.csv')
package_list <- as.character(previously_installed$Package)
package_list

install.lib <- package_list[!package_list %in% installed.packages()]   
for(lib in install.lib) install.packages(lib, dependencies = TRUE)

#------------------------------Biocmanger-----------------------

# Read your saved list
pkgs <- read.csv("previously_installed.csv", stringsAsFactors = FALSE)$Package

# Make sure helpers are installed
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")

# Reinstall CRAN packages
cran_to_install <- setdiff(pkgs, installed.packages()[,"Package"])
install.packages(cran_to_install, dependencies = TRUE)

# Reinstall Bioconductor packages
BiocManager::install(pkgs, ask = FALSE, update = TRUE)

# Optional: reinstall GitHub remotes (if you have any)
# Example: remotes::install_github("user/repo")
# (You can put your GitHub packages manually in a small vector)
gh_pkgs <- c(
  # "user/repo",
)
if (length(gh_pkgs)) remotes::install_github(gh_pkgs, dependencies = TRUE)



#---------------github--------------------
# 1) Install the required genome package from Bioconductor (hg19)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("BSgenome", "BSgenome.Hsapiens.UCSC.hg19"), ask = FALSE, update = TRUE)


options(repos = BiocManager::repositories())  # ensure Bioc repos are set
options(timeout = 36000)                       # 10h to survive slow links
BiocManager::install("BSgenome.Hsapiens.UCSC.hg19", ask = FALSE, Ncpus = parallel::detectCores())


 # 2) Then install deconstructSigs from GitHub
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("raerose01/deconstructSigs")

# 3) Load to verify
library(deconstructSigs)

#---------------------------------
install.packages("devtools")
library("devtools")

devtools::install_github("gersteinlab/siglasso")

#--------------------
library(signature.tools.lib)




# Fresh R session recommended
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Hsapiens.1000genomes.hs37d5)
library(indelsig.tools.lib)
library(signature.tools.lib)

sessionInfo()


