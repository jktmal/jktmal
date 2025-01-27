setwd("/home/jasiek/Dropbox/mcbs/project/")

options(mc.cores = parallel::detectCores())
options(mc.cores = 4)

source("/home/jasiek/Dropbox/mcbs/project/scripts/data_prep.R")
source("/home/jasiek/Dropbox/mcbs/project/scripts/Indiana_fit.R")
source("/home/jasiek/Dropbox/mcbs/project/scripts/Ohio_fit.R")
