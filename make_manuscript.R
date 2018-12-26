#system("sudo apt-get install pdftk")
setwd("~/Desktop/ff_study")
library(rmarkdown)
source("FF_figs.R")

rmarkdown::render("FF_manuscript.Rmd")
rmarkdown::render("FF_tables.Rmd")
rmarkdown::render("FF_figs.Rmd")
rmarkdown::render("supl.Rmd")



# getting today's data to add to the filename
x=date()
x=strsplit(x," ")[[1]]

system(paste("pdftk FF_manuscript.pdf FF_tables.pdf FF_figs.pdf supl.pdf cat output submission", 
             x[2],x[3],".pdf", sep = "_"))
