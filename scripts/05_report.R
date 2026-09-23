#!/usr/bin/env Rscript
# 보고서 렌더링: report/report.Rmd → report/report.html
source("R/00_setup.R")
rmarkdown::render(proj_path("report", "report.Rmd"), output_dir = proj_path("report"), quiet = TRUE)
cat("report/report.html 생성\n")
