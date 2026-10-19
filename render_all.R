# Render every case study to a GitHub-flavoured Markdown report
reports <- c("01-competing-risks/competing-risks.Rmd",
             "02-iptw-survival/iptw-survival.Rmd",
             "03-gee-vs-glmm/gee-vs-glmm.Rmd")
for (f in reports) rmarkdown::render(f, envir = new.env())
