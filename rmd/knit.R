
# figure SI ---------------------------------------------------------------

knitr::knit(input = here::here("rmd/figure_si.Rmd"), 
            output = here::here("rmd/figure_si.md")) 

rmarkdown::pandoc_convert(input = here::here("rmd/figure_si.md"),
                          to = "latex",
                          output = here::here("tex/figure_si.tex"),
                          options = c("--variable=graphics"))

## remove "\\pandocbounded"
txt <- readLines("tex/figure_si.tex")
txt <- gsub("\\\\pandocbounded\\{([^}]*)\\}", "\\1", txt)
writeLines(txt, "tex/figure_si.tex")