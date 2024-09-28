
# figure SI ---------------------------------------------------------------

knitr::knit(input = here::here("rmd/figure_si.Rmd"), 
            output = here::here("rmd/figure_si.md")) 

rmarkdown::pandoc_convert(input = here::here("rmd/figure_si.md"),
                          to = "latex",
                          output = here::here("rmd/figure_si.tex"))

# table -------------------------------------------------------------------

knitr::knit(input = here::here("rmd/table.Rmd"), 
            output = here::here("rmd/table.md")) 

rmarkdown::pandoc_convert(input = here::here("rmd/table.md"),
                          to = "latex",
                          output = here::here("rmd/table.tex"))

