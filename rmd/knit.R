knitr::knit(input = "rmd/table.Rmd", 
            output = "rmd/table.md") 

rmarkdown::pandoc_convert(input = here::here("rmd/table.md"),
                          to = "latex",
                          output = here::here("rmd/table.tex"))
