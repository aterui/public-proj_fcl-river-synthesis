#' DESCRIPTION:
#' Set figure theme

default_theme <- theme_bw() +
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = alpha("steelblue", 0.05)),
        axis.title = element_text(size = 14),
        axis.text =  element_text(size = 10))
