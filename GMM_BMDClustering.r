library(mixtools)
library(tidyverse)

deg.full_2 <- deg.full_1 %>%
  filter(set == "troph1h")


mix <- normalmixEM(deg.full_2$BMD.zSD, k = 2)
  
ggplot(deg.full_2, aes(x = BMD.zSD)) +
  geom_histogram(binwidth = 4, fill = "grey") +
  mapply(
    function(mean, sd, lambda, n, binwidth) {
      stat_function(
        fun = function(x) {
          (dnorm(x, mean = mean, sd = sd)) * 
n * binwidth * lambda
        }
      )
    },
    mean = mix[["mu"]], #mean
    sd = mix[["sigma"]], #standard deviation
    lambda = mix[["lambda"]], #amplitude
    n = length(deg.full_2$BMD.zSD), #sample size
    binwidth = 4 #binwidth used for histogram
  ) +
  labs(title = "GMM: 1h Trophs", y = "n (Genes)", x = "BMD (nM DHA)") +  theme_cowplot(font_size = 20)
