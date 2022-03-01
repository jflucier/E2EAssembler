rm(list=ls());gc()
#knitr::opts_chunk$set(fig.path='figs/')
knitr::opts_chunk$set(cache = FALSE, warning = FALSE,
                      message = FALSE, cache.lazy = FALSE)
my_plot_hook <- function(x, options)
  paste("\n", knitr::hook_plot_tex(x, options), "\n")
knitr::knit_hooks$set(plot = my_plot_hook)

library(tidyverse)
library(ggsci)
library(zoo)
library(fGarch)

args <- commandArgs(trailingOnly = TRUE)
print(args[1])
genotype1 <- read.csv(args[1], sep = "\t")

list_genotype <- unique(genotype1[c("genotype")])
for (geno in list_genotype$genotype){
    print(geno)
    telomere_length_wt <- subset(genotype1,  genotype1$genotype == geno & genotype1$telomere != "ML" & genotype1$telomere != "MR")
    wt <- telomere_length_wt
    print(head(wt))

    AllTelo <- wt %>% group_by(telomere) %>% summarise(Average = mean(length), Sd = sd(length), Count = n())
    Resume <- wt %>% group_by(telomere) %>% summarise(Average = mean(length), Sd = sd(length), Count = n(), Min = min(length), Max = max(length), Median = quantile(length, c(0.5)), QuantileBas = quantile(length, c(0.2)), Quantilehaut = quantile(length, c(0.8)))
    d <- wt %>% group_by(telomere) %>% group_map(~ density(.x$length)$x[which.max(density(.x$length)$y)])
    Resume$DensityMax <- c(as.numeric(d))
    print(Resume, width = Inf)

    m <- wt %>% summarise(Mean = mean(length))
    Overall <- Resume %>% summarise(MeanofMean = mean(Average), MeanDensity = mean(DensityMax), MeanQB = mean(QuantileBas), MeanQH = mean(Quantilehaut))
    Overall$mean <- c(m$Mean)
    print(Overall, width = Inf)

    histogram_individualtelowt <- ggplot(data=wt, aes(x=length)) +
      geom_histogram(aes(y = ..density..), color="black", alpha = 0.5)+
      scale_fill_manual(values=c("chocolate1", "darkcyan"))+
      geom_density(lwd=1, alpha = 0.2) +
      scale_color_manual(values=c("chocolate1", "darkcyan")) +
      facet_wrap(~telomere) +
      geom_vline(data=Overall, aes(xintercept= MeanofMean), color="black", lwd=1)+
      geom_vline(data=Overall, aes(xintercept= mean), color="black", lwd=1, linetype="dotted")+
      scale_x_continuous(limits = c(100, 800)) +
      theme_classic()

    pdf(
        file = paste("report/",geno,".pdf",sep = ""),
        width = 8.5,
        height = 11
    )
    print(histogram_individualtelowt)
    dev.off()

}
