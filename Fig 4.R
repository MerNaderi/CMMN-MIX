library(GGally)
library(gclus)
data("wine")
Data = wine

Data[, 1] = as.character(Data[, 1])


##### Mine

GG = ggpairs(Data,          # Data frame
        columns = 14:1,        # Columns
        aes(color = c("Class 1", "Class 2", "Class 3")[class],  # Color by group (cat. variable)
            alpha = 0.8), upper = list(continuous = wrap("cor", size = 2.5)))
ggsave(GG, file="name.png", width = 13,
  height = 13, device="png")





        