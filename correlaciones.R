rm(list=ls())

df_conj_filled_sc <- readRDS("C:/Users/jumar/OneDrive/Escritorio/TFM/Datos/df_conj_filled_sc.rds")
Y <- readRDS("C:/Users/jumar/OneDrive/Escritorio/TFM/Datos/Y.rds")

ind <- which(Y$station==230)

df_madrid<-df_conj_filled_sc[ind,]

library(corrplot)

corrplot(cor(df_madrid[2:46]),
         order = 'hclust',
         method = 'color')

# Annual trend of a geopotencial
var <- 'g300'

plot(df_madrid[, 't'], df_madrid[, var], pch = 19)

library(dplyr)

for (i in 2:46){
  var <- colnames(df_madrid)[i]
  
  media <- df_madrid %>%
    group_by(t) %>%
    summarise(media = mean(.data[[var]], na.rm = T)) %>%
    as.data.frame()
  
  plot(media[, 't'], media[, 'media'], type = 'b', pch = 19, main = var)
}

