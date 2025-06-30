rm(list=ls())

df_conj_filled_sc <- readRDS("C:/Users/jumar/OneDrive/Escritorio/TFM/Datos/df_conj_filled_sc.rds")

library(corrplot)

corrplot(cor(df_conj_filled_sc[2:46]),
         order = 'hclust',
         method = 'color')

