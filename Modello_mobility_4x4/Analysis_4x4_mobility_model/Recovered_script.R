#Library loading
library(selectr)
library(xml2)
library(rvest)
library(stringr)
library(janitor)
library(readr)
library(tidyr)
library(dplyr)


df_2 <-data.frame()

for (i in 0:99) {
   df <- tail(read.csv(paste0("/home/claudia/Scrivania/Tesi NYC pandemic /Modello_mobility_prova4x4/OUTputMOBILE1/outputMATRICES",i,".txt")),3)[1,]
   df <- as.data.frame(df)
   df <- df %>% separate(df, c("Manhattan", "Bronx","Brooklyn", "Queens","extra"),sep = " ") %>% select(Manhattan,Bronx,Brooklyn,Queens)
   df_2 <- rbind(df_2,df)
   }

df_2$Manhattan <- as.numeric(df_2$Manhattan)
df_2$Bronx <- as.numeric(df_2$Bronx)
df_2$Brooklyn <- as.numeric(df_2$Brooklyn)
df_2$Queens <- as.numeric(df_2$Queens)

summary(df_2)

write.csv(df_2,"/home/claudia/Scrivania/Tesi NYC pandemic /Modello_mobility_prova4x4/Analysis_4x4_mobility_model/novarcomm_2.5.csv",fileEncoding = "UTF-8", row.names = FALSE)


boxplot(df_2)
boxplot(df_2$Manhattan)
boxplot(df_2$Bronx)
   
