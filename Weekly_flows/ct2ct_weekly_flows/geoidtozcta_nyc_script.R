#Library loading
library(selectr)
library(xml2)
library(rvest)
library(stringr)
library(janitor)
library(readr)
library(tidyr)
library(dplyr)


#ZCTA id table
df_nyc <- read.csv("/home/claudia/Scrivania/Tesi NYC pandemic /Weekly_flows/ct2ct_weekly_flows/nyc_zcta_geoid36.csv")

unique_zcta <- data.frame(unique(df_nyc$ZCTA5))
colnames(unique_zcta) <- c("zcta")
unique_zcta <- unique_zcta%>%arrange(zcta)
unique_zcta$zcta_id <- seq.int(nrow(unique_zcta))-1

df_nyc <- merge(df_nyc,unique_zcta[ , c("zcta","zcta_id")], by.x = "ZCTA5", by.y = "zcta",all.x = TRUE)
summary(df_nyc)

#Loop per caricare i dati da geoid a zcta
df_week <- data.frame()
system.time(
for (i in c(0:19)) {
  df <- read.csv(paste0("/home/claudia/Scrivania/Tesi NYC pandemic /Weekly_flows/ct2ct_weekly_flows/settimana1/w1_file",i,".csv"))
  
  df <- df%>%filter(substr(df$geoid_o,1,2)==36 & substr(df$geoid_d,1,2)==36)

  df$zcta_o <- df_nyc$zcta_id[match(df$geoid_o,df_nyc$GEOID)] #how to do index-match in R like excel
  df$zcta_d <- df_nyc$zcta_id[match(df$geoid_d,df_nyc$GEOID)]
  
  df <- na.omit(df)

  df_week <-rbind(df_week,df)
}
)
summary(df_week)

#write.csv(df_nyc, "/home/claudia/Scrivania/Tesi NYC pandemic /Weekly_flows/ct2ct_weekly_flows/nyc_zcta_geoid_list.csv")
#write.csv(df_week,"/home/claudia/Scrivania/Tesi NYC pandemic /Weekly_flows/ct2ct_weekly_flows/week1_ct2ct_geoid.csv")

#Sum of pop flows with data
traffic <- df_week%>%filter(zcta_o!=zcta_d) %>%group_by(zcta_o,zcta_d) %>%summarise(sum_pop = round(sum(pop_flows)))
summary(traffic)


#now let's create the complete traffic file with all possible combinations of zcta 
traffic_complete <- read.csv("/home/claudia/Scrivania/Tesi NYC pandemic /Weekly_flows/ct2ct_weekly_flows/combinations_zcta_nyc.csv")
summary(traffic_complete)

# df = pd.DataFrame(columns=['zcta_o', 'zcta_d'], index=range(31152)) #non salvo le combinazioni 0-0, 1-1, etc
# t=0;
# for i in range(177):
#   for j in range(177):
#   if(i!=j):
#   df.iloc[t] = [i,j]
# t = t+1

#grow another column "pop flows" with double index matching
traffic_complete$pop_flows <- traffic$sum_pop[match(traffic_complete$zcta_o, traffic$zcta_o) & match(traffic_complete$zcta_d, traffic$zcta_d)]
traffic_complete$pop_flows <- replace_na(traffic_complete$pop_flows,0)

#traffic should be divided by seven since its a weekly average
traffic_complete$pop_flows <- round(traffic_complete$pop_flows/7)
traffic_complete$X <- NULL


#open file of tot population per zcta and create a control file for the traffic out from each zcta
df_pop_US <- read.csv("/home/claudia/Scrivania/pop_US.csv")
test <- df_pop_US #assign the dataframe of whole population in the US to "test"
test <- filter(test,sex==''& race=='') #keep all sexes and races together
test <- as.data.frame(test %>% group_by(zipCode,minAge) %>% summarise(population=sum(population))) #sum of population by zipcode and agegroup
test$minAge <- as.factor(test$minAge) #put minage as factor in order to keep the minage=NA
test <- filter(test,is.na(minAge)) #keep only rows where the age is missing to inlude all age groups
colnames(test)[1] <- "zcta"
test <- test %>% select(zcta,population) #remove the min age variable

group_pop <- merge(x=test, y=unique_zcta, by="zcta", all.x=TRUE)
control_file <- drop_na(group_pop)

#control that for each week pop_flow out from a zcta is lower then total population
control_file$pop_flows <- traffic_complete$pop_flows[match(control_file$zcta_id, traffic_complete$zcta_o)]

control_file$population_2020 <- round(control_file$population*1.02)
control_file$check <- if_else(control_file$pop_flows>control_file$population_2020,1,0)

summary(control_file)

total_population_zcta <- control_file%>%select(zcta_id,population_2020)

#from https://www1.nyc.gov/assets/planning/download/pdf/planning-level/nyc-population/new-population/current-populatiion-estimattes.pdf
#we redefine the total population per zcta updated with the general increasing of 2% of the NYC's population
write.csv(total_population_zcta, "/home/claudia/Scrivania/Tesi NYC pandemic /Weekly_flows/ct2ct_weekly_flows/total_population_zcta.csv")

write.csv(traffic_complete, "/home/claudia/Scrivania/Tesi NYC pandemic /Weekly_flows/ct2ct_weekly_flows/commuting_w1.csv")

write.table(traffic_complete, file="commuting_w1", append = FALSE, sep = " ", dec = ".", row.names = FALSE, col.names = FALSE)
