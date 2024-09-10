library(portalr)
library(data.table)
library(lubridate)
library(dplyr)


download_observations("path/to/store/data")

rodent_abundance_by_plot <- abundance(path = "path/to/store/data", time = "date", level = "plot")

data<-rodent_abundance_by_plot%>%
  mutate(censusdate=as.Date(censusdate))%>%
  mutate(Year=year(censusdate),
         Month=month(censusdate),
         Day=day(censusdate))%>%
  filter(Year>=2010)%>%
  mutate(Sampling_date_order = dense_rank(censusdate))%>%
  select(Year, Month, Day, Sampling_date_order, everything(), -censusdate, -treatment)%>%
  rename(Plot=plot)

write.csv(data, "~/.julia/dev/MetaCommunityMetrics/data/rodent_abundance_data.csv", row.names=FALSE)



