# benchmarks/result.R
#This script is used to calculate the speedup and combine the results from julia and R.
#Load necessary libraries
library(dplyr)
#Read in the result
full_result_julia<-read.csv("result/benchmark_result_full_df_julia.csv")
full_result_r<-read.csv("result/benchmark_result_full_df_r.csv")

medium_result_julia<-read.csv("result/benchmark_result_medium_df_julia.csv")
medium_result_r<-read.csv("result/benchmark_result_medium_df_r.csv")

small_result_julia<-read.csv("result/benchmark_result_small_df_julia.csv")
small_result_r<-read.csv("result/benchmark_result_small_df_r.csv")

#Reorganize the data
full_result<-full_result_julia%>%
  select(TestCase, Time_median)%>%
  rename(Time_median_julia=Time_median)%>%
  left_join(full_result_r)%>%
  select(TestCase, Time_median_julia, Time_median)%>%
  rename(Time_median_r=Time_median)%>%
  filter(!is.na(Time_median_r))%>%
  mutate(Speedup=Time_median_r/Time_median_julia)

medium_result<-medium_result_julia%>%
  select(TestCase, Time_median)%>%
  rename(Time_median_julia=Time_median)%>%
  left_join(medium_result_r)%>%
  select(TestCase, Time_median_julia, Time_median)%>%
  rename(Time_median_r=Time_median)%>%
  filter(!is.na(Time_median_r))%>%
  mutate(Speedup=Time_median_r/Time_median_julia)


small_result<-small_result_julia%>%
  select(TestCase, Time_median)%>%
  rename(Time_median_julia=Time_median)%>%
  left_join(small_result_r)%>%
  select(TestCase, Time_median_julia, Time_median)%>%
  rename(Time_median_r=Time_median)%>%
  filter(!is.na(Time_median_r))%>%
  mutate(Speedup=Time_median_r/Time_median_julia)

write.csv(small_result, "result/small_df_speedup.csv")
write.csv(medium_result, "result/medium_df_speedup.csv")
write.csv(full_result, "result/full_df_speedup.csv")







