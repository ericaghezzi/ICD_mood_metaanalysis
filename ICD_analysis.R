#### SET-UP ####
#load packages
library(readr)
library(tidyverse)
library(metafor)
library(ggthemes)
library(estmeansd)
library(grid)

#grand SD calculation 
grand.sd   <- function(S, M, N) {sqrt(weighted.mean(S^2 + M^2, N) -
                                        weighted.mean(M, N)^2)}
#### GATHER DATA ####

#set working directory
setwd("R:/GitHub/ICD_mood/")

#import data
data <- read_csv("data_forupload.csv", col_types = cols(days = col_number(),
                                              grant_timepoint = col_number(),
                                              percentage_with_condition = col_number(),
                                              N_with_condition = col_number(),
                                              N_total = col_number(),
                                              average = col_number(),
                                              variance = col_number(),
                                              range_ll = col_number(),
                                              range_ul = col_number(),
                                              N_total_icd2 = col_number(),
                                              average_icd2 = col_number(),
                                              variance_icd2 = col_number(),
                                              variance_range_ul_icd2 = col_number(),
                                              variance_range_ll_icd2 = col_number(),
                                              N_total_icd3 = col_number(),
                                              average_icd3 = col_number(),
                                              variance_icd3 = col_number(),
                                              N_total_icd4 = col_number(),
                                              average_icd4 = col_number(),
                                              variance_icd4 = col_number(),
                                              time_average = col_number(),
                                              time_variance = col_number(),
                                              time_variance_ll = col_number(),
                                              time_variance_ul = col_number(),
                                              time_average2 = col_number(),
                                              time_variance2 = col_number(),
                                              time_variance_ll2 = col_number(),
                                              time_variance_ul2 = col_number(), 
                                              time_average3 = col_number(),
                                              time_variance3 = col_number(),
                                              time_variance_ll3 = col_number(),
                                              time_variance_ul3 = col_number(),
                                              time_group1_N = col_number(),
                                              time_group2_N = col_number(),
                                              time_group3_N = col_number(),
                                              N_total_comparison = col_number(),
                                              average_comparison = col_number(),
                                              variance_comparison = col_number(),
                                              range_ll_comparison = col_number(),
                                              range_ul_comparison = col_number(),
                                              yi = col_number(),
                                              vi = col_number(),
                                              timepoint = col_number())) %>%
  filter(!is.na(covidence_id))

#create id variable (unique id for each row)
data <- data %>%
  mutate(id = row_number())

#### WRANGLE: TIMEPOINT DATA ####

#transform median to mean/SD (group 1,2,3 ICD timepoint data)
for (i in 1:length(data$id)) {
  if(data[i, "time_average_type"] == "median" && !is.na(data[i, "time_average_type"])) {
    if(data[i, "time_variance_range_type"] == "range" && !is.na(data[i, "time_variance_range_type"])) { #s1 (min, med, max)
      # mean and sd calc
      #choose the correct N for group 1
      if(!is.na(data[i, "time_group1_N"])) {
        n <- data[[i, "time_group1_N"]]
      }
      else {
        n <- data[[i, "N_total"]]
      }
      
      mean_calc <- qe.mean.sd(min.val = data[[i, "time_variance_ll"]], 
                              med.val = data[[i, "time_average"]], 
                              max.val = data[[i, "time_variance_ul"]], 
                              n = n)
      
      data[i, "time_average"] <- as.numeric(mean_calc["est.mean"])
      data[i, "time_average_type"] <- as.character("mean")
      data[i, "time_variance"] <- as.numeric(mean_calc["est.sd"])
      data[i, "time_variance_type"] <- as.character("SD")
    }
    if(data[i, "time_variance_range_type"] == "IQR" && !is.na(data[i, "time_variance_range_type"])) { #s2 (q1,q3,med)
      ##del mean and sd calc
      #choose the correct N for group 1
      if(!is.na(data[i, "time_group1_N"])) {
        n <- data[[i, "time_group1_N"]]
      }
      else {
        n <- data[[i, "N_total"]]
      }
      
      mean_calc <- qe.mean.sd(q1.val = data[[i, "time_variance_ll"]], 
                              med.val = data[[i, "time_average"]], 
                              q3.val = data[[i, "time_variance_ul"]], 
                              n = n)
      
      data[i, "time_average"] <- as.numeric(mean_calc["est.mean"])
      data[i, "time_average_type"] <- as.character("mean")
      data[i, "time_variance"] <- as.numeric(mean_calc["est.sd"])
      data[i, "time_variance_type"] <- as.character("SD")
    }
  }
  if(data[i, "time_average_type"] == "median" && !is.na(data[i, "time_group2"])) { #do the same for group 2 timepoint
    if(data[i, "time_average_type"] == "median" && !is.na(data[i, "time_average_type"])) {
      if(data[i, "time_variance_range_type"] == "range" && !is.na(data[i, "time_variance_range_type"])) { #s1 (min, med, max)
        # mean and sd calc
        mean_calc <- qe.mean.sd(min.val = data[[i, "time_variance_ll2"]], 
                                med.val = data[[i, "time_average2"]], 
                                max.val = data[[i, "time_variance_ul2"]], 
                                n = data[[i, "time_group2_N"]])
        
        data[i, "time_average2"] <- as.numeric(mean_calc["est.mean"])
        data[i, "time_variance2"] <- as.numeric(mean_calc["est.sd"])
      }
      if(data[i, "time_variance_range_type"] == "IQR" && !is.na(data[i, "time_variance_range_type"])) { #s2 (q1,q3,med)
        ##del mean and sd calc
        mean_calc <- qe.mean.sd(q1.val = data[[i, "time_variance_ll2"]], 
                                med.val = data[[i, "time_average2"]], 
                                q3.val = data[[i, "time_variance_ul2"]], 
                                n = data[[i, "time_group2_N"]])
        
        data[i, "time_average2"] <- as.numeric(mean_calc["est.mean"])
        data[i, "time_variance2"] <- as.numeric(mean_calc["est.sd"])
      }
    }
  }
  if(data[i, "time_average_type"] == "median" && !is.na(data[i, "time_group3"])) { #do the same for group 3 timepoint
    if(data[i, "time_average_type"] == "median" && !is.na(data[i, "time_average_type"])) {
      if(data[i, "time_variance_range_type"] == "range" && !is.na(data[i, "time_variance_range_type"])) { #s1 (min, med, max)
        # mean and sd calc
        mean_calc <- qe.mean.sd(min.val = data[[i, "time_variance_ll3"]], 
                                med.val = data[[i, "time_average3"]], 
                                max.val = data[[i, "time_variance_ul3"]], 
                                n = data[[i, "time_group3_N"]])
        
        data[i, "time_average3"] <- as.numeric(mean_calc["est.mean"])
        data[i, "time_variance3"] <- as.numeric(mean_calc["est.sd"])
      }
      if(data[i, "time_variance_range_type"] == "IQR" && !is.na(data[i, "time_variance_range_type"])) { #s2 (q1,q3,med)
        ##del mean and sd calc
        mean_calc <- qe.mean.sd(q1.val = data[[i, "time_variance_ll3"]], 
                                med.val = data[[i, "time_average3"]], 
                                q3.val = data[[i, "time_variance_ul3"]], 
                                n = data[[i, "time_group3_N"]])
        
        data[i, "time_average3"] <- as.numeric(mean_calc["est.mean"])
        data[i, "time_variance3"] <- as.numeric(mean_calc["est.sd"])
      }
    }
  }
}

#pooled mean/sd for timepoint data reported across multiple groups    

for (i in 1:length(data$id)) {
  if(data[i, "time_average_type"] == "mean" && !is.na(data[i, "time_average_type"])) { 
    if(!is.na(data[i, "time_group3"])) {
      mean <- c(data[[i, "time_average"]], data[[i, "time_average2"]], data[[i, "time_average3"]])
      n <- c(data[[i, "time_group1_N"]], data[[i, "time_group2_N"]], data[[i, "time_group2_N"]])
      sd <- c(data[[i, "time_variance"]], data[[i, "time_variance2"]], data[[i, "time_variance3"]])
      
      data[i, "time_average"] <- weighted.mean(mean, n)
      data[i, "time_variance"] <- grand.sd(sd,mean,n)
    }
    if(!is.na(data[i, "time_group2"]) && is.na(data[i, "time_group3"])) {
      mean <- c(data[[i, "time_average"]], data[[i, "time_average2"]])
      n <- c(data[[i, "time_group1_N"]], data[[i, "time_group2_N"]])
      sd <- c(data[[i, "time_variance"]], data[[i, "time_variance2"]])
      
      data[i, "time_average"] <- weighted.mean(mean, n)
      data[i, "time_variance"] <- grand.sd(sd,mean,n)
    }
  }
}

#calculate days represented by 'broad' timepoint data
for(i in 1:length(data$id)) {
  if(data[i, "time_broad_units"] == "day" && !is.na(data[i, "time_broad_units"])) {
    data[i,"days"] <- data[i, "time_broad_numerical"]
  }
  if(data[i, "time_broad_units"] == "week" && !is.na(data[i, "time_broad_units"])) {
    data[i,"days"] <- (data[i, "time_broad_numerical"]*7)
  }
  if(data[i, "time_broad_units"] == "month" && !is.na(data[i, "time_broad_units"])) {
    data[i,"days"] <- (data[i, "time_broad_numerical"]*30.417)
  }
  if(data[i, "time_broad_units"] == "year" && !is.na(data[i, "time_broad_units"])) {
    data[i,"days"] <- (data[i, "time_broad_numerical"]*365)
  }
}

#calculate days represented by 'average' timepoint data
for(i in 1:length(data$id)) {
  if(data[i, "time_average_units"] == "day" && !is.na(data[i, "time_average_units"])) {
    data[i,"days"] <- data[i, "time_average"]
  }
  if(data[i, "time_average_units"] == "week" && !is.na(data[i, "time_average_units"])) {
    data[i,"days"] <- (data[i, "time_average"]*7)
  }
  if(data[i, "time_average_units"] == "month" && !is.na(data[i, "time_average_units"])) {
    data[i,"days"] <- (data[i, "time_average"]*30.417)
  }
  if(data[i, "time_average_units"] == "year" && !is.na(data[i, "time_average_units"])) {
    data[i,"days"] <- (data[i, "time_average"]*365)
  }
}

#use day data to allocate into timepoint groups (groups without broad/average data already manually completed)
#(1 = pre (<discharge), 2 = discharge-5 months, 3 = 6-12months, 4 = 12months+, 5 = unspecified (already classified))
for(i in 1:length(data$id)) {
  if(is.na(data[i, "timepoint"])) {
    if(data[i, "days"] <0 && !is.na(data[i, "days"])) {
      data[i, "timepoint"] <- 1
    }
    if(data[i, "days"] >0 && data[i, "days"] < 182.5 && !is.na(data[i, "days"])) {
      data[i, "timepoint"] <- 2
    }
    if(data[i, "days"] >= 182.5 && data[i, "days"] < 365 && !is.na(data[i, "days"])) {
      data[i, "timepoint"] <- 3
    }
    if(data[i, "days"] >= 365 && !is.na(data[i, "days"])) {
      data[i, "timepoint"] <- 4
    }
  }
}

#create second timepoint variable with slightly different coding (groups without broad/average data already manually completed)
#(1 = pre (<0days), 2 = 0-5 months, 3 = 6-12months, 4 = 12months+, 5 = unspecified (already classified))
for(i in 1:length(data$id)) {
  if(is.na(data[i, "timepoint2"])) {
    if(data[i, "days"] <0 && !is.na(data[i, "days"])) {
      data[i, "timepoint2"] <- 1
    }
    if(data[i, "days"] >0 && data[i, "days"] < 182.5 && !is.na(data[i, "days"])) {
      data[i, "timepoint2"] <- 2
    }
    if(data[i, "days"] >= 182.5 && data[i, "days"] < 365 && !is.na(data[i, "days"])) {
      data[i, "timepoint2"] <- 3
    }
    if(data[i, "days"] >= 365 && !is.na(data[i, "days"])) {
      data[i, "timepoint2"] <- 4
    }
  }
}

#### WRANGLE: CONTINUOUS DATA ####

#transform median/mean 95%CI/mean SEM to mean/SD (group 1,2 average data)
for (i in 1:length(data$id)) {
  if(data[i, "average_type"] == "median" && !is.na(data[i, "average_type"])) {
    if(data[i, "range_type"] == "range" && !is.na(data[i, "range_type"])) { #s1 (min, med, max)
      # mean and sd calc
      mean_calc <- qe.mean.sd(min.val = data[[i, "range_ll"]], 
                              med.val = data[[i, "average"]], 
                              max.val = data[[i, "range_ul"]], 
                              n = data[[i, "N_total"]])
      
      data[i, "average"] <- as.numeric(mean_calc["est.mean"])
      data[i, "average_type"] <- as.character("mean")
      data[i, "variance"] <- as.numeric(mean_calc["est.sd"])
      data[i, "variance_type"] <- as.character("SD")
    }
    if(data[i, "range_type"] == "IQR" && !is.na(data[i, "range_type"])) { #s2 (q1,q3,med)
      ##del mean and sd calc
      mean_calc <- qe.mean.sd(q1.val = data[[i, "range_ll"]], 
                              med.val = data[[i, "average"]], 
                              q3.val = data[[i, "range_ul"]], 
                              n = data[[i, "N_total"]])
      
      data[i, "average"] <- as.numeric(mean_calc["est.mean"])
      data[i, "average_type"] <- as.character("mean")
      data[i, "variance"] <- as.numeric(mean_calc["est.sd"])
      data[i, "variance_type"] <- as.character("SD")
    }
  }
  if(data[i, "average_type_icd2"] == "median" && !is.na(data[i, "average_type_icd2"])) {
    if(data[i, "variance_range_type_icd2"] == "range" && !is.na(data[i, "variance_range_type_icd2"])) { #s1 (min, med, max)
      # mean and sd calc
      mean_calc <- qe.mean.sd(min.val = data[[i, "variance_range_ll_icd2"]], 
                              med.val = data[[i, "average_icd2"]], 
                              max.val = data[[i, "variance_range_ul_icd2"]], 
                              n = data[[i, "N_total_icd2"]])
      
      data[i, "average_icd2"] <- as.numeric(mean_calc["est.mean"])
      data[i, "average_type_icd2"] <- as.character("mean")
      data[i, "variance_icd2"] <- as.numeric(mean_calc["est.sd"])
      data[i, "variance_range_type_icd2"] <- as.character("SD")
    }
    if(data[i, "variance_range_type_icd2"] == "IQR" && !is.na(data[i, "variance_range_type_icd2"])) { #s2 (q1,q3,med)
      ##del mean and sd calc
      mean_calc <- qe.mean.sd(q1.val = data[[i, "variance_range_ll_icd2"]], 
                              med.val = data[[i, "average_icd2"]], 
                              q3.val = data[[i, "variance_range_ul_icd2"]], 
                              n = data[[i, "N_total_icd2"]])
      
      data[i, "average_icd2"] <- as.numeric(mean_calc["est.mean"])
      data[i, "average_type_icd2"] <- as.character("mean")
      data[i, "variance_icd2"] <- as.numeric(mean_calc["est.sd"])
      data[i, "variance_range_type_icd2"] <- as.character("SD")
    }
  }
  if(data[i, "variance_type_icd2"] == "SEM" && !is.na(data[i, "variance_type_icd2"])) { #calculate standard deviation from SEM data
    data[i, "variance_icd2"] <- data[[i, "variance"]]*sqrt(data[[i, "N_total_icd2"]])
    data[i, "variance_type_icd2"] <- as.character("SD")
  }
  if(data[i, "variance_type"] == "SEM" && !is.na(data[i, "variance_type"])) { #calculate standard deviation from SEM data
    data[i, "variance"] <- data[[i, "variance"]]*sqrt(data[[i, "N_total"]])
    data[i, "variance_type"] <- as.character("SD")
  }
  if(data[i, "variance_type"] == "SE" && !is.na(data[i, "variance_type"])) { #calculate standard deviation from SEM data
    data[i, "variance"] <- data[[i, "variance"]]*sqrt(data[[i, "N_total"]])
    data[i, "variance_type"] <- as.character("SD")
  }
}

#pooled mean/sd for average data reported across multiple groups   

for (i in 1:length(data$id)) {
  if(data[i, "average_type"] == "mean" && !is.na(data[i, "average_type"])) { #only run for mean data
    if(data[i, "icd4"] == "yes" && !is.na(data[i, "icd4"])) {
    mean <- c(data[[i, "average"]], data[[i, "average_icd2"]], data[[i, "average_icd3"]], data[[i, "average_icd4"]])
    n <- c(data[[i, "N_total"]], data[[i, "N_total_icd2"]], data[[i, "N_total_icd3"]], data[[i, "N_total_icd4"]])
    sd <- c(data[[i, "variance"]], data[[i, "variance_icd2"]], data[[i, "variance_icd3"]], data[[i, "variance_icd4"]])
    
    data[i, "average"] <- weighted.mean(mean, n)
    data[i, "variance"] <- grand.sd(sd,mean,n)
    data[i, "N_total"] <- sum(data[[i, "N_total"]], data[[i, "N_total_icd2"]], data[[i, "N_total_icd3"]], data[[i, "N_total_icd4"]])
    }
    if(data[i, "icd3"] == "yes" && !is.na(data[i, "icd3"]) && is.na(data[i, "icd4"])) {
    mean <- c(data[[i, "average"]], data[[i, "average_icd2"]], data[[i, "average_icd3"]])
    n <- c(data[[i, "N_total"]], data[[i, "N_total_icd2"]], data[[i, "N_total_icd3"]])
    sd <- c(data[[i, "variance"]], data[[i, "variance_icd2"]], data[[i, "variance_icd3"]])
    
    data[i, "average"] <- weighted.mean(mean, n)
    data[i, "variance"] <- grand.sd(sd,mean,n)
    data[i, "N_total"] <- sum(data[[i, "N_total"]], data[[i, "N_total_icd2"]], data[[i, "N_total_icd3"]])
    }
    if(data[i, "icd2"] == "yes" && !is.na(data[i, "icd2"]) && is.na(data[i, "icd3"])) {
    mean <- c(data[[i, "average"]], data[[i, "average_icd2"]])
    n <- c(data[[i, "N_total"]], data[[i, "N_total_icd2"]])
    sd <- c(data[[i, "variance"]], data[[i, "variance_icd2"]])
    
    data[i, "average"] <- weighted.mean(mean, n)
    data[i, "variance"] <- grand.sd(sd,mean,n)
    data[i, "N_total"] <- sum(data[[i, "N_total"]], data[[i, "N_total_icd2"]])
    }
  }
}
    
#### WRANGLE: CATEGORICAL DATA ####
#calculate N with condition for all categorical data (if % or N without only reported)
for (i in 1:length(data$id)) {
  if(data[i, "data_type"] == "cat" && !is.na(data[i, "data_type"])) {
    if(is.na(data[i, "N_with_condition"])) {
      if(!is.na(data[i, "percentage_with_condition"])) {
        data[i, "N_with_condition"] <- (((data[i, "percentage_with_condition"])/100) * data[i, "N_total"])
      }
      if(!is.na(data[i, "N_without_condition"])) {
        data[i, "N_with_condition"] <- (data[i, "N_total"] - data[i, "N_without_condition"])
      }
    }
  }
}

#### WRANGLE: COMPARISON DATA ####

#transform median to mean/SD (group 1,2 comparison group data)
#not currently transforming M/SEM data (as there is none)

for (i in 1:length(data$id)) {
  if(data[i, "average_type_comparison"] == "median" && !is.na(data[i, "average_type_comparison"])) {
    if(data[i, "range_type_comparison"] == "range" && !is.na(data[i, "range_type_comparison"])) { #s1 (min, med, max)
      # mean and sd calc
      mean_calc <- qe.mean.sd(min.val = data[[i, "range_ll_comparison"]], 
                              med.val = data[[i, "average_comparison"]], 
                              max.val = data[[i, "range_ul_comparison"]], 
                              n = data[[i, "N_total_comparison"]])
      
      data[i, "average_comparison"] <- as.numeric(mean_calc["est.mean"])
      data[i, "average_type_comparison"] <- as.character("mean")
      data[i, "variance_comparison"] <- as.numeric(mean_calc["est.sd"])
      data[i, "variance_type_comparison"] <- as.character("SD")
    }
    if(data[i, "range_type_comparison"] == "IQR" && !is.na(data[i, "range_type_comparison"])) { #s2 (q1,q3,med)
      ##del mean and sd calc
      mean_calc <- qe.mean.sd(q1.val = data[[i, "range_ll_comparison"]], 
                              med.val = data[[i, "average_comparison"]], 
                              q3.val = data[[i, "range_ul_comparison"]], 
                              n = data[[i, "N_total_comparison"]])
      
      data[i, "average_comparison"] <- as.numeric(mean_calc["est.mean"])
      data[i, "variance_comparison"] <- as.numeric(mean_calc["est.sd"])
    }
  }
  if(data[i, "average_type_comparison"] == "median" && !is.na(data[i, "comparison2_details"])) { #do the same for group 2 timepoint
    if(data[i, "time_average_type"] == "median" && !is.na(data[i, "time_average_type"])) {
      if(data[i, "range_type_comparison"] == "range" && !is.na(data[i, "range_type_comparison"])) { #s1 (min, med, max)
        # mean and sd calc
        mean_calc <- qe.mean.sd(min.val = data[[i, "range_ll_comparison2"]], 
                                med.val = data[[i, "average_comparison2"]], 
                                max.val = data[[i, "range_ul_comparison2"]], 
                                n = data[[i, "N_total_comparison2"]])
        
        data[i, "average_comparison2"] <- as.numeric(mean_calc["est.mean"])
        data[i, "variance_comparison2"] <- as.numeric(mean_calc["est.sd"])
      }
      if(data[i, "range_type_comparison"] == "IQR" && !is.na(data[i, "range_type_comparison"])) { #s2 (q1,q3,med)
        ##del mean and sd calc
        mean_calc <- qe.mean.sd(q1.val = data[[i, "range_ll_comparison2"]], 
                                med.val = data[[i, "average_comparison2"]], 
                                q3.val = data[[i, "range_ul_comparison2"]], 
                                n = data[[i, "N_total_comparison2"]])
        
        data[i, "average_comparison2"] <- as.numeric(mean_calc["est.mean"])
        data[i, "variance_comparison2"] <- as.numeric(mean_calc["est.sd"])
      }
    }
  }
}

#pooled mean/sd for timepoint data reported across multiple groups 
for (i in 1:length(data$id)) {
  if(data[i, "average_type_comparison"] == "mean" && !is.na(data[i, "average_type_comparison"])) { 
    if(!is.na(data[i, "comparison2_details"])) {
      mean <- c(data[[i, "average_comparison"]], data[[i, "average_comparison2"]])
      n <- c(data[[i, "N_total_comparison"]], data[[i, "N_total_comparison2"]])
      sd <- c(data[[i, "variance_comparison"]], data[[i, "variance_comparison2"]])
      
      data[i, "average_comparison"] <- weighted.mean(mean, n)
      data[i, "variance_comparison"] <- grand.sd(sd,mean,n)
      data[i, "N_total_comparison"] <- sum(data[[i, "N_total_comparison"]], data[[i, "N_total_comparison2"]])
    }
  }
}

#calculate N with condition for all categorical data (if % or N without only reported)
for (i in 1:length(data$id)) {
  if(data[i, "data_type"] == "cat" && !is.na(data[i, "data_type"])) {
    if(data[i, "comparison_group_yesno"] == "yes" && !is.na(data[i, "comparison_group_yesno"])) 
    if(is.na(data[i, "N_with_condition_comparison"])) {
      if(!is.na(data[i, "percentage_with_condition_comparison"])) {
        data[i, "N_with_condition_comparison"] <- (((data[i, "percentage_with_condition_comparison"])/100) * data[i, "N_total_comparison"])
      }
      if(!is.na(data[i, "N_without_condition_comparison"])) {
        data[i, "N_with_condition_comparison"] <- (data[i, "N_total_comparison"] - data[i, "N_without_condition_comparison"])
      }
    }
  }
}



#### CREATE ANALYSIS GROUPS ####

### MEASURE ###
#create data frame which includes all mood measures
any <- data %>%
  mutate(measure_n = 4)

#create measure variable which categorises mood measures
#1 = anxiety, 2 = depression, 3 = ptsd, 0 = any others (ad - they will be encompassed by measure_n = 4 above)
data <- data %>%
  mutate(measure_n = ifelse(measure_short == "a", 1, ifelse(measure_short == "d", 2, ifelse(measure_short == "p", 3, 0))))

#bind the two data frames together so all can be grouped by measure_n as 1,2,3 and then all are included in 4
data <- rbind(data,any)

### TIMEPOINTS ###
#1 = pre-discharge, 2 = discharge-5m, 3 = 6-12m, 4 = 12m+, 5= unspecified post, 6 = pre-ICD, 7 = post-ICD (including unspecified post)

#create data frame with all post-ICD timepoints (using timepoint2 to categorise as timepoint2 = 1 represents all pre ICD data)
post <- data %>%
  filter(timepoint2 > 1) %>%
  mutate(timepoint = 7)

#create data frame with all pre-ICD timepoints
pre <- data %>%
  filter(timepoint2 == 1) %>%
  mutate(timepoint = 6)

#bind data, pre and post df together so can be grouped based on timepoint
data <- rbind(data,pre,post)

### DATA TYPE/SEVERITY ###
#1 = cont, 2 = cat(clinically sig), 3 = cont/cat(clinically sig), 4 = cat(mild+), 5 = cont/cat(mild+), 6 = all cat/cont

#create cont only df (include cohen's d data)
cont <- data %>%
  filter(data_type == "cont" | data_type == "d")

#create cat only df
cat <- data %>%
  filter(data_type == "cat")

##make df for cat(clinically significant cutoff/diagnosis only)
#remove all mild/mild+
no_mild <- data %>%
  filter(data_type == "cat")%>%
  filter(mild != 1 | is.na(mild)) %>%
  filter(mild_plus != 1 | is.na(mild_plus))
#filter this data for only those which are clinically significant
clin_data <- no_mild %>%
  filter(clin_sig == 1)
#filter for those which are mood diagnoses
diag_data <- no_mild %>%
  filter(cat_measure == "diagnosis")
#cat data(clinically significant/diagnosis)
clin <- rbind(clin_data, diag_data) 

#create df of all at least mild cat data
mild <- data %>%
  filter(data_type == "cat") %>%
  filter(!is.na(mild_plus))

#create df of all at least mod cat data
mod <- data %>%
  filter(data_type == "cat") %>%
  filter(!is.na(mod_plus))

#create df for each data type/severity
type1 <- cont %>%
  mutate(type = 1)

type2 <- rbind(clin_data, diag_data) %>%
  mutate(type = 2)

type3 <- rbind(cont, clin_data, diag_data) %>%
  mutate(type = 3)

type4 <- mild %>%
  mutate(type = 4)

type5 <- rbind(cont, mild) %>%
  mutate(type = 5)

type6 <- data %>%
  mutate(type = 6)

type7 <- mod %>%
  mutate(type = 7)

#bind all types together so can be grouped based on type
data <- rbind(type1,type2,type3,type4,type5,type6,type7)

###CREATE SUBGROUPS DF ####

sub_data <- data %>%
  filter(!is.na(subgroup_level))

data <- data %>%
  filter(is.na(subgroup_level))

#### BETWEEN: ISOLATE DATA ####

#create df with just comparison group data
comp_data <- data %>%
  filter(timepoint != is.na(T)) %>%
  filter(comparison_group_yesno == "yes") %>%
  mutate(g = as.numeric(NA), #create variables for effect size data
         g_vi = as.numeric(NA),
         g_sei = as.numeric(NA),
         or = as.numeric(NA),
         or_vi = as.numeric(NA),
         or_sei = as.numeric(NA))


#### BETWEEN: EFFECT SIZE ####

#effect size calculations
for (i in 1:length(comp_data$id)) {
  if(comp_data[i, "average_type"] == "mean" && comp_data[i, "variance_type"] == "SD" && !is.na(comp_data[i, "average_type"]) && !is.na(comp_data[i, "variance_type"])) { #take all mean/SD data
    #run effect size (Hedges g) calc (use mean, sd, N for both groups)
    es <- escalc(measure = "SMD", 
                 m1i = comp_data[[i, "average"]], 
                 sd1i = comp_data[[i, "variance"]], 
                 n1i = comp_data[[i, "N_total"]], 
                 m2i = comp_data[[i, "average_comparison"]], 
                 sd2i = comp_data[[i, "variance_comparison"]], 
                 n2i = comp_data[[i, "N_total_comparison"]])
    es_summary <- summary.escalc(es)
    
    #copy across effect size into data frame
    comp_data[i, "g"] <- as.numeric(es$yi[1])
    comp_data[i, "g_vi"] <- as.numeric(es$vi[1])
    comp_data[i, "g_sei"] <- as.numeric(es_summary$sei)
  }
  if(comp_data[i,"data_type"]=="d" && !is.na(comp_data[i,"data_type"])){
    #convert cohen's d to effect size
    comp_data[i,"g"] <- as.numeric((1 - 3/(4*(comp_data[[i,"group1_cohensd_N"]]+comp_data[[i,"group2_cohensd_N"]]-2) - 1)) * comp_data[[i,"cohen_d_comparison"]])
    comp_data[i,"g_vi"] <- as.numeric(1/comp_data[[i,"group1_cohensd_N"]]+ 1/comp_data[[i,"group2_cohensd_N"]] + comp_data[[i,"g"]]^2/(2*(comp_data[[i,"group1_cohensd_N"]]+comp_data[[i,"group2_cohensd_N"]])))
    
    #add transformation for if comparison group is higher (*-1)
    if(comp_data[i,"higherscore_cohensd"]=="comparison" && !is.na(comp_data[i,"higherscore_cohensd"])){
      comp_data[i,"g"] <- comp_data[[i,"g"]]*-1
    }
  }
  if(comp_data[i, "data_type"] == "cat" && !is.na(comp_data[i, "data_type"])) { #take all categorical data
    #run effect size (OR) calc (assign N with condition, N without, total N for both groups)
    or <- escalc(measure = "OR", 
                 ai = comp_data[[i, "N_with_condition"]], 
                 bi = comp_data[[i, "N_total"]] - comp_data[[i, "N_with_condition"]], 
                 n1i = comp_data[[i, "N_total"]], 
                 ci = comp_data[[i, "N_with_condition_comparison"]], 
                 di = comp_data[[i, "N_total_comparison"]] - comp_data[[i, "N_with_condition_comparison"]], 
                 n2i = comp_data[[i, "N_total_comparison"]])
    or_summary <- summary.escalc(or) #use this to calculate sei (which is also able to be calculated for raw OR data and so both can be pooled together)
    
    #run effect size (Hedges g) calc (assign N with condition, N without, total N for both groups)
    g_or <- escalc(measure = "OR2DN", 
                   ai = comp_data[[i, "N_with_condition"]], 
                   bi = comp_data[[i, "N_total"]] - comp_data[[i, "N_with_condition"]], 
                   n1i = comp_data[[i, "N_total"]], 
                   ci = comp_data[[i, "N_with_condition_comparison"]], 
                   di = comp_data[[i, "N_total_comparison"]] - comp_data[[i, "N_with_condition_comparison"]], 
                   n2i = comp_data[[i, "N_total_comparison"]])
    g_or_summary <- summary.escalc(g_or)
    
    #copy across effect size into data frame
    comp_data[i, "or"] <- as.numeric(or$yi[1])
    comp_data[i, "or_vi"] <- as.numeric(or$vi[1]) 
    comp_data[i, "or_sei"] <- as.numeric(or_summary$sei)
    
    comp_data[i, "g"] <- as.numeric(g_or$yi[1]) 
    comp_data[i, "g_vi"] <- as.numeric(g_or$vi[1])
    comp_data[i, "g_sei"] <- as.numeric(g_or_summary$sei)
  }
  if(comp_data[i, "type"] == 2 | comp_data[i, "type"] == 4 | comp_data[i, "type"] == 7 && !is.na(comp_data[i, "data_type"])) { #take all solo categorical data types and rerun calculation but make g/g_vi OR/OR_vi
    #run effect size (OR) calc (assign N with condition, N without, total N for both groups)
    or <- escalc(measure = "OR",
                 ai = comp_data[[i, "N_with_condition"]],
                 bi = comp_data[[i, "N_total"]] - comp_data[[i, "N_with_condition"]],
                 n1i = comp_data[[i, "N_total"]],
                 ci = comp_data[[i, "N_with_condition_comparison"]],
                 di = comp_data[[i, "N_total_comparison"]] - comp_data[[i, "N_with_condition_comparison"]],
                 n2i = comp_data[[i, "N_total_comparison"]])
    or_summary <- summary.escalc(or) #use this to calculate sei (which is also able to be calculated for raw OR data and so both can be pooled together)

    #copy across effect size into data frame
    comp_data[i, "g"] <- as.numeric(or$yi[1])
    comp_data[i, "g_vi"] <- as.numeric(or$vi[1])
    comp_data[i, "g_sei"] <- as.numeric(or_summary$sei)
  }
}

outliers <- comp_data %>%
  filter(g>=2 || g<=-2)
#no outliers

#### BETWEEN: ANALYSIS GROUPS ####
#1 = all, 2 = partners, 3 = cardiac, 4 = general

#create df for each comparison group type
partners <- comp_data %>%
  filter(comparison_group == "partners") %>%
  mutate(comp_group = 2)

cardiac <- comp_data %>%
  filter(comparison_group == "cardiac") %>%
  mutate(comp_group = 3)

general <- comp_data %>%
  filter(comparison_group == "general") %>%
  mutate(comp_group = 4)

comp_data <- comp_data %>%
  mutate(comp_group = 1)

#bind the 3 df together so can be grouped in analyses
comp_data <- rbind(comp_data,partners,cardiac,general)

#### BETWEEN: ANALYSIS DF ####

#average across necessary variables to make complete analysis file (each study represented up to once per analysis)
b_analysis <- comp_data %>%
  group_by(covidence_id, measure_n, timepoint, type, comp_group) %>%
  mutate(g = mean(g, na.rm = T),
         g_vi = mean(g_vi, na.rm = T)) %>% #obtain average g and st err within each study (average across domains), removing NAs from calculation
  filter(row_number()==1) %>% #make each study only appear on one row (per measure, timepoint, type, group)
  ungroup() %>%
  filter(!is.na(g))

b_total_studies <- b_analysis %>%
  filter(comp_group !=1) %>%
  group_by(covidence_id, comp_group) %>%
  filter(row_number()==1)

b_total_partners <- b_total_studies %>%
  filter(comp_group ==2)

b_total_cardiac <- b_total_studies %>%
  filter(comp_group ==3)

b_total_general <- b_total_studies %>%
  filter(comp_group ==4)

#### PREVALENCE: ISOLATE  DATA ####

#create df with just prevalence data (type2 [clin/diagnosis cat data] + type4 [mild+ cat data] + type7 [mod+ cat data])
prev_data <- rbind(type2
                   ,type4
                   #,type7
                   ) %>%
  filter(timepoint != is.na(T)) %>%
  filter(is.na(comparison_group_yesno)) %>%
  mutate(g = as.numeric(NA), #create variables for effect size data
         g_vi = as.numeric(NA),
         g_sei = as.numeric(NA),
         or = as.numeric(NA),
         or_vi = as.numeric(NA),
         or_sei = as.numeric(NA))

#### PREVALENCE: EFFECT SIZE ####

#effect size calculations
for(i in 1:length(prev_data$id)) {
  es <- escalc(xi = prev_data[[i, "N_with_condition"]], 
               ni = prev_data[[i, "N_total"]], 
               measure ="PR") 
  prev_data[i, "g"] <- as.numeric(es$yi[1])
  prev_data[i, "g_vi"] <- as.numeric(es$vi[1])
}

outliers <- prev_data %>%
  filter(g>=2 || g<=-2)
#no outliers

#### PREVALENCE: ANALYSIS DF ####

#average across necessary variables to make complete analysis file (each study represented up to once per analysis)
p_analysis <- prev_data %>%
  group_by(covidence_id, measure_n, timepoint, type) %>%
  mutate(g = mean(g, na.rm = T),
         g_vi = mean(g_vi, na.rm = T)) %>% #obtain average g and st err within each study (average across domains), removing NAs from calculation
  filter(row_number()==1) %>% #make each study only appear on one row (per measure, timepoint, type)
  ungroup() %>%
  filter(!is.na(g))

#### C PREVALENCE: ISOLATE DATA ####

#create df with just comparison group prevalence data (type2 [clin/diagnosis cat data] + type4 [mild+ cat data])
compprev_data <- rbind(type2
                       ,type4
                       #,type7
                       ) %>%
  filter(timepoint != is.na(T)) %>%
  filter(comparison_group_yesno == "yes") %>%
  mutate(g = as.numeric(NA), #create variables for effect size data
         g_vi = as.numeric(NA),
         g_sei = as.numeric(NA),
         or = as.numeric(NA),
         or_vi = as.numeric(NA),
         or_sei = as.numeric(NA))

#### C PREVALENCE: ANALYSIS GROUPS ####
#1 = all, 2 = partners, 3 = cardiac, 4 = general

#create df for each comparison group type
partners <- compprev_data %>%
  filter(comparison_group == "partners") %>%
  mutate(comp_group = 2)

cardiac <- compprev_data %>%
  filter(comparison_group == "cardiac") %>%
  mutate(comp_group = 3)

general <- compprev_data %>%
  filter(comparison_group == "general") %>%
  mutate(comp_group = 4)

compprev_data <- compprev_data %>%
  mutate(comp_group = 1)

#bind the 3 df together so can be grouped in analyses
compprev_data <- rbind(compprev_data,partners,cardiac,general)

#### C PREVALENCE: EFFECT SIZE  ####

#effect size calculations
for(i in 1:length(compprev_data$id)) {
  es <- escalc(xi = compprev_data[[i, "N_with_condition_comparison"]], 
               ni = compprev_data[[i, "N_total_comparison"]], 
               measure ="PR") 
  compprev_data[i, "g"] <- as.numeric(es$yi[1])
  compprev_data[i, "g_vi"] <- as.numeric(es$vi[1])
}

outliers <- compprev_data %>%
  filter(g>=2 || g<=-2)
#no outliers

#### C PREVALENCE: ANALYSIS DF ####

#average across necessary variables to make complete analysis file (each study represented up to once per analysis)
cp_analysis <- compprev_data %>%
  group_by(covidence_id, measure_n, timepoint, comp_group, type) %>%
  mutate(g = mean(g, na.rm = T),
         g_vi = mean(g_vi, na.rm = T)) %>% #obtain average g and st err within each study (average across domains), removing NAs from calculation
  filter(row_number()==1) %>% #make each study only appear on one row (per measure, timepoint, type)
  ungroup() %>%
  filter(!is.na(g))


#### WITHIN: ISOLATE DATA ####

#get one row per study per test/subscale per timepoint
with_data <- rbind(type1,type2) %>% #change this to data when want to do all at once
  filter(is.na(subgroup_level)) %>%
  group_by(covidence_id, type,measure_n,test_used, test_subscale, timepoint) %>%
  mutate(mean = weighted.mean(average, variance),
         sd = grand.sd(variance, average, N_total),
         n_with = mean(N_with_condition),
         n = min(N_total)) %>% #take minimum N value to be conservative. Should all be the same N but could be cases where not equal
  select(covidence_id,type,timepoint,measure_n,test_used, test_subscale,mean,sd,n_with,n) %>%
  arrange(covidence_id,timepoint,measure_n,test_used, test_subscale) %>%
  filter(row_number()==1) #one row per measure (a,d,p,ad), per test/sub scale per timepoint per study

#### WITHIN: ANALYSIS GROUPS ####

#create each timepoint comparison group df
# 1 = 1+2, 2 = 2+3, 3 = 3+4, 4 = 6+7

#1+2
w_1_2 <- with_data %>%
  filter(timepoint == 1 | timepoint == 2) %>%
  group_by(covidence_id, type, measure_n, test_used, test_subscale) %>%
  mutate(nrow = n()) %>%
  filter(nrow>1) %>% #this ensures for each combination of the above that there is 2 rows (1 for each timepoint)
  mutate(w_id = cur_group_id()) %>%#this creates a unique group ID. Thought I would need this to join by but could have just entered my grouping variables into the by=c("") argument
  select(covidence_id, type, w_id,measure_n, test_used, test_subscale, timepoint,mean,sd,n_with,n)

w_1 <- w_1_2 %>%
  filter(timepoint == 1)

w_2 <- w_1_2 %>%
  filter(timepoint == 2)

w_comp_1 <- full_join(w_1,w_2,by=c("w_id","covidence_id","type","measure_n","test_used", "test_subscale")) %>%
  mutate(within = 1)

#2+3
w_2_3 <- with_data %>%
  filter(timepoint == 2 | timepoint == 3) %>%
  group_by(covidence_id, type, measure_n, test_used, test_subscale) %>%
  mutate(nrow = n()) %>%
  filter(nrow>1) %>% #this ensures for each combination of the above that there is 2 rows (1 for each timepoint)
  mutate(w_id = cur_group_id()) %>%#this creates a unique group ID. Thought I would need this to join by but could have just entered my grouping variables into the by=c("") argument
  select(covidence_id, type, w_id,measure_n, test_used, test_subscale, timepoint,mean,sd,n_with,n)

w_2 <- w_2_3 %>%
  filter(timepoint == 2)

w_3 <- w_2_3 %>%
  filter(timepoint == 3)

w_comp_2 <- full_join(w_2,w_3,by=c("w_id","covidence_id","type","measure_n","test_used", "test_subscale")) %>%
  mutate(within = 2)

#3+4
w_3_4 <- with_data %>%
  filter(timepoint == 3 | timepoint == 4) %>%
  group_by(covidence_id, type, measure_n, test_used, test_subscale) %>%
  mutate(nrow = n()) %>%
  filter(nrow>1) %>% #this ensures for each combination of the above that there is 2 rows (1 for each timepoint)
  mutate(w_id = cur_group_id()) %>%#this creates a unique group ID. Thought I would need this to join by but could have just entered my grouping variables into the by=c("") argument
  select(covidence_id, type, w_id,measure_n, test_used, test_subscale, timepoint,mean,sd,n_with,n)

w_3 <- w_3_4 %>%
  filter(timepoint == 3)

w_4 <- w_3_4 %>%
  filter(timepoint == 4)

w_comp_3 <- full_join(w_3,w_4,by=c("w_id","covidence_id","type","measure_n","test_used", "test_subscale")) %>%
  mutate(within = 3)

#6+7
#create timepoint 6-7 comparison df
w_5_6 <- with_data %>%
  filter(timepoint == 6 | timepoint == 7) %>%
  group_by(covidence_id, type, measure_n, test_used, test_subscale) %>%
  mutate(nrow = n()) %>%
  filter(nrow>1) %>% #this ensures for each combination of the above that there is 2 rows (1 for each timepoint)
  mutate(w_id = cur_group_id()) %>%#this creates a unique group ID. Thought I would need this to join by but could have just entered my grouping variables into the by=c("") argument
  select(covidence_id, type, w_id,measure_n, test_used, test_subscale, timepoint,mean,sd,n_with,n)

w_5 <- w_5_6 %>%
  filter(timepoint == 6)

w_6 <- w_5_6 %>%
  filter(timepoint == 7)

w_comp_4 <- full_join(w_5,w_6,by=c("w_id","covidence_id","type","measure_n","test_used", "test_subscale")) %>%
  mutate(within = 4)

#JOIN ALL COMPARISON TIME POINT DFs
with_data_join <- rbind(w_comp_1,w_comp_2,w_comp_3,w_comp_4) %>%
  mutate(r = 0.6) #actually impute correct r values

#### WITHIN: EFFECT SIZE ####

#get the cont data so es for cat/cont can be calculated separately
with_data_type1 <- with_data_join %>%
  filter(type == 1)

with_data_type2 <- with_data_join %>%
  filter(type == 2)

#calc SMCR es for all cont data
with_data_type1 <- escalc(measure = "SMCR", m1i = mean.x, m2i = mean.y, sd1i = sd.x, ni = n.x, ri = r, data = with_data_type1)

with_data <- with_data_type1

outliers <- with_data %>%
  filter(yi>=2 || yi<=-2)
#no outliers

#### WITHIN: ANALYSIS DF ####

#AVERAGE ACROSS EACH MEASURE TYPE WITHIN EACH COVIDENCE ID, WITHIN GROUP, TYPE (create w_analysis file)
w_analysis <- with_data %>%
  group_by(covidence_id, type, measure_n, within) %>% #add in type here if later add analyses based on data type
  mutate(g = yi,
         g_vi = vi) %>%
  mutate(g = mean(g, na.rm = T),
         g_vi = mean(g_vi, na.rm = T)) %>% #obtain average g and st err within each study (average across domains), removing NAs from calculation
  filter(row_number()==1) %>% #make each study only appear on one row (per measure, timepoint, type)
  ungroup() %>%
  filter(!is.na(g))

#### SUBGROUP: ISOLATE DATA ####

#get one row per study per test/subscale per timepoint (by imposed definition)
sub_data_edit <- sub_data %>% 
  ###***** THE BELOW ROW OF DATA makes it so all timepoints are included, not just pre/post
  mutate(timepoint3 = timepoint) %>%
  filter(type == 1 | type == 2) %>% #cat and cont data
  group_by(covidence_id, type,measure_n,test_used, test_subscale, timepoint3,subgroup_level) %>%
  mutate(mean = weighted.mean(average, variance),
         sd = grand.sd(variance, average, N_total),
         n_with = mean(N_with_condition),
         n = min(N_total)) %>% #take minimum N value to be conservative. Should all be the same N but could be cases where not equal
  select(covidence_id,type,timepoint3,subgroup_level,measure_n,test_used, test_subscale,mean,sd,n_with,n) %>%
  arrange(covidence_id,timepoint3,measure_n,test_used, test_subscale,subgroup_level) %>%
  filter(row_number()==1) #one row per measure (a,d,p,ad), per test/sub scale per timepoint per study


#### SUBGROUP: ANALYSIS GROUPS ####

#create subgroup ID variable for the below subgroup dfs
#shock/noshock = 1; sex = 2; indication = 3

# create shock df with each shock v no shock on a new row of data
shock_sub <- sub_data_edit %>%
  filter(subgroup_level == "shock" | subgroup_level == "no shock") %>%
  group_by(covidence_id,timepoint3, type, measure_n, test_used, test_subscale) %>% #group by study timepoint AND timepoint initially (there will be multiple timepoints within timepoint if just group by that)
  mutate(nrow = n())%>%
  filter(nrow>1) %>% #this ensures for each combination of the above that there is 2 rows (1 for each subgroup level)
  select(covidence_id, type,subgroup_level,measure_n, test_used, test_subscale,timepoint3,mean,sd,n_with,n)

shock <- shock_sub %>%
  filter(subgroup_level == "shock")

noshock <- shock_sub %>%
  filter(subgroup_level == "no shock")

shock_dat <- full_join(shock,noshock,by=c("covidence_id","type","measure_n","test_used", "test_subscale","timepoint3")) %>%
  mutate(subgroup = 1)

# create sex df with each female v male on a new row of data
sex_sub <- sub_data_edit %>%
  filter(subgroup_level == "female" | subgroup_level == "male") %>%
  group_by(covidence_id,timepoint3, type, measure_n, test_used, test_subscale) %>% #group by study timepoint AND timepoint initially (there will be multiple timepoints within timepoint if just group by that)
  mutate(nrow = n())%>%
  filter(nrow>1) %>% #this ensures for each combination of the above that there is 2 rows (1 for each subgroup level)
  select(covidence_id, type,subgroup_level,measure_n, test_used, test_subscale,timepoint3,mean,sd,n_with,n)

female <- sex_sub %>%
  filter(subgroup_level == "female")

male <- sex_sub %>%
  filter(subgroup_level == "male")

sex_dat <- full_join(female,male,by=c("covidence_id","type","measure_n","test_used", "test_subscale","timepoint3")) %>%
  mutate(subgroup = 2)

# create indication df with each primary v secondary on a new row of data
ind_sub <- sub_data_edit %>%
  filter(subgroup_level == "primary" | subgroup_level == "secondary") %>%
  group_by(covidence_id,timepoint3, type, measure_n, test_used, test_subscale) %>% #group by study timepoint AND timepoint initially (there will be multiple timepoints within timepoint if just group by that)
  mutate(nrow = n())%>%
  filter(nrow>1) %>% #this ensures for each combination of the above that there is 2 rows (1 for each subgroup level)
  select(covidence_id, type,subgroup_level,measure_n, test_used, test_subscale,timepoint3,mean,sd,n_with,n)

primary <- ind_sub %>%
  filter(subgroup_level == "primary")

secondary <- ind_sub %>%
  filter(subgroup_level == "secondary")

ind_dat <- full_join(primary,secondary,by=c("covidence_id","type","measure_n","test_used", "test_subscale","timepoint3")) %>%
  mutate(subgroup = 3)

#join all subgroups together for ease of ES calculation
subgroup_all <- rbind(shock_dat,sex_dat,ind_dat)

#### SUBGROUP: EFFECT SIZE ####

#ES calculation for type 1
subgroup_type1 <- subgroup_all %>%
  filter(type == 1)

sub_es1 <- escalc(measure = "SMD",
                 m1i = mean.x,
                 sd1i = sd.x,
                 n1i = n.x,
                 m2i = mean.y,
                 sd2i = sd.y,
                 n2i = n.y,
                 data = subgroup_type1)

subgroup_type2 <- subgroup_all %>%
  filter(type == 2)

sub_es2 <- escalc(measure = "OR",
                  ai = n_with.x,
                  bi = n.x - n_with.x,
                  n1i = n.x,
                  ci = n_with.y,
                  di = n.y - n_with.y,
                  n2i = n.y,
                  data = subgroup_type2)

sub_es <- rbind(sub_es1,sub_es2)

#add in cohen's d for sex
sex_d <- sub_data %>%
  filter(subgroup_level == "both")%>% 
  mutate(timepoint3 = timepoint) %>%
  filter(type == 1 | type == 2) %>% #cat and cont data
  mutate(yi = as.numeric((1 - 3/(4*(group1_cohensd_N+group2_cohensd_N-2) - 1)) * cohen_d_comparison)) %>%
  mutate(vi = as.numeric(1/group1_cohensd_N+ 1/group2_cohensd_N + yi^2/(2*(group1_cohensd_N+group2_cohensd_N)))) %>%
 mutate(n.x = group1_cohensd_N,
         n.y = group2_cohensd_N) %>% 
  mutate(subgroup_level.x = "female",
         subgroup_level.y = "male",
         mean.x = as.numeric(NA),
         mean.y = as.numeric(NA),
         sd.x = as.numeric(NA),
         sd.y = as.numeric(NA),
         n_with.x = as.numeric(NA),
         n_with.y = as.numeric(NA),
         subgroup = 2) %>%
  select(covidence_id, type,subgroup_level.x,measure_n, test_used, test_subscale,timepoint3,mean.x,sd.x,n_with.x,n.x,subgroup_level.y,mean.y,sd.y,n_with.y,n.y,subgroup,yi,vi)

#join cohen's d df to main sub_es df

sub_es <- rbind(sub_es,sex_d)


outliers <- sub_es %>%
  filter(yi>=2 || yi<=-2)
#no outliers

#### SUBGROUP: ANALYSIS DF ####

#average across measure/timepoint/subgroup/type so one ES per study per analysis
sub_es_mean <- sub_es %>%
  group_by(covidence_id,type,subgroup,measure_n,timepoint3) %>%
  mutate(yi = mean(yi,na.rm = T),
         vi = mean(vi,na.rm = T),
         n.x = min(n.x),
         n.y = min(n.y)) %>%
  filter(row_number()==1)%>%
  ungroup() %>%
  select(covidence_id,type,subgroup,timepoint3,measure_n,subgroup_level.x,n.x,subgroup_level.y,n.y,yi,vi)

#### META-ANALYSIS SET-UP ####

#create df for analysis output
output <- data.frame(matrix(ncol=24,nrow=0, dimnames=list(NULL, c("time", "measure", "analysis", "type", "comp_group","comp_time","subgroup", "N_comp", "k", "N", "b", "ci.lb", "ci.ub", "p","I2", "tau2", "se.tau2","QE","QEp",
                                                                  "eggers_int","eggers_p","trim_k0","trim_b","trim_se"))))

#assign variable types as needed
output$time <- as.numeric(output$time)
output$measure <- as.numeric(output$measure)
output$analysis <- as.character(output$analysis)
output$type <- as.numeric(output$type)
output$comp_group <- as.numeric(output$comp_group)
output$comp_time <- as.numeric(output$comp_time)
output$subgroup <- as.numeric(output$subgroup)
output$N_comp <- as.numeric(output$N_comp)
output$k <- as.numeric(output$k)
output$N <- as.numeric(output$N)
output$b <- as.numeric(output$b)
output$ci.lb <- as.numeric(output$ci.lb)
output$ci.ub <- as.numeric(output$ci.ub)
output$p <- as.numeric(output$p)
output$I2 <- as.numeric(output$I2)
output$tau2 <- as.numeric(output$tau2)
output$se.tau2 <- as.numeric(output$se.tau2)
output$QE <- as.numeric(output$QE)
output$QEp <- as.numeric(output$QEp)
output$eggers_int <- as.numeric(output$eggers_int)
output$eggers_p <- as.numeric(output$eggers_p)
output$trim_k0 <- as.numeric(output$trim_k0)
output$trim_b <- as.numeric(output$trim_b)
output$trim_se <- as.numeric(output$trim_se)

#### BETWEEN: META-ANALYSIS ####
#get df with studies included in the meta-analysis
b_studies <- b_analysis %>%
  group_by(measure_n,timepoint,type,comp_group) %>%
  mutate(nrow = n()) %>%
  filter(nrow > 2) %>%
  ungroup() %>%
  group_by(covidence_id) %>%
  filter(row_number()==1) %>%
  select(covidence_id)

#create vectors for each variable to loop over
time <- c(1,2,3,4,6,7)
measure <- c(1,2,3,4)
#type <- c(1,2,3,4,5,6)
type <- c(1,2)
group <- c(1,2,3,4)

#can't seem to add loop for group. There is no output when this 4th loop is added. Just run separately but keep adding to output df for now

##GROUP 1 (all comparison groups)
for(m in measure) {
  for(i in time) {
      for(t in type) {
        df <- b_analysis %>%
          filter(measure_n == m) %>%
          filter(timepoint == i) %>%
          filter(type == t) %>%
          filter(comp_group == 1)
        
        if(nrow(df) > 2) {
          ma <- rma(yi = g, vi = g_vi, data = df, method = "PM", test = "knha")
          
          output <- output %>%
            add_row(time=i,
                    measure=m,
                    analysis="between",
                    type=t,
                    comp_group=1,
                    comp_time=NA,
                    subgroup = NA,
                    N_comp=sum(df$N_total_comparison),
                    k=ma$k,
                    N=sum(df$N_total),
                    b=ma$b,
                    ci.lb=ma$ci.lb,
                    ci.ub=ma$ci.ub,
                    p=ma$pval,
                    I2=ma$I2,
                    tau2=ma$tau2,
                    se.tau2=ma$se.tau2,
                    QE=ma$QE,
                    QEp=ma$QEp,
                    eggers_int = NA,
                    eggers_p=NA,
                    trim_k0 = NA,
                    trim_b = NA,
                    trim_se = NA)
          }
      }
      }
}

##GROUP 2
for(m in measure) {
  for(i in time) {
    for(t in type) {
      df <- b_analysis %>%
        filter(measure_n == m) %>%
        filter(timepoint == i) %>%
        filter(type == t) %>%
        filter(comp_group == 2)
      
      if(nrow(df) > 2) {
        ma <- rma(yi = g, vi = g_vi, data = df, method = "PM", test = "knha")
        
        output <- output %>%
          add_row(time=i,
                  measure=m,
                  analysis="between",
                  type=t,
                  comp_group=2,
                  comp_time=NA,
                  N_comp=sum(df$N_total_comparison),
                  k=ma$k,
                  N=sum(df$N_total),
                  b=ma$b,
                  ci.lb=ma$ci.lb,
                  ci.ub=ma$ci.ub,
                  p=ma$pval,
                  I2=ma$I2,
                  tau2=ma$tau2,
                  se.tau2=ma$se.tau2,
                  QE=ma$QE,
                  QEp=ma$QEp,
                  eggers_int = NA,
                  eggers_p=NA,
                  trim_k0 = NA,
                  trim_b = NA,
                  trim_se = NA)
      }
    }
  }
}

##GROUP 3
for(m in measure) {
  for(i in time) {
    for(t in type) {
      df <- b_analysis %>%
        filter(measure_n == m) %>%
        filter(timepoint == i) %>%
        filter(type == t) %>%
        filter(comp_group == 3)
      
      if(nrow(df) > 2) {
        ma <- rma(yi = g, vi = g_vi, data = df, method = "PM", test = "knha")
        
        output <- output %>%
          add_row(time=i,
                  measure=m,
                  analysis="between",
                  type=t,
                  comp_group=3,
                  comp_time=NA,
                  N_comp=sum(df$N_total_comparison),
                  k=ma$k,
                  N=sum(df$N_total),
                  b=ma$b,
                  ci.lb=ma$ci.lb,
                  ci.ub=ma$ci.ub,
                  p=ma$pval,
                  I2=ma$I2,
                  tau2=ma$tau2,
                  se.tau2=ma$se.tau2,
                  QE=ma$QE,
                  QEp=ma$QEp,
                  eggers_int = NA,
                  eggers_p=NA,
                  trim_k0 = NA,
                  trim_b = NA,
                  trim_se = NA)
      }
    }
  }
}

##GROUP 4
for(m in measure) {
  for(i in time) {
    for(t in type) {
      df <- b_analysis %>%
        filter(measure_n == m) %>%
        filter(timepoint == i) %>%
        filter(type == t) %>%
        filter(comp_group == 4)
      
      if(nrow(df) > 2) {
        ma <- rma(yi = g, vi = g_vi, data = df, method = "PM", test = "knha")
        
        output <- output %>%
          add_row(time=i,
                  measure=m,
                  analysis="between",
                  type=t,
                  comp_group=4,
                  comp_time=NA,
                  N_comp=sum(df$N_total_comparison),
                  k=ma$k,
                  N=sum(df$N_total),
                  b=ma$b,
                  ci.lb=ma$ci.lb,
                  ci.ub=ma$ci.ub,
                  p=ma$pval,
                  I2=ma$I2,
                  tau2=ma$tau2,
                  se.tau2=ma$se.tau2,
                  QE=ma$QE,
                  QEp=ma$QEp,
                  eggers_int = NA,
                  eggers_p=NA,
                  trim_k0 = NA,
                  trim_b = NA,
                  trim_se = NA)
      }
    }
  }
}


#### PUBLICATION BIAS ####

id <- 1:nrow(output)

#create funnel plots
for(i in id) {
  if(output[i,"analysis"]=="between") {
    if(output[i,"type"]==1) {
    time <- output[i,"time"]
    measure <- as.numeric(output[i,"measure"])
    type <- output[i,"type"]
    group <- output[i,"comp_group"]
    funnel_name <- paste0("between_type",type,"_time",time,"_measure",measure,"_comp",group,".tiff")
    
    ti <- output[[i,"time"]]
    ty <- output[[i,"type"]]
    me<- output[[i,"measure"]]
    gr<- output[[i,"comp_group"]]
    
    df <- b_analysis %>%
      filter(type == ty) %>%
      filter(timepoint == ti) %>%
      filter(measure_n == me) %>%
      filter(comp_group == gr)
    
    ma <- rma(yi = g, vi = g_vi, data = df, method = "PM", test = "knha")
    
    tiff(filename = funnel_name, width = 8, height = 6, units = "in", res = 300, compression = "lzw")
    funnel(ma, trans = exp, xlab = "Hedges' g")
    dev.off()
    }
      if(output[i,"type"]==2) {
        time <- output[i,"time"]
        measure <- as.numeric(output[i,"measure"])
        type <- as.numeric(output[i,"type"])
        group <- output[i,"comp_group"]
        funnel_name <- paste0("between_type",type,"_time",time,"_measure",measure,"_comp",group,".tiff")
        
        ti <- output[[i,"time"]]
        ty <- output[[i,"type"]]
        me<- output[[i,"measure"]]
        gr<- output[[i,"comp_group"]]
        
        df <- b_analysis %>%
          filter(type == ty) %>%
          filter(timepoint == ti) %>%
          filter(measure_n == me) %>%
          filter(comp_group == gr)
        
        ma <- rma(yi = g, vi = g_vi, data = df, method = "PM", test = "knha")
        
        tiff(filename = funnel_name, width = 8, height = 6, units = "in", res = 300, compression = "lzw")
        funnel(ma, trans = exp, xlab = "Hedges' g")
        dev.off()
      }
    }
  }

#publication bias assessment (for analyses with k>9)
for(i in id) {
  if(output[i,"k"]>9){
    if(output[i,"analysis"]=="between") {
      
      ti <- output[[i,"time"]]
      ty <- output[[i,"type"]]
      me<- output[[i,"measure"]]
      gr<- output[[i,"comp_group"]]
      
      df <- b_analysis %>%
        filter(type == ty) %>%
        filter(timepoint == ti) %>%
        filter(measure_n == me) %>%
        filter(comp_group == gr)
    
    ma <- rma(yi = g, vi = g_vi, data = df, method = "PM", test = "knha")
    
    reg <- regtest(ma, model = "rma")
    
    output[i,"eggers_int"] <- reg$zval
    output[i,"eggers_p"]<- reg$pval
    
    if(output[i,"eggers_p"]<0.1){
      trim <- trimfill(ma)
      
      output[i,"trim_k0"] <- trim$k0
      output[i,"trim_b"] <- trim$b
      output[i,"trim_se"] <- trim$se
      
    }
  }
  }
}

#### PREVALENCE: META-ANALYSIS ####
#get df with studies included in the meta-analysis
p_studies <- p_analysis %>%
  group_by(measure_n,timepoint,type) %>%
  mutate(nrow = n()) %>%
  filter(nrow > 2) %>%
  ungroup() %>%
  group_by(covidence_id) %>%
  filter(row_number()==1) %>%
  select(covidence_id)

#create vectors for each variable to loop over
time <- c(1,2,3,4,5,6,7)
measure <- c(1,2,3,4)
type <- c(2,4,7)

for(m in measure) {
  for(i in time) {
    for(t in type) {
      df <- p_analysis %>%
        filter(measure_n == m) %>%
        filter(timepoint == i) %>%
        filter(type == t)
      
      if(nrow(df) > 2) {
        ma <- rma(yi = g, vi = g_vi, data = df, method = "PM", test = "knha")
        
        output <- output %>%
          add_row(time=i,
                  measure=m,
                  analysis="prevalence",
                  type=t,
                  comp_group=NA,
                  comp_time=NA,
                  subgroup = NA,
                  N_comp=NA,
                  k=ma$k,
                  N=sum(df$N_total),
                  b=ma$b*100,
                  ci.lb=ma$ci.lb*100,
                  ci.ub=ma$ci.ub*100,
                  p=ma$pval,
                  I2=ma$I2,
                  tau2=ma$tau2,
                  se.tau2=ma$se.tau2,
                  QE=ma$QE,
                  QEp=ma$QEp,
                  eggers_int = NA,
                  eggers_p=NA,
                  trim_k0 = NA,
                  trim_b = NA,
                  trim_se = NA)
      }
    }
  }
}

#### PUBLICATION BIAS ####

id <- 1:nrow(output)

#create funnel plots
for(i in id) {
  if(output[i,"analysis"]=="prevalence") {
    time <- output[i,"time"]
    measure <- as.numeric(output[i,"measure"])
    type <- output[i,"type"]
    funnel_name <- paste0("prevalence_type",type,"_measure",measure,"_time",time,".tiff")
    
    ti <- output[[i,"time"]]
    ty <- output[[i,"type"]]
    me<- output[[i,"measure"]]
    
    df <- p_analysis %>%
      filter(type == ty) %>%
      filter(timepoint == ti) %>%
      filter(measure_n == me)
    
    ma <- rma(yi = g, vi = g_vi, data = df, method = "PM", test = "knha")
    
    tiff(filename = funnel_name, width = 8, height = 6, units = "in", res = 300, compression = "lzw")
    funnel(ma, trans = exp, xlab = "%")
    dev.off()
  }
}

#publication bias assessment (for analyses with k>9)
for(i in id) {
  if(output[i,"k"]>9){
    if(output[i,"analysis"]=="prevalence") {
      
      ti <- output[[i,"time"]]
      ty <- output[[i,"type"]]
      me<- output[[i,"measure"]]
      
      df <- p_analysis %>%
        filter(type == ty) %>%
        filter(timepoint == ti) %>%
        filter(measure_n == me)
      
      ma <- rma(yi = g, vi = g_vi, data = df, method = "PM", test = "knha")
      
      reg <- regtest(ma, model = "rma")
      
      output[i,"eggers_int"] <- reg$zval
      output[i,"eggers_p"]<- reg$pval
      
      if(output[i,"eggers_p"]<0.1){
        trim <- trimfill(ma)
        
        output[i,"trim_k0"] <- trim$k0
        output[i,"trim_b"] <- trim$b
        output[i,"trim_se"] <- trim$se
        
      }
    }
  }
}

#### C PREVALENCE: META-ANALYSIS ####

#get df with studies included in the meta-analysis
cp_studies <- cp_analysis %>%
  group_by(measure_n,timepoint,type,comp_group) %>%
  mutate(nrow = n()) %>%
  filter(nrow > 2) %>%
  ungroup() %>%
  group_by(covidence_id) %>%
  filter(row_number()==1) %>%
  select(covidence_id)

#create vectors for each variable to loop over
time <- c(1,2,3,4,5,6,7)
measure <- c(1,2,3,4)
type <- c(2,4,7)

#group 1
for(m in measure) {
  for(i in time) {
    for(t in type) {
      df <- cp_analysis %>%
        filter(measure_n == m) %>%
        filter(timepoint == i) %>%
        filter(type == t) %>%
        filter(comp_group == 1)
      
      if(nrow(df) > 2) {
        ma <- rma(yi = g, vi = g_vi, data = df, method = "PM", test = "knha")
        
        output <- output %>%
          add_row(time=i,
                  measure=m,
                  analysis="comp_prevalence",
                  type=t,
                  comp_group=1,
                  comp_time=NA,
                  subgroup=NA,
                  N_comp=NA,
                  k=ma$k,
                  N=sum(df$N_total_comparison),
                  b=ma$b*100,
                  ci.lb=ma$ci.lb*100,
                  ci.ub=ma$ci.ub*100,
                  p=ma$pval,
                  I2=ma$I2,
                  tau2=ma$tau2,
                  se.tau2=ma$se.tau2,
                  QE=ma$QE,
                  QEp=ma$QEp,
                  eggers_int = NA,
                  eggers_p=NA,
                  trim_k0 = NA,
                  trim_b = NA,
                  trim_se = NA)
      }
    }
  }
}

#group 2
for(m in measure) {
  for(i in time) {
    for(t in type) {
      df <- cp_analysis %>%
        filter(measure_n == m) %>%
        filter(timepoint == i) %>%
        filter(type == t) %>%
        filter(comp_group == 2)
      
      if(nrow(df) > 2) {
        ma <- rma(yi = g, vi = g_vi, data = df, method = "PM", test = "knha")
        
        output <- output %>%
          add_row(time=i,
                  measure=m,
                  analysis="comp_prevalence",
                  type=t,
                  comp_group=2,
                  comp_time=NA,
                  subgroup = NA,
                  N_comp=NA,
                  k=ma$k,
                  N=sum(df$N_total_comparison),
                  b=ma$b*100,
                  ci.lb=ma$ci.lb*100,
                  ci.ub=ma$ci.ub*100,
                  p=ma$pval,
                  I2=ma$I2,
                  tau2=ma$tau2,
                  se.tau2=ma$se.tau2,
                  QE=ma$QE,
                  QEp=ma$QEp,
                  eggers_int = NA,
                  eggers_p=NA,
                  trim_k0 = NA,
                  trim_b = NA,
                  trim_se = NA)
      }
    }
  }
}

#group 3
for(m in measure) {
  for(i in time) {
    for(t in type) {
      df <- cp_analysis %>%
        filter(measure_n == m) %>%
        filter(timepoint == i) %>%
        filter(type == t) %>%
        filter(comp_group == 3)
      
      if(nrow(df) > 2) {
        ma <- rma(yi = g, vi = g_vi, data = df, method = "PM", test = "knha")
        
        output <- output %>%
          add_row(time=i,
                  measure=m,
                  analysis="comp_prevalence",
                  type=t,
                  comp_group=3,
                  comp_time=NA,
                  subgroup = NA,
                  N_comp=NA,
                  k=ma$k,
                  N=sum(df$N_total_comparison),
                  b=ma$b*100,
                  ci.lb=ma$ci.lb*100,
                  ci.ub=ma$ci.ub*100,
                  p=ma$pval,
                  I2=ma$I2,
                  tau2=ma$tau2,
                  se.tau2=ma$se.tau2,
                  QE=ma$QE,
                  QEp=ma$QEp,
                  eggers_int = NA,
                  eggers_p=NA,
                  trim_k0 = NA,
                  trim_b = NA,
                  trim_se = NA)
      }
    }
  }
}

#group 3
for(m in measure) {
  for(i in time) {
    for(t in type) {
      df <- cp_analysis %>%
        filter(measure_n == m) %>%
        filter(timepoint == i) %>%
        filter(type == t) %>%
        filter(comp_group == 4)
      
      if(nrow(df) > 2) {
        ma <- rma(yi = g, vi = g_vi, data = df, method = "PM", test = "knha")
        
        output <- output %>%
          add_row(time=i,
                  measure=m,
                  analysis="comp_prevalence",
                  type=t,
                  comp_group=4,
                  comp_time=NA,
                  subgroup = NA,
                  N_comp=NA,
                  k=ma$k,
                  N=sum(df$N_total_comparison),
                  b=ma$b*100,
                  ci.lb=ma$ci.lb*100,
                  ci.ub=ma$ci.ub*100,
                  p=ma$pval,
                  I2=ma$I2,
                  tau2=ma$tau2,
                  se.tau2=ma$se.tau2,
                  QE=ma$QE,
                  QEp=ma$QEp,
                  eggers_int = NA,
                  eggers_p=NA,
                  trim_k0 = NA,
                  trim_b = NA,
                  trim_se = NA)
      }
    }
  }
}

#### PUBLICATION BIAS ####

id <- 1:nrow(output)

#create funnel plots
for(i in id) {
  if(output[i,"analysis"]=="comp_prevalence") {
    time <- output[i,"time"]
    measure <- as.numeric(output[i,"measure"])
    type <- output[i,"type"]
    funnel_name <- paste0("comp_prevalence_type",type,"_measure",measure,"_time",time,".tiff")
    
    ti <- output[[i,"time"]]
    ty <- output[[i,"type"]]
    me<- output[[i,"measure"]]
    
    df <- cp_analysis %>%
      filter(type == ty) %>%
      filter(timepoint == ti) %>%
      filter(measure_n == me)
    
    ma <- rma(yi = g, vi = g_vi, data = df, method = "PM", test = "knha")
    
    tiff(filename = funnel_name, width = 8, height = 6, units = "in", res = 300, compression = "lzw")
    funnel(ma, trans = exp, xlab = "%")
    dev.off()
  }
}

#publication bias assessment (for analyses with k>9)
for(i in id) {
  if(output[i,"k"]>9){
    if(output[i,"analysis"]=="comp_prevalence") {
      
      ti <- output[[i,"time"]]
      ty <- output[[i,"type"]]
      me<- output[[i,"measure"]]
      
      df <- cp_analysis %>%
        filter(type == ty) %>%
        filter(timepoint == ti) %>%
        filter(measure_n == me)
      
      ma <- rma(yi = g, vi = g_vi, data = df, method = "PM", test = "knha")
      
      reg <- regtest(ma, model = "rma")
      
      output[i,"eggers_int"] <- reg$zval
      output[i,"eggers_p"]<- reg$pval
      
      if(output[i,"eggers_p"]<0.1){
        trim <- trimfill(ma)
        
        output[i,"trim_k0"] <- trim$k0
        output[i,"trim_b"] <- trim$b
        output[i,"trim_se"] <- trim$se
        
      }
    }
  }
}

#### WITHIN: META-ANALYSIS ####
#get df with studies included in the meta-analysis
w_studies <- w_analysis %>%
  filter(type == 1) %>%
  group_by(type,measure_n,within) %>%
  mutate(nrow = n()) %>%
  filter(nrow > 2) %>%
  ungroup() %>%
  group_by(covidence_id) %>%
  filter(row_number()==1) %>%
  select(covidence_id)

#get df with studies that could potentially be in data analysis
w_total_studies <- w_analysis %>%
  group_by(covidence_id) %>%
  filter(row_number()==1)

#add loop in for "within" variable - this is the timepoint comparison
#create vectors for each variable to loop over
within <- c(1,2,3,4)
measure <- c(1,2,3,4)
type <- c(1,2)
#group <- c(1,2,3)

#analysis
for(t in type) {
for(m in measure) {
  for(w in within) {
      df <- w_analysis %>%
        filter(type == t) %>%
        filter(measure_n == m) %>%
        filter(within == w)
      
      if(nrow(df) > 2) {
        ma <- rma(yi = g, vi = g_vi, data = df, method = "PM", test = "knha")
        
        output <- output %>%
          add_row(time=NA,
                  measure=m,
                  analysis="within",
                  type=t, 
                  comp_group=NA,
                  comp_time=w,
                  subgroup=NA,
                  N_comp=sum(df$n.y), #N at timepoint 2
                  k=ma$k,
                  N=sum(df$n.x), #N at timepoint 1
                  b=ma$b,
                  ci.lb=ma$ci.lb,
                  ci.ub=ma$ci.ub,
                  p=ma$pval,
                  I2=ma$I2,
                  tau2=ma$tau2,
                  se.tau2=ma$se.tau2,
                  QE=ma$QE,
                  QEp=ma$QEp,
                  eggers_int = NA,
                  eggers_p=NA,
                  trim_k0 = NA,
                  trim_b = NA,
                  trim_se = NA)
      }
  }
}
}

w_clin <- output %>%
  filter(type == 2) %>%
  filter(analysis == "within")

sig <- output %>%
  filter(p<0.05)

#### PUBLICATION BIAS ####

id <- 1:nrow(output)

#create funnel plots
for(i in id) {
  if(output[i,"analysis"]=="within") {
    time <- output[i,"comp_time"]
    measure <- as.numeric(output[i,"measure"])
    type <- output[i,"type"]
    funnel_name <- paste0("within_type",type,"_measure",measure,"_comptime",time,".tiff")
    
    ti <- output[[i,"comp_time"]]
    ty <- output[[i,"type"]]
    me<- output[[i,"measure"]]
    
    df <- w_analysis %>%
      filter(type == ty) %>%
      filter(within == ti) %>%
      filter(measure_n == me)
    
    ma <- rma(yi = g, vi = g_vi, data = df, method = "PM", test = "knha")
    
    tiff(filename = funnel_name, width = 8, height = 6, units = "in", res = 300, compression = "lzw")
    funnel(ma, trans = exp, xlab = "SMCR")
    dev.off()
  }
}

#publication bias assessment (for analyses with k>9)
for(i in id) {
  if(output[i,"k"]>9){
    if(output[i,"analysis"]=="within") {
      
      ti <- output[[i,"comp_time"]]
      ty <- output[[i,"type"]]
      me<- output[[i,"measure"]]
      
      df <- w_analysis %>%
        filter(type == ty) %>%
        filter(within == ti) %>%
        filter(measure_n == me)
      
      ma <- rma(yi = g, vi = g_vi, data = df, method = "PM", test = "knha")
      
      reg <- regtest(ma, model = "rma")
      
      output[i,"eggers_int"] <- reg$zval
      output[i,"eggers_p"]<- reg$pval
      
      if(output[i,"eggers_p"]<0.1){
        trim <- trimfill(ma)
        
        output[i,"trim_k0"] <- trim$k0
        output[i,"trim_b"] <- trim$b
        output[i,"trim_se"] <- trim$se
        
      }
    }
  }
}

#### SUBGROUP: META-ANALYSIS ####

s_studies <- sub_es_mean %>%
  group_by(measure_n,timepoint3,type,subgroup) %>%
  mutate(nrow = n()) %>%
  filter(nrow > 2) %>%
  ungroup() %>%
  group_by(covidence_id) %>%
  filter(row_number()==1) %>%
  select(covidence_id)

s_studies_sub <- sub_es_mean %>%
  group_by(measure_n,timepoint3,type,subgroup) %>%
  mutate(nrow = n()) %>%
  filter(nrow > 2) %>%
  ungroup() %>%
  group_by(covidence_id,subgroup) %>%
  filter(row_number()==1) %>%
  select(covidence_id)

#create vectors for each variable to loop over
time <- c(1,2,3,4,6,7)
measure <- c(1,2,3,4)
type <- c(1,2)
sub <- c(1,2,3)

##TYPE 1 (all comparison groups)
for(m in measure) {
  for(i in time) {
    for(s in sub) {
      df <- sub_es_mean %>%
        filter(measure_n == m) %>%
        filter(timepoint3 == i) %>%
        filter(subgroup == s) %>%
        filter(type == 1)
      
      if(nrow(df) > 2) {
        ma <- rma(yi = yi, vi = vi, data = df, method = "PM", test = "knha")
        
        output <- output %>%
          add_row(time=i,
                  measure=m,
                  analysis="subgroup",
                  type=1,
                  comp_group=NA,
                  comp_time=NA,
                  subgroup = s,
                  N_comp=sum(df$n.x),
                  k=ma$k,
                  N=sum(df$n.y),
                  b=ma$b,
                  ci.lb=ma$ci.lb,
                  ci.ub=ma$ci.ub,
                  p=ma$pval,
                  I2=ma$I2,
                  tau2=ma$tau2,
                  se.tau2=ma$se.tau2,
                  QE=ma$QE,
                  QEp=ma$QEp,
                  eggers_int = NA,
                  eggers_p=NA,
                  trim_k0 = NA,
                  trim_b = NA,
                  trim_se = NA)
      }
    }
  }
}

##TYPE 1 (all comparison groups)
for(m in measure) {
  for(i in time) {
    for(s in sub) {
      df <- sub_es_mean %>%
        filter(measure_n == m) %>%
        filter(timepoint3 == i) %>%
        filter(subgroup == s) %>%
        filter(type == 2)
      
      if(nrow(df) > 2) {
        ma <- rma(yi = yi, vi = vi, data = df, method = "PM", test = "knha")
        
        output <- output %>%
          add_row(time=i,
                  measure=m,
                  analysis="subgroup",
                  type=2,
                  comp_group=NA,
                  comp_time=NA,
                  subgroup = s,
                  N_comp=sum(df$n.x),
                  k=ma$k,
                  N=sum(df$n.y),
                  b=ma$b,
                  ci.lb=ma$ci.lb,
                  ci.ub=ma$ci.ub,
                  p=ma$pval,
                  I2=ma$I2,
                  tau2=ma$tau2,
                  se.tau2=ma$se.tau2,
                  QE=ma$QE,
                  QEp=ma$QEp,
                  eggers_int = NA,
                  eggers_p=NA,
                  trim_k0 = NA,
                  trim_b = NA,
                  trim_se = NA)
      }
    }
  }
}

#### PUBLICATION BIAS ####

id <- 1:nrow(output)

#create funnel plots
for(i in id) {
  if(output[i,"analysis"]=="subgroup") {
    if(output[i,"type"]==1) {
      time <- output[i,"time"]
      measure <- as.numeric(output[i,"measure"])
      type <- output[i,"type"]
      group <- output[i,"subgroup"]
      funnel_name <- paste0("subgroup_type",type,"_time",time,"_measure",measure,"_sub",group,".tiff")
      
      ti <- output[[i,"time"]]
      ty <- output[[i,"type"]]
      me<- output[[i,"measure"]]
      gr<- output[[i,"subgroup"]]
      
      df <- sub_es_mean %>%
        filter(type == ty) %>%
        filter(timepoint3 == ti) %>%
        filter(measure_n == me) %>%
        filter(subgroup == gr)
      
      ma <- rma(yi = yi, vi = vi, data = df, method = "PM", test = "knha")
      
      tiff(filename = funnel_name, width = 8, height = 6, units = "in", res = 300, compression = "lzw")
      funnel(ma, trans = exp, xlab = "Hedges' g")
      dev.off()
    }
    if(output[i,"type"]==2) {
      time <- output[i,"time"]
      measure <- as.numeric(output[i,"measure"])
      type <- output[i,"type"]
      group <- output[i,"subgroup"]
      funnel_name <- paste0("subgroup_type",type,"_time",time,"_measure",measure,"_sub",group,".tiff")
      
      ti <- output[[i,"time"]]
      ty <- output[[i,"type"]]
      me<- output[[i,"measure"]]
      gr<- output[[i,"subgroup"]]
      
      df <- sub_es_mean %>%
        filter(type == ty) %>%
        filter(timepoint3 == ti) %>%
        filter(measure_n == me) %>%
        filter(subgroup == gr)
      
      ma <- rma(yi = yi, vi = vi, data = df, method = "PM", test = "knha")
      
      tiff(filename = funnel_name, width = 8, height = 6, units = "in", res = 300, compression = "lzw")
      funnel(ma, trans = exp, xlab = "Hedges' g")
      dev.off()
    }
  }
}

#publication bias assessment (for analyses with k>9)
for(i in id) {
  if(output[i,"k"]>9){
    if(output[i,"analysis"]=="subgroup") {
      
      ti <- output[[i,"time"]]
      ty <- output[[i,"type"]]
      me<- output[[i,"measure"]]
      gr<- output[[i,"subgroup"]]
      
      df <- sub_es_mean %>%
        filter(type == ty) %>%
        filter(timepoint3 == ti) %>%
        filter(measure_n == me) %>%
        filter(subgroup == gr)
      
      ma <- rma(yi = yi, vi = vi, data = df, method = "PM", test = "knha")
      
      reg <- regtest(ma, model = "rma")
      
      output[i,"eggers_int"] <- reg$zval
      output[i,"eggers_p"]<- reg$pval
      
      if(output[i,"eggers_p"]<0.1){
        trim <- trimfill(ma)
        
        output[i,"trim_k0"] <- trim$k0
        output[i,"trim_b"] <- trim$b
        output[i,"trim_se"] <- trim$se
        
      }
    }
  }
}


#### CLEAN OUTPUT FILE ####

#remove unspecified post-ICD timepoint
output <- output %>%
  filter(time != 5 | is.na(time))


#add some labels to output doc
type <- c(1,2,3,4,5,6,7)
type_lab <- c("Continuous","Clinically significant/diagnosis","Continuous/clinically significant","At least mild symptoms","Continuous/mild+","All data","Moderate+")
type_labs <- data.frame(type,type_lab)

output <- left_join(output,type_labs,by="type")

measure <- c(1:4)
meas_lab <- c("Anxiety","Depression","PTSD","Any mood disorder")
meas_labs <- data.frame(measure,meas_lab)

output <- left_join(output,meas_labs,by="measure")

#change time numbers to match with the following definition
output <- output %>%
  mutate(time = ifelse(time == 1,4,
                       ifelse(time == 2,3,
                              ifelse(time == 3,2,
                                     ifelse(time == 4,1,
                                            ifelse(time == 5,5,
                                                   ifelse(time == 6,7,
                                                          ifelse(time == 7,6,NA))))))))

time <- c(1:7)
time_labs <- c("12 months+","6 months - 12 months","Discharge - 5 months", "Pre-discharge","Unspecified post-ICD", "Post-ICD", "Pre-ICD")
time_labs <- data.frame(time, time_labs)

output <- left_join(output,time_labs,by="time")

subgroup <- c(1:3)
sub_labs <- c("Shock vs no shock","Female vs male","Primary vs secondary indication")
sub_labs <- data.frame(subgroup,sub_labs)

output <- left_join(output,sub_labs,by="subgroup")

#add comp group labels
comp_group <- c(1:4)
comp_lab <- c("Any comparison group","Partners","Cardiac","General")
comp_labs <- data.frame(comp_group,comp_lab)

output <- left_join(output,comp_labs,by="comp_group")

comp_time <- c(1:4)
comp_time_labs <- c("Pre-discharge vs discharge-5months","Discharge-5months vs 6-12months","6-12months vs 12months+","Pre-ICD vs post-ICD")
comp_time_labs <- data.frame(comp_time,comp_time_labs)

output <- left_join(output,comp_time_labs,by="comp_time")


sig <- output%>%
  filter(p<.05)

#### PREVALENCE OUTPUT ####

##prevalence graph
prevalence <- output %>%
  filter(analysis == "prevalence")%>%
  arrange(desc(time))

prevalence$time_labs <- factor(prevalence$time_labs, ordered = T)

prevalence$time_labs <- fct_reorder(prevalence$time_labs, prevalence$time)

#add facet for mood disorder variables
prevalence$meas_lab <- fct_reorder(prevalence$meas_lab,prevalence$measure)

tiff(filename = "prevalence.tiff",width=10,height=11,units = "in",compression="lzw",res=600)
prevalence %>%
  ggplot(mapping = aes(x = time_labs, y = b, ymin = ci.lb, ymax = ci.ub))+
  geom_pointrange(aes(col = type_lab), size = 1, fatten = 2, position = position_dodge2(width = 0.5))+
  scale_color_manual(values = c(rgb(135,135,135,maxColorValue = 255),
                                rgb(174,16,34,maxColorValue = 255),
                                rgb(30,169,145,maxColorValue = 255))) +
  geom_vline(mapping = aes(xintercept = 4.5)) +
  scale_x_discrete(name = "Timepoint") +
  ylab("%")+
  coord_flip(ylim=c(0,50)) +
  theme_clean()+
  ggforce::facet_col(~ meas_lab, scales = "free", space = "free") +
  labs(colour="Mood cut-off")+
  theme(strip.text.x = element_text(face="bold"),
        plot.background = element_blank())
dev.off()

#prev graph for graphical abstract
g_abstract <- prevalence %>%
  filter(type_lab == "Clinically significant/diagnosis") %>%
  filter(time_labs == "Post-ICD") %>%
  filter(meas_lab != "Any mood disorder")

png(filename = "prev_gabstract.png",width=7,height=4,units = "in",res=600,bg = "transparent")
g_abstract %>%
  ggplot(mapping = aes(x = measure, y = b, ymin = ci.lb, ymax = ci.ub))+
  geom_pointrange(size = 2, fatten = 2, position = position_dodge2(width = 0.5), color = rgb(0,86,148,maxColorValue = 255))+
  scale_x_discrete(name = "Mood", limits = c("Anxiety","Depression","PTSD"))+
  ylab("Prevalence (%)")+
  coord_flip(ylim= c(0,30)) +
  theme_clean()+
  theme(strip.text.x = element_text(face="bold"),
        plot.background = element_blank())
dev.off()

##create prevalence output table
prev_output <- prevalence %>%
  mutate(b = round(b, digits =2)) %>%
  mutate(ci.lb = round(ci.lb, digits =2)) %>%
  mutate(ci.ub = round(ci.ub, digits =2)) %>%
  mutate(p = round(p, digits =3)) %>%
  mutate(I2 = round(I2, digits =2)) %>%
  mutate(tau2 = round(tau2, digits =2)) %>%
  mutate(se.tau2 = round(se.tau2, digits =2)) %>%
  mutate(QE = round(QE, digits =2)) %>%
  mutate(QEp = round(QEp, digits =2)) %>%
  select(meas_lab,time_labs,type_lab,k:trim_se)

write.csv(prev_output, file = "prev_output.csv")

#### C PREVALENCE OUTPUT ####

##comp prevalence graph

#1 = all, 2 = partners, 3 = cardiac, 4 = general

comp_prevalence_partner <- output %>%
  filter(analysis == "comp_prevalence") %>%
  filter(comp_group == 2) %>%
  arrange(desc(time))

comp_prevalence_cardiac <- output %>%
  filter(analysis == "comp_prevalence") %>%
  filter(comp_group == 3) %>%
  arrange(desc(time))

comp_prevalence_general <- output %>%
  filter(analysis == "comp_prevalence") %>%
  filter(comp_group == 4) %>%
  arrange(desc(time))

comp_prevalence_partner$time_labs <- factor(comp_prevalence_partner$time_labs, ordered = T)

comp_prevalence_partner$time_labs <- fct_reorder(comp_prevalence_partner$time_labs, comp_prevalence_partner$time)

#add facet for mood disorder variables
comp_prevalence_partner$meas_lab <- fct_reorder(comp_prevalence_partner$meas_lab,comp_prevalence_partner$measure)

##create comp prevalence output table
compprev_partner_output <- comp_prevalence_partner %>%
  mutate(b = round(b, digits =2)) %>%
  mutate(ci.lb = round(ci.lb, digits =2)) %>%
  mutate(ci.ub = round(ci.ub, digits =2)) %>%
  mutate(p = round(p, digits =3)) %>%
  mutate(I2 = round(I2, digits =2)) %>%
  mutate(tau2 = round(tau2, digits =2)) %>%
  mutate(se.tau2 = round(se.tau2, digits =2)) %>%
  mutate(QE = round(QE, digits =2)) %>%
  mutate(QEp = round(QEp, digits =2)) %>%
  select(meas_lab,time_labs,type_lab,k:trim_se)

write.csv(compprev_partner_output, file = "compprev_partner_output.csv")

#### BETWEEN OUTPUT ####

#between groups graph
between_cont <- output %>%
  filter(analysis == "between") %>%
  filter(type == 1) %>%
  arrange(desc(time))

text_low <- textGrob("Less mood symptoms in ICD group",gp=gpar(fontsize=7),just = "right")
text_high <- textGrob("More mood symptoms in ICD group",gp=gpar(fontsize=7),just="left")

between_cont$time_labs <- fct_inorder(between_cont$time_labs)

#add facet for mood disorder variables
between_cont$meas_lab <- fct_reorder(between_cont$meas_lab,between_cont$measure)

tiff(filename = "between_continuous.tiff",width=10,height=8,units = "in",compression="lzw",res=600)
between_cont %>%
  ggplot(mapping = aes(x = time_labs, y = b, ymin = ci.lb, ymax = ci.ub))+
  geom_hline(aes(yintercept = 0), linetype="dashed", color="#666666") +
  geom_pointrange(aes(col = comp_lab), size = 1, fatten = 2, position = position_dodge2(width = 0.5))+
  scale_color_manual(values = c(rgb(135,135,135,maxColorValue = 255),
                                rgb(174,16,34,maxColorValue = 255),
                                rgb(0,86,148,maxColorValue = 255))) +
  geom_vline(aes(xintercept = 4.5)) +
  scale_x_discrete(name = "Timepoint") +
  ylab("Hedges' g")+
  coord_flip(ylim = c(-0.4,0.4)) +
  theme_clean()+
  annotation_custom(text_low,xmin=0.6,xmax=0.6,ymin=-0.05,ymax=-0.05)+
  annotation_custom(text_high,xmin=0.6,xmax=0.6,ymin=0.05,ymax=0.05)+
  ggforce::facet_col(~ meas_lab, scales = "free", space = "free") +
  labs(colour="Comparison group")+
  theme(strip.text.x = element_text(face="bold"),
        plot.background = element_blank())
dev.off()

between_clin <- output %>%
  filter(analysis == "between") %>%
  filter(type == 2) %>%
  arrange(desc(time))

#CHANGE TO OR (NOT LOG ODDS)
between_clin <- between_clin %>%
  mutate(exp_or = exp(b),
         exp_cilb = exp(ci.lb),
         exp_ciub = exp(ci.ub))

text_low <- textGrob("Lower odds of mood disorder in ICD",gp=gpar(fontsize=7),just = "right")
text_high <- textGrob("Higher odds of mood disorder in ICD",gp=gpar(fontsize=7),just="left")

between_clin$time_labs <- fct_inorder(between_clin$time_labs)

#add facet for mood disorder variables
between_clin$meas_lab <- fct_reorder(between_clin$meas_lab,between_clin$measure)

tiff(filename = "between_clinical.tiff",width=10,height=7,units = "in",compression="lzw",res=600)
between_clin %>%
  ggplot(mapping = aes(x = time_labs, y = exp_or, ymin = exp_cilb, ymax = exp_ciub))+
  geom_hline(aes(yintercept = 1), linetype="dashed", color="#666666") +
  geom_pointrange(aes(col = comp_lab), size = 1, fatten = 2, position = position_dodge2(width = 0.5))+
  scale_color_manual(values = c(rgb(135,135,135,maxColorValue = 255),
                                rgb(0,86,148,maxColorValue = 255))) +
  geom_vline(aes(xintercept = 4.5)) +
  scale_x_discrete(name = "Timepoint")+
  ylab("OR")+
  coord_flip(ylim = c(0,4)) +
  theme_clean()+
  annotation_custom(text_high,xmin=0.6,xmax=0.6,ymin=1.05,ymax=1.05)+
  ggforce::facet_col(~ meas_lab, scales = "free", space = "free") +
  labs(colour="Comparison group")+
  theme(strip.text.x = element_text(face="bold"),
        plot.background = element_blank())
dev.off()

#set up between groups output file
between_clin_output <- between_clin %>%
  mutate(b = exp(b)) %>% #this has already been changed above?
  mutate(ci.lb = exp(ci.lb)) %>%
  mutate(ci.ub = exp(ci.ub)) %>%
  mutate(b = round(b, digits =2)) %>%
  mutate(ci.lb = round(ci.lb, digits =2)) %>%
  mutate(ci.ub = round(ci.ub, digits =2)) %>%
  mutate(p = round(p, digits =3)) %>%
  mutate(I2 = round(I2, digits =2)) %>%
  mutate(tau2 = round(tau2, digits =2)) %>%
  mutate(se.tau2 = round(se.tau2, digits =2)) %>%
  mutate(QE = round(QE, digits =2)) %>%
  mutate(QEp = round(QEp, digits =2)) %>%
  select(meas_lab,time_labs,comp_lab,type_lab,k:trim_se) %>%
  arrange(time_labs,meas_lab,comp_lab)

write.csv(between_clin_output, file = "between_clin_output.csv")

#cont output file
between_cont_output <- between_cont%>%
  mutate(b = round(b, digits =2)) %>%
  mutate(ci.lb = round(ci.lb, digits =2)) %>%
  mutate(ci.ub = round(ci.ub, digits =2)) %>%
  mutate(p = round(p, digits =3)) %>%
  mutate(I2 = round(I2, digits =2)) %>%
  mutate(tau2 = round(tau2, digits =2)) %>%
  mutate(se.tau2 = round(se.tau2, digits =2)) %>%
  mutate(QE = round(QE, digits =2)) %>%
  mutate(QEp = round(QEp, digits =2)) %>%
  select(meas_lab,time_labs,comp_lab,type_lab,k:trim_se) %>%
  arrange(time_labs,meas_lab,comp_lab)

write.csv(between_cont_output, file = "between_cont_output.csv")

#### WITHIN OUTPUT ####

#within groups graph
within_cont <- output %>%
  filter(analysis == "within") %>%
  filter(type == 1)

#cont output file
#no need to add time for this analysis

within_cont_output <- within_cont %>%
  mutate(b = round(b, digits =2)) %>%
  mutate(ci.lb = round(ci.lb, digits =2)) %>%
  mutate(ci.ub = round(ci.ub, digits =2)) %>%
  mutate(p = round(p, digits =3)) %>%
  mutate(I2 = round(I2, digits =2)) %>%
  mutate(tau2 = round(tau2, digits =2)) %>%
  mutate(se.tau2 = round(se.tau2, digits =2)) %>%
  mutate(QE = round(QE, digits =2)) %>%
  mutate(QEp = round(QEp, digits =2)) %>%
  select(meas_lab,comp_time_labs,type_lab,k:trim_se) %>%
  arrange(meas_lab,comp_time_labs)

write.csv(within_cont_output, file = "within_cont_output.csv")

#set up annotations
text_low <- textGrob("Less mood symptoms at T1",gp=gpar(fontsize=7),just = "right")
text_high <- textGrob("More mood symptoms at T1",gp=gpar(fontsize=7),just = "left")

#add facet for mood disorder variables
within_cont$meas_lab <- fct_reorder(within_cont$meas_lab,within_cont$measure)

tiff(filename = "within_continuous.tiff",width=8,height=8,units = "in",compression="lzw",res=600)
within_cont %>%
  ggplot(mapping = aes(x = comp_time, y = b, ymin = ci.lb, ymax = ci.ub))+
  geom_hline(aes(yintercept = 0), linetype="dashed", color="#666666") +
  geom_pointrange(size = 1, fatten = 2, position = position_dodge2(width = 0.5))+
  scale_x_discrete(name = "Timepoint comparison", limits = c("Pre-discharge vs discharge-5months","Discharge-5months vs 6-12months","6-12months vs 12months+","Pre-ICD vs post-ICD"))+
  ylab("SMCR")+
  coord_flip(ylim= c(-0.4,0.3)) +
  theme_clean()+
  annotation_custom(text_low,xmin=0.6,xmax=0.6,ymin=-0.05,ymax=-0.05)+
  annotation_custom(text_high,xmin=0.6,xmax=0.6,ymin=0.05,ymax=0.05)+
  ggforce::facet_col(~ meas_lab, scales = "free", space = "free")+
  theme(strip.text.x = element_text(face="bold"),
        plot.background = element_blank())
dev.off()

#### SUBGROUP OUTPUT ####

#between groups graph
subgroup_cont <- output %>%
  filter(analysis == "subgroup") %>%
  filter(type == 1) %>%
  arrange(desc(time))

text_low <- textGrob("Less mood symptoms in subgroup 1",gp=gpar(fontsize=7),just = "right")
text_high <- textGrob("More mood symptoms in subgroup 1",gp=gpar(fontsize=7),just="left")

subgroup_cont$time_labs <- fct_inorder(subgroup_cont$time_labs)

#add facet for mood disorder variables
subgroup_cont$meas_lab <- fct_reorder(subgroup_cont$meas_lab,subgroup_cont$measure)

tiff(filename = "subgroup_continuous.tiff",width=10,height=8,units = "in",compression="lzw",res=600)
subgroup_cont %>%
  ggplot(mapping = aes(x = time_labs, y = b, ymin = ci.lb, ymax = ci.ub))+
  geom_hline(aes(yintercept = 0), linetype="dashed", color="#666666") +
  geom_pointrange(aes(
    col = sub_labs
    ), size = 1, fatten = 2, position = position_dodge2(width = 0.5))+
  scale_color_manual(values = c(rgb(135,135,135,maxColorValue = 255),
                                rgb(174,16,34,maxColorValue = 255),
                                rgb(0,86,148,maxColorValue = 255))) +
  geom_vline(aes(xintercept = 4.5)) +
  scale_x_discrete(name = "Timepoint") +
  ylab("Hedges' g")+
  coord_flip(ylim = c(-0.8,0.8)) +
  theme_clean()+
  ggforce::facet_col(~ meas_lab, scales = "free", space = "free") +
  labs(colour="Comparison group")+
  theme(strip.text.x = element_text(face="bold"),
        plot.background = element_blank())
dev.off()

#cont output file
subgroup_cont_output <- subgroup_cont %>%
  mutate(b = round(b, digits =2)) %>%
  mutate(ci.lb = round(ci.lb, digits =2)) %>%
  mutate(ci.ub = round(ci.ub, digits =2)) %>%
  mutate(p = round(p, digits =3)) %>%
  mutate(I2 = round(I2, digits =2)) %>%
  mutate(tau2 = round(tau2, digits =2)) %>%
  mutate(se.tau2 = round(se.tau2, digits =2)) %>%
  mutate(QE = round(QE, digits =2)) %>%
  mutate(QEp = round(QEp, digits =2)) %>%
  select(meas_lab,time_labs,sub_labs,type_lab,k:trim_se) %>%
  arrange(sub_labs,time_labs,meas_lab)

write.csv(subgroup_cont_output, file = "subgroup_cont_output.csv")

subgroup_clin <- output %>%
  filter(analysis == "subgroup") %>%
  filter(type == 2) %>%
  arrange(desc(time))

subgroup_clin <- left_join(subgroup_clin,comp_labs,by="comp_group")

#CHANGE TO OR (NOT LOG ODDS)
subgroup_clin <- subgroup_clin %>%
  mutate(exp_or = exp(b),
         exp_cilb = exp(ci.lb),
         exp_ciub = exp(ci.ub))


text_low <- textGrob("Lower odds of mood disorder in Subgroup 1",gp=gpar(fontsize=7),just = "right")
text_high <- textGrob("Higher odds of mood disorder in Subgroup 1",gp=gpar(fontsize=7),just="left")

subgroup_clin$time_labs <- fct_inorder(subgroup_clin$time_labs)

#add facet for mood disorder variables
subgroup_clin$meas_lab <- fct_reorder(subgroup_clin$meas_lab,subgroup_clin$measure)

tiff(filename = "subgroup_clinical.tiff",width=10,height=7,units = "in",compression="lzw",res=600)
subgroup_clin %>%
  ggplot(mapping = aes(x = time_labs, y = exp_or, ymin = exp_cilb, ymax = exp_ciub))+
  geom_hline(aes(yintercept = 1), linetype="dashed", color="#666666") +
  geom_pointrange(aes(col = sub_labs), size = 1, fatten = 2, position = position_dodge2(width = 0.5))+
  scale_color_manual(values = c(rgb(135,135,135,maxColorValue = 255),
                                #rgb(174,16,34,maxColorValue = 255),
                                rgb(0,86,148,maxColorValue = 255))) +
  geom_vline(aes(xintercept = 4.5)) +
  scale_x_discrete(name = "Timepoint")+
  ylab("OR")+
  coord_flip(ylim = c(0,4)) +
  theme_clean()+
  ggforce::facet_col(~ meas_lab, scales = "free", space = "free") +
  labs(colour="Comparison group")+
  theme(strip.text.x = element_text(face="bold"),
        plot.background = element_blank())
dev.off()

#set up between groups output file
subgroup_clin_output <- subgroup_clin %>%
  mutate(b = exp(b)) %>%
  mutate(ci.lb = exp(ci.lb)) %>%
  mutate(ci.ub = exp(ci.ub)) %>%
  mutate(b = round(b, digits =2)) %>%
  mutate(ci.lb = round(ci.lb, digits =2)) %>%
  mutate(ci.ub = round(ci.ub, digits =2)) %>%
  mutate(p = round(p, digits =3)) %>%
  mutate(I2 = round(I2, digits =2)) %>%
  mutate(tau2 = round(tau2, digits =2)) %>%
  mutate(se.tau2 = round(se.tau2, digits =2)) %>%
  mutate(QE = round(QE, digits =2)) %>%
  mutate(QEp = round(QEp, digits =2)) %>%
  select(meas_lab,time_labs,sub_labs,type_lab,k:trim_se) %>%
  arrange(sub_labs,time_labs,meas_lab)

write.csv(subgroup_clin_output, file = "subgroup_clin_output.csv")