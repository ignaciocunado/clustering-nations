library(haven)
library(dplyr)

df <- read.csv("pisa2018_dataset.csv")
data=read.csv("OECD_data.csv")
data = data %>% group_by(CNT) %>% filter(n()>3) %>% ungroup()
data = as.data.frame(data)
data$CNT <- as.numeric(as.factor(data$CNT))

scuole <- unique(df$CNTSCHID)
df2 <- data.frame(CNTSCHID = scuole)
df2 <- df %>%
  group_by(CNTSCHID) %>%
  summarise(PV1MATH = mean(PV1MATH, na.rm = TRUE))  # na.rm = TRUE rimuove i NA

data <- data %>%
  left_join(df2, by = "CNTSCHID")

dataOriginal <- read.csv("OECD_data.csv")

