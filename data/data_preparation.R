## 1. DATA IMPORT 

#__________________
## Libraries import
library(foreign)
library(dplyr)
library(ggplot2)
library(naniar)
library(finalfit)
library(MissMech)


#_____________________________________________
## 1.A) Schools questionnaire data file import

#_________________
## Raw Data import
## downloaded from  https://webfs.oecd.org/pisa2018/SPSS_SCH_QQQ.zip
SCH = read.spss('data/schools/raw_data_schools.sav', reencode='utf-8')

sch0 = data.frame('CNT' = SCH[["CNT"]],
                 'CNTSCHID' = SCH[["CNTSCHID"]],
                 'PRIVATESCH' = SCH[["PRIVATESCH"]],
                 'STRATIO' = as.numeric(as.character(SCH[["STRATIO"]])),
                 'SCHSIZE' = as.numeric(as.character(SCH[["SCHSIZE"]]))
) 

sch0$PRIVATESCH[sch0$PRIVATESCH == "PRIVATE"] <- "private"
sch0$PRIVATESCH[sch0$PRIVATESCH == "PUBLIC"] <- "public"
sch0$PRIVATESCH[sch0$PRIVATESCH == "PUBLIC "] <- "public"
sch0$PRIVATESCH[sch0$PRIVATESCH == "public "] <- "public"

sch0$PRIVATESCH[sch0$PRIVATESCH == "invalid"] <- NA
sch0$PRIVATESCH[sch0$PRIVATESCH == "       "] <- NA
sch0$PRIVATESCH[sch0$PRIVATESCH == "missing"] <- NA

sch = sch0

#_____________________________________________
## 1.B) Students questionnaire data file import

#_____________
## Raw Data import
## downloaded from https://webfs.oecd.org/pisa2018/SPSS_STU_QQQ.zip
STU = read.spss('data/students/raw_data_students.sav', reencode='utf-8')

stu0 = data.frame('CNT' = STU[["CNT"]],
                 'CNTSCHID' = STU[["CNTSCHID"]],
                 'PV1MATH' = as.numeric(as.character(STU[["PV1MATH"]])),
                 'ESCS' = as.numeric(as.character(STU[["ESCS"]])) )


# https://www.oecd.org/pisa/pisa-for-development/pisafordevelopment2018technicalreport/PISA-D%20TR%20Chapter%2015%20-%20Proficiency%20Scale%20Construction%20-%20final.pdf
stu0$MATH1below = ifelse(stu0$PV1MATH <= 482.38, 1, 0) # low achieving students definition

by_sch = stu0 %>% group_by(CNTSCHID)

schools = by_sch %>% summarise(
  SCH_TESTED = n(),
  sum_MATH1below = sum(MATH1below), 
  mean_ESCS = mean(ESCS) 
)



#___________________________________________
## 2) CREATION OF THE PREPROCESSED DATAFRAME

df0 <- merge(x = sch, y=schools, 
            by = 'CNTSCHID', all.x=TRUE)

df = df0 

df = df[df$CNT != 'Vietnam' & df$CNT != 'Singapore' & df$CNT != 'United Kingdom', ] 

# restrict to schools with more than 10 students and strictly less than 20
df = df %>% filter(SCH_TESTED >= 10, SCH_TESTED < 20) 
# df = df %>% filter(SCHSIZE >= 10)

unique(df$CNT)
length(unique(df$CNTSCHID))

# scale mean_ESCS with respect to the country CNT
boxplot(df$mean_ESCS ~ df$CNT)
df <- within(df, mean_ESCS_std <- ave(mean_ESCS, CNT, FUN=function(x) (scale(x))))
boxplot(df$mean_ESCS_std ~ df$CNT)
# boxplot(df$mean_ESCS ~ df$PRIVATESCH)

df$Y_MATH1 = df$sum_MATH1below
df$Y_MATH1_rate = df$sum_MATH1below/df$SCH_TESTED


missing_plot(df, dependent = "Y_MATH1", explanatory = c("mean_ESCS_std", "SCHSIZE", 'CNT'))

missing_pattern(df, dependent = "Y_MATH1", explanatory = c("mean_ESCS_std", "SCHSIZE", "CNT"))

missing_pairs(df, dependent = "Y_MATH1", explanatory = c("mean_ESCS_std", "SCHSIZE", "PRIVATESCH"), position = "fill")

# mcar_test(data.frame(y = df$Y_MATH1, z = df$mean_ESCS_std, W = df$SCHSIZE))

# DOING LISTWISE DELETION?
df = na.omit(df) 
df

unique(df$CNT)

# compute Y_BIN_MATH
df$Y_BIN_MATH1 = ifelse(df$Y_MATH1 > 10, 1, 0)
table(df$Y_BIN_MATH1)

# save as a .csv
write.csv(df,"/data/csv/OECD_data.csv", row.names = FALSE)



#_________________
## A) POISSON response

# fit the model without fixed intercept
modPoi = glmer(Y_MATH1 ~ -1 + (1 | CNT) + scale(mean_ESCS_std) + scale(SCHSIZE), 
               data = df, offset = SCH_TESTED,
               family = poisson(link = "log"))

summary(modPoi)

df$pred_GLMM_Poi <- predict(modPoi, df, type="response")

for(i in sort(as.vector(ranef(modPoi)$CNT$`(Intercept)`))){
  print(i)
}


rr1 <- ranef(modPoi)
dd <- as.data.frame(rr1)
if (require(ggplot2)) {
  ggplot(dd, aes(y=grp,x=condval)) +
    geom_point() + facet_wrap(~term,scales="free_x") +
    geom_errorbarh(aes(xmin=condval -2*condsd,
                       xmax=condval +2*condsd), height=0) + 
    theme_bw() + xlab("") + ylab('') +
    ggtitle("Poisson response")
}



# DATA EXPLORATION

ggplot(df, aes(x = SCHSIZE, y = sum_MATH1below)) +
  geom_point(alpha = 0.6) +
  #geom_smooth(method = "loess", color = "blue", se = TRUE) +
  labs(
    x = "School Size (Number of Students)",
    y = "Number of Low-Achieving Math Students",
    title = "Relationship Between School Size and Low-Achieving Students"
  ) +
  theme_minimal()


ggplot(df, aes(x = CNT, y = Y_MATH1_rate, fill = CNT)) +
  geom_boxplot() +
  labs(
    x = "Country",
    y = "Rate of Low-Achieving Math Students",
    title = "Distribution of Rate of Low-Achieving Math Students by Country"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  # Rotate x-axis labels for readability


ggplot(df, aes(x = mean_ESCS_std, y = sum_MATH1below)) +
  geom_point(alpha = 0.6, color = "red") +
  #geom_smooth(method = "lm", color = "darkgreen", se = TRUE) +
  labs(
    x = "Standardized Mean ESCS",
    y = "Number of Low-Achieving Math Students",
    title = "Relationship Between Mean ESCS and Low-Achieving Students"
  ) +
  theme_minimal()

ggplot(df, aes(x = CNT, y = mean_ESCS_std, fill = CNT)) +
  geom_boxplot() +
  labs(
    x = "Country",
    y = "Standardized Mean ESCS",
    title = "Distribution of Standardized ESCS by Country"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  # Rotate x-axis labels for readability


ggplot(df, aes(x = PRIVATESCH, y = SCHSIZE, fill = PRIVATESCH)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white") +
  labs(
    x = "School Type",
    y = "School Size (Number of Students)",
    title = "Distribution of School Size by School Type"
  ) +
  theme_minimal()





