library(dplyr); library(lmerTest); library(psych)

# read in data

df <- as.data.frame(read.csv("df.csv"))

df <- df[which(df$group !="REL" | is.na(df$group)),]
df$group<-droplevels(df$group)

# change variable types

headArray <- function(first,last, df = df) match(first, names(df)):match(last, names(df))

df[headArray("ID", "ihh_race", df)] <- data.frame(lapply(df[headArray("ID", "ihh_race", df)], function(x) as.factor(x)))
df[headArray("age", "cgnq24b_4_score", df)] <- data.frame(lapply(df[headArray("age", "cgnq24b_4_score", df)], function(x) as.numeric(as.character(x))))

# fix ksoq varialbes

df$ksoq2[df$ksoq2 == 8] <- NA
df$ksoq4[df$ksoq4 == 8] <- NA
df$ksoq6[df$ksoq6 == 8] <- NA
df[headArray("ksoq1", "ksoq6", df)] <- df[headArray("ksoq1", "ksoq6", df)] - 1
df[df$sex == "f", headArray("ksoq1", "ksoq6", df)] <- 6 - df[df$sex == "f", headArray("ksoq1", "ksoq6", df)]

describe(df[df$sex == "m", headArray("ksoq1", "ksoq9", df)])
describe(df[df$sex == "f", headArray("ksoq1", "ksoq9", df)])

# scale scores for each question (seperately by sex)

df[(ncol(df) + 1):(ncol(df) + 31) ] <- df[headArray("cgnq1_score", "cgnq24b_4_score", df)]
df <- df %>%
  group_by(sex) %>%
  mutate_at(names(df[(ncol(df) - 30):ncol(df)]), scale)

names(df) <- c(names(df[1:(ncol(df) - 31)]), paste("cgnq", 1:23, "_scaled", sep = ""),
               paste("cgnq24a", 1:4, "_scaled", sep = ""), paste("cgnq24b", 1:4, "_scaled", sep = ""))

describe(df[df$sex == "m", headArray("cgnq1_scaled", "cgnq3_scaled", df)])
describe(df[df$sex == "f", headArray("cgnq1_scaled", "cgnq3_scaled", df)])

##############

# remove raw responses, unscaled scores; calculate  composites

df <- df[-c(headArray("cgnq1_raw", "cgnq24b_4_score", df))] 

df <- df[-c(match("cgnq24a2_scaled", names(df)), match("cgnq24b2_scaled", names(df)))]
df$cgn_comp <- rowMeans(df[headArray("cgnq1_scaled", "cgnq24b4_scaled", df)], na.rm = TRUE)
df$so <- rowMeans(df[c("Kinsey_Attr", "Kinsey_Fantasy", "ksoq5")], na.rm =TRUE)
df$mean_T <- rowMeans(df[c("AM_T","PM_T")], na.rm = TRUE)
df$all_T <- rowMeans(df[c("mean_T","ihh_T")], na.rm = TRUE)

# this probs isn't necessary for you but I need it; turns NaNs in NAs

is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
df[is.nan(df)] <- NA

# describe

describe(df[df$study=="ihh", c("cgn_comp", "sdi_mean", "all_T", "so",  "soi_behav", "soi_attit", "soi_des")])
describe(df[df$study=="msu", c("cgn_comp", "all_T", "so", "age")])

# average cgn for each sex

mean(df$cgn_comp[df$sex=="m"], na.rm = TRUE)
mean(df$cgn_comp[df$sex=="f"], na.rm = TRUE)
mean(df$all_T[df$sex=="m" & df$study == "msu"], na.rm = TRUE)
mean(df$all_T[df$sex=="f"], na.rm = TRUE)

################## Race/ethnicity for both
df$race_combined<- 5 #other
df$race_combined[df$ihh_race == 5 & df$ihh_ethnicity == "n"] <- 1 #white non hispanic
df$race_combined[df$ihh_race == 2 & df$ihh_ethnicity == "n"] <- 2 #asian NH
df$race_combined[df$ihh_race == 3 & df$ihh_ethnicity == "n"] <- 3 #black NH
df$race_combined[df$ihh_ethnicity == "h"] <- 4

df$race_combined[df$msu_ethnicity == 6] <- 1 #white non hispanic
df$race_combined[df$msu_ethnicity == 2] <- 2
df$race_combined[df$msu_ethnicity == 3] <- 3
df$race_combined[df$msu_ethnicity == 4] <- 4

# calculate SDI mean z scores for men and women
df <- df %>%
  group_by(sex) %>%
  mutate(SDI_Z = scale(sdi_mean))


cgnT_ihh <- lm(cgn_comp ~ all_T, data = subset(df, study == "msu" & group != "REL" & sex == "m")); summary(cgnT_ihh); hist(residuals(cgnT_ihh))
cgnT_msu <- lmer(cgn_comp ~ log(all_T) + (1|sibID), data = subset(df, study == "msu" & sex == "m")); summary(cgnT_msu); hist(residuals(cgnT_msu))


soT_ihh <- lm(ksoq9 ~ log(all_T), data = subset(df, study == "ihh" & group == "CON" & sex == "f")); summary(soT_ihh); hist(residuals(soT_ihh))
soT_msu <- lmer(cgn_comp ~ all_T + (1|sibID), data = subset(df, study == "msu" & sex == "f")); summary(soT_msu); hist(residuals(soT_msu))

soCGN_ihh <- lm(ksoq9 ~ cgn_comp, data = subset(df, study == "ihh" & group == "CON" & sex == "f")); summary(soCGN_ihh); hist(residuals(soCGN_ihh))

soCGN_msu <- lm(so ~ cgn_comp, data = subset(df, study == "msu" & sex == "m")); summary(soCGN_msu); hist(residuals(soCGN_msu))

amt_so <- lm(formula = scale(log(am_t)) ~ so); summary(amt_so)
pmt_cgn <- lm(formula = scale(log(pm_t)) ~ cgn); summary(pmt_cgn)
pmt_so <- lm(formula = scale(log(pm_t)) ~ so); summary(pmt_so)
cgn_so <- lm(formula = cgn ~ so); summary(cgn_so)

scatter.smooth(scale(log(am_t)),cgn)

hist(MSUMEN_AMT, breaks = 20, freq = FALSE)
cor.test(MSUMEN_AMT,MSUMEN_CGNQ, method="spearman")
shapiro.test(MSUMEN_AMT)
qqnorm(MSUMEN_AMT)
hist(MSUMEN_AMT, breaks = 20, freq = TRUE)
kurtosis(MSUMEN_AMT_RAW)
cor.test(MSUMEN_AMT_RAW,MSUMEN_CGNQ, method="spearman")
shapiro.test(MSUMEN_AMT_RAW)
qqnorm(MSUMEN_AMT_RAW)
hist(MSUMEN_AMT_RAW, breaks = 20, freq = TRUE)
scatter.smooth(MSUMEN_AMT_RAW,MSUMEN_PMT_RAW)
scatter(MSUMEN_AMT,MSUMEN_PMT)
hist(MSUMEN_CGNQ, breaks = 20, freq = TRUE)
mean(MSUMEN_CGNQ)
mean(MSUDATA$average_24_cgnq)
cor.test(MSUMEN_AMT, MSUMEN_CGNQ)
rma(-0.1579)
plot(MSUMEN_AMT,MSUMEN_CGNQ, cex=1.3, pch=19, col="brown2", main = "Testosterone plotted against CGN", xlab="Testosterone (logged and standardized)", ylab="CGN composite score", bty="n", xlim=c(-3,3), ylim=c(-1.5,1))
lineMSU <- lm(MSUMEN_AMT ~ MSUMEN_CGNQ)
abline(lineMSU, col="brown2", lwd=3)
par(new=T)
#
#
#
#
MSUDATA<- as.data.frame(read.csv("~/Documents/Penn State/CHH/MSU/MSUCSV.csv"))
MSUMEN<- MSUDATA[(1:191),]
# MSUMEN<- MSUMEN[-c(81, 108, 125, 126, 127),]
#127 has CGNQ 5 SD lower than mean; 108, 125, 126 have AM way higher than PMtest; 81 possible blood contamination
MSUMEN_CGNQ <- as.numeric(as.character((MSUMEN[,9])))
MSUMEN_AMT_RAW <- as.numeric(as.character((MSUMEN[,10])))
MSUMEN_PMT_RAW <- as.numeric(as.character((MSUMEN[,11])))
MSUMEN_AMT <- as.numeric(as.character((MSUMEN[,14])))
MSUMEN_PMT <- as.numeric(as.character((MSUMEN[,15])))
MSUMEN_ORIENT <- as.numeric(as.character((MSUMEN[,16])))
MSUMEN_VAR <- data.frame(MSUMEN_CGNQ, MSUMEN_AMT, MSUMEN_AMT_RAW, MSUMEN_PMT, MSUMEN_PMT_RAW, MSUMEN_ORIENT)
colnames(MSUMEN_VAR) <- c("CGNQ", "AMTEST", "AMTESTRAW", "PMTEST", "PMTESTRAW", "ORIENT")
MSUMEN_VAR
MSUMEN_AMTxCGNQ <- lm(formula = AMTEST ~ CGNQ, data = MSUMEN_VAR)
MSUMEN_AMTxORIENT <- lm(formula = AMTEST ~ ORIENT, data = MSUMEN_VAR)
MSUMEN_PMTxCGNQ <- lm(formula = PMTEST ~ CGNQ, data = MSUMEN_VAR)
MSUMEN_PMTxORIENT <- lm(formula = PMTEST ~ ORIENT, data = MSUMEN_VAR)
MSUMEN_AMTxCGNQxORIENT <- lm(CGNQ ~ MSUMEN_AMT + ORIENT, data = MSUMEN_VAR)
summary(MSUMEN_AMTxCGNQ)
summary(MSUMEN_AMTxORIENT)
summary(MSUMEN_CGNQxORIENT)
summary(MSUMEN_PMTxCGNQ)
summary(MSUMEN_PMTxORIENT)
summary(MSUMEN_AMTxCGNQxORIENT)
scatter.smooth(MSUMEN_AMT,MSUMEN_CGNQ)
hist(MSUMEN_AMT, breaks = 20, freq = FALSE)
cor.test(MSUMEN_AMT,MSUMEN_CGNQ, method="spearman")
shapiro.test(MSUMEN_AMT)
qqnorm(MSUMEN_AMT)
hist(MSUMEN_AMT, breaks = 20, freq = TRUE)
kurtosis(MSUMEN_AMT_RAW)
cor.test(MSUMEN_AMT_RAW,MSUMEN_CGNQ, method="spearman")
shapiro.test(MSUMEN_AMT_RAW)
qqnorm(MSUMEN_AMT_RAW)
hist(MSUMEN_AMT_RAW, breaks = 20, freq = TRUE)
scatter.smooth(MSUMEN_AMT_RAW,MSUMEN_PMT_RAW)
scatter(MSUMEN_AMT,MSUMEN_PMT)
hist(MSUMEN_CGNQ, breaks = 20, freq = TRUE)
mean(MSUMEN_CGNQ)
mean(MSUDATA$average_24_cgnq)
cor.test(MSUMEN_AMT, MSUMEN_CGNQ)
rma(-0.1579)
plot(MSUMEN_AMT,MSUMEN_CGNQ, cex=1.3, pch=19, col="brown2", main = "Testosterone plotted against gender conformity", xlab="Testosterone Z-scores", cex.lab=1.5, ylab="Gender conformity", bty="n", xlim=c(-3,2.6), ylim=c(-1,1), mgp = c(2.3, 1, 0))
lineMSU <- lm(MSUMEN_AMT ~ MSUMEN_CGNQ)
abline(lineMSU, col="brown2", lwd=3)
par(new=T)

ggplot(MSUMEN, aes(x = MSUMEN_AMT, y = MSUMEN_CGNQ)) + 
  geom_point(color="red", size=3) +
  stat_smooth(method = "lm", color = "red") + 
  xlim(-3, 3.2) +
  ylim(-1, 1) +
  ggtitle("Testosterone plotted against gender conformity", subtitle = NULL) +
  theme(plot.title = element_text(face="bold", size=18, hjust=.4)) +
  xlab("Testosterone Z-scores") +
  ylab("Gender conformity") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold", vjust=.1)) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
par(new=T)

ggplot() + points(MSUMEN_VAR, aes(x = AMTEST, y = CGNQ))
length(MSUMEN_VAR$CGNQ)
length(MSUMEN_VAR$AMTEST)
CHH_FRAME <- data.frame(c(ALL_VAR$CGNQ,0,0), c(ALL_VAR$TEST,0,0))
ggplot(MSUMEN, aes(x = MSUMEN_AMT, y = MSUMEN_CGNQ)) +
         geom_point() +
         geom_point(data = CHH_FRAME, color = "red")
comp7xT_MSU <- lm(MSUMEN_COMP7 ~ MSUMEN_AMT)
summary(comp7xT_MSU)
cor.test(MSUMEN_COMP7, MSUMEN_AMT)
mean(MSUMEN_CGNQ)
sd(MSUMEN_CGNQ)
mean(MSUMEN_AMT)
sd(MSUMEN_AMT)
