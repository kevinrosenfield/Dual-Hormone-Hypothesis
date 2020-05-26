library(pacman); pacman::p_load(tidyverse, dplyr, lmerTest, psych, sjPlot, sjmisc, sjlabelled, data.table, reshape2, ggcorrplot, stargazer, gridExtra)

###############################################################################################################################
################################################                               ################################################
################################################        DATA PREPERATION       ################################################
################################################                               ################################################
###############################################################################################################################

# read in data

df <- as.data.frame(read.csv("df.csv"))

df <- df[which(df$group =="CON" | is.na(df$group)),]
df$group<-droplevels(df$group)
df$session[is.na(df$session)] <- 1

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
df[df$sex == "f", headArray("ksoq1", "ksoq6", df)] <- 6 - df[df$sex == "f", headArray("ksoq1", "ksoq6", df)]
df[df$sex == "m", "ksoq8"] <- 6 - df[df$sex == "m", "ksoq8"]
df[df$sex == "f", "ksoq7"] <- 6 - df[df$sex == "f", "ksoq7"] 

# scale scores for each question (seperately by sex)

df[(ncol(df) + 1):(ncol(df) + 31) ] <- df[headArray("cgnq1_score", "cgnq24b_4_score", df)]
df <- df %>%
  group_by(sex) %>%
  mutate_at(names(df[(ncol(df) - 30):ncol(df)]), scale)

names(df) <- c(names(df[1:(ncol(df) - 31)]), paste("cgnq", 1:23, "_scaled", sep = ""),
               paste("cgnq24a", 1:4, "_scaled", sep = ""), paste("cgnq24b", 1:4, "_scaled", sep = ""))

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

# composite SOI psychology scores

df$soi_psych <- rowMeans(df[c("soi_des", "soi_attit")], na.rm = TRUE)
df$soi_mean <- rowMeans(df[c("soit1", "soit9")], na.rm = TRUE)

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

# calculate mean SDI z scores for men and women
df <- df %>%
  group_by(sex) %>%
  mutate(sdi_mean.z = scale(sdi_mean))

df <- df %>%
  group_by(subID) %>%
  mutate(ihh_C.cm = mean(ihh_C, na.rm =TRUE),
         ihh_T.cm = mean(ihh_T, na.rm =TRUE),
         sdi_mean.cm = mean(sdi_mean, na.rm = TRUE),
         sdi_sol.cm = mean(sdi_sol, na.rm = TRUE),
         sdi_dyad.cm = mean(sdi_dyad, na.rm = TRUE),
         soi_mean.cm = mean(soi_mean, na.rm = TRUE),
         sdi_des.cm = mean(soi_des, na.rm = TRUE),
         sdi_attit.cm = mean(soi_attit, na.rm = TRUE),
         sdi_behav.cm = mean(soi_behav, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(ihh_C.cwc = ihh_C - ihh_C.cm,
         ihh_T.cwc = ihh_T - ihh_T.cm,
         ihh_T.ccm = ihh_T.cm - mean(ihh_T.cm, na.rm = TRUE),
         ihh_C.ccm = ihh_C.cm - mean(ihh_C.cm, na.rm = TRUE))

df[is.nan(df)] <- NA


###############################################################################################################################
################################################                               ################################################
################################################     DESCRIPTIVE STATISTICS    ################################################
################################################                               ################################################
###############################################################################################################################

df[df$study == "ihh",] %>%
  group_by(sex, n_sessions) %>%
  summarise(
    n = n()
  )

samples <- df %>%
  pivot_wider(names_from = session,
              values_from = colnames(df[12:length(colnames(df))])) %>%
  filter(study == "ihh")

samples %>%
  filter((!is.na(sdi_mean_1) & !is.na(ihh_T_1) & !is.na(ihh_C_1) & (is.na(sdi_mean_2) | is.na(ihh_T_2) | is.na(ihh_C_2))) |
           ((is.na(sdi_mean_1) | is.na(ihh_T_1) | is.na(ihh_C_1)) & !is.na(sdi_mean_2) & !is.na(ihh_T_2) & !is.na(ihh_C_2))) %>%
  group_by(sex) %>%
  summarize(n()) %>%
  rename(sdi_b = `n()`) %>%
  bind_cols(samples %>%
          filter(!is.na(sdi_mean_1) & !is.na(ihh_T_1) & !is.na(ihh_C_1)) %>%
          filter(!is.na(sdi_mean_2) & !is.na(ihh_T_2) & !is.na(ihh_C_2)) %>%
          group_by(sex) %>%
          summarize(n()) %>%
          select(-sex) %>%
            rename(sdi_w = `n()`)) %>%
  bind_cols(samples %>%
              filter((!is.na(soi_mean_1) & !is.na(ihh_T_1) & !is.na(ihh_C_1) & (is.na(soi_mean_2) | is.na(ihh_T_2) | is.na(ihh_C_2))) |
                       ((is.na(soi_mean_1) | is.na(ihh_T_1) | is.na(ihh_C_1)) & !is.na(soi_mean_2) & !is.na(ihh_T_2) & !is.na(ihh_C_2))) %>%
              group_by(sex) %>%
              summarize(n()) %>%
              rename(soi_b = `n()`)) %>%
              bind_cols(samples %>%
                          filter(!is.na(soi_mean_1) & !is.na(ihh_T_1) & !is.na(ihh_C_1)) %>%
                          filter(!is.na(soi_mean_2) & !is.na(ihh_T_2) & !is.na(ihh_C_2)) %>%
                          group_by(sex) %>%
                          summarize(n()) %>%
                          select(-sex)) %>%
  select(-sex1)

df %>%
  group_by(subID) %>%
  summarise(
    sessions = n_distinct(session)
  ) %>%
  group_by(sessions) %>%
  summarise(
    n = n()
  )

df %>%
  filter(sex == "m", study =="ihh") %>%
  group_by(subID) %>%
  select(cgn_comp, sdi_mean, ihh_T, so, soi_mean, age) %>%
  summarise_all(funs(mean)) %>%
  describe

df %>%
  filter(sex == "f", study =="ihh") %>%
  group_by(subID) %>%
  select(cgn_comp, sdi_mean, ihh_T, so, soi_mean, age) %>%
  summarise_all(funs(mean)) %>%
  describe


describe(df[df$sex == "m", headArray("ksoq1", "ksoq9", df)])
describe(df[df$sex == "f", headArray("ksoq1", "ksoq9", df)])

ksoq_pca_m <- prcomp(na.omit(df[df$sex == "m", headArray("ksoq1", "ksoq9", df)]))
ksoq_pca_f <- prcomp(na.omit(df[df$sex == "f", headArray("ksoq1", "ksoq9", df)]))

ggbiplot(ksoq_pca_f)
ggbiplot(ksoq_pca_m)

detach(package:ggbiplot); detach(package:plyr)

describe(df[df$study=="ihh", c("cgn_comp", "sdi_mean", "all_T", "so",  "soi_behav", "soi_attit", "soi_des", "age")])
describe(df[df$study=="msu", c("cgn_comp", "all_T", "so", "age")])

mean(df$cgn_comp[df$sex=="m"], na.rm = TRUE)
mean(df$cgn_comp[df$sex=="f"], na.rm = TRUE)
mean(df$all_T[df$sex=="m" & df$study == "msu"], na.rm = TRUE)
mean(df$all_T[df$sex=="f"], na.rm = TRUE)

describe(df[df$sex == "m", headArray("cgnq1_scaled", "cgnq3_scaled", df)])
describe(df[df$sex == "f", headArray("cgnq1_scaled", "cgnq3_scaled", df)])

df.m <- df %>% filter(sex == "m", study == "ihh")
df.f <- df %>% filter(sex == "f", study == "ihh")

pairs.panels(
  df.m %>% select(sdi_mean, sdi_sol, sdi_dyad), stars= TRUE, method = "pearson", hist.col = "gray", cex.labels=1.9, ellipses = FALSE)
  
pairs.panels(
  df.f %>% select(sdi_mean, sdi_sol, sdi_dyad), stars= TRUE, method = "pearson", hist.col = "gray", cex.labels=1.9, ellipses = FALSE)
  
pairs.panels(
  df.m %>% select(soi_mean, soi_attit, soi_behav, soi_des), stars= TRUE, method = "pearson", hist.col = "gray", cex.labels=1.9, ellipses = FALSE)

pairs.panels(
  df.f %>% select(soi_mean, soi_attit, soi_behav, soi_des), stars= TRUE, method = "pearson", hist.col = "gray", cex.labels=1.9, ellipses = FALSE)

###############################################################################################################################
##############################################                                     ############################################
##############################################    DUAL-HORMONE HYPOTHESIS TESTS    ############################################
##############################################                                     ############################################
###############################################################################################################################


### HORMONES UNCENTERED: Within- and between-individual effects ARE NOT isolated

# men

m.d <- lmer(sdi_mean ~ ihh_T*ihh_C + (1 | subID), data = df.m) # sdi-2

m.sol <- lmer(sdi_sol ~ ihh_T*ihh_C + (1 | subID), data = df.m) # sdi-2[solitary]

m.dyad <- lmer(sdi_dyad ~ ihh_T*ihh_C + (1 | subID), data = df.m) # sdi-2[dyadic]

m.s <- lmer(soi_mean ~ ihh_T*ihh_C + (1 | subID), data = df.m) # soiR

m.att <- lmer(soi_attit ~ ihh_T*ihh_C + (1 | subID), data = df.m) # soi-R[attitude]

m.beh <- lmer(soi_behav ~ ihh_T*ihh_C + (1 | subID), data = df.m) # soi-R[behavior]

m.des <- lmer(soi_des ~ ihh_T*ihh_C + (1 | subID), data = df.m) # soi-R[desire]

# women

w.d <- lmer(sdi_mean ~ ihh_T*ihh_C + (1 | subID), data = df.f) # sdi-2

w.sol <- lmer(sdi_sol ~ ihh_T*ihh_C + (1 | subID), data = df.f) # sdi-2[solitary]

w.dyad <- lmer(sdi_dyad ~ ihh_T*ihh_C + (1 | subID), data = df.f) # sdi-2[dyadic]

w.s <- lmer(soi_mean ~ ihh_T*ihh_C + (1 | subID), data = df.f) # soiR

w.att <- lmer(soi_attit ~ ihh_T*ihh_C + (1 | subID), data = df.f) # soi-R[attitude]

w.beh <- lmer(soi_behav ~ ihh_T*ihh_C + (1 | subID), data = df.f) # soi-R[behavior]

w.des <- lmer(soi_des ~ ihh_T*ihh_C + (1 | subID), data = df.f) # soi-R[desire]


### HORMONES ARE CENTERED: Within(cwc)- and between(cm)-individual effects ARE isolated

# men

m.d_w <- lmer(sdi_mean ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + (1 | subID), data = df.m) # sdi-2

m.sol_w <- lmer(sdi_sol ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + (1 | subID), data = df.m) # sdi-2[solitary]

m.dyad_w <- lmer(sdi_dyad ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + (1 | subID), data = df.m) # sdi-2[dyadic]

m.s_w <- lmer(soi_mean ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + (1 | subID), data = df.m) # soiR

m.att_w <- lmer(soi_attit ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + (1 | subID), data = df.m) # soi-R[attitude]

m.beh_w <- lmer(soi_behav ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + (1 | subID), data = df.m) # soi-R[behavior]

m.des_w <- lmer(soi_des ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + (1 | subID), data = df.m) # soi-R[desire]

# women

w.d_w <- lmer(sdi_mean ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + (1 | subID), data = df.f) # sdi-2

w.sol_w <- lmer(sdi_sol ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + (1 | subID), data = df.f) # sdi-2[solitary]

w.dyad_w <- lmer(sdi_dyad ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + (1 | subID), data = df.f) # sdi-2[dyadic]

w.s_w <- lmer(soi_mean ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + (1 | subID), data = df.f) # soiR

w.att_w <- lmer(soi_attit ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + (1 | subID), data = df.f) # soi-R[attitude]

w.beh_w <- lmer(soi_behav ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + (1 | subID), data = df.f) # soi-R[behavior]

w.des_w <- lmer(soi_des ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + (1 | subID), data = df.f) # soi-R[desire]

# Have to do this to use stargazer

class(w.d) <- "lmerMod"
class(w.sol) <- "lmerMod"
class(w.dyad) <- "lmerMod"
class(w.s) <- "lmerMod"
class(w.att) <- "lmerMod"
class(w.beh) <- "lmerMod"
class(w.des) <- "lmerMod"

class(m.d) <- "lmerMod"
class(m.sol) <- "lmerMod"
class(m.dyad) <- "lmerMod"
class(m.s) <- "lmerMod"
class(m.att) <- "lmerMod"
class(m.beh) <- "lmerMod"
class(m.des) <- "lmerMod"

class(w.d_w) <- "lmerMod"
class(w.sol_w) <- "lmerMod"
class(w.dyad_w) <- "lmerMod"
class(w.s_w) <- "lmerMod"
class(w.att_w) <- "lmerMod"
class(w.beh_w) <- "lmerMod"
class(w.des_w) <- "lmerMod"

class(m.d_w) <- "lmerMod"
class(m.sol_w) <- "lmerMod"
class(m.dyad_w) <- "lmerMod"
class(m.s_w) <- "lmerMod"
class(m.att_w) <- "lmerMod"
class(m.beh_w) <- "lmerMod"
class(m.des_w) <- "lmerMod"

# FOR BETWEEN SUBJECT MODELS, WOMEN
stargazer(w.d, w.sol, w.dyad, w.s, w.att, w.beh, w.des, single.row = F, align = F, intercept.bottom = T, no.space = T, 
          dep.var.caption = "", type = "text", out = "women_between.html", report = "vc*p", 
          star.char = c("*", "**"), digits = 3, star.cutoffs = c(0.05, 0.01), 
          notes = "*p<0.05;**p<0.01", notes.append = F,
          column.labels = c("Total desire", "Solitary desire", "Dyadic Desire", "Uncommited Sex Total",
                            "Uncommited Sex Attitude", "Uncommited Sex Behavior", "Uncommited Sex Desire" ), omit.stat = "all")

# FOR BETWEEN SUBJECT MODELS, MEN
stargazer(m.d, m.sol, m.dyad, m.s, m.att, m.beh, m.des, single.row = F, align = F, intercept.bottom = T, no.space = T, 
          dep.var.caption = "", type = "text", out = "men_between.html", report = "vc*p", 
          star.char = c("*", "**"), digits = 3, star.cutoffs = c(0.05, 0.01), 
          notes = "*p<0.05;**p<0.01", notes.append = F,
          column.labels = c("Total desire", "Solitary desire", "Dyadic Desire", "Uncommited Sex Total",
                            "Uncommited Sex Attitude", "Uncommited Sex Behavior", "Uncommited Sex Desire" ), omit.stat = "all")

# FOR WITHIN SUBJECT MODELS, WOMEN
stargazer(w.d_w, w.sol_w, w.dyad_w, w.s_w, w.att_w, w.beh_w, w.des_w, single.row = F, align = F, intercept.bottom = T, no.space = T, 
          dep.var.caption = "", type = "text", out = "women_within.html", report = "vc*p", 
          star.char = c("*", "**"), digits = 3, star.cutoffs = c(0.05, 0.01), 
          notes = "*p<0.05;**p<0.01", notes.append = F,
          column.labels = c("Total desire", "Solitary desire", "Dyadic Desire", "Uncommited Sex Total",
                            "Uncommited Sex Attitude", "Uncommited Sex Behavior", "Uncommited Sex Desire" ), omit.stat = "all")

# FOR WITHIN SUBJECT MODELS, MEN
stargazer(m.d_w, m.sol_w, m.dyad_w, m.s_w, m.att_w, m.beh_w, m.des_w, single.row = F, align = F, intercept.bottom = T, no.space = T, 
          dep.var.caption = "", type = "text", out = "men_within.html", report = "vc*p", 
          star.char = c("*", "**"), digits = 3, star.cutoffs = c(0.05, 0.01), 
          notes = "*p<0.05;**p<0.01", notes.append = F,
          column.labels = c("Total desire", "Solitary desire", "Dyadic Desire", "Uncommited Sex Total",
                            "Uncommited Sex Attitude", "Uncommited Sex Behavior", "Uncommited Sex Desire" ), omit.stat = "all")













# RCGN predict current sexual orientation in women?

soCGN_ihh_f <- lmer(ksoq9 ~ cgn_comp + (1|subID), data = subset(df, study == "ihh" & group == "CON" & sex == "f")); summary(soCGN_ihh_f); hist(residuals(soCGN_ihh_f))
soCGN_msu_f <- lmer(so ~ cgn_comp + (1|sibID), data = subset(df, study == "msu" & sex == "f")); summary(soCGN_msu_f); hist(residuals(soCGN_msu_f))

# RCGN predict current sexual orientation in men?

soCGN_ihh_m <- lm(ksoq9 ~ cgn_comp, data = subset(df, study == "ihh" & group == "CON" & sex == "f")); summary(soCGN_ihh_m); hist(residuals(soCGN_ihh_m))
soCGN_msu_m <- lm(so ~ cgn_comp, data = subset(df, study == "msu" & sex == "f")); summary(soCGN_msu_m); hist(residuals(soCGN_msu_m))

# Does current T predict RCGN in women?

cgnT_ihh_f <- lm(cgn_comp ~ all_T, data = subset(df, study == "msu" & group != "REL" & sex == "f")); summary(cgnT_ihh_f); hist(residuals(cgnT_ihh_f))
cgnT_msu_f <- lmer(cgn_comp ~ log(all_T) + (1|sibID), data = subset(df, study == "msu" & sex == "f")); summary(cgnT_msu_f); hist(residuals(cgnT_msu_f))

# Does current T predict RCGN in men?

cgnT_ihh_m <- lm(cgn_comp ~ all_T, data = subset(df, study == "msu" & group != "REL" & sex == "m")); summary(cgnT_ihh_m); hist(residuals(cgnT_ihh_m))
cgnT_msu_m <- lmer(cgn_comp ~ log(all_T) + (1|sibID), data = subset(df, study == "msu" & sex == "m")); summary(cgnT_msu_m); hist(residuals(cgnT_msu_m))

tab_model(cgnT_msu_f)

# Does current T predict current sexual orientation in women?

soT_ihh_f <- lm(ksoq9 ~ log(all_T), data = subset(df, study == "ihh" & group == "CON" & sex == "f")); summary(soT_ihh_f); hist(residuals(soT_ihh_f))
soT_msu_f <- lmer(cgn_comp ~ all_T + (1|sibID), data = subset(df, study == "msu" & sex == "f")); summary(soT_msu_f); hist(residuals(soT_msu_f))

# Does current T predict current sexual orientation in men?

soT_ihh_m <- lm(ksoq9 ~ log(all_T), data = subset(df, study == "ihh" & group == "CON" & sex == "m")); summary(soT_ihh_m); hist(residuals(soT_ihh_m))
soT_msu_m <- lmer(cgn_comp ~ all_T + (1|sibID), data = subset(df, study == "msu" & sex == "m")); summary(soT_msu_m); hist(residuals(soT_msu_m))

# Does current T predict current sexual desire in women?

sdiT_ihh_f <- lm(sdi_mean ~ log(all_T), data = subset(df, study == "ihh" & group == "CON" & sex == "f")); summary(sdiT_ihh_f); hist(residuals(sdiT_ihh_f))

# Does current T predict current sexual desire in men?

sdiT_ihh_m <- lm(sdi_sol ~ sdi_dyad, data = subset(df, study == "ihh" & group == "CON" & sex == "m")); summary(sdiT_ihh_m); hist(residuals(sdiT_ihh_m))

# Does current T predict current interest in casual sex in women?

soibT_ihh_f <- lm(soi_behav ~ log(all_T), data = subset(df, study == "ihh" & group == "CON" & sex == "f")); summary(soibT_ihh_f); hist(residuals(soibT_ihh_f))
soidT_ihh_f <- lm(log(all_T) ~ soi_psych, data = subset(df, study == "ihh" & group == "CON" & sex == "f")); summary(soidT_ihh_f); hist(residuals(soidT_ihh_f))
soiaT_ihh_f <- lm(soi_attit ~ log(all_T), data = subset(df, study == "ihh" & group == "CON" & sex == "f")); summary(soiaT_ihh_f); hist(residuals(soiaT_ihh_f))

# Does current T predict current sexual casual sex in men?

soibT_ihh_m <- lm(log(all_T) ~ soi_behav + soi_psych , data = subset(df, study == "ihh" & group == "CON" & sex == "m")); summary(soibT_ihh_m); hist(residuals(soibT_ihh_m))
soidT_ihh_m <- lm(soi_des ~ log(all_T), data = subset(df, study == "ihh" & group == "CON" & sex == "m")); summary(soidT_ihh_m); hist(residuals(soidT_ihh_m))
soiaT_ihh_m <- lm(soi_attit ~ log(all_T), data = subset(df, study == "ihh" & group == "CON" & sex == "m")); summary(soiaT_ihh_m); hist(residuals(soiaT_ihh_m))
