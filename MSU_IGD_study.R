library(pacman); pacman::p_load(tidyverse, lmerTest, psych, sjPlot, sjmisc, sjlabelled,
                                data.table, reshape2, ggcorrplot, stargazer, gridExtra, htmlTable)

###############################################################################################################################
################################################                               ################################################
################################################        DATA PREPERATION       ################################################
################################################                               ################################################
###############################################################################################################################

# read in data

df <- as.data.frame(read.csv("df.csv"))

# change variable types

df <- df %>%
  mutate_at(vars(ID:ihh_race), as.factor) %>%
  mutate_at(vars(age:cgnq24b4_score), as.character) %>%
  mutate_at(vars(age:cgnq24b4_score), as.numeric)

# fix ksoq varialbes

df$ksoq2[df$ksoq2 == 8 | df$ksoq4 == 8 | df$ksoq6 == 8] <- NA

df <- df %>%
  mutate_at(vars(starts_with("ksoq")), function(.) .  - 1)

df.f <- df %>%
  filter(sex == "f") %>%
  mutate_at(vars(ksoq1:ksoq7), function (.) 6 - .)

df <- df %>%
  filter(sex == "m") %>%
  mutate_at(vars(ksoq8), function (.) 6 - .) %>%
  bind_rows(df.f)

remove(df.f)

# scale cgnq scores for each question (seperately by sex)

df[(ncol(df) + 1):(ncol(df) + 31)] <- df %>%
  select(c(sex, cgnq1_score:cgnq24b4_score)) %>%
  group_by(sex) %>%
  mutate_at(vars(cgnq1_score:cgnq24b4_score), scale) %>%
  rename_at(vars(ends_with("score")), funs(str_replace(., "score", "scaled"))) %>%
  ungroup() %>%
  select(-sex)

# remove CGNQ raw responses and unscaled scores; calculate  composites

df <- df %>%
  select(-c(cgnq1_raw:cgnq24b4_score, cgnq24a2_scaled, cgnq24b2_scaled)) %>%
  mutate(cgn_comp = rowMeans(select(., cgnq1_scaled:cgnq24b4_scaled), na.rm = TRUE),
         so = rowMeans(select(., c(Kinsey_Attr, Kinsey_Fantasy, ksoq5)), na.rm = TRUE),
         mean_T = rowMeans(select(., c(AM_T, PM_T)), na.rm = TRUE),
         soi_mean = rowMeans(select(., c(soit1:soit9)), na.rm = TRUE),
         soi_psych = rowMeans(select(., c(soi_des, soi_attit)), na.rm = TRUE))

# homogonize race/ethnicity variables for both studies

df$race_combined <- 5 #other
df$race_combined[(df$ihh_race == 5 & df$ihh_ethnicity == "n") | df$msu_ethnicity == 6] <- 1 #white non hispanic
df$race_combined[(df$ihh_race == 2 & df$ihh_ethnicity == "n") | df$msu_ethnicity == 2] <- 2 #asian NH
df$race_combined[(df$ihh_race == 3 & df$ihh_ethnicity == "n") | df$msu_ethnicity == 3] <- 3 #black NH
df$race_combined[df$ihh_ethnicity == "h" | df$msu_ethnicity == 4] <- 4

# calculate mean SDI z scores for men and women
df <- df %>%
  group_by(sex) %>%
  mutate(sdi_mean.z = scale(sdi_mean))

# calculate subject means and subject mean-center hormone variables; center subject-means; scale subject means and subject-centered data within each sex

df <- df %>%
  group_by(subID) %>%
  mutate(ihh_C.cm = mean(ihh_C, na.rm =TRUE),
         ihh_T.cm = mean(ihh_T, na.rm =TRUE),
         ihh_P.cm = mean(ihh_P, na.rm =TRUE),
         ihh_E.cm = mean(ihh_E, na.rm =TRUE),
         ihh_C.cwc = ihh_C - mean(ihh_C, na.rm =TRUE),
         ihh_T.cwc = ihh_T - mean(ihh_T, na.rm =TRUE),
         ihh_P.cwc = ihh_P - mean(ihh_P, na.rm =TRUE),
         ihh_E.cwc = ihh_E - mean(ihh_E, na.rm =TRUE),
         binary_T = ifelse(ihh_T == min(ihh_T), "low", "high"),
         binary_C = ifelse(ihh_C == min(ihh_C), "low", "high")) %>%
  ungroup() %>%
  group_by(sex) %>%
  mutate(ihh_C.cwc = ihh_C.cwc / sd(ihh_C.cwc, na.rm = TRUE),
         ihh_T.cwc = ihh_T.cwc / sd(ihh_T.cwc, na.rm = TRUE),
         ihh_P.cwc = ihh_P.cwc / sd(ihh_P.cwc, na.rm = TRUE),
         ihh_E.cwc = ihh_E.cwc / sd(ihh_E.cwc, na.rm = TRUE),
         ihh_C.cm = (ihh_C.cm - mean(ihh_C.cm, na.rm = TRUE)) / sd(ihh_C.cm, na.rm = TRUE),
         ihh_T.cm = (ihh_T.cm - mean(ihh_T.cm, na.rm = TRUE)) / sd(ihh_T.cm, na.rm = TRUE),
         ihh_P.cm = (ihh_P.cm - mean(ihh_P.cm, na.rm = TRUE)) / sd(ihh_P.cm, na.rm = TRUE),
         ihh_E.cm = (ihh_E.cm - mean(ihh_E.cm, na.rm = TRUE)) / sd(ihh_E.cm, na.rm = TRUE)) %>%
  ungroup()

# break out datasets by sex/study

df_ihh.m <- df %>% filter(sex == "m", study == "ihh")
df_ihh.f <- df %>% filter(sex == "f", study == "ihh")
df_msu.m <- df %>% filter(sex == "m", study == "msu")
df_msu.f <- df %>% filter(sex == "f", study == "msu")

###############################################################################################################################
################################################                               ################################################
################################################     DESCRIPTIVE STATISTICS    ################################################
################################################                               ################################################
###############################################################################################################################

# count men and women in PSU study with 1 or 2 sessions

df %>%
  filter(study == "ihh", session == 1) %>%
  select(n_sessions, sex) %>%
  group_by(n_sessions, sex) %>%
  summarise(n = n())

# Table 1

tableRows = c("SOI-R ", "SOI: Attitudes", "SOI: Behavior", "SOI: Desire", "SDI-2",
              "SDI: Solitary", "SDI: Dyadic")

tableCols = c("n", "Mean", "SD")

table <- cbind(df_ihh.m %>%
  group_by(subID) %>%
  select(soi_mean, soi_attit, soi_behav, soi_des, sdi_mean, sdi_sol, sdi_dyad) %>%
  summarise_all(funs(mean)) %>%
  describe, as.numeric(rep("",8)))

table <- cbind(table, df_ihh.f %>%
  group_by(subID) %>%
  select(soi_mean, soi_attit, soi_behav, soi_des, sdi_mean, sdi_sol, sdi_dyad) %>%
  summarise_all(funs(mean)) %>%
  describe)

options(table_counter = TRUE)

htmlTable(round(table[2:8,c(2:4, 14, 16:18)], 2),
          align = "c",
          caption = "Responses to SOI-R, SDI-2, and subscales",
          rnames = tableRows, 
          header = c(tableCols, "", tableCols),
          col.rgroup = c("none", "#F7F7F7"),
          cgroup = c("Men", "", "Women"), n.cgroup = c(3, 1, 3))

# Figure 1a-d: Correlations of scales and their subscales; each scale/sex combination presented separately

pairs.panels(
  df_ihh.m %>% select(sdi_mean, sdi_sol, sdi_dyad), stars= TRUE, method = "pearson", hist.col = "gray", cex.labels=1.9, cex.cor = 2, ellipses = FALSE)
  
pairs.panels(
  df_ihh.f %>% select(sdi_mean, sdi_sol, sdi_dyad), stars= TRUE, method = "pearson", hist.col = "gray", cex.labels=1.9, cex.cor = 2, ellipses = FALSE)
  
pairs.panels(
  df_ihh.m %>% select(soi_mean, soi_attit, soi_behav, soi_des), stars= TRUE, method = "pearson", hist.col = "gray", cex.labels=1.9, cex.cor = 2, ellipses = FALSE)

pairs.panels(
  df_ihh.f %>% select(soi_mean, soi_attit, soi_behav, soi_des), stars= TRUE, method = "pearson", hist.col = "gray", cex.labels=1.9, cex.cor = 2, ellipses = FALSE)

# Supplementary Figure 1: Correlations of ALL scales and subscales, seperately by sex

pairs.panels(
  df_ihh.m %>% select(soi_mean, soi_attit, soi_behav, soi_des, sdi_mean, sdi_sol, sdi_dyad), stars= TRUE, method = "pearson", hist.col = "gray", cex.labels=1.9, cex.cor = 2, ellipses = FALSE)

pairs.panels(
  df_ihh.f %>% select(soi_mean, soi_attit, soi_behav, soi_des, sdi_mean, sdi_sol, sdi_dyad), stars= TRUE, method = "pearson", hist.col = "gray", cex.labels=1.9, cex.cor = 2, ellipses = FALSE)

# sex differences on SOI-R, SDI-2, and subject-mean hormones

ihh.sex <- df %>%
  filter(study == "ihh", !is.na(soi_mean), !is.nan(soi_mean), !is.na(sex)) %>%
  select(subID, soi_mean, sex, sdi_mean, ihh_C, ihh_T) %>%
  group_by(subID) %>%
  mutate(soi_mean = mean(soi_mean),
         sdi_mean = mean(sdi_mean),
         C_mean = mean(soi_mean),
         T_mean = mean(sdi_mean))

wilcox.test(ihh.sex$soi_mean ~ ihh.sex$sex)
t.test(ihh.sex$sdi_mean ~ ihh.sex$sex)
t.test(ihh.sex$T_mean ~ ihh.sex$sex)
t.test(ihh.sex$C_mean ~ ihh.sex$sex)

# Within-subject changes in scales, subscales, and hormones, broken out by hormonal contraception use

ihh.hc <- df %>%
  select(ID, sex, n_sessions, session, HC, ihh_C, ihh_T, ihh_P, ihh_E, soi_behav,
         soi_attit, soi_des, soi_mean, sdi_dyad, sdi_sol, sdi_mean) %>%
  filter(!is.na(HC), sex == 'f', n_sessions == 2) %>%
  pivot_longer(cols = ID)

for (k in c(0, 1)){
  for (j in colnames(ihh.hc[5:15])) {
    number = 0
    total = 0
    df_dummy = ihh.hc %>% filter(!is.nan(j), HC == k)
    for (i in unique(df_dummy$value)) {
      if (sum(df_dummy$value == i) == 2) {
        var = abs(df_dummy[[j]][df_dummy$value == i & df_dummy$session == 1] -
                    df_dummy[[j]][df_dummy$value == i & df_dummy$session == 2])
        if (!is.na(var)) {
          number = number + 1
          total = total + var
        }
      }
    }
    hc <- ifelse(k == 1, "hc", "nc")
    print(paste(j, hc, total / number))
  }
}

# Raw hormone data plots; data for both sexes

par(mfrow=c(3,2))

plot(df$ihh_C[!is.na(df$ihh_C)], ylab = "Cortisol")
plot(df$ihh_T[!is.na(df$ihh_T)], ylab = "Testosterone")
plot(df$ihh_P[!is.na(df$ihh_P)], ylab = "Progesterone")
plot(df$ihh_E[!is.na(df$ihh_E)], ylab = "Estradiol")

plot(df$ihh_C[!is.na(df$ihh_C) & df$ihh_C < 40], ylab = "Cortisol - 1 Outlier") # cortisol without 20 sd outlier
plot(df$ihh_T[!is.na(df$ihh_T) & df$ihh_T < 145], ylab = "Testosterone - 3 Outliers") # testosterone without ~ 20 sd outlier

# Hormone data: Each point is a subject mean, global-mean centered on zero, scaled (to sd = 1) within sex; data for both sexes

plot(df$ihh_C.cm[!is.na(df$ihh_C.cm)], ylab = "Cortisol")
plot(df$ihh_T.cm[!is.na(df$ihh_T.cm)], ylab = "Testosterone")
plot(df$ihh_P.cm[!is.na(df$ihh_P.cm)], ylab = "Progesterone")
plot(df$ihh_E.cm[!is.na(df$ihh_E.cm)], ylab = "Estradiol")

plot(df$ihh_C.cm[!is.na(df$ihh_C.cm) & df$ihh_C.cm < 19.5], ylab = "Cortisol") # cortisol without ~ 20 sd outlier
plot(df$ihh_T.cm[!is.na(df$ihh_T.cm) & df$ihh_T.cm < 9.5], ylab = "Testosterone") # testosterone without ~ 20 sd outlier

# Hormone data: Each point is a session, subject-mean centered on zero, scaled (to sd = 1) within sex; data for both sexes

par(mfrow=c(2,2))

plot(df$ihh_C.cwc[!is.na(df$ihh_C.cwc)])
plot(df$ihh_T.cwc[!is.na(df$ihh_T.cwc)])
plot(df$ihh_P.cwc[!is.na(df$ihh_P.cwc)])
plot(df$ihh_E.cwc[!is.na(df$ihh_E.cwc)])

###############################################################################################################################
##############################################                                     ############################################
##############################################    DUAL-HORMONE HYPOTHESIS TESTS    ############################################
##############################################                                     ############################################
###############################################################################################################################

# are SOI and/or SDI higher in higher T sessions in women?

soi_binary_T_f <- df_ihh.f %>%
  filter(n_sessions == 2, !is.na(soi_mean) & !is.na(binary_T)) %>%
  select(subID, binary_T, binary_C, soi_mean) %>%
  pivot_wider(id_cols = subID, names_from = c(binary_T), values_from = c(soi_mean, binary_C))

sdi_binary_T_f <- df_ihh.f %>%
  filter(n_sessions == 2, !is.na(sdi_mean) & !is.na(binary_T)) %>%
  select(subID, binary_T, binary_C, sdi_mean) %>%
  pivot_wider(id_cols = subID, names_from = c(binary_T), values_from = c(sdi_mean, binary_C))

t.test(soi_binary_T_f$soi_mean_low, soi_binary_T_f$soi_mean_high, paired = TRUE)
t.test(sdi_binary_T_f$sdi_mean_low, sdi_binary_T_f$sdi_mean_high, paired = TRUE)

# are SOI and/or SDI higher in higher T sessions in men?

soi_binary_T_m <- df_ihh.m %>%
  filter(n_sessions == 2, !is.na(soi_mean) & !is.na(binary_T)) %>%
  select(subID, binary_T, binary_C, soi_mean) %>%
  pivot_wider(id_cols = subID, names_from = c(binary_T), values_from = c(soi_mean, binary_C))

sdi_binary_T_m <- df_ihh.m %>%
  filter(n_sessions == 2, !is.na(sdi_mean) & !is.na(binary_T)) %>%
  select(subID, binary_T, binary_C, sdi_mean) %>%
  pivot_wider(id_cols = subID, names_from = c(binary_T), values_from = c(sdi_mean, binary_C))

t.test(soi_binary_T_m$soi_mean_low, soi_binary_T_m$soi_mean_high, paired = TRUE)
t.test(sdi_binary_T_m$sdi_mean_low, sdi_binary_T_m$sdi_mean_high, paired = TRUE)

# are SOI and/or SDI higher in higher T sessions that are also higher/lower C sessions (additive/antagonistic effects) in women?

soi_add_C_f <- soi_binary_T_f %>%
  filter(binary_C_low == "low")

soi_ant_C_f <- soi_binary_T_f %>%
  filter(binary_C_low == "high")

sdi_add_C_f <- sdi_binary_T_f %>%
  filter(binary_C_low == "low")

sdi_ant_C_f <- sdi_binary_T_f %>%
  filter(binary_C_low == "high")

t.test(soi_add_C_f$soi_mean_low, soi_add_C_f$soi_mean_high, paired = TRUE)
t.test(soi_ant_C_f$soi_mean_low, soi_ant_C_f$soi_mean_high, paired = TRUE)

t.test(sdi_add_C_f$sdi_mean_low, sdi_add_C_f$sdi_mean_high, paired = TRUE)
t.test(sdi_ant_C_f$sdi_mean_low, sdi_ant_C_f$sdi_mean_high, paired = TRUE)

# are SOI and/or SDI higher in higher T sessions that are also higher/lower C sessions (additive/antagonistic effects) in men?

soi_add_C_m <- soi_binary_T_m %>%
  filter(binary_C_low == "low")

soi_ant_C_m <- soi_binary_T_m %>%
  filter(binary_C_low == "high")

sdi_add_C_m <- sdi_binary_T_m %>%
  filter(binary_C_low == "low")

sdi_ant_C_m <- sdi_binary_T_m %>%
  filter(binary_C_low == "high")

t.test(soi_add_C_m$soi_mean_low, soi_add_C_m$soi_mean_high, paired = TRUE)
t.test(soi_ant_C_m$soi_mean_low, soi_ant_C_m$soi_mean_high, paired = TRUE)

t.test(sdi_add_C_m$sdi_mean_low, sdi_add_C_m$sdi_mean_high, paired = TRUE)
t.test(sdi_ant_C_m$sdi_mean_low, sdi_ant_C_m$sdi_mean_high, paired = TRUE)

### Preliminarry analysis: HORMONES UNCENTERED: Within- and between-individual effects ARE NOT isolated

# men

plot(df_ihh.m$ihh_T[df_ihh.m$ihh_T < 150], df_ihh.m$sdi_mean[df_ihh.m$ihh_T < 150])

df_ihh.m.bayes <- df_ihh.m[!is.na(df_ihh.m$ihh_T) & !is.na(df_ihh.m$sdi_mean) & df_ihh.m$ihh_T < 150,]

m.bayes <-rethinking::map(
  alist(
    sdi_mean ~ dnorm(mu, sigma),
    mu <- a + b*ihh_T,
    a ~ dnorm(2, 1),
    b ~ dnorm(0.1, 0.01),
    sigma ~ dunif(0, 75)
  ),
  data = df_ihh.m.bayes)

rethinking::precis(m.bayes)

m.samples <- extract.samples(m.bayes)

dens(head(m.samples$b))

dens(df$ihh_T)

T.seq <- seq( from=20 , to=120 , by=1 )

mu <- link(m.bayes, data=data.frame(ihh_T=T.seq))

plot( sdi_mean ~ ihh_T , df_ihh.m.bayes , type='n')

for (i in 1:500) {
  points( T.seq , mu[i,] , pch=16 , col=col.alpha(rangi2,0.1) )
}

mu.mean <- apply(mu, 2, mean)
mu.HPDI <- apply(mu, 2, HPDI, prob = 0.89)

plot( sdi_mean ~ ihh_T , df_ihh.m.bayes, col=col.alpha(rangi2,0.5) )
lines(T.seq, mu.mean)
shade(mu.HPDI, T.seq)


sim.SDI <- sim( m.bayes , data=list(ihh_T=T.seq) )

HDPI.SDI <- apply(sim.SDI, 2, HPDI, 0.89)

shade(HDPI.SDI, T.seq)

m.d <- lmer(sdi_mean ~ ihh_T*ihh_C + (1 | subID), data = df_ihh.m) # sdi-2

m.sol <- lmer(sdi_sol ~ ihh_T*ihh_C + (1 | subID), data = df_ihh.m) # sdi-2[solitary]

m.dyad <- lmer(sdi_dyad ~ ihh_T*ihh_C + (1 | subID), data = df_ihh.m) # sdi-2[dyadic]

m.s <- lmer(soi_mean ~ ihh_T*ihh_C + (1 | subID), data = df_ihh.m) # soiR

m.att <- lmer(soi_attit ~ ihh_T*ihh_C + (1 | subID), data = df_ihh.m) # soi-R[attitude]

m.beh <- lmer(soi_behav ~ ihh_T*ihh_C + (1 | subID), data = df_ihh.m) # soi-R[behavior]

m.des <- lmer(soi_des ~ ihh_T*ihh_C + (1 | subID), data = df_ihh.m) # soi-R[desire]

# women

w.d <- lmer(sdi_mean ~ ihh_T*ihh_C*HC + (1 | subID), data = df_ihh.f) # sdi-2

w.sol <- lmer(sdi_sol ~ ihh_T*ihh_C*HC + (1 | subID), data = df_ihh.f) # sdi-2[solitary]

w.dyad <- lmer(sdi_dyad ~ ihh_T*ihh_C*HC + (1 | subID), data = df_ihh.f) # sdi-2[dyadic]

w.s <- lmer(soi_mean ~ ihh_T*ihh_C*HC + (1 | subID), data = df_ihh.f) # soiR

w.att <- lmer(soi_attit ~ ihh_T*ihh_C*HC + (1 | subID), data = df_ihh.f) # soi-R[attitude]

w.beh <- lmer(soi_behav ~ ihh_T*ihh_C*HC + (1 | subID), data = df_ihh.f) # soi-R[behavior]

w.des <- lmer(soi_des ~ ihh_T*ihh_C*HC + (1 | subID), data = df_ihh.f) # soi-R[desire]


### HORMONES ARE CENTERED: Within(cwc)- and between(cm)-individual effects ARE isolated

# men

m.d_w <- lmer(sdi_mean ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + (1 | subID), data = df_ihh.m) # sdi-2

m.sol_w <- lmer(sdi_sol ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + (1 | subID), data = df_ihh.m) # sdi-2[solitary]

m.dyad_w <- lmer(sdi_dyad ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + (1 | subID), data = df_ihh.m) # sdi-2[dyadic]

m.s_w <- lmer(soi_mean ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + (1 | subID), data = df_ihh.m) # soiR

m.att_w <- lmer(soi_attit ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + (1 | subID), data = df_ihh.m) # soi-R[attitude]

m.beh_w <- lmer(soi_behav ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + (1 | subID), data = df_ihh.m) # soi-R[behavior]

m.des_w <- lmer(soi_des ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + (1 | subID), data = df_ihh.m) # soi-R[desire]

# women

w.d_w <- lmer(sdi_mean ~ ihh_T.cwc*ihh_C.cwc*HC + ihh_T.cm*ihh_C.cm*HC*HC + (1 | subID), data = df_ihh.f) # sdi-2

w.sol_w <- lmer(sdi_sol ~ ihh_T.cwc*ihh_C.cwc*HC + ihh_T.cm*ihh_C.cm*HC + (1 | subID), data = df_ihh.f) # sdi-2[solitary]

w.dyad_w <- lmer(sdi_dyad ~ ihh_T.cwc*ihh_C.cwc*HC + ihh_T.cm*ihh_C.cm*HC + (1 | subID), data = df_ihh.f) # sdi-2[dyadic]

w.s_w <- lmer(soi_mean ~ ihh_T.cwc*ihh_C.cwc*HC + ihh_T.cm*ihh_C.cm*HC + (1 | subID), data = df_ihh.f) # soiR

w.att_w <- lmer(soi_attit ~ ihh_T.cwc*ihh_C.cwc*HC + ihh_T.cm*ihh_C.cm*HC + (1 | subID), data = df_ihh.f) # soi-R[attitude]

w.beh_w <- lmer(soi_behav ~ ihh_T.cwc*ihh_C.cwc*HC + ihh_T.cm*ihh_C.cm*HC + (1 | subID), data = df_ihh.f) # soi-R[behavior]

w.des_w <- lmer(soi_des ~ ihh_T.cwc*ihh_C.cwc*HC + ihh_T.cm*ihh_C.cm*HC + (1 | subID), data = df_ihh.f) # soi-R[desire]

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

# Figure

figure1 <- df %>%
  filter(study == "ihh", n_sessions == 2) %>%
  pivot_wider(id_cols = subID, names_from = session, values_from = c(ihh_C, ihh_T, sdi_sol)) %>%
  mutate(C_diff = ihh_C_2 - ihh_C_1,
         T_diff = ihh_T_2 - ihh_T_1 ,
         sdi_sol_diff = sdi_sol_2 - sdi_sol_1) %>%
  select(C_diff, T_diff, sdi_sol_diff) %>%
#  filter(T_diff < 100) # view figure without one extreme T difference outlier

ggplot(figure1, aes(T_diff, sdi_sol_diff, size= C_diff)) + geom_point()


## Robustness tests

# women taking hormonal contraception

df_ihh.hc <- df_ihh.f[df_ihh.f$HC == 1,]

hc.d_w <- lmer(sdi_mean ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + (1 | subID), data = df_ihh.hc) # sdi-2

hc.sol_w <- lmer(sdi_sol ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + (1 | subID), data = df_ihh.hc) # sdi-2[solitary]

hc.dyad_w <- lmer(sdi_dyad ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + (1 | subID), data = df_ihh.hc) # sdi-2[dyadic]

hc.s_w <- lmer(soi_mean ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + (1 | subID), data = df_ihh.hc) # soiR

hc.att_w <- lmer(soi_attit ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + (1 | subID), data = df_ihh.hc) # soi-R[attitude]

hc.beh_w <- lmer(soi_behav ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + (1 | subID), data = df_ihh.hc) # soi-R[behavior]

hc.des_w <- lmer(soi_des ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + (1 | subID), data = df_ihh.hc) # soi-R[desire]

# women NOT taking hormonal contraception

df_ihh.nc <- df_ihh.f[df_ihh.f$HC == 0,]

nc.d_w <- lmer(sdi_mean ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + (1 | subID), data = df_ihh.nc) # sdi-2

nc.sol_w <- lmer(sdi_sol ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + (1 | subID), data = df_ihh.nc) # sdi-2[solitary]

nc.dyad_w <- lmer(sdi_dyad ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + (1 | subID), data = df_ihh.nc) # sdi-2[dyadic]

nc.s_w <- lmer(soi_mean ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + (1 | subID), data = df_ihh.nc) # soiR

nc.att_w <- lmer(soi_attit ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + (1 | subID), data = df_ihh.nc) # soi-R[attitude]

nc.beh_w <- lmer(soi_behav ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + (1 | subID), data = df_ihh.nc) # soi-R[behavior]

nc.des_w <- lmer(soi_des ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + (1 | subID), data = df_ihh.nc) # soi-R[desire]

class(hc.d_w) <- "lmerMod"
class(hc.sol_w) <- "lmerMod"
class(hc.dyad_w) <- "lmerMod"
class(hc.s_w) <- "lmerMod"
class(hc.att_w) <- "lmerMod"
class(hc.beh_w) <- "lmerMod"
class(hc.des_w) <- "lmerMod"

class(nc.d_w) <- "lmerMod"
class(nc.sol_w) <- "lmerMod"
class(nc.dyad_w) <- "lmerMod"
class(nc.s_w) <- "lmerMod"
class(nc.att_w) <- "lmerMod"
class(nc.beh_w) <- "lmerMod"
class(nc.des_w) <- "lmerMod"

# FOR WITHIN SUBJECT MODELS, WOMEN taking hormonal contraception

stargazer(hc.d_w, hc.sol_w, hc.dyad_w, hc.s_w, hc.att_w, hc.beh_w, hc.des_w, single.row = F, align = F, intercept.bottom = T, no.space = T, 
          dep.var.caption = "", type = "text", out = "hc_within.html", report = "vc*p", 
          star.char = c("*", "**"), digits = 3, star.cutoffs = c(0.05, 0.01), 
          notes = "*p<0.05;**p<0.01", notes.append = F,
          column.labels = c("Total desire", "Solitary desire", "Dyadic Desire", "Uncommited Sex Total",
                            "Uncommited Sex Attitude", "Uncommited Sex Behavior", "Uncommited Sex Desire" ), omit.stat = "all")

# FOR WITHIN SUBJECT MODELS, WOMEN not taking hormonal contraception
stargazer(nc.d_w, nc.sol_w, nc.dyad_w, nc.s_w, nc.att_w, nc.beh_w, nc.des_w, single.row = F, align = F, intercept.bottom = T, no.space = T, 
          dep.var.caption = "", type = "text", out = "nc_within.html", report = "vc*p", 
          star.char = c("*", "**"), digits = 3, star.cutoffs = c(0.05, 0.01), 
          notes = "*p<0.05;**p<0.01", notes.append = F,
          column.labels = c("Total desire", "Solitary desire", "Dyadic Desire", "Uncommited Sex Total",
                            "Uncommited Sex Attitude", "Uncommited Sex Behavior", "Uncommited Sex Desire" ), omit.stat = "all")


# do effects of C and T remain after accounting for effects of ovarian hormones on SDI and/or SOI in non-contracepting women?

ov.d_w <- lmer(sdi_mean ~ ihh_P.cwc*ihh_E.cwc + ihh_P.cm*ihh_E.cm + ihh_C.cwc*ihh_T.cwc +
                 ihh_C.cm*ihh_T.cm + (1 | subID), data = df_ihh.nc) # sdi-2

ov.sol_w <- lmer(sdi_sol ~ ihh_P.cwc*ihh_E.cwc + ihh_P.cm*ihh_E.cm + ihh_C.cwc*ihh_T.cwc +
                   ihh_C.cm*ihh_T.cm + (1 | subID), data = df_ihh.nc) # sdi-2[solitary]

ov.dyad_w <- lmer(sdi_dyad ~ ihh_P.cwc*ihh_E.cwc + ihh_P.cm*ihh_E.cm + ihh_C.cwc*ihh_T.cwc +
                    ihh_C.cm*ihh_T.cm + (1 | subID), data = df_ihh.nc) # sdi-2[dyadic]

ov.s_w <- lmer(soi_mean ~ ihh_P.cwc*ihh_E.cwc + ihh_P.cm*ihh_E.cm + ihh_C.cwc*ihh_T.cwc +
                 ihh_C.cm*ihh_T.cm + (1 | subID), data = df_ihh.nc) # soiR

ov.att_w <- lmer(soi_attit ~ ihh_P.cwc*ihh_E.cwc + ihh_P.cm*ihh_E.cm + ihh_C.cwc*ihh_T.cwc +
                   ihh_C.cm*ihh_T.cm + (1 | subID), data = df_ihh.nc) # soi-R[attitude]

ov.beh_w <- lmer(soi_behav ~ ihh_P.cwc*ihh_E.cwc + ihh_P.cm*ihh_E.cm + ihh_C.cwc*ihh_T.cwc +
                   ihh_C.cm*ihh_T.cm + (1 | subID), data = df_ihh.nc) # soi-R[behavior]

ov.des_w <- lmer(soi_des ~ ihh_P.cwc*ihh_E.cwc + ihh_P.cm*ihh_E.cm + ihh_C.cwc*ihh_T.cwc +
                   ihh_C.cm*ihh_T.cm + (1 | subID), data = df_ihh.nc) # soi-R[desire]

class(ov.d_w) <- "lmerMod"
class(ov.sol_w) <- "lmerMod"
class(ov.dyad_w) <- "lmerMod"
class(ov.s_w) <- "lmerMod"
class(ov.att_w) <- "lmerMod"
class(ov.beh_w) <- "lmerMod"
class(ov.des_w) <- "lmerMod"

# FOR WITHIN SUBJECT MODELS, WOMEN not taking hormonal contraception
stargazer(ov.d_w, ov.sol_w, ov.dyad_w, ov.s_w, ov.att_w, ov.beh_w, ov.des_w, single.row = F, align = F, intercept.bottom = T, no.space = T, 
          dep.var.caption = "", type = "text", out = "ov_within.html", report = "vc*p", 
          star.char = c("*", "**"), digits = 3, star.cutoffs = c(0.05, 0.01), 
          notes = "*p<0.05;**p<0.01", notes.append = F,
          column.labels = c("Total desire", "Solitary desire", "Dyadic Desire", "Uncommited Sex Total",
                            "Uncommited Sex Attitude", "Uncommited Sex Behavior", "Uncommited Sex Desire" ), omit.stat = "all")

lapply(c("tidyverse", "dplyr", "lmerTest", "psych", "sjPlot", "sjmisc", "sjlabelled",
"data.table", "reshape2", "ggcorrplot", "stargazer", "gridExtra", "htmlTable"), FUN = citation)

# Main analyses with predictor hormones log-transformed before centering and scaling

df <- df %>%
  group_by(subID) %>%
  mutate(ihh_C.log = log(ihh_C),
         ihh_T.log = log(ihh_T),
         ihh_P.log = log(ihh_P),
         ihh_E.log = log(ihh_E)) %>%
  mutate(ihh_C.log.cm = mean(ihh_C.log, na.rm =TRUE),
         ihh_T.log.cm = mean(ihh_T.log, na.rm =TRUE),
         ihh_P.log.cm = mean(ihh_P.log, na.rm =TRUE),
         ihh_E.log.cm = mean(ihh_E.log, na.rm =TRUE),
         ihh_C.log.cwc = ihh_C.log - mean(ihh_C.log, na.rm =TRUE),
         ihh_T.log.cwc = ihh_T.log - mean(ihh_T.log, na.rm =TRUE),
         ihh_P.log.cwc = ihh_P.log - mean(ihh_P.log, na.rm =TRUE),
         ihh_E.log.cwc = ihh_E.log - mean(ihh_E.log, na.rm =TRUE)) %>%
  ungroup() %>%
  group_by(sex) %>%
  mutate(ihh_C.log.cwc = ihh_C.log.cwc / sd(ihh_C.log.cwc, na.rm = TRUE),
         ihh_T.log.cwc = ihh_T.log.cwc / sd(ihh_T.log.cwc, na.rm = TRUE),
         ihh_P.log.cwc = ihh_P.log.cwc / sd(ihh_P.log.cwc, na.rm = TRUE),
         ihh_E.log.cwc = ihh_E.log.cwc / sd(ihh_E.log.cwc, na.rm = TRUE),
         ihh_C.log.cm = (ihh_C.log.cm - mean(ihh_C.log.cm, na.rm = TRUE)) / sd(ihh_C.log.cm, na.rm = TRUE),
         ihh_T.log.cm = (ihh_T.log.cm - mean(ihh_T.log.cm, na.rm = TRUE)) / sd(ihh_T.log.cm, na.rm = TRUE),
         ihh_P.log.cm = (ihh_P.log.cm - mean(ihh_P.log.cm, na.rm = TRUE)) / sd(ihh_P.log.cm, na.rm = TRUE),
         ihh_E.log.cm = (ihh_E.log.cm - mean(ihh_E.log.cm, na.rm = TRUE)) / sd(ihh_E.log.cm, na.rm = TRUE))


# Log-transformed raw hormone data plots; data for both sexes

par(mfrow=c(2,2))

plot(df$ihh_C.log[!is.na(df$ihh_C.log)], ylab = "Log Cortisol")
plot(df$ihh_T.log[!is.na(df$ihh_T.log)], ylab = "Log Testosterone")
plot(df$ihh_P.log[!is.na(df$ihh_P.log)], ylab = "Log Progesterone")
plot(df$ihh_E.log[!is.na(df$ihh_E.log)], ylab = "Log Estradiol")

# Log-transformed hormone data: Each point is a subject mean, global-mean centered on zero, scaled (to sd = 1) within sex; data for both sexes

plot(df$ihh_C.log.cm[!is.na(df$ihh_C.log.cm)], ylab = "Log Cortisol")
plot(df$ihh_T.log.cm[!is.na(df$ihh_T.log.cm)], ylab = "Log Testosterone")
plot(df$ihh_P.log.cm[!is.na(df$ihh_P.log.cm)], ylab = "Log Progesterone")
plot(df$ihh_E.log.cm[!is.na(df$ihh_E.log.cm)], ylab = "Log Estradiol")

# Log-transformed hormone data: Each point is a session, subject-mean centered on zero, scaled (to sd = 1) w/in sex; data for both sexes

plot(df$ihh_C.log.cwc[!is.na(df$ihh_C.log.cwc)], ylab = "Log Cortisol")
plot(df$ihh_T.log.cwc[!is.na(df$ihh_T.log.cwc)], ylab = "Log Testosterone")
plot(df$ihh_P.log.cwc[!is.na(df$ihh_P.log.cwc)], ylab = "Log Progesterone")
plot(df$ihh_E.log.cwc[!is.na(df$ihh_E.log.cwc)], ylab = "Log Estradiol")

# Models

### HORMONES ARE CENTERED: Within(cwc)- and between(cm)-individual effects ARE isolated; hormones are log-transformed

# men

m.d.log <- lmer(sdi_mean ~ ihh_T.log.cwc*ihh_C.log.cwc + ihh_T.log.cm*ihh_C.log.cm + (1 | subID), data = df_ihh.m) # sdi-2

m.sol.log <- lmer(sdi_sol ~ ihh_T.log.cwc*ihh_C.log.cwc + ihh_T.log.cm*ihh_C.log.cm + (1 | subID), data = df_ihh.m) # sdi-2[solitary]

m.dyad.log <- lmer(sdi_dyad ~ ihh_T.log.cwc*ihh_C.log.cwc + ihh_T.log.cm*ihh_C.log.cm + (1 | subID), data = df_ihh.m) # sdi-2[dyadic]

m.s.log <- lmer(soi_mean ~ ihh_T.log.cwc*ihh_C.log.cwc + ihh_T.log.cm*ihh_C.log.cm + (1 | subID), data = df_ihh.m) # soiR

m.att.log <- lmer(soi_attit ~ ihh_T.log.cwc*ihh_C.log.cwc + ihh_T.log.cm*ihh_C.log.cm + (1 | subID), data = df_ihh.m) # soi-R[attitude]

m.beh.log <- lmer(soi_behav ~ ihh_T.log.cwc*ihh_C.log.cwc + ihh_T.log.cm*ihh_C.log.cm + (1 | subID), data = df_ihh.m) # soi-R[behavior]

m.des.log <- lmer(soi_des ~ ihh_T.log.cwc*ihh_C.log.cwc + ihh_T.log.cm*ihh_C.log.cm + (1 | subID), data = df_ihh.m) # soi-R[desire]

# women

w.d.log <- lmer(sdi_mean ~ ihh_T.log.cwc*ihh_C.log.cwc*HC + ihh_T.log.cm*ihh_C.log.cm*HC + (1 | subID), data = df_ihh.f) # sdi-2

w.sol.log <- lmer(sdi_sol ~ ihh_T.log.cwc*ihh_C.log.cwc*HC + ihh_T.log.cm*ihh_C.log.cm*HC + (1 | subID), data = df_ihh.f) # sdi solitary

w.dyad.log <- lmer(sdi_dyad ~ ihh_T.log.cwc*ihh_C.log.cwc*HC + ihh_T.log.cm*ihh_C.log.cm*HC + (1 | subID), data = df_ihh.f) # sdi dyadic

w.s.log <- lmer(soi_mean ~ ihh_T.log.cwc*ihh_C.log.cwc*HC + ihh_T.log.cm*ihh_C.log.cm*HC + (1 | subID), data = df_ihh.f) # soiR

w.att.log <- lmer(soi_attit ~ ihh_T.log.cwc*ihh_C.log.cwc*HC + ihh_T.log.cm*ihh_C.log.cm*HC + (1 | subID), data = df_ihh.f) # attitude

w.beh.log <- lmer(soi_behav ~ ihh_T.log.cwc*ihh_C.log.cwc*HC + ihh_T.log.cm*ihh_C.log.cm*HC + (1 | subID), data = df_ihh.f) # behavior

w.des.log <- lmer(soi_des ~ ihh_T.log.cwc*ihh_C.log.cwc*HC + ihh_T.log.cm*ihh_C.log.cm*HC + (1 | subID), data = df_ihh.f) # soi desire

# Have to do this to use stargazer

class(w.d.log) <- "lmerMod"
class(w.sol.log) <- "lmerMod"
class(w.dyad.log) <- "lmerMod"
class(w.s.log) <- "lmerMod"
class(w.att.log) <- "lmerMod"
class(w.beh.log) <- "lmerMod"
class(w.des.log) <- "lmerMod"

class(m.d.log) <- "lmerMod"
class(m.sol.log) <- "lmerMod"
class(m.dyad.log) <- "lmerMod"
class(m.s.log) <- "lmerMod"
class(m.att.log) <- "lmerMod"
class(m.beh.log) <- "lmerMod"
class(m.des.log) <- "lmerMod"

# FOR WITHIN SUBJECT MODELS, WOMEN
stargazer(w.d.log, w.sol.log, w.dyad.log, w.s.log, w.att.log, w.beh.log, w.des.log, single.row = F, align = F, intercept.bottom = T,
          no.space = T, dep.var.caption = "", type = "text", out = "women_within.html", report = "vc*p", 
          star.char = c("*", "**"), digits = 3, star.cutoffs = c(0.05, 0.01), 
          notes = "*p<0.05;**p<0.01", notes.append = F,
          column.labels = c("Total desire", "Solitary desire", "Dyadic Desire", "Uncommited Sex Total",
                            "Uncommited Sex Attitude", "Uncommited Sex Behavior", "Uncommited Sex Desire" ), omit.stat = "all")

# FOR WITHIN SUBJECT MODELS, MEN
stargazer(m.d.log, m.sol.log, m.dyad.log, m.s.log, m.att.log, m.beh.log, m.des.log, single.row = F, align = F, intercept.bottom = T,
          no.space = T, dep.var.caption = "", type = "text", out = "menithin.html", report = "vc*p", 
          star.char = c("*", "**"), digits = 3, star.cutoffs = c(0.05, 0.01), 
          notes = "*p<0.05;**p<0.01", notes.append = F,
          column.labels = c("Total desire", "Solitary desire", "Dyadic Desire", "Uncommited Sex Total",
                            "Uncommited Sex Attitude", "Uncommited Sex Behavior", "Uncommited Sex Desire" ), omit.stat = "all")

