library(pacman); pacman::p_load(lmerTest, psych, sjPlot, sjmisc, sjlabelled,
                                data.table, reshape2, ggcorrplot, stargazer, gridExtra, htmlTable, tidyverse, dplyr)

###############################################################################################################################
################################################                               ################################################
################################################        DATA PREPERATION       ################################################
################################################                               ################################################
###############################################################################################################################

# read in data

df <- as.data.frame(read.csv("df.csv"))

# change variable types

df <- df %>%
  filter(study == "ihh") %>%
  #filter(age <= 45) %>% # uncomment this line to exclude subject age 45+
  select(-c(study, ihh_race, msu_ethnicity, ihh_ethnicity, AM_T, PM_T, Kinsey_Attr:cgnq24b4_score)) %>%
  mutate_at(vars(Date), as.Date) %>%
  mutate_at(vars(ID:IUD), as.factor) %>%
  mutate_at(vars(time, age:sdi_mean), as.character) %>%
  mutate_at(vars(time, age:sdi_mean), as.numeric)

# calculate mean soi and soi psychology scores

df <- df %>%
  mutate(soi_mean = rowMeans(select(., c(soit1:soit9)), na.rm = TRUE),
         soi_psych = rowMeans(select(., c(soi_des, soi_attit)), na.rm = TRUE))

# Calculate days between samples

dfDate <- df %>%
  filter(n_sessions == 2, subID !="CON858") %>%
  select(subID, session, Date) %>%
  pivot_wider(names_from = session, values_from = Date) %>%
  `colnames<-`(c("subID", "one", "two"))

dateDiff = vector()
for (i in 1:nrow(dfDate)) {
  dateDiff = append(dateDiff, abs(julian(as.Date(dfDate$one[i]), as.Date(dfDate$two)[i])[1]))
}

dfDate <- bind_cols("subID" = dfDate["subID"], "diff" = dateDiff)
df <- df %>% left_join(dfDate)
df %>% filter(session == 2) %>% select(diff) %>% describe()

# calculate subject means and subject mean-center hormone variables; center subject-means; scale subject means and subject-centered data within each sex

df <- df %>%
  group_by(subID) %>%
  mutate(ihh_C.cm = mean(ihh_C, na.rm =TRUE),
         ihh_T.cm = mean(ihh_T, na.rm =TRUE),
         ihh_P.cm = mean(ihh_P, na.rm =TRUE),
         ihh_E.cm = mean(ihh_P, na.rm =TRUE),
         soi.cm = mean(soi_mean, na.rm =TRUE),
         soiB.cm = mean(soi_behav, na.rm =TRUE),
         soiA.cm = mean(soi_attit, na.rm =TRUE),
         soiD.cm = mean(soi_des, na.rm =TRUE),
         sdi.cm = mean(sdi_mean, na.rm =TRUE),
         sdiD.cm = mean(sdi_dyad, na.rm =TRUE),
         sdiS.cm = mean(sdi_sol, na.rm =TRUE),
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

df_2session <- df %>%
  filter(n_sessions == 2) %>%
  select(-c(ID, group, sibID, IUD, ihh_P, ihh_E, soi_mean, soi_psych, soit1:sdi14, ihh_C.cm:binary_C)) %>%
  group_by(subID) %>%
  mutate(C_diff = ihh_C[session == 1] - ihh_C[session == 2],
         T_diff = ihh_T[session == 1] - ihh_T[session == 2],
         sdi_diff = sdi_mean[session == 1] - sdi_mean[session == 2],
         sol_diff = sdi_sol[session == 1] - sdi_sol[session == 2],
         dyad_diff = sdi_dyad[session == 1] - sdi_dyad[session == 2])

df_2session.men <- df_2session %>%
  filter(sex == "m") %>%
  mutate()
df_2session.women <- df_2session %>% filter(sex == "f")

# log-transformed

df.log <- df %>%
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


# break out datasets by sex/study/hormonal contraception use

df.m <- df %>% filter(sex == "m")
df.f <- df %>% filter(sex == "f")
df.hc <- df.f[df.f$HC == 1,]
df.nc <- df.f[df.f$HC == 0,]
df.log.m <- df.log %>% filter(sex == "m")
df.log.f <- df.log %>% filter(sex == "f")

# create function for easier viewing of model output

table.model <- function(model) {
  model$coefficients %>%
    as.data.frame() %>%
    rownames_to_column(var = "Effect") %>%
    mutate_if(is.numeric, round, 3) %>%
    rename(p = `Pr(>|t|)`)
}

###############################################################################################################################
################################################                               ################################################
################################################     DESCRIPTIVE STATISTICS    ################################################
################################################                               ################################################
###############################################################################################################################

# inter-item reliability tests (cronbach's alpha for SOI and SDI questions)

df_sdi <- df %>%
  dplyr::select(soit1:soit9)

df_sdi <- df %>%
  dplyr::select(sdi1:sdi14)

library(ltm)
cronbach.alpha(df_soi, na.rm = T)
cronbach.alpha(df_sdi, na.rm = T)
detach(package:ltm)
remove(df_sdi)
remove(df_soi)

df %>% filter(sex == "f", session == 1) %>% count()
df %>% filter(sex == "m", session == 1) %>% count()

# differences between subjects who attended 1 or 2 sessions

compareNumberSessions <- df %>%
  dplyr::select(subID, n_sessions, sex, session, ihh_C.cm, ihh_T.cm, soi.cm,  soiA.cm, soiB.cm, soiD.cm, sdi.cm, sdiD.cm, sdiS.cm) %>%
  filter(session == 1) %>%
  pivot_wider(names_from = n_sessions, values_from = ihh_C.cm:sdiS.cm)

compareNumberSessionsM <- compareNumberSessions %>% filter(sex == "m")
compareNumberSessionsF <- compareNumberSessions %>% filter(sex == "f")

for (i in seq(4,20,2)) {
  print(t.test(compareNumberSessionsM[i], compareNumberSessionsM[i + 1])[3])
}

for (i in seq(4,20,2)) {
  print(t.test(compareNumberSessionsF[i], compareNumberSessionsF[i + 1])[3])
}

# n of subjects who attended 1 or 2 sessions

df %>% filter(sex == "f", n_sessions == 1) %>% count()
df %>% filter(sex == "m", n_sessions == 1) %>% count()


df %>% filter(sex == "f", n_sessions == 2) %>% count()/2
df %>% filter(sex == "m", n_sessions == 2) %>% count()/2

# Table 1 - age, total n, and n of subs with at least one set of T, C, and either SDI or SOI samples

options(table_counter = TRUE)

df %>%
  filter(session == 1) %>%
  group_by(sex, HC) %>%
  summarise("Age (mean)" = round(mean(age), 2), "SD" = round(sd(age), 2), n = n()) %>%
  filter(!is.na(HC) | sex == "m") %>%
  ungroup() %>%
  dplyr::select(-c(HC, sex)) %>%
  htmlTable(align = "c",
            caption = "Age statistics broken out by subgroup",
            rnames = c("Women (NC)", "Women (HC)", "Men"))

# Table 1 cont'd - SOI-R and SDI-2 Statistics

table3Rows = c("SOI-R ", "SOI: Attitudes", "SOI: Behavior", "SOI: Desire", "SDI-2",
               "SDI: Solitary", "SDI: Dyadic")
table3Cols = c("n", "Mean", "SD")

table3 <- df.m %>%
  group_by(subID) %>%
  dplyr::select(soi_mean, soi_attit, soi_behav, soi_des, sdi_mean, sdi_sol, sdi_dyad, ihh_T, ihh_C) %>%
  filter(!is.na(ihh_T), !is.na(ihh_C)) %>%
  summarise_all(funs(mean)) %>%
  describe %>%
  bind_cols(as.numeric(rep("",10)), df.nc %>%
              group_by(subID) %>%
              dplyr::select(soi_mean, soi_attit, soi_behav, soi_des, sdi_mean, sdi_sol, sdi_dyad, ihh_T, ihh_C) %>%
              filter(!is.na(ihh_T), !is.na(ihh_C)) %>%
              summarise_all(funs(mean)) %>%
              describe %>%
              bind_cols(as.numeric(rep("",10)), df.hc %>%
                          group_by(subID) %>%
                          dplyr::select(soi_mean, soi_attit, soi_behav, soi_des, sdi_mean, sdi_sol, sdi_dyad, ihh_T, ihh_C) %>%
                          filter(!is.na(ihh_T), !is.na(ihh_C)) %>%
                          summarise_all(funs(mean)) %>%
                          describe))

htmlTable(round(table3[2:8,c(2:4, 14, 16:18, 28, 30:32)], 2),
          align = "c",
          caption = "Responses to SOI-R, SDI-2, and subscales",
          rnames = table3Rows, 
          header = c(table3Cols, "", table3Cols, "", table3Cols),
          col.rgroup = c("none", "#F7F7F7"),
          cgroup = c("Men", "", "Women (NC)", "", "Women (HC)"), n.cgroup = c(3, 1, 3, 1, 3))

# Table 1 continued - Raw hormone average levels

table4Rows = c("T", "C", "P", "E")
table4Cols = c("n", "Mean", "SD")

table4 <- df.m %>%
  group_by(subID) %>%
  dplyr::select(ihh_T, ihh_C, ihh_P, ihh_E) %>%
  summarise_all(funs(mean)) %>%
  describe %>%
  bind_cols(as.numeric(rep("",5)), df.hc %>%
              group_by(subID) %>%
              dplyr::select(ihh_T, ihh_C, ihh_P, ihh_E) %>%
              summarise_all(funs(mean)) %>%
              describe, as.numeric(rep("",5)), df.nc %>%
              group_by(subID) %>%
              dplyr::select(ihh_T, ihh_C, ihh_P, ihh_E) %>%
              summarise_all(funs(mean)) %>%
              describe)
table4 <- round(table4, 2)
table4[4:5, 2:4] <- "NA"

htmlTable(table4[2:5,c(14, 2:4, 14, 16:18, 28, 30:32)],
          align = "c",
          caption = "Raw hormone descriptive statistics",
          rnames = table4Rows, 
          header = c("", table4Cols, "", table4Cols, "", table4Cols),
          col.rgroup = c("none", "#F7F7F7"),
          cgroup = c("", "Men", "", "Women (HC)", "", "Women (NC)"), n.cgroup = c(1, 3, 1, 3, 1, 3))

# ESM Table 1 - n of subs with one or twos set of T, C, and either SDI or SOI samples

df %>%
  filter(session == 1) %>%
  dplyr::select(n_sessions, sex, HC, ihh_T, ihh_C, sdi_mean, soi_mean) %>%
  group_by(n_sessions, sex, HC) %>%
  summarise(n = n(),
            "n (SOI)" = sum(!is.na(ihh_T) & !is.na(ihh_C) & !is.na(soi_mean)),
            "n (SDI)" = sum(!is.na(ihh_T) & !is.na(ihh_C) & !is.na(sdi_mean))) %>%
  filter(!is.na(HC) | sex == "m") %>%
  ungroup() %>%
  dplyr::select(-c(HC, sex, n_sessions)) %>%
  htmlTable(align = "c",
            col.rgroup = c("none", "#F7F7F7"),
            tspanner = c("Single session","2 Sessions"), n.tspanner = c(3,3),
            caption = "Sample sizes broken out by number of sessions",
            rnames = c("Women (NC)", "Women (HC)", "Men",
                       "Women (NC)", "Women (HC)", "Men"))

# sex and hormonal contraception group differences on SOI-R, SDI-2, and subject-mean hormones

ihh.sex <- df %>%
  filter(!is.na(soi_mean), !is.nan(soi_mean), !is.na(sex)) %>%
  dplyr::select(subID, HC, sex, soi_mean, soi_attit, soi_behav, soi_des, sdi_mean, sdi_dyad, sdi_sol, ihh_C, ihh_T, ihh_P, ihh_E) %>%
  group_by(subID) %>%
  mutate(soi_mean = mean(soi_mean),
         sdi_mean = mean(sdi_mean),
         C_mean = mean(soi_mean),
         T_mean = mean(sdi_mean))

ihh.women <- ihh.sex %>% filter(sex == "f")

for (i in 4:12){
  p <- t.test(unlist(ihh.sex[i])~ unlist(ihh.sex$sex))[3][[1]]
  print(paste(names(ihh.women[i]), round(p,4)))
}

for (i in 4:14){
  p <- t.test(unlist(ihh.women[i])~ unlist(ihh.women$HC))[3][[1]]
  print(paste(names(ihh.women[i]), round(p,4)))
}

remove(ihh.sex)
remove(ihh.women)

###############################################################################################################################
##############################################                                     ############################################
##############################################          Data visualization         ############################################
##############################################                                     ############################################
###############################################################################################################################


# Raw hormone data plots; data for both sexes

par(mfrow=c(3,2))

plot(df$ihh_C[!is.na(df$ihh_C)], ylab = "Cortisol")
plot(df$ihh_T[!is.na(df$ihh_T)], ylab = "Testosterone")
plot(df$ihh_P[!is.na(df$ihh_P)], ylab = "Progesterone")
plot(df$ihh_E[!is.na(df$ihh_E)], ylab = "Estradiol")

plot(df$ihh_C[!is.na(df$ihh_C) & df$ihh_C < 40], ylab = "Cortisol - 1 Outlier") # cortisol w/out 20 sd outlier
plot(df$ihh_T[!is.na(df$ihh_T) & df$ihh_T < 145], ylab = "Testosterone - 3 Outliers") # testosterone w/out 20 sd outlier

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

# Supplementary Figures 1-4 : Correlations of scales and their subscales; each scale/sex combination presented separately

pairs.panels(
  df.m %>% select(sdi_mean, sdi_sol, sdi_dyad), stars= TRUE, method = "pearson", hist.col = "gray", cex.labels=1.9, cex.cor = 2, ellipses = FALSE)

pairs.panels(
  df.f %>% select(sdi_mean, sdi_sol, sdi_dyad), stars= TRUE, method = "pearson", hist.col = "gray", cex.labels=1.9, cex.cor = 2, ellipses = FALSE)

pairs.panels(
  df.m %>% select(soi_mean, soi_attit, soi_behav, soi_des), stars= TRUE, method = "pearson", hist.col = "gray", cex.labels=1.9, cex.cor = 2, ellipses = FALSE)

pairs.panels(
  df.f %>% select(soi_mean, soi_attit, soi_behav, soi_des), stars= TRUE, method = "pearson", hist.col = "gray", cex.labels=1.9, cex.cor = 2, ellipses = FALSE)

# Supplementary Figures 5-6: Correlations of ALL scales and subscales, seperately by sex

pairs.panels(
  df.m %>% select(soi_mean, soi_attit, soi_behav, soi_des, sdi_mean, sdi_sol, sdi_dyad), stars= TRUE, method = "pearson", hist.col = "gray", cex.labels=1.9, cex.cor = 2, ellipses = FALSE)

pairs.panels(
  df.f %>% select(soi_mean, soi_attit, soi_behav, soi_des, sdi_mean, sdi_sol, sdi_dyad), stars= TRUE, method = "pearson", hist.col = "gray", cex.labels=1.9, cex.cor = 2, ellipses = FALSE)

# Supplementary Figure 7: Log-transformed raw hormone data plots; data for both sexes

par(mfrow=c(2,2))

plot(df.log$ihh_C.log[!is.na(df.log$ihh_C.log)], ylab = "Log Cortisol")
plot(df.log$ihh_T.log[!is.na(df.log$ihh_T.log)], ylab = "Log Testosterone")
plot(df.log$ihh_P.log[!is.na(df.log$ihh_P.log)], ylab = "Log Progesterone")
plot(df.log$ihh_E.log[!is.na(df.log$ihh_E.log)], ylab = "Log Estradiol")

# Supplementary Figure 8: Log-transformed hormone data: Each point is a subject mean, global-mean centered on zero, scaled (to sd = 1) within sex; data for both sexes

plot(df.log$ihh_C.log.cm[!is.na(df.log$ihh_C.log.cm)], ylab = "Log Cortisol")
plot(df.log$ihh_T.log.cm[!is.na(df.log$ihh_T.log.cm)], ylab = "Log Testosterone")
plot(df.log$ihh_P.log.cm[!is.na(df.log$ihh_P.log.cm)], ylab = "Log Progesterone")
plot(df.log$ihh_E.log.cm[!is.na(df.log$ihh_E.log.cm)], ylab = "Log Estradiol")

# Supplementary Figure 9: Log-transformed hormone data: Each point is a session, subject-mean centered on zero, scaled (to sd = 1) w/in sex; data for both sexes

plot(df.log$ihh_C.log.cwc[!is.na(df.log$ihh_C.log.cwc)], ylab = "Log Cortisol")
plot(df.log$ihh_T.log.cwc[!is.na(df.log$ihh_T.log.cwc)], ylab = "Log Testosterone")
plot(df.log$ihh_P.log.cwc[!is.na(df.log$ihh_P.log.cwc)], ylab = "Log Progesterone")
plot(df.log$ihh_E.log.cwc[!is.na(df.log$ihh_E.log.cwc)], ylab = "Log Estradiol")


women <- df %>% filter(sex == 'f', !is.na(ihh_C), !is.na(ihh_T), ihh_T.cwc != 0)

quantiles_women <-df.log %>% filter(sex == 'f', !is.na(ihh_C), !is.na(ihh_T), ihh_T.cwc != 0) %>%
  as.data.frame() %>%
  select(ihh_C) %>%
  quantile(unlist(.), probs=c(0.25, 0.5, 0.75, 1.0)) %>%
  findInterval(women$ihh_C,.,rightmost.closed=TRUE) %>%
  as.factor()

women_plot <- ggplot(women, aes(ihh_T.cwc, sdi_dyad, col = quantiles_women)) +
  geom_smooth(method = lm, aes(group = quantiles_women), se = FALSE) +
  geom_point(method = lm, aes(group = quantiles_women)) +
  scale_colour_manual(name="Cortisol", 
                      values=c("#BDBDBD", "#848484", "#585858","#1C1C1C"),
                      breaks=c(0, 1,2,3),
                      labels = c("Bottom 25%", "25-50%", "50-75%", "Top 25%")) +
  xlim(-5,5) + ylim(0,10) + theme(legend.position = "none", axis.title.x=element_text(size=15), axis.title.y=element_text(size=15)) +
  xlab("Subject mean-centered testosterone") + ylab("Dyadic desire (Women)") +
  #theme(plot.margin=margin(30,30,30,30), axis.title=element_text(size=17,face="bold"), legend.position = c(0.8, 0.8)) +
  annotate(geom="text", x=-3.75, y=1.3, label="A", fontface=2, size=13)

men <- df %>% filter(sex == 'm', !is.na(ihh_C), !is.na(ihh_T), ihh_T.cwc != 0)

quantiles_men <-df %>% filter(sex == 'm', !is.na(ihh_C), !is.na(ihh_T), ihh_T.cwc != 0) %>%
  select(ihh_C) %>%
  as.data.frame() %>%
  quantile(unlist(.), probs=c(0.25, 0.5, 0.75, 1.0)) %>%
  findInterval(men$ihh_C,.,rightmost.closed=TRUE) %>%
  as.factor()

men_plot <- ggplot(men, aes(ihh_T.cwc, sdi_sol, col = quantiles_men)) +
  geom_smooth(method = lm, aes(group = quantiles_men), se = FALSE) +
  geom_point(method = lm, aes(group = quantiles_men)) +
  scale_colour_manual(name="Cortisol", 
                      values=c("#BDBDBD", "#848484", "#585858","#1C1C1C"),
                      breaks=c(0, 1,2,3),
                      labels = c("Bottom 25%", "25-50%", "50-75%", "Top 25%")) +
  xlab("Subject mean-centered testosterone") + ylab("Solitary desire (Men)") +
  xlim(-5,5) + ylim(0,10) + theme(legend.position= "none", axis.title.x=element_text(size=15), axis.title.y=element_text(size=15)) +
  #theme(plot.margin=margin(30,30,30,30), axis.title=element_text(size=17,face="bold"), legend.position = c(0.8, 0.8)) +
  annotate(geom="text", x=-3.75, y=1.3, label="B", fontface=2, size=13)

grid.arrange(women_plot, men_plot, nrow = 1)


###############################################################################################################################
##############################################                                     ############################################
##############################################        Preliminary Analyses         ############################################
##############################################                                     ############################################
###############################################################################################################################

# are SOI and/or SDI higher in higher T sessions in women?

soi_binary_T_f <- df.f %>%
  filter(n_sessions == 2, !is.na(soi_mean) & !is.na(binary_T)) %>%
  select(subID, binary_T, binary_C, soi_mean) %>%
  pivot_wider(id_cols = subID, names_from = c(binary_T), values_from = c(soi_mean, binary_C))

sdi_binary_T_f <- df.f %>%
  filter(n_sessions == 2, !is.na(sdi_mean) & !is.na(binary_T)) %>%
  select(subID, binary_T, binary_C, sdi_mean) %>%
  pivot_wider(id_cols = subID, names_from = c(binary_T), values_from = c(sdi_mean, binary_C))

t.test(soi_binary_T_f$soi_mean_low, soi_binary_T_f$soi_mean_high, paired = TRUE)
t.test(sdi_binary_T_f$sdi_mean_low, sdi_binary_T_f$sdi_mean_high, paired = TRUE)

# are SOI and/or SDI higher in higher T sessions in men?

soi_binary_T_m <- df.m %>%
  filter(n_sessions == 2, !is.na(soi_mean) & !is.na(binary_T)) %>%
  select(subID, binary_T, binary_C, soi_mean) %>%
  pivot_wider(id_cols = subID, names_from = c(binary_T), values_from = c(soi_mean, binary_C))

sdi_binary_T_m <- df.m %>%
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

### Preliminary analysis: HORMONES UNCENTERED: Within- and between-individual effects ARE NOT isolated

# men

m.d <- table.model(summary(lmer(sdi_mean ~ ihh_T*ihh_C + (1 | subID), data = df.m))) # sdi-2

m.sol <- table.model(summary(lmer(sdi_sol ~ ihh_T*ihh_C + (1 | subID), data = df.m))) # sdi-2[solitary]

m.dyad <- table.model(summary(lmer(sdi_dyad ~ ihh_T*ihh_C + (1 | subID), data = df.m))) # sdi-2[dyadic]

m.s <- table.model(summary(lmer(soi_mean ~ ihh_T*ihh_C + (1 | subID), data = df.m))) # soiR

m.att <- table.model(summary(lmer(soi_attit ~ ihh_T*ihh_C + (1 | subID), data = df.m))) # soi-R[attitude]

m.beh <- table.model(summary(lmer(soi_behav ~ ihh_T*ihh_C + (1 | subID), data = df.m))) # soi-R[behavior]

m.des <- table.model(summary(lmer(soi_des ~ ihh_T*ihh_C + (1 | subID), data = df.m))) # soi-R[desire]


# women

w.d <- table.model(summary(lmer(sdi_mean ~ ihh_T*ihh_C*HC + (1 | subID), data = df.f))) # sdi-2

w.sol <- table.model(summary(lmer(sdi_sol ~ ihh_T*ihh_C*HC + (1 | subID), data = df.f))) # sdi-2[solitary]

w.dyad <- table.model(summary(lmer(sdi_dyad ~ ihh_T*ihh_C*HC + (1 | subID), data = df.f))) # sdi-2[dyadic]

w.s <- table.model(summary(lmer(soi_mean ~ ihh_T*ihh_C*HC + (1 | subID), data = df.f))) # soiR

w.att <- table.model(summary(lmer(soi_attit ~ ihh_T*ihh_C*HC + (1 | subID), data = df.f))) # soi-R[attitude]

w.beh <- table.model(summary(lmer(soi_behav ~ ihh_T*ihh_C*HC + (1 | subID), data = df.f))) # soi-R[behavior]

w.des <- table.model(summary(lmer(soi_des ~ ihh_T*ihh_C*HC + (1 | subID), data = df.f))) # soi-R[desire]


###############################################################################################################################
##############################################                                     ############################################
##############################################    DUAL-HORMONE HYPOTHESIS TESTS    ############################################
##############################################                                     ############################################
###############################################################################################################################

### Main Analyses: hormones centered; Within(cwc)- and between(cm)-individual effects ARE isolated

# men

m.d_w <- table.model(summary(lmer(sdi_mean ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + (1 | subID), data = df.m))) # sdi-2

m.sol_w <- table.model(summary(lmer(sdi_sol ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + (1 | subID), data = df.m))) # sdi-2[solitary]

m.dyad_w <- table.model(summary(lmer(sdi_dyad ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + (1 | subID), data = df.m))) # sdi-2[dyadic]

m.s_w <- table.model(summary(lmer(soi_mean ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + (1 | subID), data = df.m))) # soiR

m.att_w <- table.model(summary(lmer(soi_attit ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + (1 | subID), data = df.m))) # soi-R[attitude]

m.beh_w <- table.model(summary(lmer(soi_behav ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + (1 | subID), data = df.m))) # soi-R[behavior]

m.des_w <- table.model(summary(lmer(soi_des ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + (1 | subID), data = df.m))) # soi-R[desire]

# women

w.d_w <- table.model(summary(lmer(sdi_mean ~ ihh_T.cwc*ihh_C.cwc*HC + ihh_T.cm*ihh_C.cm*HC + (1 | subID), data = df.f))) # sdi-2

w.sol_w <- table.model(summary(lmer(sdi_sol ~ ihh_T.cwc*ihh_C.cwc*HC + ihh_T.cm*ihh_C.cm*HC + (1 | subID), data = df.f))) # sdi-2[solitary]

w.dyad_w <- table.model(summary(lmer(sdi_dyad ~ ihh_T.cwc*ihh_C.cwc*HC + ihh_T.cm*ihh_C.cm*HC + (1 | subID), data = df.f))) # sdi-2[dyadic]

w.s_w <- table.model(summary(lmer(soi_mean ~ ihh_T.cwc*ihh_C.cwc*HC + ihh_T.cm*ihh_C.cm*HC + (1 | subID), data = df.f))) # soiR

w.att_w <- table.model(summary(lmer(soi_attit ~ ihh_T.cwc*ihh_C.cwc*HC + ihh_T.cm*ihh_C.cm*HC + (1 | subID), data = df.f))) # soi-R[attitude]

w.beh_w <- table.model(summary(lmer(soi_behav ~ ihh_T.cwc*ihh_C.cwc*HC + ihh_T.cm*ihh_C.cm*HC + (1 | subID), data = df.f))) # soi-R[behavior]

w.des_w <- table.model(summary(lmer(soi_des ~ ihh_T.cwc*ihh_C.cwc*HC + ihh_T.cm*ihh_C.cm*HC + (1 | subID), data = df.f))) # soi-R[desire]

##### Table 2.1 - Sociosexual orientation (women NOT broken out by hormonal contraception use)

table2_f.1cols <- c("Estimate", "df", "t", "p")
table2_f.1rows <- c("Intercept", "T (within)", "C (within)", "Contraception use", "T (between)",
                "C (between)", "T x C (within)", "T (within) x Contraception use",
                "C (within) x Contraception use", "T x C (between)", "C (between) x Contraception use",
                "C (between) x Contraception use", "T x C (within) x Contraception use", "T x C (between) x Contraception use")
table2_f.1 <- round(bind_rows(bind_cols(w.s_w[-c(1,3)], w.att_w[-c(1,3)], w.beh_w[-c(1,3)], w.des_w[-c(1,3)])),3)
                          

htmlTable(table2_f.1,
          align = "c",
          caption = "SOI-R Model results",
          rnames = c(table2_f.1rows, table2_f.1rows, table2_f.1rows),
          header = c(table2_f.1cols, table2_f.1cols, table2_f.1cols, table2_f.1cols),
          col.rgroup = c("none", "#F7F7F7"),
          cgroup = list(c("SOI-R"), c("General", "Attitudes", "Behavior", "Desire")),
          n.cgroup = list(c(4), c(4, 4, 4, 4)),
          tspanner = c("Women"), n.tspanner = c(14),
          css.tspanner = "padding-top: 0.5em; font-weight:bold",
          css.cell = "padding-right: 0.5em")

##### Table 2.2 - Sexual desire (women NOT broken out by hormonal contraception use)

table2_f.2cols <- c("Estimate", "df", "t", "p")
table2_f.2rows <- c("Intercept", "T (within)", "C (within)", "Contraception use", "T (between)",
                    "C (between)", "T x C (within)", "T (within) x Contraception use",
                    "C (within) x Contraception use", "T x C (between)", "C (between) x Contraception use",
                    "C (between) x Contraception use", "T x C (within) x Contraception use", "T x C (between) x Contraception use")
table2_f.2 <- round(bind_rows(bind_cols(w.d_w[-c(1,3)], w.sol_w[-c(1,3)], w.dyad_w[-c(1,3)])),2)


htmlTable(table2_f.2,
          align = "c",
          caption = "SDI-2 Model results",
          rnames = c(table2_f.2rows, table2_f.2rows, table2_f.2rows),
          header = c(table2_f.2cols, table2_f.2cols, table2_f.2cols),
          col.rgroup = c("none", "#F7F7F7"),
          cgroup = list(c("SDI-2"), c("General", "Solitary", "Dyadic")),
          n.cgroup = list(c(12), c(4, 4, 4)),
          tspanner = c("Women"), n.tspanner = c(14),
          css.tspanner = "padding-top: 0.5em; font-weight:bold",
          css.cell = "padding-right: 0.5em")

##### Table 2.3 - Sociosexual orientation (men)

table2_m.1cols <- c("Estimate", "df", "t", "p")
table2_m.1rows <- c("Intercept", "T (within)", "C (within)", "T (between)",
                "C (between)", "T x C (within)", "T x C (between)")
table2_m.1 <- round(bind_rows(bind_cols(m.s_w[-c(1,3)], m.att_w[-c(1,3)], m.beh_w[-c(1,3)], m.des_w[-c(1,3)])),2)

htmlTable(table2_m.1,
          align = "c",
          caption = "SOI-R Model results",
          rnames = c(table2_m.1rows, table2_m.1rows, table2_m.1rows),
          header = c(table2_m.1cols, table2_m.1cols, table2_m.1cols, table2_m.1cols),
          col.rgroup = c("none", "#F7F7F7"),
          cgroup = list(c("SOI-R"), c("General", "Attitudes", "Behavior", "Desire")),
          n.cgroup = list(c(16), c(4, 4, 4,4)),
          tspanner = c("Men"), n.tspanner = c(7),
          css.tspanner = "padding-top: 0.5em; font-weight:bold",
          css.cell = "padding-right: 0.5em")

##### Table 2.4 - Sexual desire (men)

table2_m.2cols <- c("Estimate", "df", "t", "p")
table2_m.2rows <- c("Intercept", "T (within)", "C (within)", "T (between)",
                "C (between)", "T x C (within)", "T x C (between)")
table2_m.2 <- round(bind_rows(bind_cols(m.d_w[-c(1,3)], m.sol_w[-c(1,3)], m.dyad_w[-c(1,3)])),2)

htmlTable(table2_m.2,
          align = "c",
          caption = "SDI-2 Model results",
          rnames = c(table2_m.2rows, table2_m.2rows, table2_m.2rows),
          header = c(table2_m.2cols, table2_m.2cols, table2_m.2cols),
          col.rgroup = c("none", "#F7F7F7"),
          cgroup = list(c("SDI-2"), c("General", "Solitary", "Dyadic")),
          n.cgroup = list(c(12), c(4, 4, 4)),
          tspanner = c("Men"), n.tspanner = c(7),
          css.tspanner = "padding-top: 0.5em; font-weight:bold",
          css.cell = "padding-right: 0.5em")


###############################################################################################################################
##############################################                                     ############################################
##############################################         Robustness Analyses         ############################################
##############################################                                     ############################################
###############################################################################################################################

# women taking hormonal contraception

hc.d_w <- table.model(summary(lmer(sdi_mean ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + (1 | subID), data = df.hc))) # sdi-2

hc.sol_w <- table.model(summary(lmer(sdi_sol ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + (1 | subID), data = df.hc))) # sdi-2[solitary]

hc.dyad_w <- table.model(summary(lmer(sdi_dyad ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + (1 | subID), data = df.hc))) # sdi-2[dyadic]

hc.s_w <- table.model(summary(lmer(soi_mean ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + (1 | subID), data = df.hc))) # soiR

hc.att_w <- table.model(summary(lmer(soi_attit ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + (1 | subID), data = df.hc))) # soi-R[attitude]

hc.beh_w <- table.model(summary(lmer(soi_behav ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + (1 | subID), data = df.hc))) # soi-R[behavior]

hc.des_w <- table.model(summary(lmer(soi_des ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + (1 | subID), data = df.hc))) # soi-R[desire]

# women NOT taking hormonal contraception

nc.d_w <- table.model(summary(lmer(sdi_mean ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + (1 | subID), data = df.nc))) # sdi-2

nc.sol_w <- table.model(summary(lmer(sdi_sol ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + (1 | subID), data = df.nc))) # sdi-2[solitary]

nc.dyad_w <- table.model(summary(lmer(sdi_dyad ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + (1 | subID), data = df.nc))) # sdi-2[dyadic]

nc.s_w <- table.model(summary(lmer(soi_mean ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + (1 | subID), data = df.nc))) # soiR

nc.att_w <- table.model(summary(lmer(soi_attit ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + (1 | subID), data = df.nc))) # soi-R[attitude]

nc.beh_w <- table.model(summary(lmer(soi_behav ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + (1 | subID), data = df.nc))) # soi-R[behavior]

nc.des_w <- table.model(summary(lmer(soi_des ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + (1 | subID), data = df.nc))) # soi-R[desire]

##### ESM Table 1 - Sociosexual orientation (women broken out by hormonal contraception use)

table3cols <- c("Estimate", "df", "t", "p")
table3rows <- c("Intercept", "T (within)", "C (within)", "T (between)",
                "C (between)", "T x C (within)", "T x C (between)")
table3 <- round(bind_rows(bind_cols(m.s_w[-c(1,3)], m.att_w[-c(1,3)], m.beh_w[-c(1,3)], m.des_w[-c(1,3)]),
                          bind_cols(nc.s_w[-c(1,3)], nc.att_w[-c(1,3)], nc.beh_w[-c(1,3)], nc.des_w[-c(1,3)]),
                          bind_cols(hc.s_w[-c(1,3)], hc.att_w[-c(1,3)], hc.beh_w[-c(1,3)], hc.des_w[-c(1,3)])),2)

htmlTable(table3,
          align = "c",
          caption = "SOI-R Model results",
          rnames = c(table3rows, table3rows, table3rows),
          header = c(table3cols, table3cols, table3cols, table3cols),
          col.rgroup = c("none", "#F7F7F7"),
          cgroup = list(c("SOI-R"), c("General", "Attitudes", "Behavior", "Desire")),
          n.cgroup = list(c(4), c(4, 4, 4, 4)),
          tspanner = c("Men","Women not usng hormonal contraception", "Women using hormonal contraception"), n.tspanner = c(7,7),
          css.tspanner = "padding-top: 0.5em; font-weight:bold",
          css.cell = "padding-right: 0.5em")

##### ESM Table 2 - Sociosexuality (women broken out by hormonal contraception use)

table4cols <- c("Estimate", "df", "t", "p")
table4rows <- c("Intercept", "T (within)", "C (within)", "T (between)",
                "C (between)", "T x C (within)", "T x C (between)")
table4 <- round(bind_rows(bind_cols(m.d_w[-c(1,3)], m.sol_w[-c(1,3)], m.dyad_w[-c(1,3)]),
                          bind_cols(nc.d_w[-c(1,3)], nc.sol_w[-c(1,3)], nc.dyad_w[-c(1,3)]),
                          bind_cols(hc.d_w[-c(1,3)], hc.sol_w[-c(1,3)], hc.dyad_w[-c(1,3)])),2)

htmlTable(table4,
          align = "c",
          caption = "SDI-2 Model results",
          rnames = c(table4rows, table4rows, table4rows),
          header = c(table4cols, table4cols, table4cols),
          col.rgroup = c("none", "#F7F7F7"),
          cgroup = list(c("SDI-2"), c("General", "Solitary", "Dyadic")),
          n.cgroup = list(c(3), c(4, 4, 4)),
          tspanner = c("Men","Women not usng hormonal contraception", "Women using hormonal contraception"), n.tspanner = c(7,7),
          css.tspanner = "padding-top: 0.5em; font-weight:bold",
          css.cell = "padding-right: 0.5em")

#####

# do effects of C and T remain after accounting for effects of ovarian hormones on SDI and/or SOI in non-contracepting women?

ov.d_w <- table.model(summary(lmer(sdi_mean ~ ihh_P.cwc*ihh_E.cwc + ihh_P.cm*ihh_E.cm + ihh_C.cwc*ihh_T.cwc +
                 ihh_C.cm*ihh_T.cm + (1 | subID), data = df.nc))) # sdi-2

ov.sol_w <- table.model(summary(lmer(sdi_sol ~ ihh_P.cwc*ihh_E.cwc + ihh_P.cm*ihh_E.cm + ihh_C.cwc*ihh_T.cwc +
                   ihh_C.cm*ihh_T.cm + (1 | subID), data = df.nc))) # sdi-2[solitary]

ov.dyad_w <- table.model(summary(lmer(sdi_dyad ~ ihh_P.cwc*ihh_E.cwc + ihh_P.cm*ihh_E.cm + ihh_C.cwc*ihh_T.cwc +
                    ihh_C.cm*ihh_T.cm + (1 | subID), data = df.nc))) # sdi-2[dyadic]

ov.s_w <- table.model(summary(lmer(soi_mean ~ ihh_P.cwc*ihh_E.cwc + ihh_P.cm*ihh_E.cm + ihh_C.cwc*ihh_T.cwc +
                 ihh_C.cm*ihh_T.cm + (1 | subID), data = df.nc))) # soiR

ov.att_w <- table.model(summary(lmer(soi_attit ~ ihh_P.cwc*ihh_E.cwc + ihh_P.cm*ihh_E.cm + ihh_C.cwc*ihh_T.cwc +
                   ihh_C.cm*ihh_T.cm + (1 | subID), data = df.nc))) # soi-R[attitude]

ov.beh_w <- table.model(summary(lmer(soi_behav ~ ihh_P.cwc*ihh_E.cwc + ihh_P.cm*ihh_E.cm + ihh_C.cwc*ihh_T.cwc +
                   ihh_C.cm*ihh_T.cm + (1 | subID), data = df.nc))) # soi-R[behavior]

ov.des_w <- table.model(summary(lmer(soi_des ~ ihh_P.cwc*ihh_E.cwc + ihh_P.cm*ihh_E.cm + ihh_C.cwc*ihh_T.cwc +
                   ihh_C.cm*ihh_T.cm + (1 | subID), data = df.nc))) # soi-R[desire]

# With predictor hormones log-transformed before centering and scaling

# men

m.d.log <- table.model(summary(lmer(sdi_mean ~ ihh_T.log.cwc*ihh_C.log.cwc +
                                      ihh_T.log.cm*ihh_C.log.cm + (1 | subID), data = df.log.m))) # sdi-2

m.sol.log <- table.model(summary(lmer(sdi_sol ~ ihh_T.log.cwc*ihh_C.log.cwc +
                                        ihh_T.log.cm*ihh_C.log.cm + (1 | subID), data = df.log.m))) # sdi-2[solitary]

m.dyad.log <- table.model(summary(lmer(sdi_dyad ~ ihh_T.log.cwc*ihh_C.log.cwc +
                                         ihh_T.log.cm*ihh_C.log.cm + (1 | subID), data = df.log.m))) # sdi-2[dyadic]

m.s.log <- table.model(summary(lmer(soi_mean ~ ihh_T.log.cwc*ihh_C.log.cwc +
                                      ihh_T.log.cm*ihh_C.log.cm + (1 | subID), data = df.log.m))) # soiR

m.att.log <- table.model(summary(lmer(soi_attit ~ ihh_T.log.cwc*ihh_C.log.cwc +
                                        ihh_T.log.cm*ihh_C.log.cm + (1 | subID), data = df.log.m))) # soi-R[attitude]

m.beh.log <- table.model(summary(lmer(soi_behav ~ ihh_T.log.cwc*ihh_C.log.cwc +
                                        ihh_T.log.cm*ihh_C.log.cm + (1 | subID), data = df.log.m))) # soi-R[behavior]

m.des.log <- table.model(summary(lmer(soi_des ~ ihh_T.log.cwc*ihh_C.log.cwc +
                                        ihh_T.log.cm*ihh_C.log.cm + (1 | subID), data = df.log.m))) # soi-R[desire]

# women

w.d.log <- table.model(summary(lmer(sdi_mean ~ ihh_T.log.cwc*ihh_C.log.cwc*HC +
                                      ihh_T.log.cm*ihh_C.log.cm*HC + (1 | subID), data = df.log.f))) # sdi-2

w.sol.log <- table.model(summary(lmer(sdi_sol ~ ihh_T.log.cwc*ihh_C.log.cwc*HC +
                                        ihh_T.log.cm*ihh_C.log.cm*HC + (1 | subID), data = df.log.f))) # sdi solitary

w.dyad.log <- table.model(summary(lmer(sdi_dyad ~ ihh_T.log.cwc*ihh_C.log.cwc*HC + 
                                         ihh_T.log.cm*ihh_C.log.cm*HC + (1 | subID), data = df.log.f))) # sdi dyadic

w.s.log <- table.model(summary(lmer(soi_mean ~ ihh_T.log.cwc*ihh_C.log.cwc*HC +
                                      ihh_T.log.cm*ihh_C.log.cm*HC + (1 | subID), data = df.log.f))) # soiR

w.att.log <- table.model(summary(lmer(soi_attit ~ ihh_T.log.cwc*ihh_C.log.cwc*HC +
                                        ihh_T.log.cm*ihh_C.log.cm*HC + (1 | subID), data = df.log.f))) # attitude

w.beh.log <- table.model(summary(lmer(soi_behav ~ ihh_T.log.cwc*ihh_C.log.cwc*HC +
                                        ihh_T.log.cm*ihh_C.log.cm*HC + (1 | subID), data = df.log.f))) # behavior

w.des.log <- table.model(summary(lmer(soi_des ~ ihh_T.log.cwc*ihh_C.log.cwc*HC +
                                        ihh_T.log.cm*ihh_C.log.cm*HC + (1 | subID), data = df.log.f))) # soi desire
# age

# men

m.d_w_age <- table.model(summary(lmer(sdi_mean ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + age + (1 | subID), data = df.m))) # sdi-2

m.sol_w_age <- table.model(summary(lmer(sdi_sol ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + age + (1 | subID), data = df.m))) # sdi-2[solitary]

m.dyad_w_age <- table.model(summary(lmer(sdi_dyad ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + age + (1 | subID), data = df.m))) # sdi-2[dyadic]

m.s_w_age <- table.model(summary(lmer(soi_mean ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + age + (1 | subID), data = df.m))) # soiR

m.att_w_age <- table.model(summary(lmer(soi_attit ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + age + (1 | subID), data = df.m))) # soi-R[attitude]

m.beh_w_age <- table.model(summary(lmer(soi_behav ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + age + (1 | subID), data = df.m))) # soi-R[behavior]

m.des_w_age <- table.model(summary(lmer(soi_des ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + age + (1 | subID), data = df.m))) # soi-R[desire]

# women

w.d_w_age <- table.model(summary(lmer(sdi_mean ~ ihh_T.cwc*ihh_C.cwc*HC + ihh_T.cm*ihh_C.cm*HC + age + (1 | subID), data = df.f))) # sdi-2

w.sol_w_age <- table.model(summary(lmer(sdi_sol ~ ihh_T.cwc*ihh_C.cwc*HC + ihh_T.cm*ihh_C.cm*HC + age + (1 | subID), data = df.f))) # sdi-2[solitary]

w.dyad_w_age <- table.model(summary(lmer(sdi_dyad ~ ihh_T.cwc*ihh_C.cwc*HC + ihh_T.cm*ihh_C.cm*HC + age + (1 | subID), data = df.f))) # sdi-2[dyadic]

w.s_w_age <- table.model(summary(lmer(soi_mean ~ ihh_T.cwc*ihh_C.cwc*HC + ihh_T.cm*ihh_C.cm*HC + age + (1 | subID), data = df.f))) # soiR

w.att_w_age <- table.model(summary(lmer(soi_attit ~ ihh_T.cwc*ihh_C.cwc*HC + ihh_T.cm*ihh_C.cm*HC + age + (1 | subID), data = df.f))) # soi-R[attitude]

w.beh_w_age <- table.model(summary(lmer(soi_behav ~ ihh_T.cwc*ihh_C.cwc*HC + ihh_T.cm*ihh_C.cm*HC + age + (1 | subID), data = df.f))) # soi-R[behavior]

w.des_w_age <- table.model(summary(lmer(soi_des ~ ihh_T.cwc*ihh_C.cwc*HC + ihh_T.cm*ihh_C.cm*HC + age + (1 | subID), data = df.f))) # soi-R[desire]

# time

# men

m.d_w_time <- table.model(summary(lmer(sdi_mean ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + time + (1 | subID), data = df.m))) # sdi-2

m.sol_w_time <- table.model(summary(lmer(sdi_sol ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + time + (1 | subID), data = df.m))) # sdi-2[solitary]

m.dyad_w_time <- table.model(summary(lmer(sdi_dyad ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + time + (1 | subID), data = df.m))) # sdi-2[dyadic]

m.s_w_time <- table.model(summary(lmer(soi_mean ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + time + (1 | subID), data = df.m))) # soiR

m.att_w_time <- table.model(summary(lmer(soi_attit ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + time + (1 | subID), data = df.m))) # soi-R[attitude]

m.beh_w_time <- table.model(summary(lmer(soi_behav ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + time + (1 | subID), data = df.m))) # soi-R[behavior]

m.des_w_time <- table.model(summary(lmer(soi_des ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + time + (1 | subID), data = df.m))) # soi-R[desire]

# women

w.d_w_time <- table.model(summary(lmer(sdi_mean ~ ihh_T.cwc*ihh_C.cwc*HC + ihh_T.cm*ihh_C.cm*HC + time + (1 | subID), data = df.f))) # sdi-2

w.sol_w_time <- table.model(summary(lmer(sdi_sol ~ ihh_T.cwc*ihh_C.cwc*HC + ihh_T.cm*ihh_C.cm*HC + time + (1 | subID), data = df.f))) # sdi-2[solitary]

w.dyad_w_time <- table.model(summary(lmer(sdi_dyad ~ ihh_T.cwc*ihh_C.cwc*HC + ihh_T.cm*ihh_C.cm*HC + time + (1 | subID), data = df.f))) # sdi-2[dyadic]

w.s_w_time <- table.model(summary(lmer(soi_mean ~ ihh_T.cwc*ihh_C.cwc*HC + ihh_T.cm*ihh_C.cm*HC + time + (1 | subID), data = df.f))) # soiR

w.att_w_time <- table.model(summary(lmer(soi_attit ~ ihh_T.cwc*ihh_C.cwc*HC + ihh_T.cm*ihh_C.cm*HC + time + (1 | subID), data = df.f))) # soi-R[attitude]

w.beh_w_time <- table.model(summary(lmer(soi_behav ~ ihh_T.cwc*ihh_C.cwc*HC + ihh_T.cm*ihh_C.cm*HC + time + (1 | subID), data = df.f))) # soi-R[behavior]

w.des_w_time <- table.model(summary(lmer(soi_des ~ ihh_T.cwc*ihh_C.cwc*HC + ihh_T.cm*ihh_C.cm*HC + time + (1 | subID), data = df.f))) # soi-R[desire]

# package citations

lapply(c("tidyverse", "dplyr", "lmerTest", "psych", "sjPlot", "sjmisc", "sjlabelled",
         "data.table", "reshape2", "ggcorrplot", "stargazer", "gridExtra", "htmlTable"), FUN = citation)

