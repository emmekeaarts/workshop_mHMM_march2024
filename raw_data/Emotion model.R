
library(mHMMbayes)
emotion_data <- readRDS("./raw_data/data_Rowland2020.rds")
emotion_data$beep2 <- (emotion_data$dayno - 1) * 6 + emotion_data$beep

ggplot_emotion <- rbind(data.frame(ID = emotion_data$subj_id, 
                                 beep2 = emotion_data$beep2, 
                                 outcome = emotion_data$happy,
                                 emotion = "happy"),
                      data.frame(ID = emotion_data$subj_id, 
                                 beep2 = emotion_data$beep2, 
                                 outcome = emotion_data$excited,
                                 emotion = "excited"),
                      data.frame(ID = emotion_data$subj_id, 
                                 beep2 = emotion_data$beep2, 
                                 outcome = emotion_data$relaxed,
                                 emotion = "relaxed"),
                      data.frame(ID = emotion_data$subj_id, 
                                 beep2 = emotion_data$beep2, 
                                 outcome = emotion_data$satisfied,
                                 emotion = "satisfied"),
                      data.frame(ID = emotion_data$subj_id, 
                                 beep2 = emotion_data$beep2, 
                                 outcome = emotion_data$angry,
                                 emotion = "angry"),
                      data.frame(ID = emotion_data$subj_id, 
                                 beep2 = emotion_data$beep2, 
                                 outcome = emotion_data$anxious,
                                 emotion = "anxious"),
                      data.frame(ID = emotion_data$subj_id, 
                                 beep2 = emotion_data$beep2, 
                                 outcome = emotion_data$depressed,
                                 emotion = "depressed"),
                      data.frame(ID = emotion_data$subj_id, 
                                 beep2 = emotion_data$beep2, 
                                 outcome = emotion_data$sad,
                                 emotion = "sad"))

ggplot_emotion$emotion <- as.factor(ggplot_emotion$emotion)

ggplot(ggplot_emotion[ggplot_emotion$ID %in% c(1:5),]) +
  geom_line( mapping = aes(x = beep2, y = outcome, group = emotion, color = emotion)) +
  facet_wrap(vars(ID), ncol = 1)


emotion_mHMM <- data.frame(subj_ID = emotion_data$subj_id,
                           happy = emotion_data$happy,
                           excited = emotion_data$excited,
                           relaxed = emotion_data$relaxed,
                           satisfied = emotion_data$satisfied,
                           angry = emotion_data$angry,
                           anxious = emotion_data$anxious,
                           depressed = emotion_data$depressed,
                           sad = emotion_data$sad)


#####################
## 2 state model ####
#####################
m <- 2
n_dep <- 8

start_gamma <- matrix(c(0.7, 0.3, 
                        0.3, 0.7), byrow = TRUE, ncol = m, nrow = m)

start_emiss <- list(matrix(c(70, 15,                                     # happy
                             30, 15), byrow = TRUE, ncol = 2, nrow = m), 
                    matrix(c(70, 15,                                     # excited
                             30, 15), byrow = TRUE, ncol = 2, nrow = m), 
                    matrix(c(70, 15,                                     # relaxed
                             30, 15), byrow = TRUE, ncol = 2, nrow = m), 
                    matrix(c(70, 15,                                     # satisfied
                             30, 15), byrow = TRUE, ncol = 2, nrow = m), 
                    matrix(c(15, 10,                                     # angry
                             40, 15), byrow = TRUE, ncol = 2, nrow = m), 
                    matrix(c(15, 10,                                     # anxious
                             40, 15), byrow = TRUE, ncol = 2, nrow = m), 
                    matrix(c(15, 10,                                     # depressed
                             40, 15), byrow = TRUE, ncol = 2, nrow = m), 
                    matrix(c(15, 10,                                     # sad
                             40, 15), byrow = TRUE, ncol = 2, nrow = m))

emotion_prior_emiss <- prior_emiss_cont(
  gen = list(m = m, n_dep = n_dep),
  emiss_mu0 = list(matrix(c(70, 30), nrow = 1),  # happy
                   matrix(c(70, 30), nrow = 1),  # excited
                   matrix(c(70, 30), nrow = 1),  # relaxed
                   matrix(c(70, 30), nrow = 1),  # satisfied
                   matrix(c(15, 40), nrow = 1),  # angry
                   matrix(c(15, 40), nrow = 1),  # anxious
                   matrix(c(15, 40), nrow = 1),  # depressed
                   matrix(c(15, 40), nrow = 1)), # sad 
  emiss_K0 = rep(list(1), n_dep),
  emiss_V = rep(list(rep(5^2, m)), n_dep),
  emiss_nu = rep(list(1), n_dep),
  emiss_a0 = rep(list(rep(1.5, m)), n_dep),
  emiss_b0 = rep(list(rep(20, m)), n_dep),
)


out_2st_emotion <- mHMM(s_data = emotion_mHMM,
                                  data_distr = "continuous",
                                  gen = list(m = m, n_dep = n_dep),
                                  start_val = c(list(start_gamma), start_emiss),
                                  emiss_hyp_prior = emotion_prior_emiss,
                                  mcmc = list(J = 500, burn_in = 200))

saveRDS(out_2st_emotion, file = "out_2st_emotion.rds")

out_2st_emotion
summary(out_2st_emotion)
gamma_group <- obtain_gamma(out_2st_emotion) 

plot(gamma_group)

emiss_group <- obtain_emiss(out_2st_emotion) 
plot(emiss_group)

plot(out_2st_emotion, component = "emiss", dep = 8)

#####################
## 3 state model ####
#####################
## general model properties
m <- 3
n_dep <- 8

## starting values for gamma
start_gamma_3st <- matrix(c(0.8, 0.1, 0.1,
                            0.1, 0.8, 0.1,
                            0.1, 0.1, 0.8), byrow = TRUE, ncol = m, nrow = m)

## starting values for the emission distribution 
start_emiss_3st <- list(matrix(c(80, 10,                                     # happy
                                 50, 10,
                                 20, 10), byrow = TRUE, ncol = 2, nrow = m), 
                        matrix(c(80, 10,                                     # excited
                                 50, 10,
                                 20, 10), byrow = TRUE, ncol = 2, nrow = m), 
                        matrix(c(80, 10,                                     # relaxed
                                 50, 10,
                                 20, 10), byrow = TRUE, ncol = 2, nrow = m), 
                        matrix(c(80, 10,                                     # satisfied
                                 50, 10,
                                 20, 10), byrow = TRUE, ncol = 2, nrow = m), 
                        matrix(c(10,  5,                                     # angry
                                 30, 10,
                                 50, 15), byrow = TRUE, ncol = 2, nrow = m), 
                        matrix(c(10,  5,                                     # anxious
                                 30, 10,
                                 50, 15), byrow = TRUE, ncol = 2, nrow = m), 
                        matrix(c(10,  5,                                     # depressed
                                 30, 10,
                                 50, 15), byrow = TRUE, ncol = 2, nrow = m), 
                        matrix(c(10,  5,                                     # sad
                                 30, 10,
                                 50, 15), byrow = TRUE, ncol = 2, nrow = m))

## specifying weakly informative prior for continuous emission distributions
emotion_prior_emiss_3st <- prior_emiss_cont(
  gen = list(m = m, n_dep = n_dep),
  emiss_mu0 = list(matrix(c(80, 50, 20), nrow = 1),  # happy
                   matrix(c(80, 50, 20), nrow = 1),  # excited
                   matrix(c(80, 50, 20), nrow = 1),  # relaxed
                   matrix(c(80, 50, 20), nrow = 1),  # satisfied
                   matrix(c(10, 30, 50), nrow = 1),  # angry
                   matrix(c(10, 30, 50), nrow = 1),  # anxious
                   matrix(c(10, 30, 50), nrow = 1),  # depressed
                   matrix(c(10, 30, 50), nrow = 1)), # sad 
  emiss_K0 = rep(list(1), n_dep),
  emiss_V = rep(list(rep(5^2, m)), n_dep),
  emiss_nu = rep(list(1), n_dep),
  emiss_a0 = rep(list(rep(1.5, m)), n_dep),
  emiss_b0 = rep(list(rep(20, m)), n_dep),
)

out_3st_emotion <- mHMM(s_data = emotion_mHMM,
                        data_distr = "continuous",
                        gen = list(m = m, n_dep = n_dep),
                        start_val = c(list(start_gamma_3st), start_emiss_3st),
                        emiss_hyp_prior = emotion_prior_emiss_3st,
                        mcmc = list(J = 500, burn_in = 200))

saveRDS(out_3st_emotion, file = "out_3st_emotion.rds")

summary(out_3st_emotion)
out_3st_emotion

#####################
## 4 state model ####
#####################
## general model properties
m <- 4
n_dep <- 8

## starting values for gamma
start_gamma_4st <- matrix(c(0.7, 0.1, 0.1, 0.1,
                            0.1, 0.7, 0.1, 0.1,
                            0.1, 0.1, 0.7, 0.1,
                            0.1, 0.1, 0.1, 0.7), byrow = TRUE, ncol = m, nrow = m)

## starting values for the emission distribution 
start_emiss_4st <- list(matrix(c(80, 10,                                     # happy
                                 60, 10,
                                 40, 10,
                                 20, 10), byrow = TRUE, ncol = 2, nrow = m), 
                        matrix(c(80, 10,                                     # excited
                                 60, 10,
                                 40, 10,
                                 20, 10), byrow = TRUE, ncol = 2, nrow = m), 
                        matrix(c(80, 10,                                     # relaxed
                                 60, 10,
                                 40, 10,
                                 20, 10), byrow = TRUE, ncol = 2, nrow = m), 
                        matrix(c(80, 10,                                     # satisfied
                                 60, 10,
                                 40, 10,
                                 20, 10), byrow = TRUE, ncol = 2, nrow = m), 
                        matrix(c(10,  5,                                     # angry
                                 20, 10,
                                 40, 10,
                                 60, 10), byrow = TRUE, ncol = 2, nrow = m), 
                        matrix(c(10,  5,                                     # anxious
                                 320, 10,
                                 40, 10,
                                 60, 10), byrow = TRUE, ncol = 2, nrow = m), 
                        matrix(c(10,  5,                                     # depressed
                                 20, 10,
                                 40, 10,
                                 60, 10), byrow = TRUE, ncol = 2, nrow = m), 
                        matrix(c(10,  5,                                     # sad
                                 20, 10,
                                 40, 10,
                                 60, 10), byrow = TRUE, ncol = 2, nrow = m))

## specifying weakly informative prior for continuous emission distributions
emotion_prior_emiss_4st <- prior_emiss_cont(
  gen = list(m = m, n_dep = n_dep),
  emiss_mu0 = list(matrix(c(80, 60, 40, 20), nrow = 1),  # happy
                   matrix(c(80, 60, 40, 20), nrow = 1),  # excited
                   matrix(c(80, 60, 40, 20), nrow = 1),  # relaxed
                   matrix(c(80, 60, 40, 20), nrow = 1),  # satisfied
                   matrix(c(10, 20, 40, 60), nrow = 1),  # angry
                   matrix(c(10, 20, 40, 60), nrow = 1),  # anxious
                   matrix(c(10, 20, 40, 60), nrow = 1),  # depressed
                   matrix(c(10, 20, 40, 60), nrow = 1)), # sad 
  emiss_K0 = rep(list(1), n_dep),
  emiss_V = rep(list(rep(5^2, m)), n_dep),
  emiss_nu = rep(list(1), n_dep),
  emiss_a0 = rep(list(rep(1.5, m)), n_dep),
  emiss_b0 = rep(list(rep(20, m)), n_dep),
)

out_4st_emotion <- mHMM(s_data = emotion_mHMM,
                        data_distr = "continuous",
                        gen = list(m = m, n_dep = n_dep),
                        start_val = c(list(start_gamma_4st), start_emiss_4st),
                        emiss_hyp_prior = emotion_prior_emiss_4st,
                        mcmc = list(J = 500, burn_in = 200))

saveRDS(out_4st_emotion, file = "out_4st_emotion.rds")

out_4st_emotion
summary(out_4st_emotion)

## checking traceplots ####
data1 <- out_2st_emotion
J <- 500
burn_in <- 200
j <- 1 
i <- 1
  for(j in 1:n_dep){
    par(mfrow = c(m, 2))
    for(i in 1:m){
      plot(data1$emiss_mu_bar[[j]][, i], type = "l", 
           main = paste("Mu dep", j, "in state", i), 
           ylim = c(1,100))
      abline(h = mean(data1$emiss_mu_bar[[j]][((burn_in + 1): J), i]), col = "blue")
      plot(data1$emiss_sd_bar[[j]][, i], type = "l", 
           main = paste("Var dep", j, "in state", i), 
           ylim = c(0,30))
      abline(h = mean(data1$emiss_sd_bar[[j]][((burn_in + 1): J), i]), col = "blue")
    }
  }

n_subj <- 20

for(s in 1:n_subj){
  par(mfrow = c(2, m))
  for(j in 3:4){
    for(i in 1:m){
      plot(data1$PD_subj[[s]]$cont_emiss[, (j-1) * m + i], type = "l", 
           main = paste("Mu dep", j, "in state", i), 
           ylim = c(1,100))
      abline(h = mean(data1$PD_subj[[s]]$cont_emiss[((burn_in + 1): J), (j-1) * m + i]), col = "blue")
    }
  }
}
