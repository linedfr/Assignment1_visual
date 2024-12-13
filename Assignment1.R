library(tidyverse);
library(lme4)
library(MASS)
library(plotly)
library(gridExtra)
library(pbkrtest)
library(lattice)

# Part I ----------------------------------------------

## Problem 5 -----------------------------------------

rho = 0.25

n_simul = 100000

Sigma = matrix(c(1, rho, 1/2,
                 rho, 1, 1/2,
                 1/2, 1/2, 1),
               byrow = T, nrow = 3) 

eps = c(0,0,0)

set.seed(888) #kesi

dist_join <- mvrnorm(n = n_simul, eps, Sigma)

dist_adf <- data.frame(dist_join)

Sigma12 <- cor(dist_adf$X1, dist_adf$X2)
Sigma13 <- cor(dist_adf$X1, dist_adf$X3)
Sigma23 <- cor(dist_adf$X2, dist_adf$X3)

Sigma11 <- var(dist_adf$X1)
Sigma22 <- var(dist_adf$X2)
Sigma33 <- var(dist_adf$X3)

cond_sigma <- matrix(c(Sigma11 - Sigma13/Sigma33*Sigma13, Sigma12 - Sigma13/Sigma33*Sigma23,
                       Sigma12 - Sigma23/Sigma33*Sigma13, Sigma22-Sigma23/Sigma33*Sigma23),
                     nrow = 2, byrow = T)

cond_sigma_theo <- matrix(c(3/4, rho-1/4,
                            rho-1/4, 3/4),
                          nrow = 2, byrow = T)

cond_xi <- c(NA, NA)

cond_xi[1] <- mean(Sigma13 * Sigma33^(-1) * dist_adf$X3)
cond_xi[2] <- mean(Sigma23 * Sigma33^(-1) * dist_adf$X3)

cond_xi_theo <- c(0,0)

cond_xi_theo
cond_xi

cond_sigma
cond_sigma_theo

# Part II ----------------------------------------------------------

load("assignment2024-1.Rdata")

func_hist <- function(variable) {
  ggplot(data = likingdata, aes(x = !!sym(variable))) +
    geom_bar(width = 0.5, fill = "skyblue") +
    theme_bw()
}

p1 <- func_hist("Class")
p2 <- func_hist("Liking")
p3 <- func_hist("School")
p4 <- func_hist("Grade")
p5 <- func_hist("Age")
p6 <- func_hist("Gender")

grid.arrange(p1, p2, p3, p4, p5, p6, nrow = 2)

## Problem 2 -----------------------------------------------------

fit1 <- lmer(Liking ~ Product + (1|Participant) + (1|Class), data=likingdata)
fit2 <- lmer(Liking ~ ProdVersion + ProdType + (1|Participant) + (1|Class), data=likingdata)

var_corr <- data.frame(VarCorr(fit1))

# Note that for i = 1,...,450, var(Y_i) = tau_1^2 + tau_2^2 + sigma^2. Kan udledes på samme måde som i MM3

# Observations from the same participant. Note that observations from the same participants also are observations from the same class

var_corr |>
  transmute(correlation = (sdcor[grp == "Participant"]^2 + sdcor[grp == "Class"]^2) / 
              (sdcor[grp == "Participant"]^2 + sdcor[grp == "Class"]^2 + sdcor[grp == "Residual"]^2)) |>
  distinct()

corparticipant <- (0.86638^2 + 0.45160^2)/(0.86638^2 + 0.45160^2 + 1.57387^2)
corparticipant

# Observations from two different participants from the same class

var_corr |>
  transmute(correlation = (sdcor[grp == "Class"]^2) / (sdcor[grp == "Participant"]^2 + sdcor[grp == "Class"]^2 + sdcor[grp == "Residual"]^2)) |>
  distinct()

corclass <- 0.45160^2/(0.86638^2 + 0.45160^2 + 1.57387^2)
corclass

## Problem 3 ------------------------------------------------------

anova(fit1, fit2)

## Problem 4 ----------------------------------------------------

version <- ggplot(likingdata, aes(x = Liking)) +
  geom_bar(position = "dodge", color = "black", fill = "lightblue") +
  labs(
    x = "Liking",
    y = "Frequency"
  ) +
  scale_x_continuous(breaks = 1:7) +  
  theme_minimal() +
  facet_wrap(~ProdVersion)
theme(
  legend.position = "right"
)

type <- ggplot(likingdata, aes(x = Liking, fill = ProdType)) +
  geom_bar(position = "dodge", color = "black") +
  labs(
    x = "Liking",
    y = "Frequency"
  ) +
  scale_fill_manual(values = c("lightgoldenrod1", "lightpink")) +
  scale_x_continuous(breaks = 1:7) +  
  theme_minimal() +
  theme(
    legend.position = "right"
  )

gridExtra::grid.arrange(version, type, ncol = 1)

drop1(fit2, test = "Chisq")

## Problem 5 -----------------------------------------------------------

X <- fit2 %>% model.matrix()
Z <- getME(fit2, "Z")
beta <- as.vector(fixef(fit2))
tauP <- as.data.frame(VarCorr(fit2))$sdcor[1] 
tauC <- as.data.frame(VarCorr(fit2))$sdcor[2]
sigma <- as.data.frame(VarCorr(fit2))$sdcor[3]

set.seed(1)

M <- 2000 #number of simulations

deltasim <- matrix(NA,M,3)

for (i in 1:M){
  B1 <- rnorm(75, mean = 0, sd = tauP)
  B2 <- rnorm(5, mean = 0, sd = tauC)
  eps <- rnorm(450, mean = 0, sd = sigma)
  B <- c(B1,B2)
  y <- X %*% beta + Z %*% B + eps
  y <- y %>% as.numeric() # NB. This seems to be necessary
  
  data <- cbind(likingdata, y)
  
  lmm2 <- lmer(y ~ ProdVersion + ProdType + (1|Participant) + (1|Class), data=data)
  deltasim[i,1] <- fixef(lmm2)[4]
  deltasim[i,2:3] <- (lmm2 %>% confint(method="Wald"))[7,]
}

deltasim <- deltasim %>% data.frame()
names(deltasim) <- c("est","lower","upper")

delta <- fixef(fit2)[4]
delta_c <- (fit2 %>% confint(method="Wald"))[7,]

ggplot(deltasim, aes(est)) +
  geom_histogram(bins = 60, fill = "lightblue", color = "black") +
  geom_vline(xintercept = delta, color = "blue") +
  geom_vline(xintercept = delta_c[1], linetype = "dashed") +
  geom_vline(xintercept = delta_c[2], linetype = "dashed") +
  labs(
    x = "delta from simulations",  # New x-axis label
    y = "Count",  # New y-axis label
  ) +
  theme_bw()

delta_c <- (fit2 %>% confint(method="Wald"))[7,]

delta_df <- cbind(delta, delta_c)

data_CI_vs2 <- deltasim |>
  transmute(est,
            lower = delta_c[1],
            upper = delta_c[2],
            in_CI = ifelse(est >= lower & est <=upper,1,0))

nominal_ci_2 <- sum(data_CI_vs2$in_CI)/M
nominal_ci_2

bias_n <- delta - mean(deltasim$est)

## Problem 6 ---------------------------------------------------------

### Student t ----------------------------------------

X <- fit2 %>% model.matrix()
Z <- getME(fit2, "Z")
beta <- as.vector(fixef(fit2))
tauP <- as.data.frame(VarCorr(fit2))$sdcor[1] #skal man kalde via X her?
tauC <- as.data.frame(VarCorr(fit2))$sdcor[2]
sigma <- as.data.frame(VarCorr(fit2))$sdcor[3]

df <- 3

set.seed(1708)

M <- 2000 #number of simulations

deltasim_t <- matrix(NA,M,3)

for (i in 1:M){
  B1 <- rt(75, 3)
  B1_scaled <- (B1 - mean(0)) / sqrt(3/(3-2)) * sqrt(tauP^2) + 0
  print(var(B1_scaled))
  print(mean(B1_scaled))
  B2 <- rt(5, 3)
  B2_scaled <- (B2 - mean(0)) / sqrt(3/(3-2)) * sqrt(tauC^2) + 0
  print(var(B2_scaled))
  print(mean(B2_scaled))
  eps <- rt(450, 3)
  eps_scaled <- (eps - mean(0)) / sqrt(3/(3-2)) * sqrt(sigma^2) + 0
  
  B <- c(B1_scaled,B2_scaled)
  y <- X %*% beta + Z %*% B + eps_scaled
  y <- y %>% as.numeric() # NB. This seems to be necessary
  
  data <- cbind(likingdata, y)
  
  lmm2 <- lmer(y ~ ProdVersion + ProdType + (1|Participant) + (1|Class), data=data)
  
  deltasim_t[i,1] <- fixef(lmm2)[4]
  deltasim_t[i,2:3] <- (lmm2 %>% confint(method="Wald"))[7,]
}

deltasim_t <- deltasim_t %>% data.frame()
names(deltasim_t) <- c("est","lower","upper")

data_CI_t <- deltasim_t |>
  transmute(est,
            lower = delta_c[1],
            upper = delta_c[2],
            in_CI = ifelse(est >= lower & est <=upper,1,0))

nominal_ci_t <- sum(data_CI_t$in_CI)/M
nominal_ci_t

bias_t <- delta - mean(deltasim_t$est)

ggplot(deltasim_t, aes(est)) +
  geom_histogram(bins = 60, fill = "lightblue", color = "black") +
  geom_vline(xintercept = delta, color = "blue") +
  geom_vline(xintercept = delta_c[1], linetype = "dashed") +
  geom_vline(xintercept = delta_c[2], linetype = "dashed") +
  labs(
    x = "delta from simulations",  # New x-axis label
    y = "Count",  # New y-axis label
  ) +
  theme_bw()

### Exponential --------------------------------------------------

X <- fit2 %>% model.matrix()
Z <- getME(fit2, "Z")
beta <- as.vector(fixef(fit2))
tauP <- as.data.frame(VarCorr(fit2))$sdcor[1] #skal man kalde via X her?
tauC <- as.data.frame(VarCorr(fit2))$sdcor[2]
sigma <- as.data.frame(VarCorr(fit2))$sdcor[3]

set.seed(888)

M <- 2000 #number of simulations

deltasim_exp <- matrix(NA,M,3)

for (i in 1:M){
  B1 <- rexp(75, 1/tauP) 
  B1_scaled <- (B1 - tauP)
  
  B2 <- rexp(5, 1/tauC)
  B2_scaled <- (B2 - tauC)
  
  eps <- rexp(450, 1/sigma)
  eps_scaled <- (eps - sigma)
  
  B <- c(B1_scaled,B2_scaled)
  y <- X %*% beta + Z %*% B + eps_scaled
  y <- y %>% as.numeric() # NB. This seems to be necessary
  
  data <- cbind(likingdata, y)
  
  lmm2 <- lmer(y ~ ProdVersion + ProdType + (1|Participant) + (1|Class), data=data)
  deltasim_exp[i,1] <- fixef(lmm2)[4]
  deltasim_exp[i,2:3] <- (lmm2 %>% confint(method="Wald"))[7,]
}

deltasim_exp <- deltasim_exp %>% data.frame()
names(deltasim_exp) <- c("est","lower","upper")

data_CI_exp <- deltasim_exp |>
  transmute(est,
            lower = delta_c[1],
            upper = delta_c[2],
            in_CI = ifelse(est >= lower & est <=upper,1,0))

data_CI_exp <- sum(data_CI_exp$in_CI)/M
data_CI_exp

bias_exp <- delta - mean(deltasim_exp$est)

ggplot(deltasim_exp, aes(est)) +
  geom_histogram(bins = 60, fill = "lightblue", color = "black") +
  geom_vline(xintercept = delta, color = "blue") +
  geom_vline(xintercept = delta_c[1], linetype = "dashed") +
  geom_vline(xintercept = delta_c[2], linetype = "dashed") +
  labs(
    x = "delta from simulations",  # New x-axis label
    y = "Count",  # New y-axis label
  ) +
  theme_bw()

## Problem 7 ------------------------------------------------

### Residuals -----------------------------------------

residuals <- resid(fit2)
fitted_values <- fitted(fit2)

plot_data <- data.frame(
  Fitted = fitted_values,
  Residuals = residuals
)

qqplot <- ggplot(data.frame(residuals = residuals), aes(sample = residuals)) +
  stat_qq(size = 0.7) +
  stat_qq_line(color = "blue", size = 0.7) +
  labs(title = "QQ Plot of Residuals",
       x = "Theoretical Quantiles",
       y = "Sample Quantiles") +
  theme_bw() +
  theme(plot.title = element_text(size = 9)) +
  coord_fixed(ratio = 0.8)

residplot <- ggplot(plot_data, aes(x = Fitted, y = Residuals)) +
  geom_point(color = "black", alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "blue") +
  labs(
    title = "Residuals vs Fitted Values",
    x = "Fitted Values",
    y = "Residuals"
  ) +
  theme_bw() +
  theme(plot.title = element_text(size = 9)) +
  coord_fixed(ratio = 0.9)

histplot <- ggplot(plot_data, aes(residuals)) +
  # Histogram scaled to density
  geom_histogram(aes(y = ..density..), bins = 25, color = "black", fill = "lightblue") +
  # Overlay the normal density curve
  stat_function(fun = dnorm, args = list(mean = 0, sd = sigma), color = "blue", size = 0.5) +
  # Add theme
  theme_bw() +
  # Add labels
  labs(title = "Histogram with Density N(0,sigma^2)",
       x = "Residuals",
       y = "Density") +
  theme(plot.title = element_text(size = 8)) +
  coord_fixed(ratio = 38)

grid.arrange(qqplot, residplot, histplot, nrow = 1)

### Random effects --------------------------------------------

fit2 %>% ranef() %>% dotplot()

random_fit <- ranef(fit2)

qq_participant <- ggplot(data.frame(residuals = random_fit$Participant[,1]), aes(sample = random_fit$Participant[,1])) +
  stat_qq() +
  stat_qq_line(color = "blue") +
  labs(title = "QQ Plot of B_participant",
       x = "Theoretical Quantiles",
       y = "Sample Quantiles") +
  theme_bw()

ggplot(data.frame(residuals = random_fit$Class[,1]), aes(sample = random_fit$Class[,1])) +
  stat_qq() +
  stat_qq_line(color = "blue") +
  labs(title = "QQ Plot of B_Class",
       x = "Theoretical Quantiles",
       y = "Sample Quantiles") +
  theme_bw()

hist_participant <- ggplot(data.frame(data = random_fit$Participant[,1]), aes(data)) +
  # Histogram scaled to density
  geom_histogram(aes(y = ..density..), bins = 20, color = "black", fill = "lightblue", alpha = 1) +
  # Overlay the normal density curve
  stat_function(fun = dnorm, args = list(mean = 0, sd = tauP), color = "blue", size = 0.7) +
  # Add theme
  theme_bw() +
  # Add labels
  labs(title = "Histogram with density N(0, tau_P^2)",
       x = "Data",
       y = "Density")

grid.arrange(qq_participant, hist_participant, nrow = 1)

## problem 8 ---------------------------------------

### Discrete ----------------

X <- fit2 %>% model.matrix()
Z <- getME(fit2, "Z")
beta <- as.vector(fixef(fit2))
tauP <- as.data.frame(VarCorr(fit2))$sdcor[1] #skal man kalde via X her?
tauC <- as.data.frame(VarCorr(fit2))$sdcor[2]
sigma <- as.data.frame(VarCorr(fit2))$sdcor[3]

M <- 2000 #number of simulations

set.seed(200)

deltasim_discrete <- matrix(NA,M,3)

for (i in 1:M){
  B1 <- rnorm(75, mean = 0, sd = tauP)
  B1_discrete <- round(B1,0)
  B2 <- rnorm(5, mean = 0, sd = tauC)
  eps <- rnorm(450, mean = 0, sd = sigma)
  B <- c(B1,B2)
  y <- X %*% beta + Z %*% B + eps
  y <- as.numeric(y) # NB. This seems to be necessary
  
  y <- round(y,0)
  
  data <- cbind(likingdata, y)
  
  lmm2 <- lmer(y ~ ProdVersion + ProdType + (1|Participant) + (1|Class), data=data)
  
  deltasim_discrete[i,1] <- fixef(lmm2)[4]
  deltasim_discrete[i,2:3] <- (lmm2 %>% confint(method="Wald"))[7,]
}

deltasim_discrete <- deltasim_discrete %>% data.frame()
names(deltasim_discrete) <- c("est","lower","upper")

data_CI_disc <- deltasim_discrete |>
  transmute(est,
            lower = delta_c[1],
            upper = delta_c[2],
            in_CI = ifelse(est >= lower & est <=upper,1,0))

data_CI_disc <- sum(data_CI_disc$in_CI)/M
data_CI_disc

bias_discrete <- delta - mean(deltasim_discrete$est)
bias_discrete

ggplot(deltasim_discrete, aes(est)) +
  geom_histogram(bins = 60, fill = "lightblue", color = "black") +
  geom_vline(xintercept = delta, color = "blue") +
  geom_vline(xintercept = delta_c[1], linetype = "dashed") +
  geom_vline(xintercept = delta_c[2], linetype = "dashed") +
  labs(
    x = "delta from simulations",  # New x-axis label
    y = "Count",  # New y-axis label
  ) +
  theme_bw()

### Discrete and bound --------------------------------

X <- fit2 %>% model.matrix()
Z <- getME(fit2, "Z")
beta <- as.vector(fixef(fit2))
tauP <- as.data.frame(VarCorr(fit2))$sdcor[1] #skal man kalde via X her?
tauC <- as.data.frame(VarCorr(fit2))$sdcor[2]
sigma <- as.data.frame(VarCorr(fit2))$sdcor[3]

M <- 2000 #number of simulations

set.seed(541)

deltasim_db <- matrix(NA,M,3)

for (i in 1:M){
  B1 <- rnorm(75, mean = 0, sd = tauP)
  B1_discrete <- round(B1,0)
  B2 <- rnorm(5, mean = 0, sd = tauC)
  eps <- rnorm(450, mean = 0, sd = sigma)
  B <- c(B1,B2)
  y <- X %*% beta + Z %*% B + eps
  y <- as.numeric(y) # NB. This seems to be necessary
  
  y <- case_when(y < 1 ~ 1,
                 y > 7 ~ 7,
                 TRUE ~ round(y,0))
  
  data <- cbind(likingdata, y)
  
  lmm2 <- lmer(y ~ ProdVersion + ProdType + (1|Participant) + (1|Class), data=data)
  
  deltasim_db[i,1] <- fixef(lmm2)[4]
  deltasim_db[i,2:3] <- (lmm2 %>% confint(method="Wald"))[7,]
}

deltasim_db <- deltasim_db %>% data.frame()
names(deltasim_db) <- c("est","lower","upper")

data_CI_db <- deltasim_db |>
  transmute(est,
            lower = delta_c[1],
            upper = delta_c[2],
            in_CI = ifelse(est >= lower & est <=upper,1,0))

data_CI_db <- sum(data_CI_db$in_CI)/M
data_CI_db

bias_db <- delta - mean(deltasim_db$est)
bias_db

ggplot(deltasim_db, aes(est)) +
  geom_histogram(bins = 60, fill = "lightblue", color = "black") +
  geom_vline(xintercept = delta, color = "blue") +
  geom_vline(xintercept = delta_c[1], linetype = "dashed") +
  geom_vline(xintercept = delta_c[2], linetype = "dashed") +
  labs(
    x = "delta from simulations",  # New x-axis label
    y = "Count",  # New y-axis label
  ) +
  theme_bw()


