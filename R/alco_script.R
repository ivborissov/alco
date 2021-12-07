library(dMod)
library(dplyr)
library(ggplot2)

# Define reactions
reactions <- NULL
reactions <- addReaction(reactions, "q_a", "C", rate = "k_a * q_a", description = "vabs")
reactions <- addReaction(reactions, "C", "", rate = "Vmax * C / (Km + (C / blood))", description = "vel")

# Define events
#myevents <- addEvent(NULL, var = "q_a", time = 0, value = "dose_first * body_weight", method = "add")
myevents <- addEvent(NULL, var = "q_a", time = 0.333, value = "dose_first * body_weight", method = "add")
myevents <- addEvent(myevents, var = "q_a", time = 0.667, value = "dose_first * body_weight", method = "add")
myevents <- addEvent(myevents, var = "q_a", time = 1.5, value = "dose_last * body_weight", method = "add")

# Translate into ODE model
alco_model <- odemodel(reactions, modelname = "Alcomodel", events = myevents)

# Prediction function
x <- Xs(alco_model, names = "C")

# Define observations
observables <- eqnvec(
  BrAC = "C / blood" 
  #sigma1 = "sqrt(sigma_add^2 + C^2 * sigma_prop^2)"
)

# Generate observation function
g <- Y(observables, reactions, compile = TRUE, modelname = "obsfn", attach.input = FALSE)

# Identity transformation
trafo <- list(scn1 = c(q_a = 13.02, C = 0, k_a = "k_a", Vmax = "Vmax",
Km = "Km", Vd = "Vd", dose_first = 0.186,
dose_last = 0.112, body_weight = 70, blood = "Vd*70",
sigma_add = 0.1, sigma_prop = 0))


# Generate parameter transformation
p <- P(trafo)

# Experimental Data
rawdata <- read.csv("D:/INSYSBIO/repos/alco/R/data.csv")
df <- data.frame(
  name = "BrAC",
  time = rawdata$t,
  value = rawdata$measurement,
  sigma = rawdata$prob.sigma,
  condition = "scn1"
)
data <- as.datalist(df, split.by = "condition")

times <- df$time

pars <- c(k_a=0.062, Vmax=0.136, Km=0.09601, Vd=0.457)

out <- (g*x*p)(times, pars)

counter <- 0

obj <- normL2(data, g*x*p)

fitpars <- c(Vd = 0.5474208905879729, Vmax = 0.25500177168711313, Km = 0.12558367707084078, k_a = 4.518139485708639)

myprofile <- profile(obj2, pars = fitpars, whichPar = names(fitpars))

jpeg(file="./R/k_a.jpeg")
plotProfile(myprofile)
dev.off()

# Prediction Bands for BrAC
out_df <- as.data.frame(out$scn1)


prediction_band <- do.call(rbind, lapply(df$time, function(t) {
  
  cat("Computing prediction profile for t =", t, "\n")
  
  sigma <- data$scn1$sigma[data$scn1$time == t]
  obj.validation <- datapointL2(name = "BrAC",
                                time = t,
                                value = "d1",
                                sigma = 1e-3,
                                condition = "scn1")
  
  myfit <- trust(obj + obj.validation, 
        parinit = c(d1 = 1, fitpars), 
        rinit = 1, rmax = 10)
  
  #d1_val <- out_df$BrAC[out_df$time == t]
  counter <<- 0
  obj2 <- function(...)
  {
    counter <<- counter + 1
    return((obj + obj.validation)(...))
  }

  class(obj2) <- c("objfn", "fn")

  profile_prediction <- profile(obj2,
                                myfit$argument, "d1")
  
  d1 <- confint(profile_prediction, val.column = "value")
  
  # Output
  data.frame(time = t, condition = "scn1", name = "BrAC",  d1[-1], iter = counter)
}))

jpeg(file="./R/r_conf_band.jpeg")
plot((g*x*p)(times, fitpars), data) + 
  geom_ribbon(aes(x = time, ymin = lower, ymax = upper), 
              data = prediction_band,
              lty = 0, alpha = .3)
dev.off()
