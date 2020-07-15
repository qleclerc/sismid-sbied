
## SETUP ########

library(doParallel)
library(doRNG)
library(tidyverse)
library(pomp)

registerDoParallel()
registerDoRNG(2488820)

#load data
read_csv(paste0("https://kingaa.github.io/sbied/stochsim/",
                "Measles_Consett_1948.csv")) %>%
  select(week,reports=cases) -> meas
meas %>% as.data.frame() %>% head()

meas %>%
  ggplot(aes(x=week,y=reports))+
  geom_line()+
  geom_point()



## POMP WITH R FUNCTIONS ########

#define transitions
sir_step <- function (S, I, R, H, N, Beta, mu_IR, delta.t, ...)
{
  dN_SI <- rbinom(n=1,size=S,prob=1-exp(-Beta*I/N*delta.t))
  dN_IR <- rbinom(n=1,size=I,prob=1-exp(-mu_IR*delta.t))
  S <- S - dN_SI
  I <- I + dN_SI - dN_IR
  R <- R + dN_IR
  H <- H + dN_IR;
  c(S = S, I = I, R = R, H = H)
}

#define initial values
sir_rinit <- function (N, eta, ...) {
  c(S = round(N*eta), I = 1, R = round(N*(1-eta)), H = 0)
}

#create pomp object
meas %>%
  pomp(times="week",t0=0,
       rprocess=euler(sir_step,delta.t=1/7),
       rinit=sir_rinit,
       accumvars="H",
       statenames=c("S","I","R","H"),
       paramnames=c("Beta","mu_IR","N","eta","rho")
  ) -> measSIR


#define density for reporting
sir_dmeas <- function (reports, H, rho, log, ...) {
  dbinom(x=reports, size=H, prob=rho, log=log)
}

#define sampling for reporting
sir_rmeas <- function (H, rho, ...) {
  c(reports=rbinom(n=1, size=H, prob=rho))
}

#update pomp object
measSIR %>%
  pomp(
    rmeasure=sir_rmeas,
    dmeasure=sir_dmeas
  ) -> measSIR

#run 20 simulations
measSIR %>%
  simulate(
    params=c(Beta=7.5,mu_IR=0.5,rho=0.5,eta=0.03,N=38000),
    nsim=20,format="data.frame",include.data=TRUE
  ) -> sims

sims %>%
  ggplot(aes(x=week,y=reports,group=.id,color=.id=="data"))+
  geom_line()+
  guides(color=FALSE)



## POMP WITH C SNIPPETS ########


sir_step <- Csnippet("
  double dN_SI = rbinom(S,1-exp(-Beta*I/N*dt));
  double dN_IR = rbinom(I,1-exp(-mu_IR*dt));
  S -= dN_SI;
  I += dN_SI - dN_IR;
  R += dN_IR;
  H += dN_IR;
  ")

sir_rinit <- Csnippet("
  S = nearbyint(eta*N);
  I = 1;
  R = nearbyint((1-eta)*N);
  H = 0;
  ")

sir_dmeas <- Csnippet("
  lik = dbinom(reports,H,rho,give_log);
  ")

sir_rmeas <- Csnippet("
  reports = rbinom(H,rho);
  ")

meas %>%
  pomp(times="week",t0=0,
       rprocess=euler(sir_step,delta.t=1/7),
       rinit=sir_rinit,
       rmeasure=sir_rmeas,
       dmeasure=sir_dmeas,
       accumvars="H",
       statenames=c("S","I","R","H"),
       paramnames=c("Beta","mu_IR","N","eta","rho")
  ) -> measSIR

measSIR %>%
  simulate(
    params=c(Beta=7.5,mu_IR=0.5,rho=0.5,eta=0.03,N=38000),
    nsim=20,format="data.frame",include.data=TRUE
  ) -> sims

sims %>%
  ggplot(aes(x=week,y=reports,group=.id,color=.id=="data"))+
  geom_line()+
  guides(color=FALSE)



## TWEAKING PARAMETERS ########

#default
measSIR@params <- c(Beta=7.5,mu_IR=0.5,rho=0.5,eta=0.03,N=38000)

#to play around with
custom_reporting = 0.35
custom_eta = 521/custom_reporting/38000
measSIR@params <- c(Beta=18,mu_IR=0.2,rho=custom_reporting,eta=custom_eta,N=38000)

measSIR %>%
  simulate(
    nsim=100,format="data.frame",include.data=TRUE
  ) -> sims

sims %>%
  ggplot(aes(x=week,y=reports,group=.id,color=.id=="data"))+
  geom_line()+
  guides(color=FALSE)



## SEIR POMP ########


seir_step <- Csnippet("
  double dN_SE = rbinom(S,1-exp(-Beta*I/N*dt));
  double dN_EI = rbinom(E,1-exp(-mu_EI*dt));
  double dN_IR = rbinom(I,1-exp(-mu_IR*dt));
  S -= dN_SE;
  E += dN_SE - dN_EI;
  I += dN_EI - dN_IR;
  R += dN_IR;
  H += dN_IR;
  ")

seir_rinit <- Csnippet("
  S = nearbyint(eta*N);
  E = 0;
  I = 1;
  R = nearbyint((1-eta)*N);
  H = 0;
  ")

seir_dmeas <- Csnippet("
  lik = dbinom(reports,H,rho,give_log);
  ")

seir_rmeas <- Csnippet("
  reports = rbinom(H,rho);
  ")

meas %>%
  pomp(times="week",t0=0,
       rprocess=euler(seir_step,delta.t=1/7),
       rinit=seir_rinit,
       rmeasure=seir_rmeas,
       dmeasure=seir_dmeas,
       accumvars="H",
       statenames=c("S","E","I","R","H"),
       paramnames=c("Beta","mu_EI","mu_IR","N","eta","rho")
  ) -> measSEIR

measSEIR %>%
  simulate(
    params=c(Beta=7.5,mu_EI=0.5,mu_IR=0.5,rho=0.5,eta=0.03,N=38000),
    nsim=20,format="data.frame",include.data=TRUE
  ) -> sims

sims %>%
  ggplot(aes(x=week,y=reports,group=.id,color=.id=="data"))+
  geom_line()+
  guides(color=FALSE)


#default
measSEIR@params <- c(Beta=7.5,mu_EI=0.5,mu_IR=0.5,rho=0.5,eta=0.03,N=38000)

#to play around with
measSEIR@params <- c(Beta=20,mu_EI=1/1,mu_IR=1/3,rho=0.5,eta=0.05,N=38000)

measSEIR %>%
  simulate(
    nsim=100,format="data.frame",include.data=TRUE
  ) -> sims

sims %>%
  ggplot(aes(x=week,y=reports,group=.id,color=.id=="data"))+
  geom_line()+
  guides(color=FALSE)

