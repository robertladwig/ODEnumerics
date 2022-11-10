# Burchard et al. (2003) "A high-order conservative Patankar-type discretisation
# for stiff systems of production–destruction equations
library(tidyverse)

#### PROBLEM SET, TWO COUPLED ODES
# dc1/dt = c2 - a c1
# dc2/dt = a c1 - c2

dc1dt <- function(c1, c2, a){
  return(c2 - a * c1)
}

dc2dt <- function(c1, c2, a){
  return(a * c1 - c2)
}

#### PARAMETERS
a = 5
c10 = 0.9
c20 = 0.1

c1inf = (c10 + c20) / (a + 1)

dt = 0.25

t = seq(0, 1.8, dt)

c = c10/c1inf - 1

#### ANALYTICAL SOLUTION
c1 = (1 + c *  exp(-(a+1)* t)) * c1inf
c2 = c10 + c20 - c1

#### NUMERICAL EXPLICIT FORWARD EULER SOLUTION
c1euler = rep(NA, length(t))
c2euler = rep(NA, length(t))
c1euler[1] = c10
c2euler[1] = c20

for (i in t[2:length(t)]){
  c1euler[match(i, t)] = c1euler[match(i, t) -1] +
    dt * (c2euler[match(i, t) -1] - c1euler[match(i, t) -1] * a)
  c2euler[match(i, t)] = c2euler[match(i, t) -1] +
    dt * (a * c1euler[match(i, t) -1] - c2euler[match(i, t) -1])
}

#### NUMERICAL 4TH ORDER RUNGE KUTTA SOLUTION
c1rk = rep(NA, length(t))
c2rk = rep(NA, length(t))
c1rk[1] = c10
c2rk[1] = c20

for (i in t[2:length(t)]){
  k1 = dc1dt(c1rk[match(i, t) -1], c2rk[match(i, t) -1], a)
  l1 = dc2dt(c1rk[match(i, t) -1], c2rk[match(i, t) -1], a)
  
  k2 = dc1dt(c1rk[match(i, t) -1] + 1/2 * dt * k1, c2rk[match(i, t) -1] + 1/2 * dt * l1, a)
  l2 = dc2dt(c1rk[match(i, t) -1] + 1/2 * dt * k1, c2rk[match(i, t) -1] + 1/2 * dt * l1, a)
  
  k3 = dc1dt(c1rk[match(i, t) -1] + 1/2 * dt * k2, c2rk[match(i, t) -1] + 1/2 * dt * l2, a)
  l3 = dc2dt(c1rk[match(i, t) -1] + 1/2 * dt * k2, c2rk[match(i, t) -1] + 1/2 * dt * l2, a)
  
  k4 = dc1dt(c1rk[match(i, t) -1] + dt * k3, c2rk[match(i, t) -1] + dt * l3, a)
  l4 = dc2dt(c1rk[match(i, t) -1] + dt * k3, c2rk[match(i, t) -1] + dt * l3, a)
  
  c1rk[match(i, t)] = c1rk[match(i, t) -1] +
    dt * (1/6 * (k1 + 2 * k2 + 2 * k3 + k4))
  
  c2rk[match(i, t)] = c2rk[match(i, t) -1] +
    dt * (1/6 * (l1 + 2 * l2 + 2 * l3 + l4))
}

#### NUMERICAL EULER PATANKAR SOLUTION
c1pteuler = rep(NA, length(t))
c2pteuler = rep(NA, length(t))
c1pteuler[1] = c10
c2pteuler[1] = c20

for (i in t[2:length(t)]){
  # dc1/dt = c2 - a c1
  # dc2/dt = a c1 - c2
  p1 = c2pteuler[match(i, t) -1] 
  p2 = a * c1pteuler[match(i, t) -1]
  d1 = p2
  d2 = p1
  p = sum(p1,p2)
  d = sum(d1,d2)
  c1pteuler[match(i, t)] = (c1pteuler[match(i, t) -1] + dt * p1) /
    (1 + dt * d1 / c1pteuler[match(i, t) -1])

  c2pteuler[match(i, t)] = (c2pteuler[match(i, t) -1] + dt * p2) /
    (1 + dt * d2 / c2pteuler[match(i, t) -1])
}

#### NUMERICAL MODIFIED EULER PATANKAR SOLUTION
c1modpteuler = rep(NA, length(t))
c2modpteuler = rep(NA, length(t))
c1modpteuler[1] = c10
c2modpteuler[1] = c20

len_y0 = length(c(c10,c20))
eye = diag(len_y0)
eye = ifelse(eye == 1, TRUE, FALSE)
eyetilde = ifelse(eye == 1, FALSE, TRUE)
avec = eye * 0
r = avec[,1] * 0
for (i in t[2:length(t)]){
  # dc1/dt = c2 - a c1
  # dc2/dt = a c1 - c2
  p1 = c(0, c2modpteuler[match(i, t) -1])
  p2 = c(a * c1modpteuler[match(i, t) -1],0)
  d1 = p2
  d2 = p1
  p = matrix(c(p1,p2), nrow = len_y0, ncol = len_y0, byrow = T)
  d = matrix(c(d1,d2), nrow = len_y0, ncol = len_y0, byrow = T)
  
  ydat = c(c1modpteuler[match(i, t) -1],
           c2modpteuler[match(i, t) -1])
  
  avec[eye] = dt * c(sum(d1),
                     sum(d2)) / ydat +1
  
  c_rep = rbind(ydat, ydat)
  
  avec[eyetilde] = -dt * p[eyetilde] / c_rep[eyetilde]
  
  r =  ydat + dt * p[eye]
  c1modpteuler[match(i, t)] = solve(avec, r)[1]
  c2modpteuler[match(i, t)] = solve(avec, r)[2]
}

#### NUMERICAL MODIFIED RUNGE-KUTTA PATANKAR SOLUTION
c1modptrk = rep(NA, length(t))
c2modptrk = rep(NA, length(t))
c1modptrk[1] = c10
c2modptrk[1] = c20

len_y0 = length(c(c10,c20))
eye = diag(len_y0)
eye = ifelse(eye == 1, TRUE, FALSE)
eyetilde = ifelse(eye == 1, FALSE, TRUE)
avec = eye * 0
r = avec[,1] * 0
for (i in t[2:length(t)]){
  # dc1/dt = c2 - a c1
  # dc2/dt = a c1 - c2
  p1 = c(0, c2modpteuler[match(i, t) -1])
  p2 = c(a * c1modpteuler[match(i, t) -1],0)
  d1 = p2
  d2 = p1
  p0 = matrix(c(p1,p2), nrow = len_y0, ncol = len_y0, byrow = T)
  d0 = matrix(c(d1,d2), nrow = len_y0, ncol = len_y0, byrow = T)
  
  ydat = c(c1modpteuler[match(i, t) -1],
           c2modpteuler[match(i, t) -1])
  
  avec[eye] = dt * c(sum(d1),
                     sum(d2)) / ydat +1
  
  c_rep = rbind(ydat, ydat)
  
  avec[eyetilde] = -dt * p0[eyetilde] / c_rep[eyetilde]
  
  r =  ydat + dt * p0[eye]
  
  cproxy1 = solve(avec, r)[1]
  cproxy2 = solve(avec, r)[2]
  
  p1 = c(0, cproxy2)
  p2 = c(a * cproxy1,0)
  d1 = p2
  d2 = p1
  p = matrix(c(p1,p2), nrow = len_y0, ncol = len_y0, byrow = T)
  d = matrix(c(d1,d2), nrow = len_y0, ncol = len_y0, byrow = T)
  
  p = 0.5 * (p0 + p)
  d = 0.5 * (d0 + d)
  
  avec[eye] = dt * c(sum(d[1:2]),
                     sum(d[3:4])) / c(cproxy1, cproxy2) +1
  
  c_rep = rbind(c(cproxy1, cproxy2), c(cproxy1, cproxy2))
  
  avec[eyetilde] = -dt * p[eyetilde] / c_rep[eyetilde]
  
  r = ydat + dt * p[eye]
  
  c1modptrk[match(i, t)] = solve(avec, r)[1]
  c2modptrk[match(i, t)] = solve(avec, r)[2]
}

#### PLOTTING OF RESULTS
df = rbind(
  data.frame('time' = t, 'conc' = c1, 'type' = 'c1', 'solver' = 'Analytical'),
  data.frame('time' = t, 'conc' = c1euler, 'type' = 'c1', 'solver' = 'Forward Euler'),
  data.frame('time' = t, 'conc' = c1rk, 'type' = 'c1', 'solver' = 'Runge-Kutta 4th'),
  data.frame('time' = t, 'conc' = c1pteuler, 'type' = 'c1', 'solver' = 'Patankar-Euler'),
  data.frame('time' = t, 'conc' = c1modpteuler, 'type' = 'c1', 'solver' = 'Mod. Patankar-Euler'),
  data.frame('time' = t, 'conc' = c1modptrk, 'type' = 'c1', 'solver' = 'Mod. Patankar-RK2'),
  data.frame('time' = t, 'conc' = c2, 'type' = 'c2', 'solver' = 'Analytical'),
  data.frame('time' = t, 'conc' = c2euler, 'type' = 'c2', 'solver' = 'Forward Euler'),
  data.frame('time' = t, 'conc' = c2rk, 'type' = 'c2', 'solver' = 'Runge-Kutta 4th'),
  data.frame('time' = t, 'conc' = c2pteuler, 'type' = 'c2', 'solver' = 'Patankar-Euler'),
  data.frame('time' = t, 'conc' = c2modpteuler, 'type' = 'c2', 'solver' = 'Mod. Patankar-Euler'),
  data.frame('time' = t, 'conc' = c2modptrk, 'type' = 'c2', 'solver' = 'Mod. Patankar-RK2')
  )

ggplot(df) +
  geom_line(aes(time, conc, col = solver, linetype = type)) +
  theme_bw()
