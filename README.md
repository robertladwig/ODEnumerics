# Numerical schemes for ODEs

Non-negative and conservative discretisations are needed for water quality modeling. Most numerical schemes are neither. [Burchard, Deleersnijder and Meister](https://www.sciencedirect.com/science/article/pii/S0168927403001016) highlighted the implementation of high-order conservative Patankar-type discretisations to water quality ODEs. This repository includes code to replicate their findings in R for the simple linear model:

The package can be installed through

```{r gh-installation, eval = FALSE}
# dc1/dt = c2 - a c1
# dc2/dt = a c1 - c2
```

Code for forward Euler, Runge-Kutta 4th order, Patankar-Euler, modified Patankar-Euler and Mod. Patankar Runge-Kutta 2nd order are included.


<a href="url"><img src="figs/solver.png" width=100% height=100% ></a> 
