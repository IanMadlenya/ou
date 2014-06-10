Ornstein Uhlenbeck/Vasicek Model of Crop Prices in Alberta 
========================================================


## 1 Introduction

We use the historic data to estimate the parameters for Ornstein Uhlenbeck process of crop prices in Alberta. And then we use the estimated parameters for Monte Carlo simulation.

The Ornstein-Uhlenbeck or Vasicek process is a stochastic process which is stationary, Gaussian, and Markovian. 

Over time, the process tends to drift towards its long-term mean: such a process is called mean-reverting.

The Ornstein-Uhlenbeck or Vasicek process is the unique solution to the following stochastic differential equation:(Stochastic Differential Equation, 2008, p44.)


"The Ornstein-Uhlenbeck process is one of several approaches used to model (with modifications) interest rates, currency exchange rates, and commodity prices stochastically. The parameter $\mu$ represents the equilibrium or mean value supported by fundamentals; $\sigma$ the degree of volatility around it caused by shocks, and $\theta$ the rate by which these shocks dissipate and the variable reverts towards the mean."(Wikipedia)



## 4  Conclusion

This documents show we can use ouFit.ML function to estimate the parameter $\theta$, $\mu$, and $\sigma$.
And them use the **"yuima"** package to generate the simulative price vector/matrices. 
