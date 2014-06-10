Ornstein Uhlenbeck/Vasicek Model of Crop Prices in Alberta 
========================================================


## 1 Introduction

We use the historic data to estimate the parameters for Ornstein Uhlenbeck process of crop prices in Alberta. And then we use the estimated parameters for Monte Carlo simulation.

## 2 Calibrating Ornstein Uhlenbeck Model

### 2.1 Method

The Ornstein-Uhlenbeck or Vasicek process is   a stochastic process which is stationary, Gaussian, and Markovian. 

Over time, the process tends to drift towards its long-term mean: such a process is called mean-reverting.

The Ornstein-Uhlenbeck or Vasicek process is the unique solution to the following stochastic differential equation:(Stochastic Differential Equation, 2008, p44.)

Ornstein-Uhlenbeck model
$$dX_t =  - \theta_2 X_t dt +\theta_3dW_t,    $$

Vasicek modified it to 

$$dX_t = (\theta_1 - \theta_2X_t)dt +\theta_3dW_t,    X_0 = x_0$$

$(\theta_1 - \theta_2X_t)dt$ is **deterministic part**; $\theta_3dW_t$ is **stochastic part**.


$dW_t$ is the Brownian motion, which follows random normal distributioin $N(0,t)$.

$(\theta_1 - \theta_2X_t)$ is **drift**; $\theta_3$ is **diffusion**.


In finance, the model more often is written as:

$$ dS_t = \theta(\mu - S_t)dt + \sigma dW_t$$



where  $\sigma$ is interpreted as the  instantaneous volatility ,   $\sigma^2/(2\theta)$ is the long term variance;  $\mu$ is the long-run equilibrium value of the process, and $\theta$ is the speed of reversion.

$\theta$ increases the speed at which the system stabilizes around the long term mean $\mu$.

$\sigma$ increases the amount of  randomness entering the system.


"The Ornstein-Uhlenbeck process is one of several approaches used to model (with modifications) interest rates, currency exchange rates, and commodity prices stochastically. The parameter $\mu$ represents the equilibrium or mean value supported by fundamentals; $\sigma$ the degree of volatility around it caused by shocks, and $\theta$ the rate by which these shocks dissipate and the variable reverts towards the mean."(Wikipedia)



For $\theta_2>0$,  $X_t~$ are assumed iid $$N(\frac{\theta_1}{\theta_2},\frac{\theta_3^{2}}{2\theta_2})$$.

$$\frac{\theta_1}{\theta_2} = \mu   ~~\mbox{and}~~ \theta_2 = \theta ~~\mbox{and}~~ \theta_3 = \sigma$$.

For any $t>=0$, the density of the distribution of $X_t$ given $X_0 = x_0$, with mean and variance respectively:

$$m(t,x) = E_\theta(X_t ~|~ X_0 = x_0) = \frac{\theta_1}{\theta_2} + (x_0-\frac{\theta_1}{\theta_2})e^{-\theta_2t}$$

and

$$v(t,x) = Var_\theta(X_t ~|~ X_0 = x_0) =\frac{\theta_3^2 (1- e^{-2 \theta_2 t} )}{2 \theta_2}  $$



### 2.1.0 Data processing

We got monthly price data of crops in Alberta from **Statistics Canada**[http://www5.statcan.gc.ca/subject-sujet/result-resultat?pid=920&id=2024&lang=eng&type=ARRAY&sortType=3&pageNum=0](http://www5.statcan.gc.ca/subject-sujet/result-resultat?pid=920&id=2024&lang=eng&type=ARRAY&sortType=3&pageNum=0).




```r
setwd("E:/Dropbox/book/economics/485/projects/pricemodel/oumodel")

crop.price <- read.csv("abcropprice.csv", header = T, sep = ",")
# set the date format
crop.price[, 1] <- as.Date(crop.price[, 1], format = "%d/%m/%Y")

head(crop.price)
```

```
##         Date Wheat.excluding.durum  Oats Barley Canola Dry.peas
## 1 1985-01-01                 166.0 114.6  127.4  342.0       NA
## 2 1985-02-01                 167.2 116.2  127.4  347.3       NA
## 3 1985-03-01                 166.9 118.4  130.3  350.0       NA
## 4 1985-04-01                 164.8 115.7  127.4  363.8       NA
## 5 1985-05-01                 166.8 105.7  125.6  354.6       NA
## 6 1985-06-01                 166.6 102.2  121.4  349.6       NA
```


For the Ornstein Uhlenbeck Model, long term data has more breaks or jumps, so we only use the short term data which is from stable period and capture the current characteristics of the stochastic process. We decide to use the data after 2007.



```r
price <- crop.price[crop.price$Date > "2006-12-01", ]
head(price)
```

```
##           Date Wheat.excluding.durum  Oats Barley Canola Dry.peas
## 265 2007-01-01                 169.9 149.0  134.5  339.1      181
## 266 2007-02-01                 170.7 152.6  142.7  350.1      189
## 267 2007-03-01                 174.5 163.4  144.8  353.0      205
## 268 2007-04-01                 167.2 165.3  147.2  354.0      230
## 269 2007-05-01                 165.2 165.3  154.4  362.6      239
## 270 2007-06-01                 167.4 167.4  160.8  370.8      249
```



### 2.2  Calibration using Maximum Likelihood estimates

There are many ways to calculate the parameters. We tried two ways. First, we use package "sde.sim" (Stochastic Differential Equation, 2008); 

Recap
$$\frac{\theta_1}{\theta_2} = \mu   ~~\mbox{and}~~ \theta_2 = \theta ~~\mbox{and}~~ \theta_3 = \sigma$$.


#### 2.2.1 "sde.sim" package

Get the necessary packeages and functions

```r
# install.packages('stats4') install.packages('sde')
require(stats4)
```

```
## Loading required package: stats4
```

```r
require(sde)
```

```
## Loading required package: sde
## Loading required package: MASS
## Loading required package: fda
## Loading required package: splines
## Loading required package: Matrix
## 
## Attaching package: 'fda'
## 
## The following object is masked from 'package:graphics':
## 
##     matplot
## 
## Loading required package: zoo
## 
## Attaching package: 'zoo'
## 
## The following objects are masked from 'package:base':
## 
##     as.Date, as.Date.numeric
```

```
## 
## To check the errata corrige of the book, type vignette("sde.errata")
```


```
# Two functions from page 114

 
# Functions dcOU evaluates the conditional density of the process. It generates the density of the x which from cerntain distribution
dcOU<- function (x, t, x0 , theta , log = FALSE ){
        Ex <- theta [1] / theta [2]+( x0 - theta [1] / theta [2]) * exp (- theta [2] *t)
        Vx <- theta [3]^2 *(1- exp (-2* theta [2] *t))/(2* theta [2])
        dnorm (x, mean =Ex , sd = sqrt (Vx), log = log )
}


# Function OU.lik generates the log likelihood function of X for the MLE process. 
OU.lik <- function ( theta1 , theta2 , theta3 ){
        n <- length (X)
        # deltat is the time interval between observations
        dt <- deltat (X)-sum (dcOU (X [2: n], dt , X [1:(n -1)] , c( theta1 , theta2 , theta3 ), log = TRUE ))
}
# The function OU.lik needs as input the three parameters and assumes that sample observations of the process are available in the current R workspace in the ts object X.

```

##### 2.2.1.1 Calculate the parameters for wheat

Mle is very sensitive for the start value. We choose theta1 =20, theta2 =0.1 , theta3 =1.
We cannot get it converg. It seems not working.
(The result value which the process converges to might not be the unique solution. )
```
wheat <- as.vector(price[,2]) # convert to vector 
head(wheat)
X <- ts(data=wheat)
mle(OU.lik , start = list ( theta1 =20, theta2= 0.1 , theta3 =10) , method ="L-BFGS-B", lower =c(-Inf ,0 ,0)) -> fitwheat
summary ( fitwheat )
coef(fitwheat)[[1]]/coef(fitwheat)[[2]]
```



#### 2.2.2 Function "ouFit.ML"


Second, we manually use the function from "Calibrating the Ornstein-Uhlenbeck (Vasicek) model"[http://www.sitmo.com/article/calibrating-the-ornstein-uhlenbeck-model/](http://www.sitmo.com/article/calibrating-the-ornstein-uhlenbeck-model/)

Function "ouFit.ML" is saved in file "oufit.R".


```r
# function for Calibration using Maximum Likelihood estimates
ouFit.ML = function(spread) {
    n = length(spread)
    delta = 1  # delta 
    Sx = sum(spread[1:n - 1])
    Sy = sum(spread[2:n])
    Sxx = sum((spread[1:n - 1])^2)
    Syy = sum((spread[2:n])^2)
    Sxy = sum((spread[1:n - 1]) * (spread[2:n]))
    mu = (Sy * Sxx - Sx * Sxy)/((n - 1) * (Sxx - Sxy) - (Sx^2 - Sx * Sy))
    theta = -log((Sxy - mu * Sx - mu * Sy + (n - 1) * mu^2)/(Sxx - 2 * mu * 
        Sx + (n - 1) * mu^2))/delta
    a = exp(-theta * delta)
    sigmah2 = (Syy - 2 * a * Sxy + a^2 * Sxx - 2 * mu * (1 - a) * (Sy - a * 
        Sx) + (n - 1) * mu^2 * (1 - a)^2)/(n - 1)
    sigma = sqrt((sigmah2) * 2 * theta/(1 - a^2))
    theta = list(theta = theta, mu = mu, sigma = sigma, sigmah2 = sigmah2)
    return(theta)
}

```


##### 2.2.2.1 Calculate the parameters for wheat


```r
# source('oufit.R') # not necessary

wheat <- as.vector(price[, 2])  # convert to vector 
head(wheat)
```

```
## [1] 169.9 170.7 174.5 167.2 165.2 167.4
```

```r
plot(wheat, type = "l")
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5.png) 

```r
ouFit.ML(wheat)
```

```
## $theta
## [1] 0.1198
## 
## $mu
## [1] 249.9
## 
## $sigma
## [1] 21.99
## 
## $sigmah2
## [1] 430.1
```


##### 2.2.2.2 Calculate the parameters for oats

```r
oats <- as.vector(na.omit(price[, 3]))  # get rid of n/a and convert to vector 
head(oats)
```

```
## [1] 149.0 152.6 163.4 165.3 165.3 167.4
```

```r
plot(oats, type = "l")
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6.png) 

```r
ouFit.ML(oats)
```

```
## $theta
## [1] 0.08256
## 
## $mu
## [1] 186.2
## 
## $sigma
## [1] 9.471
## 
## $sigmah2
## [1] 82.68
```


##### 2.2.2.3 Calculate the parameters for barley

```r
barley <- as.vector(na.omit(price[, 4]))  # get rid of n/a and convert to vector 
head(barley)
```

```
## [1] 134.5 142.7 144.8 147.2 154.4 160.8
```

```r
plot(barley, type = "l")
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7.png) 

```r
ouFit.ML(barley)
```

```
## $theta
## [1] 0.04971
## 
## $mu
## [1] 188.8
## 
## $sigma
## [1] 11.29
## 
## $sigmah2
## [1] 121.3
```


##### 2.2.2.4 Calculate the parameters for canola

```r
canola <- as.vector(na.omit(price[, 5]))  # get rid of n/a and convert to vector 
head(canola)
```

```
## [1] 339.1 350.1 353.0 354.0 362.6 370.8
```

```r
plot(canola, type = "l")
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8.png) 

```r
ouFit.ML(canola)
```

```
## $theta
## [1] 0.04689
## 
## $mu
## [1] 507.7
## 
## $sigma
## [1] 18.46
## 
## $sigmah2
## [1] 325.5
```



##### 2.2.2.5 Calculate the parameters for Dry peas

```r
dry.peas <- as.vector(na.omit(price[, 6]))  # get rid of n/a and convert to vector 
head(dry.peas)
```

```
## [1] 181 189 205 230 239 249
```

```r
plot(dry.peas, type = "l")
```

![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9.png) 

```r
ouFit.ML(dry.peas)
```

```
## $theta
## [1] 0.05435
## 
## $mu
## [1] 279.5
## 
## $sigma
## [1] 14.08
## 
## $sigmah2
## [1] 187.7
```




## 3 Generate the price vectors/matrices by simulation

There are two packages we can use to generate Vasicek process. We try both and get different results. It may be because both use the random method to simulation.

We find the "yuima" is better than "sde.sim"

### 3.1 "yuima" packages

The simulation follows the method mentioned on "1 4 Vasicek Model "[https://www.youtube.com/watch?v=5BpOYPNxWsA&index=4&list=PLRj8HuwfWZEsXW2pzAwAWYC8EZbD2ieOq](https://www.youtube.com/watch?v=5BpOYPNxWsA&index=4&list=PLRj8HuwfWZEsXW2pzAwAWYC8EZbD2ieOq)

Recall "yuima" package uses this notation:

$$ dS_t = \theta(\mu - S_t)dt + \sigma dW_t$$

The delta of t  
$$dt=(1/n)=(1/100)$$

by default. We can change it by using $grid=setSampling(Terminal=1,n=1000)$. 

#### 3.1.1  Generate  1000 random samples  with 1000 time periods for wheats price using the parameters from estimate


One simulation, different time different result. And time period still is 1000 ("yuimu" default: delta = 1/100, we change it to 1/1000).

```r
require(yuima)
```

```
## Loading required package: yuima
## Loading required package: expm
## 
## Attaching package: 'expm'
## 
## The following object is masked from 'package:Matrix':
## 
##     expm
## 
## ############################################
## This is YUIMA Project package.
## Check for the latest development version at:
## http://R-Forge.R-Project.org/projects/yuima
## ############################################
## 
## Attaching package: 'yuima'
## 
## The following object is masked from 'package:stats':
## 
##     simulate
```

```r
grid = setSampling(Terminal = 1, n = 1000)
m1 = setModel(drift = "theta*(mu-x)", diffusion = "sigma", state.var = "x", 
    time.var = "t", solve.var = "x", xinit = ouFit.ML(wheat)[[2]])
Xwheat = simulate(m1, true.param = list(mu = ouFit.ML(wheat)[[2]], sigma = ouFit.ML(wheat)[[3]], 
    theta = ouFit.ML(wheat)[[1]]), sampling = grid)
plot(Xwheat)
```

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-10.png) 



Increasing the number of simulation to 1000, and time period still is 100 ("yuimu" default: delta = 1/100). Plot the mean of the 1000 silumtion result, different time , very similar result.

```r
simnum = 1000
dist = c(0.31, 0.52, 0.6, 0.7, 0.95)
newsim = function(i) {
    simulate(m1, true.param = list(mu = ouFit.ML(wheat)[[2]], sigma = ouFit.ML(wheat)[[3]], 
        theta = ouFit.ML(wheat)[[1]]))@data@original.data
}
# newsim(1) simulation 1000 times, each time there are 100 time periods
sim = sapply(1:simnum, function(x) newsim(x))
# transfor to time seires format.
m2 = t(sim)
mwheat <- apply(m2, 2, mean)

# plot the mean of the 1000 time simulation for the 100 time periods
plot(mwheat, type = "l")
```

![plot of chunk unnamed-chunk-11](figure/unnamed-chunk-11.png) 

```r

# find out the quantile to decribe the distribution
tile = sapply(1:100, function(x) quantile(m2[, x], dist))
tile
```

```
##      [,1]  [,2]  [,3]  [,4]  [,5]  [,6]  [,7]  [,8]  [,9] [,10] [,11]
## 31% 249.9 248.8 248.3 248.0 247.6 247.4 247.5 247.2 246.9 246.7 246.6
## 52% 249.9 250.0 250.0 250.0 250.1 250.1 250.3 249.9 250.1 250.6 250.2
## 60% 249.9 250.6 250.7 250.7 250.8 251.1 251.3 251.3 251.3 251.8 251.6
## 70% 249.9 251.2 251.5 251.7 251.9 252.3 252.5 252.6 253.1 253.4 253.3
## 95% 249.9 253.4 254.9 256.1 257.2 257.8 258.7 260.0 260.4 260.7 261.4
##     [,12] [,13] [,14] [,15] [,16] [,17] [,18] [,19] [,20] [,21] [,22]
## 31% 246.4 246.3 245.9 245.6 245.7 245.7 245.5 245.5 245.4 245.4 245.2
## 52% 250.2 250.1 250.1 250.4 250.5 250.6 250.4 250.5 250.5 250.4 250.6
## 60% 251.8 251.7 251.7 251.9 252.1 252.4 252.3 252.5 252.3 252.4 252.5
## 70% 253.7 253.7 254.0 254.3 254.5 254.5 254.9 254.9 255.1 255.0 255.1
## 95% 261.8 262.0 262.3 262.6 263.5 262.9 264.3 265.1 265.1 266.3 266.4
##     [,23] [,24] [,25] [,26] [,27] [,28] [,29] [,30] [,31] [,32] [,33]
## 31% 245.1 244.8 244.5 244.6 244.6 244.5 244.5 244.7 244.3 244.6 244.1
## 52% 250.4 250.6 250.7 251.0 250.7 250.7 250.6 250.6 250.6 250.4 250.7
## 60% 252.6 252.9 252.9 253.1 253.0 252.9 252.9 253.0 253.1 253.2 253.2
## 70% 255.7 256.0 255.7 255.8 255.6 256.2 255.9 256.0 256.3 256.2 256.6
## 95% 266.6 267.0 267.8 268.0 268.5 268.6 268.5 268.9 269.6 269.0 269.0
##     [,34] [,35] [,36] [,37] [,38] [,39] [,40] [,41] [,42] [,43] [,44]
## 31% 244.1 244.2 244.2 243.8 243.6 244.2 244.0 243.8 243.5 243.7 243.4
## 52% 250.6 250.6 250.8 250.9 251.2 251.2 250.7 251.0 250.7 250.7 250.8
## 60% 253.4 253.7 253.7 253.4 253.7 253.7 253.9 253.7 253.5 254.1 253.7
## 70% 257.2 257.0 256.7 256.8 257.0 256.9 256.9 257.1 257.3 257.6 257.6
## 95% 269.6 269.8 270.6 271.0 270.9 271.0 272.3 272.1 273.1 273.2 272.6
##     [,45] [,46] [,47] [,48] [,49] [,50] [,51] [,52] [,53] [,54] [,55]
## 31% 243.5 243.4 243.0 243.0 242.7 243.1 242.7 242.7 242.7 242.5 242.6
## 52% 250.6 250.7 250.4 250.8 250.5 250.4 250.4 250.3 250.4 250.1 250.2
## 60% 253.8 253.9 254.0 254.2 254.0 254.0 253.9 254.1 253.7 253.5 253.5
## 70% 257.5 257.3 257.8 257.8 257.6 257.7 258.2 258.4 258.4 258.0 258.2
## 95% 273.4 272.8 273.7 273.8 274.5 274.9 274.8 274.6 274.7 274.8 275.1
##     [,56] [,57] [,58] [,59] [,60] [,61] [,62] [,63] [,64] [,65] [,66]
## 31% 242.8 242.8 242.6 242.4 241.9 241.7 241.8 241.4 241.4 241.1 240.4
## 52% 250.1 250.3 250.6 250.8 250.9 250.9 250.7 250.9 250.5 250.5 250.0
## 60% 253.9 253.7 253.8 254.0 254.0 253.9 253.8 253.5 253.8 253.6 253.7
## 70% 258.0 258.3 258.4 258.3 257.9 258.3 258.4 258.5 258.4 258.3 258.5
## 95% 275.1 275.7 275.8 276.5 276.7 277.7 277.9 277.6 278.4 277.4 277.3
##     [,67] [,68] [,69] [,70] [,71] [,72] [,73] [,74] [,75] [,76] [,77]
## 31% 240.9 241.0 240.7 240.5 240.6 240.3 240.6 240.6 240.7 240.2 240.5
## 52% 250.0 250.1 249.9 250.3 250.0 250.3 250.3 250.2 250.5 250.9 250.4
## 60% 253.3 252.9 253.3 253.6 253.4 253.5 254.0 254.1 254.0 253.9 254.4
## 70% 258.8 258.9 258.6 258.7 259.1 259.2 259.2 258.7 259.0 258.7 259.5
## 95% 277.4 278.0 279.0 278.8 278.7 278.6 279.2 278.9 279.1 279.0 279.7
##     [,78] [,79] [,80] [,81] [,82] [,83] [,84] [,85] [,86] [,87] [,88]
## 31% 240.1 240.6 240.5 240.3 240.2 240.2 240.1 239.9 240.0 239.8 240.2
## 52% 251.0 250.5 250.9 250.4 250.9 251.0 250.8 250.8 250.4 250.8 250.6
## 60% 254.7 254.3 254.6 254.5 254.8 255.1 255.1 255.0 255.5 255.0 254.9
## 70% 259.9 259.6 259.8 259.6 260.2 259.8 259.9 259.7 259.5 260.0 259.9
## 95% 280.1 280.4 280.5 280.3 280.1 281.0 280.7 281.3 282.3 281.7 281.6
##     [,89] [,90] [,91] [,92] [,93] [,94] [,95] [,96] [,97] [,98] [,99]
## 31% 239.6 239.5 239.8 239.5 239.7 239.0 239.2 239.4 238.7 239.2 239.0
## 52% 250.8 250.6 250.9 250.9 250.9 250.6 250.6 251.3 250.7 250.8 250.9
## 60% 254.9 254.7 255.3 255.7 255.6 255.5 255.7 254.8 254.8 255.5 255.0
## 70% 259.4 260.1 260.2 260.5 260.5 260.9 260.6 260.8 261.0 261.2 260.6
## 95% 282.1 282.4 282.3 282.7 283.5 283.4 283.5 283.2 283.1 282.9 282.9
##     [,100]
## 31%  238.8
## 52%  251.0
## 60%  255.4
## 70%  260.8
## 95%  282.4
```


#### 3.1.2  Generate  1000 random samples  with 100 time periods for oats price using the parameters from estimate


One simulation, different time different result. And time period still is 100 ("yuimu" default: delta = 1/100).

```r
require(yuima)
# initial value is 188
m1 = setModel(drift = "theta*(mu-x)", diffusion = "sigma", state.var = "x", 
    time.var = "t", solve.var = "x", xinit = ouFit.ML(oats)[[2]])
Xoats = simulate(m1, true.param = list(mu = ouFit.ML(oats)[[2]], sigma = ouFit.ML(oats)[[3]], 
    theta = ouFit.ML(oats)[[1]]))
plot(Xoats)
```

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-12.png) 



Increasing the number of simulation to 1000, and time period still is 100 ("yuimu" default: delta = 1/100). Plot the mean of the 1000 silumtion result, different time , very similar result.

```r
simnum = 1000
# specific qunatile (which we can pick any another quantile)
dist = c(0.31, 0.52, 0.6, 0.7, 0.95)
newsim = function(i) {
    simulate(m1, true.param = list(mu = ouFit.ML(oats)[[2]], sigma = ouFit.ML(oats)[[3]], 
        theta = ouFit.ML(oats)[[1]]))@data@original.data
}
# newsim(1) simulation 1000 times, each time there are 100 time periods
sim = sapply(1:simnum, function(x) newsim(x))
# transfor to time seires format.
m2 = t(sim)
moats <- apply(m2, 2, mean)

# plot the mean of the 1000 time simulation for the 100 time periods
plot(moats, type = "l")
```

![plot of chunk unnamed-chunk-13](figure/unnamed-chunk-13.png) 

```r

# find out the quantile to decribe the distribution
tile = sapply(1:100, function(x) quantile(m2[, x], dist))
tile
```

```
##      [,1]  [,2]  [,3]  [,4]  [,5]  [,6]  [,7]  [,8]  [,9] [,10] [,11]
## 31% 186.2 185.8 185.7 185.5 185.4 185.4 185.3 185.1 185.1 185.0 184.9
## 52% 186.2 186.3 186.4 186.4 186.4 186.5 186.5 186.5 186.4 186.4 186.5
## 60% 186.2 186.5 186.7 186.8 186.8 186.9 187.0 187.0 186.9 187.0 187.1
## 70% 186.2 186.8 187.1 187.2 187.4 187.4 187.6 187.6 187.5 187.6 187.8
## 95% 186.2 187.8 188.6 189.0 189.5 189.8 190.1 190.5 190.8 191.0 191.1
##     [,12] [,13] [,14] [,15] [,16] [,17] [,18] [,19] [,20] [,21] [,22]
## 31% 184.7 184.7 184.7 184.7 184.6 184.6 184.4 184.4 184.3 184.3 184.4
## 52% 186.4 186.5 186.6 186.5 186.7 186.7 186.6 186.5 186.6 186.6 186.6
## 60% 187.0 187.2 187.2 187.2 187.3 187.3 187.3 187.3 187.3 187.4 187.5
## 70% 187.9 187.9 188.0 188.2 188.1 188.3 188.3 188.3 188.4 188.5 188.6
## 95% 191.5 191.5 191.9 192.3 192.3 192.4 192.3 192.7 192.6 193.0 193.3
##     [,23] [,24] [,25] [,26] [,27] [,28] [,29] [,30] [,31] [,32] [,33]
## 31% 184.3 184.2 184.2 184.0 184.3 184.2 184.3 184.3 184.2 184.1 183.9
## 52% 186.6 186.5 186.4 186.6 186.5 186.5 186.6 186.7 186.7 186.7 186.6
## 60% 187.5 187.3 187.5 187.5 187.5 187.6 187.7 187.8 187.6 187.6 187.6
## 70% 188.6 188.6 188.7 188.7 188.8 188.9 188.9 188.8 189.0 188.8 189.0
## 95% 193.5 193.7 193.8 193.6 194.1 194.1 194.2 194.2 194.6 195.0 195.4
##     [,34] [,35] [,36] [,37] [,38] [,39] [,40] [,41] [,42] [,43] [,44]
## 31% 183.9 183.8 183.8 183.6 183.5 183.5 183.4 183.3 183.2 183.2 183.2
## 52% 186.6 186.6 186.7 186.6 186.8 186.7 186.6 186.6 186.6 186.7 186.6
## 60% 187.6 187.6 187.7 187.7 187.8 187.9 187.9 188.0 187.9 187.8 188.1
## 70% 189.0 189.4 189.3 189.3 189.3 189.4 189.6 189.7 189.6 189.7 189.9
## 95% 195.2 195.6 195.6 196.0 195.7 195.9 195.9 196.4 196.4 196.5 196.8
##     [,45] [,46] [,47] [,48] [,49] [,50] [,51] [,52] [,53] [,54] [,55]
## 31% 183.1 183.1 183.0 182.9 182.9 182.9 183.0 183.0 182.7 182.9 182.9
## 52% 186.7 186.7 186.9 186.8 186.7 186.8 186.6 186.7 186.7 186.7 186.9
## 60% 188.0 188.2 188.1 188.2 188.1 188.2 188.3 188.2 188.3 188.1 188.3
## 70% 190.0 189.9 189.9 189.8 190.0 189.9 190.0 189.9 190.0 190.0 190.2
## 95% 197.0 196.9 197.3 197.4 197.4 197.7 197.4 197.7 197.8 198.0 198.2
##     [,56] [,57] [,58] [,59] [,60] [,61] [,62] [,63] [,64] [,65] [,66]
## 31% 182.9 183.0 182.9 183.0 182.9 182.9 182.5 182.7 182.8 182.7 182.7
## 52% 186.9 186.8 186.9 186.8 186.8 186.7 186.8 186.7 186.7 186.6 186.8
## 60% 188.2 188.4 188.4 188.3 188.4 188.4 188.4 188.4 188.3 188.5 188.5
## 70% 190.2 190.2 190.2 190.4 190.4 190.0 190.2 190.6 190.6 190.6 190.6
## 95% 198.7 198.8 198.8 199.0 198.6 198.7 198.8 198.9 198.7 199.0 199.1
##     [,67] [,68] [,69] [,70] [,71] [,72] [,73] [,74] [,75] [,76] [,77]
## 31% 182.4 182.4 182.5 182.6 182.6 182.6 182.4 182.4 182.3 182.3 182.4
## 52% 186.9 186.8 186.7 186.8 186.9 186.7 186.8 186.7 186.6 186.7 186.5
## 60% 188.5 188.5 188.6 188.4 188.5 188.3 188.3 188.2 188.2 188.1 188.2
## 70% 190.5 190.5 190.7 190.6 190.5 190.4 190.5 190.5 190.5 190.5 190.5
## 95% 199.3 199.1 199.1 199.1 199.0 198.9 199.3 199.1 199.1 199.3 199.4
##     [,78] [,79] [,80] [,81] [,82] [,83] [,84] [,85] [,86] [,87] [,88]
## 31% 182.2 182.2 182.2 182.1 182.3 182.3 182.2 182.2 182.2 182.3 182.1
## 52% 186.4 186.6 186.4 186.4 186.7 186.4 186.5 186.5 186.6 186.6 186.6
## 60% 188.0 188.1 188.1 188.2 188.3 188.2 188.2 188.4 188.4 188.4 188.4
## 70% 190.6 190.5 190.7 190.7 190.7 190.7 190.6 190.8 190.8 190.9 191.0
## 95% 199.8 199.8 199.7 200.2 199.6 199.9 200.0 199.9 200.0 200.2 200.5
##     [,89] [,90] [,91] [,92] [,93] [,94] [,95] [,96] [,97] [,98] [,99]
## 31% 182.2 182.3 182.2 181.9 181.9 181.9 181.8 181.8 181.9 181.9 182.0
## 52% 186.6 186.6 186.6 186.4 186.6 186.5 186.5 186.6 186.6 186.6 186.4
## 60% 188.4 188.4 188.6 188.5 188.5 188.5 188.5 188.4 188.3 188.5 188.7
## 70% 190.9 191.0 190.9 190.9 190.8 190.8 190.8 191.0 191.1 191.1 191.1
## 95% 200.5 200.2 200.1 200.2 200.0 200.3 200.4 201.0 200.8 201.0 200.6
##     [,100]
## 31%  181.8
## 52%  186.4
## 60%  188.6
## 70%  191.1
## 95%  200.9
```



#### 3.1.3  Generate  1000 random samples  with 100 time periods for barley price using the parameters from estimate


One simulation, different time different result. And time period still is 100 ("yuimu" default: delta = 1/100).

```r
require(yuima)
# initial value is 188
m1 = setModel(drift = "theta*(mu-x)", diffusion = "sigma", state.var = "x", 
    time.var = "t", solve.var = "x", xinit = ouFit.ML(barley)[[2]])
Xbarley = simulate(m1, true.param = list(mu = ouFit.ML(barley)[[2]], sigma = ouFit.ML(barley)[[3]], 
    theta = ouFit.ML(barley)[[1]]))
plot(Xbarley)
```

![plot of chunk unnamed-chunk-14](figure/unnamed-chunk-14.png) 



Increasing the number of simulation to 1000, and time period still is 100 ("yuimu" default: delta = 1/100). Plot the mean of the 1000 silumtion result, different time , very similar result.

```r
simnum = 1000
# specific qunatile (which we can pick any another quantile)
dist = c(0.31, 0.52, 0.6, 0.7, 0.95)
newsim = function(i) {
    simulate(m1, true.param = list(mu = ouFit.ML(barley)[[2]], sigma = ouFit.ML(barley)[[3]], 
        theta = ouFit.ML(barley)[[1]]))@data@original.data
}
# newsim(1) simulation 1000 times, each time there are 100 time periods
sim = sapply(1:simnum, function(x) newsim(x))
# transfor to time seires format.
m2 = t(sim)
mbarley <- apply(m2, 2, mean)

# plot the mean of the 1000 time simulation for the 100 time periods
plot(mbarley, type = "l")
```

![plot of chunk unnamed-chunk-15](figure/unnamed-chunk-15.png) 

```r

# find out the quantile to decribe the distribution
tile = sapply(1:100, function(x) quantile(m2[, x], dist))
tile
```

```
##      [,1]  [,2]  [,3]  [,4]  [,5]  [,6]  [,7]  [,8]  [,9] [,10] [,11]
## 31% 188.8 188.3 188.1 187.9 187.6 187.7 187.4 187.4 187.5 187.4 187.4
## 52% 188.8 189.0 188.9 188.9 189.0 188.9 189.0 189.1 189.2 189.3 189.2
## 60% 188.8 189.2 189.3 189.4 189.5 189.4 189.7 189.7 189.8 190.1 190.0
## 70% 188.8 189.5 189.6 189.9 190.1 190.2 190.5 190.6 190.7 190.9 191.0
## 95% 188.8 190.7 191.5 192.0 192.5 192.8 193.3 193.6 194.0 194.1 194.7
##     [,12] [,13] [,14] [,15] [,16] [,17] [,18] [,19] [,20] [,21] [,22]
## 31% 187.3 187.3 187.2 187.0 186.9 186.7 186.5 186.5 186.4 186.2 186.2
## 52% 189.2 189.1 189.3 189.4 189.3 189.5 189.5 189.3 189.1 189.1 189.1
## 60% 189.9 190.1 190.2 190.2 190.3 190.3 190.4 190.4 190.2 190.4 190.1
## 70% 191.0 191.2 191.3 191.3 191.5 191.6 191.6 191.8 191.7 191.9 191.9
## 95% 195.2 195.2 195.6 196.3 196.5 196.6 196.8 197.2 197.2 197.5 197.7
##     [,23] [,24] [,25] [,26] [,27] [,28] [,29] [,30] [,31] [,32] [,33]
## 31% 185.9 185.8 185.9 185.9 185.8 185.7 185.5 185.4 185.5 185.5 185.4
## 52% 189.0 188.9 188.7 188.6 188.7 188.7 188.7 188.7 188.8 188.9 189.0
## 60% 190.0 190.0 189.8 189.7 190.0 190.1 190.2 190.0 190.2 190.2 190.1
## 70% 191.7 191.5 191.6 191.7 191.7 191.7 191.6 191.8 191.9 192.1 192.3
## 95% 197.9 197.9 198.3 198.6 198.6 199.0 199.3 199.2 199.3 199.3 199.5
##     [,34] [,35] [,36] [,37] [,38] [,39] [,40] [,41] [,42] [,43] [,44]
## 31% 185.2 185.5 185.3 185.3 185.3 185.3 185.1 185.2 185.2 185.2 185.0
## 52% 188.9 188.7 189.0 188.8 188.9 188.8 189.1 189.1 189.0 188.9 189.1
## 60% 190.2 190.2 190.3 190.5 190.5 190.4 190.4 190.5 190.7 190.8 190.7
## 70% 192.1 192.1 192.4 192.2 192.4 192.4 192.5 192.7 192.8 193.0 192.9
## 95% 199.6 199.8 200.1 200.3 200.5 200.7 200.7 200.6 200.7 200.7 201.2
##     [,45] [,46] [,47] [,48] [,49] [,50] [,51] [,52] [,53] [,54] [,55]
## 31% 185.1 184.9 184.9 184.9 184.8 184.9 184.9 184.8 184.7 184.6 184.4
## 52% 188.9 188.9 188.9 189.0 189.3 189.0 189.0 188.9 188.8 188.8 189.1
## 60% 190.6 190.9 190.9 190.9 191.0 190.8 190.8 190.9 190.8 190.6 190.7
## 70% 192.9 193.0 192.9 192.8 192.7 192.8 193.1 192.9 193.4 193.3 193.1
## 95% 201.1 201.2 201.2 201.7 202.2 202.0 201.8 202.1 202.3 202.6 202.4
##     [,56] [,57] [,58] [,59] [,60] [,61] [,62] [,63] [,64] [,65] [,66]
## 31% 184.5 184.6 184.5 184.5 184.5 184.4 184.4 184.2 184.3 184.3 184.2
## 52% 188.9 188.9 188.9 189.0 189.2 189.3 189.0 188.9 189.0 189.1 189.2
## 60% 190.7 190.7 190.6 190.8 190.9 190.7 190.8 190.6 190.7 190.6 190.9
## 70% 193.0 192.8 193.2 193.1 193.1 192.9 193.1 193.1 193.2 193.4 193.5
## 95% 202.2 202.6 202.9 203.0 203.1 203.4 203.2 203.3 203.2 203.8 203.6
##     [,67] [,68] [,69] [,70] [,71] [,72] [,73] [,74] [,75] [,76] [,77]
## 31% 184.1 184.3 184.3 184.2 184.4 184.3 184.2 184.2 184.1 183.9 184.0
## 52% 189.3 189.1 189.2 189.3 189.3 189.3 189.0 188.9 188.8 189.0 188.9
## 60% 190.7 190.9 190.9 191.0 191.0 190.9 190.6 190.7 190.9 190.7 190.9
## 70% 193.2 193.3 193.4 193.3 193.4 193.5 193.4 193.4 193.5 193.4 193.4
## 95% 203.9 203.6 204.0 204.3 204.2 204.5 204.5 204.7 205.0 204.4 204.6
##     [,78] [,79] [,80] [,81] [,82] [,83] [,84] [,85] [,86] [,87] [,88]
## 31% 183.9 183.9 184.0 183.9 184.0 183.7 183.8 183.8 183.6 183.9 183.6
## 52% 188.7 188.6 188.8 189.0 188.8 188.7 188.9 188.9 189.0 189.2 189.0
## 60% 190.9 190.9 190.7 191.0 191.1 191.1 191.2 191.1 191.3 191.2 191.1
## 70% 193.5 193.7 193.7 193.7 193.8 193.5 193.6 193.6 193.5 193.6 193.8
## 95% 204.7 204.8 204.5 204.4 205.1 204.7 205.3 205.6 205.8 205.4 206.2
##     [,89] [,90] [,91] [,92] [,93] [,94] [,95] [,96] [,97] [,98] [,99]
## 31% 183.5 183.2 183.3 183.4 183.3 183.2 183.2 183.0 183.1 183.2 183.1
## 52% 189.0 189.2 189.1 189.2 189.2 189.3 189.1 188.9 188.8 189.0 188.6
## 60% 191.3 191.2 191.2 191.0 191.2 191.0 191.1 191.1 191.0 191.1 191.1
## 70% 193.7 194.0 194.1 193.8 193.9 194.2 194.1 194.1 194.1 194.3 194.4
## 95% 206.1 206.2 206.0 206.2 205.8 206.4 206.3 206.7 206.8 206.8 206.4
##     [,100]
## 31%  183.2
## 52%  188.7
## 60%  190.9
## 70%  194.5
## 95%  206.6
```



#### 3.1.4  Generate  1000 random samples  with 100 time periods for canola price using the parameters from estimate


One simulation, different time different result. And time period still is 100 ("yuimu" default: delta = 1/100).

```r
require(yuima)
# initial value is 188
m1 = setModel(drift = "theta*(mu-x)", diffusion = "sigma", state.var = "x", 
    time.var = "t", solve.var = "x", xinit = ouFit.ML(canola)[[2]])
Xcanola = simulate(m1, true.param = list(mu = ouFit.ML(canola)[[2]], sigma = ouFit.ML(canola)[[3]], 
    theta = ouFit.ML(canola)[[1]]))
plot(Xcanola)
```

![plot of chunk unnamed-chunk-16](figure/unnamed-chunk-16.png) 



Increasing the number of simulation to 1000, and time period still is 100 ("yuimu" default: delta = 1/100). Plot the mean of the 1000 silumtion result, different time , very similar result.

```r
simnum = 1000
# specific qunatile (which we can pick any another quantile)
dist = c(0.31, 0.52, 0.6, 0.7, 0.95)
newsim = function(i) {
    simulate(m1, true.param = list(mu = ouFit.ML(canola)[[2]], sigma = ouFit.ML(canola)[[3]], 
        theta = ouFit.ML(canola)[[1]]))@data@original.data
}
# newsim(1) simulation 1000 times, each time there are 100 time periods
sim = sapply(1:simnum, function(x) newsim(x))
# transfor to time seires format.
m2 = t(sim)
mcanola <- apply(m2, 2, mean)

# plot the mean of the 1000 time simulation for the 100 time periods
plot(mcanola, type = "l")
```

![plot of chunk unnamed-chunk-17](figure/unnamed-chunk-17.png) 

```r

# find out the quantile to decribe the distribution
tile = sapply(1:100, function(x) quantile(m2[, x], dist))
tile
```

```
##      [,1]  [,2]  [,3]  [,4]  [,5]  [,6]  [,7]  [,8]  [,9] [,10] [,11]
## 31% 507.7 506.8 506.5 506.2 505.9 505.7 505.5 505.5 505.3 505.3 505.0
## 52% 507.7 507.9 507.9 507.9 507.8 507.9 508.0 508.0 508.1 508.1 508.2
## 60% 507.7 508.2 508.4 508.6 508.5 508.6 508.7 508.9 509.0 509.3 509.4
## 70% 507.7 508.7 509.1 509.3 509.6 509.7 509.9 510.0 510.3 510.9 511.1
## 95% 507.7 510.9 512.0 513.0 513.9 514.3 515.0 515.1 515.9 516.6 516.8
##     [,12] [,13] [,14] [,15] [,16] [,17] [,18] [,19] [,20] [,21] [,22]
## 31% 505.0 504.7 504.6 504.5 504.3 504.3 504.4 504.4 504.4 504.1 503.6
## 52% 508.3 508.4 508.5 508.6 508.6 508.3 508.4 508.5 508.7 508.6 508.5
## 60% 509.4 509.6 509.9 510.0 509.9 509.8 509.8 510.2 510.2 510.3 510.5
## 70% 511.3 511.6 511.2 511.6 511.9 511.8 512.0 512.4 512.4 512.6 512.5
## 95% 517.5 517.9 518.2 518.8 518.9 519.3 519.8 520.4 520.7 521.0 521.2
##     [,23] [,24] [,25] [,26] [,27] [,28] [,29] [,30] [,31] [,32] [,33]
## 31% 503.6 503.3 503.4 503.4 503.4 503.3 503.1 503.1 502.9 502.6 502.7
## 52% 508.5 508.6 508.6 508.4 508.5 508.7 508.7 508.8 508.8 508.7 508.5
## 60% 510.4 510.4 510.3 510.2 510.3 510.4 510.8 510.8 510.6 510.6 510.8
## 70% 512.5 512.9 513.1 512.9 513.1 513.2 513.3 513.3 513.4 513.3 513.7
## 95% 521.9 522.0 522.5 523.1 523.1 523.4 523.5 524.5 525.2 524.2 524.6
##     [,34] [,35] [,36] [,37] [,38] [,39] [,40] [,41] [,42] [,43] [,44]
## 31% 502.4 502.4 502.5 502.3 502.3 502.0 502.2 502.3 502.1 502.1 501.9
## 52% 508.4 508.6 508.6 508.7 508.9 508.9 509.0 509.0 509.1 509.1 509.0
## 60% 510.6 510.9 510.8 511.1 511.1 511.4 511.2 511.2 511.3 511.2 511.4
## 70% 513.9 514.0 514.1 513.7 514.1 514.2 514.3 514.2 514.2 514.3 514.6
## 95% 525.0 525.0 525.7 525.9 526.9 527.0 527.8 527.3 527.5 528.0 528.1
##     [,45] [,46] [,47] [,48] [,49] [,50] [,51] [,52] [,53] [,54] [,55]
## 31% 501.6 501.7 501.5 501.9 501.9 502.3 502.0 501.5 501.4 501.5 501.1
## 52% 508.7 508.9 508.6 508.4 508.4 508.3 508.2 508.2 508.5 507.8 508.1
## 60% 511.5 511.3 511.3 511.5 511.1 511.1 511.3 511.5 511.5 511.8 511.5
## 70% 514.8 514.4 514.5 514.4 514.3 514.6 514.9 514.9 514.9 515.0 515.0
## 95% 528.1 528.6 528.9 529.4 530.2 530.3 529.8 530.2 530.3 530.2 530.5
##     [,56] [,57] [,58] [,59] [,60] [,61] [,62] [,63] [,64] [,65] [,66]
## 31% 501.3 501.4 501.3 500.7 500.7 500.9 500.9 500.7 500.6 500.6 500.4
## 52% 508.1 508.1 508.4 508.4 508.6 508.7 508.4 508.5 508.8 508.8 508.6
## 60% 510.9 510.9 511.3 511.2 510.9 510.9 511.1 511.3 511.4 511.1 511.1
## 70% 514.9 515.0 514.6 514.6 514.5 514.8 514.8 514.6 514.5 514.9 514.7
## 95% 530.5 531.3 531.1 531.0 532.0 532.1 532.6 533.1 533.4 533.2 533.3
##     [,67] [,68] [,69] [,70] [,71] [,72] [,73] [,74] [,75] [,76] [,77]
## 31% 500.5 500.0 500.4 500.3 500.3 499.8 499.8 500.2 499.5 499.9 499.8
## 52% 508.6 508.4 508.2 508.0 508.5 508.1 508.5 508.3 508.5 508.6 508.5
## 60% 510.8 511.1 510.7 511.0 511.2 511.0 511.4 511.7 511.4 511.8 511.6
## 70% 514.6 515.2 515.4 515.1 515.5 515.6 515.4 515.4 515.6 516.0 516.0
## 95% 533.0 533.1 533.6 533.5 533.4 533.8 534.0 535.0 534.7 534.6 535.0
##     [,78] [,79] [,80] [,81] [,82] [,83] [,84] [,85] [,86] [,87] [,88]
## 31% 499.6 499.4 499.7 499.2 499.8 500.0 499.6 499.7 499.4 499.6 499.5
## 52% 508.6 509.0 508.8 508.6 508.8 508.9 508.9 509.1 509.0 509.1 509.1
## 60% 512.1 512.1 511.9 511.9 511.9 512.1 512.2 512.2 512.2 512.3 512.2
## 70% 515.7 516.0 515.8 516.1 516.1 516.3 516.1 516.3 516.5 517.0 517.2
## 95% 535.7 536.3 536.6 537.0 536.5 536.8 536.7 536.5 536.8 537.5 536.6
##     [,89] [,90] [,91] [,92] [,93] [,94] [,95] [,96] [,97] [,98] [,99]
## 31% 499.4 499.7 499.3 498.9 499.0 499.0 498.7 499.3 499.0 499.6 499.0
## 52% 509.0 509.2 509.3 509.0 508.8 509.1 509.2 509.3 508.8 508.7 508.6
## 60% 512.7 512.6 512.3 512.5 512.3 512.2 512.4 513.0 513.3 512.9 512.7
## 70% 516.7 516.6 516.9 516.8 516.9 516.9 516.8 517.4 517.5 517.5 517.4
## 95% 537.3 536.6 536.9 537.5 537.4 537.0 537.5 537.4 537.3 537.6 537.6
##     [,100]
## 31%  499.3
## 52%  508.7
## 60%  512.4
## 70%  516.7
## 95%  537.9
```



#### 3.1.5  Generate  1000 random samples  with 100 time periods for dry peas price using the parameters from estimate


One simulation, different time different result. And time period still is 100 ("yuimu" default: delta = 1/100).

```r
require(yuima)
# initial value is 188
m1 = setModel(drift = "theta*(mu-x)", diffusion = "sigma", state.var = "x", 
    time.var = "t", solve.var = "x", xinit = ouFit.ML(dry.peas)[[2]])
Xdry.peas = simulate(m1, true.param = list(mu = ouFit.ML(dry.peas)[[2]], sigma = ouFit.ML(dry.peas)[[3]], 
    theta = ouFit.ML(dry.peas)[[1]]))
plot(Xdry.peas)
```

![plot of chunk unnamed-chunk-18](figure/unnamed-chunk-18.png) 



Increasing the number of simulation to 1000, and time period still is 100 ("yuimu" default: delta = 1/100). Plot the mean of the 1000 silumtion result, different time , very similar result.

```r
simnum = 1000
# specific qunatile (which we can pick any another quantile)
dist = c(0.31, 0.52, 0.6, 0.7, 0.95)
newsim = function(i) {
    simulate(m1, true.param = list(mu = ouFit.ML(dry.peas)[[2]], sigma = ouFit.ML(dry.peas)[[3]], 
        theta = ouFit.ML(dry.peas)[[1]]))@data@original.data
}
# newsim(1) simulation 1000 times, each time there are 100 time periods
sim = sapply(1:simnum, function(x) newsim(x))
# transfor to time seires format.
m2 = t(sim)
mdry.peas <- apply(m2, 2, mean)

# plot the mean of the 1000 time simulation for the 100 time periods
plot(mdry.peas, type = "l")
```

![plot of chunk unnamed-chunk-19](figure/unnamed-chunk-19.png) 

```r

# find out the quantile to decribe the distribution
tile = sapply(1:100, function(x) quantile(m2[, x], dist))
tile
```

```
##      [,1]  [,2]  [,3]  [,4]  [,5]  [,6]  [,7]  [,8]  [,9] [,10] [,11]
## 31% 279.5 278.8 278.5 278.3 278.1 278.0 277.8 277.6 277.3 277.3 277.3
## 52% 279.5 279.5 279.6 279.6 279.6 279.7 279.6 279.7 279.6 279.7 279.7
## 60% 279.5 279.8 280.0 280.0 280.1 280.4 280.3 280.4 280.5 280.5 280.4
## 70% 279.5 280.1 280.5 280.7 280.9 281.2 281.2 281.3 281.5 281.7 281.8
## 95% 279.5 281.7 282.8 283.2 283.9 284.3 284.9 285.4 286.2 286.4 286.7
##     [,12] [,13] [,14] [,15] [,16] [,17] [,18] [,19] [,20] [,21] [,22]
## 31% 277.3 277.2 277.1 277.0 276.9 276.8 276.5 276.6 276.7 276.4 276.2
## 52% 279.7 279.7 279.8 279.7 279.7 279.8 279.8 279.9 279.8 279.8 280.0
## 60% 280.6 280.7 280.9 280.7 281.0 281.0 280.9 281.1 281.0 281.1 281.1
## 70% 281.9 282.2 282.1 282.3 282.2 282.3 282.6 282.7 282.7 282.9 283.0
## 95% 287.0 287.5 287.8 287.6 288.4 288.5 288.9 289.4 289.3 289.5 289.5
##     [,23] [,24] [,25] [,26] [,27] [,28] [,29] [,30] [,31] [,32] [,33]
## 31% 276.1 275.9 276.0 276.0 275.8 275.7 275.5 275.5 275.4 275.6 275.2
## 52% 279.7 279.9 279.6 279.7 279.9 279.9 279.7 279.8 279.7 279.9 279.8
## 60% 281.1 281.2 281.1 281.1 281.1 281.4 281.2 281.4 281.4 281.3 281.4
## 70% 282.9 282.8 283.0 282.9 283.2 283.5 283.5 283.4 283.4 283.5 283.6
## 95% 290.2 290.2 290.3 290.3 291.3 291.1 291.5 291.1 291.6 292.1 292.3
##     [,34] [,35] [,36] [,37] [,38] [,39] [,40] [,41] [,42] [,43] [,44]
## 31% 275.3 275.4 275.2 275.0 275.1 275.0 274.8 274.9 274.7 274.6 274.6
## 52% 280.0 279.9 279.9 280.1 279.9 279.7 279.5 279.6 279.5 279.7 279.8
## 60% 281.6 281.5 281.5 281.5 281.4 281.5 281.2 281.2 281.4 281.3 281.4
## 70% 283.6 283.6 283.8 283.7 283.6 283.6 283.7 283.9 284.1 284.3 284.1
## 95% 292.7 292.3 292.4 292.7 292.8 293.1 293.4 293.8 294.0 294.3 294.1
##     [,45] [,46] [,47] [,48] [,49] [,50] [,51] [,52] [,53] [,54] [,55]
## 31% 274.5 274.5 274.3 274.4 274.0 273.8 273.9 274.0 274.0 274.0 274.1
## 52% 279.6 279.6 279.6 279.9 279.7 279.8 279.8 279.7 280.0 279.8 279.9
## 60% 281.5 281.2 281.6 281.5 281.4 281.7 281.7 281.9 281.9 281.8 281.8
## 70% 284.1 284.2 284.3 284.4 284.6 284.5 284.6 284.7 284.5 284.4 284.8
## 95% 295.3 295.5 295.2 295.5 296.1 296.2 296.3 296.7 296.7 296.7 296.3
##     [,56] [,57] [,58] [,59] [,60] [,61] [,62] [,63] [,64] [,65] [,66]
## 31% 274.1 274.1 274.0 273.9 273.9 273.9 274.0 273.8 273.9 273.9 274.3
## 52% 279.9 280.0 280.0 279.9 279.9 280.0 279.8 279.9 279.9 280.0 279.9
## 60% 282.1 282.1 282.2 282.2 282.2 282.1 282.2 282.3 282.3 282.4 282.4
## 70% 284.7 285.2 284.7 284.9 285.1 285.1 284.9 285.1 285.2 285.0 285.1
## 95% 296.2 296.7 296.8 296.7 296.8 296.8 297.2 297.1 297.0 297.6 298.1
##     [,67] [,68] [,69] [,70] [,71] [,72] [,73] [,74] [,75] [,76] [,77]
## 31% 273.8 273.7 273.8 273.9 273.7 273.3 273.3 273.6 273.5 273.6 273.5
## 52% 280.0 280.1 280.0 280.2 280.2 280.0 280.3 280.2 280.0 280.1 280.2
## 60% 282.4 282.3 282.4 282.3 282.5 282.5 282.3 282.2 282.3 282.5 282.1
## 70% 285.3 285.2 285.1 285.1 285.1 285.3 285.7 285.4 285.8 285.7 285.5
## 95% 297.8 298.2 298.4 298.1 298.6 299.0 299.1 299.0 299.4 299.9 299.5
##     [,78] [,79] [,80] [,81] [,82] [,83] [,84] [,85] [,86] [,87] [,88]
## 31% 273.5 273.2 273.2 273.3 273.2 273.2 272.9 273.3 273.4 273.1 273.1
## 52% 280.1 280.0 280.2 280.2 279.9 280.3 280.4 280.1 280.2 280.1 280.2
## 60% 282.4 282.6 282.2 282.2 282.2 282.2 282.3 282.8 282.7 282.9 282.6
## 70% 285.3 285.3 285.5 285.0 285.2 285.6 285.5 285.6 285.6 285.8 285.9
## 95% 299.7 300.1 299.6 300.0 300.0 300.2 300.5 300.8 300.7 300.3 300.2
##     [,89] [,90] [,91] [,92] [,93] [,94] [,95] [,96] [,97] [,98] [,99]
## 31% 272.9 273.0 272.9 273.0 273.0 273.0 273.2 273.0 272.5 272.4 272.5
## 52% 280.2 280.2 280.1 280.1 280.1 280.1 280.1 280.3 279.8 280.3 280.3
## 60% 282.8 282.7 282.6 282.6 282.6 282.5 282.4 282.5 282.6 282.8 282.8
## 70% 286.0 286.1 285.9 285.9 285.9 285.8 286.0 285.7 285.7 286.1 285.9
## 95% 300.1 300.0 300.3 300.1 300.8 300.3 300.7 300.8 300.7 301.2 301.2
##     [,100]
## 31%  272.7
## 52%  280.4
## 60%  282.6
## 70%  286.3
## 95%  301.0
```


### 3.2 "Sde.sim" packages

Recall "sde.sim" package uses this notation:

$$dX_t = (\theta_1 - \theta_2X_t)dt +\theta_3dW_t,    X_0 = x_0$$

For "sde.sim" packages, each price simulation uses different random seed. When we try the same random seed, all prices have the same trend.




#### 3.2.1 Generate  1 random samples with 1000 time periods for wheat price using the parameters from estimate


```r
set.seed(123)

thetaw <- c(ouFit.ML(wheat)[[1]] * ouFit.ML(wheat)[[2]], ouFit.ML(wheat)[[1]], 
    ouFit.ML(wheat)[[3]])
Xwheat <- sde.sim(X0 = ouFit.ML(wheat)[[2]], model = "OU", theta = thetaw, N = 1000, 
    delta = 1)
```

```
## 
## T set to = 1000.000000
```

```r
plot(Xwheat, main = " Ornstein - Uhlenbeck ")
```

![plot of chunk unnamed-chunk-20](figure/unnamed-chunk-20.png) 


#### 3.2.2  Generate  1 random samples  with 1000 time periods for oats price using the parameters from estimate



```r
set.seed(321)

thetao <- c(ouFit.ML(oats)[[1]] * ouFit.ML(oats)[[2]], ouFit.ML(oats)[[1]], 
    ouFit.ML(oats)[[3]])
Xoats <- sde.sim(X0 = ouFit.ML(oats)[[2]], model = "OU", theta = thetao, N = 1000, 
    delta = 1)
```

```
## 
## T set to = 1000.000000
```

```r
plot(Xoats, main = " Ornstein - Uhlenbeck ")
```

![plot of chunk unnamed-chunk-21](figure/unnamed-chunk-21.png) 


#### 3.2.3 Generate  1 random samples  with 1000 time periods for barley price using the parameters from estimate


```r
set.seed(132)

thetab <- c(ouFit.ML(barley)[[1]] * ouFit.ML(barley)[[2]], ouFit.ML(barley)[[1]], 
    ouFit.ML(barley)[[3]])
# set initial value as the long term mean mu.
Xbarley <- sde.sim(X0 = ouFit.ML(barley)[[2]], model = "OU", theta = thetab, 
    N = 1000, delta = 1)
```

```
## 
## T set to = 1000.000000
```

```r
plot(Xbarley, main = " Ornstein - Uhlenbeck ")
```

![plot of chunk unnamed-chunk-22](figure/unnamed-chunk-22.png) 

#### 3.2.4 Generate  1 random samples with 1000 time periods for canola price using the parameters from estimate


```r
set.seed(312)

thetac <- c(ouFit.ML(canola)[[1]] * ouFit.ML(canola)[[2]], ouFit.ML(canola)[[1]], 
    ouFit.ML(canola)[[3]])
Xcanola <- sde.sim(X0 = ouFit.ML(canola)[[2]], model = "OU", theta = thetac, 
    N = 1000, delta = 1)
```

```
## 
## T set to = 1000.000000
```

```r
plot(Xcanola, main = " Ornstein - Uhlenbeck ")
```

![plot of chunk unnamed-chunk-23](figure/unnamed-chunk-23.png) 

#### 3.2.5 Generate  1 random samples with 1000 time periods for dry.peas price using the parameters from estimate


```r
set.seed(231)

thetad <- c(ouFit.ML(dry.peas)[[1]] * ouFit.ML(dry.peas)[[2]], ouFit.ML(dry.peas)[[1]], 
    ouFit.ML(dry.peas)[[3]])
Xdry.peas <- sde.sim(X0 = ouFit.ML(dry.peas)[[2]], model = "OU", theta = thetad, 
    N = 1000, delta = 1)
```

```
## 
## T set to = 1000.000000
```

```r
plot(Xdry.peas, main = " Ornstein - Uhlenbeck ")
```

![plot of chunk unnamed-chunk-24](figure/unnamed-chunk-24.png) 



## 4  Conclusion

This documents show we can use ouFit.ML function to estimate the parameter $\theta$, $\mu$, and $\sigma$.
And them use the **"yuima"** package to generate the simulative price vector/matrices. 

