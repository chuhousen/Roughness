R script for simultaneous estimation of roughness length (Z0) & zero-plane displacement height (d) 
  from single-level eddy covariance measurements

 The code modifies and implements the 3 approaches used in Graf et al.,2014 Boundary-Layer Meteorol.
   FP-RE-1: from logarithmic wind profile theory: WS/Ustr = 1/k*(ln((zm-d)/Z0)-beta*((zm-d)/Lm))

   FP-RE-2: from logarithmic wind profile theory: WS = Ustr/k*(ln((zm-d)/Z0)-beta*((zm-d)/Lm))

   FV-RE-1: from flux-variance similarity theory: std.Uz/Ustr = C1*(1-C2*(zm-d)/Lm)^(1/3)
                                                  std.Uz = k*C1*WS/ln((zm-d)/Z0)                 

[Desccription]
   The goal is to estimate Z0 and d from continuous measurements of wind and turbulent statistics, i.e., WS, Ustr, std.Uz..... 
   All input variables are originally measured at a 30-min time step, but are pre-filtered according to data quality & required meteorological conditions.
   i.e., not every 30 min is used in model fitting.

   While estimating Z0 & d, the time series is re-grouped to a daily step, i.e., treating each 30 min within a day as a sampling replicate.
   Z0 and d are estimated at a daily time step, as they are both functions of canopy structures (height, leaf area...), assuming no changes within a day   
   The time series of Z0 or d is assumed as an autoregressive process, and it's assumed there's correlation (cor) between the time series of Z0 & d
   Z0[t+1] = Z0[t] + Z0.err ; Z0.err ~ N(0,sigma.Z0)
   d[t+1] = d[t] + sigma.d/sigma.Z0*cor*(Z0[t+1]-Z0[t]) + d.err ; d.err ~ N(0,sigma.d.iid)

[Input variables]  
       WS: wind speed (m/s)
       Ustr: friction velocity (m/s)
       zm: measurement height (m)
       Lm: Monin–Obukhov length (m)
       std.Uz: standard deviation of vertical wind velocity (m/s)
       hc: vegetation canopy height (m) *only used in filtering data*
       ZL: Monin–Obukhov stability (unitless) *only used in filtering data*

[Constants]      
       k: 0.4 von karman constant
       C1: 1.3
       C2: 2.0 estimates for the universal constants, (Panofsky & Dutton 1984)
       beta: 6.0

[Output variables]
       Z01[DOY]: roughness length estimated for each day (DOY) by using FP-RE-1 model
       d1[DOY]: displacement height estimated for each day (DOY) by using FP-RE-1 model
       Z02[DOY]: roughness length estimated for each day (DOY) by using FP-RE-2 model
       d2[DOY]: displacement height estimated for each day (DOY) by using FP-RE-2 model
       Z03[DOY]: roughness length estimated for each day (DOY) by using FV-RE-1 model
       d3[DOY]: displacement height estimated for each day (DOY) by using FV-RE-1 model
