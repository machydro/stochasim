
#' Estimate flow velocity
#'
#' \code{velocity} returns the mean flow velocity (m/s) for a given hydraulic
#' radius, bed surface D84, and reach average gradient using Ferguson (2007)
#'
#' @param R hydraulic radius for channel cross section (m)
#' @param S reach average channel gradient (m/m)
#' @param D84 84th percentile of bed surface in main channel (mm)
#' @export
velocity = function(R, S, D84){
  g = 9.81
  D84 = D84/1000
  a1 = 5.5
  a2 = 2.2
  Res = a1*a2*(R/D84) / (a1^2 + a2^2*(R/D84)^(5/3))^(1/2)
  U = Res * sqrt(g * R * S)
  return(U)
}


#' Calculate the mean boundary shear stress
#'
#' \code{stress} calculates the mean boundary shear stress (Pa), using the
#' hydraulic radius, reach average gradient and fluid density
#'
#' @param R hydraulic radius for channel cross section (m)
#' @param S reach average channel gradient (m/m)
#' @param rho (optional) fluid density (kg/m3), default density is water
#' @export
stress = function(R, S, rho = 1000){
  g = 9.81
  stress = g * rho * R * S
  return(stress)
}


#' Estimate the critical shear stress
#'
#' \code{crit_stress} estimates the shear stress required to entrain a grain of
#' size Di from a bed surface with a median size of D50. The equations are based
#' on the equations presented in Wilcock and Crowe (2001), assuming a critical
#' dimensionless shear stress for the D50 of 0.045 (as a default)
#'
#' @param Di vector of grain sizes to analyse (mm)
#' @param D50 median grain size for the bed surface (mm)
#' @param TC50 (optional) critical dimensionlessstress for D50 (default is 0.036)
#' @export
crit_stress = function(Di, D50, TC50 = 0.036){
  g = 9.81
  rho = 1000
  rho.sed = 2650
  Tci = numeric(length(Di))
  if(length(D50) == 1){
    Tcrit = TC50 * g * (rho.sed - rho) * D50 / 1000
    Tci[which(Di < D50)] = Tcrit * (Di[which(Di < D50)]/D50)^0.12
    Tci[which(Di >= D50)] = Tcrit * (Di[which(Di >= D50)]/D50)^0.67
    return(Tci)
  }else{
    return(NA)
  }
}


#' Calculate stable grain fractions
#'
#' \code{stab_fraction} calculates the proportion of a grain size fraction that
#' is likely to remain immobile during a flood event, based on the ratio of the
#' mean boundary shear stress to the critical shear stress to entrain the that
#' grain size fraction
#'
#' @param T mean boundary shear stress (Pa)
#' @param Tci critical shear stress(es) for size fraction(s) of interest (Pa)
#' @param Map (optional) TRUE = display function, default is FALSE
#' @export
stab_fraction = function(T, Tci, Map = FALSE){
  t =  T / Tci
  stab.frac =  1 - (1)/(1 + exp(-4.3*(t - 1.5)))

 if(Map == TRUE){
    par(mfcol = c(1,1))
    tm = seq(0, 3, 0.01)
    fr = 1 - (1)/(1 + exp(-4.3*(tm - 1.5)))
    plot(tm, fr, type = "l",
         col = "red",
         lty = 2,
         ylim = c(0,1),
         xlab = "stress / critical stress",
         ylab = "stable fraction")
    abline(v = t, col = "blue", lty = 3)
    legend("topright",
           cex = 0.9,
           bty = "n",
           col = c("red","blue"),
           lty = c(2,2),
           legend = c("stability function", "output values"))
  }
  return(stab.frac)
}


#' Simulate a log-normal grain size distribution
#'
#' \code{sim_gsd} generates a log-normal grain size distribution for a
#' user-specified D50. The user can also vary the spread of the distribution.
#' The function returns a data frame containing the upper range of the size
#' classes, the mid point of each size class, the proportion in each size class,
#' and the cumulative proportion finer than the upper size class boundaries.
#'
#' @param D50 median bed surface grain size (mm)
#' @param sp (optional) spread of the distribution (default = 1)
#' @param Map (optional) LOGICAL, when TRUE, the distribution is graphed
#' @export
sim_gsd = function(D50, sp = 1, Map = FALSE){
  dist.mean = log2(D50)
  grain.sizes = 2^rnorm(100000, mean = dist.mean, sd = sp)
  lim.upper = 2^(dist.mean + 2.5*sp)
  lim.lower = 2^(dist.mean - 2.5*sp)
  grain.sizes = grain.sizes[-which(grain.sizes >lim.upper)]
  grain.sizes = grain.sizes[-which(grain.sizes < lim.lower)]
  phi.sizes = 2^seq(dist.mean - 2.5*sp, dist.mean + 2.5*sp, 0.25*sp)
  grain.dist = hist(grain.sizes, breaks = phi.sizes, plot = F)
  p = grain.dist$counts/sum(grain.dist$counts)
  cdf = cumsum(p)
  results = data.frame(grain.dist$breaks[-1], grain.dist$mids, p, cdf)
  colnames(results) = c("size_class", "size_mid", "p", "cdf")

  if(Map == TRUE){
    par(mfcol = c(1,1))
    plot(results$size_class, results$cdf,
         type = "o",
         log = "x",
         xlim = c(lim.lower, lim.upper),
         xlab = "Grain size (mm)",
         ylab = "Proportion Finer")
  }
  return(results)
}


#' Calculate sediment transport rate
#'
#' \code{bedload} estimates the bedload transport rate per unit width based on
#' the Eaton and Church (2011) function, using the bed surface D50 and the mean
#'  boundary shear stress
#'
#' @param R hydraulic radius for channel cross section (m)
#' @param D50 median surface grain size (mm)
#' @param D84 84th percentile of bed surface in main channel (mm)
#' @param S reach average channel gradient (m/m)
#' @param TC50 (optional) critical dimensionless stress for D50 (default = 0.036)
#' @export
bedload = function(R, D50, D84, S, Tc50 = 0.036){
  g = 9.81
  Gs = 1.65
  rho = 1000
  D50 = D50 / 1000
  omega = stress(R,S) * velocity(R,S, D84)
  om_star = omega / (rho*(g*Gs*D50)^(3/2))
  dcrit = Tc50*Gs*D50 / S
  Res = velocity(dcrit, S, D84) / sqrt(g*dcrit*S)
  om_crit = Res * Tc50^(3/2)
  E_star = (0.92 - 0.25*sqrt(om_crit/om_star))^9
  qb_star = E_star * om_star
  qb = qb_star * sqrt(R*g)*(D50)^(3/2)
  return(qb)
}


#' Analyze bed stability
#'
#' \code{bed_stability} is a function that takes the grain size distribution
#' data output from \code{simulate_gsd}, and an estimate of the mean boundary
#' shear stress to estimate the fraction of the bed that is stable, accounting
#' for partial mobility for stresses between 1 and 2 times the critical value.
#' Estimates of the stable proportion in each size class is added to the data
#' frame, as are estimates of the CDF. The total areal fraction of the bed that
#' is stable can be estimated by summing the variable p.stab in the output
#' dataframe.
#'
#' @param gsd data frame containing GSD information (as per \code{simulate_gsd})
#' @param stress estimated mean boundary shear stress acting on the bed (Pa)
#' @param TC50 (optional) critical dimensionless stress for D50 (default is 0.045)
#' @param Map (optional) LOGICAL, when TRUE, the distributions are graphed
#' @export
bed_stability = function(gsd, stress, TC50 = 0.036, Map = F){
  D50 = as.numeric(approx(x = gsd$cdf, y = gsd$size_class, xout = 0.5))
  grain.tci = crit_stress(gsd$size_mid, D50[2], TC50)
  gsd$p.stab = round(stab_fraction(stress, grain.tci) * gsd$p, digits = 3)
  gsd$p.mob = (1 - stab_fraction(stress, grain.tci)) * gsd$p
  index = gsd$p.mob / gsd$p < 0.02
  gsd$p.mob[index] = 0
  gsd$cdf.stab = cumsum(gsd$p.stab) / sum(gsd$p.stab)
  gsd$cdf.mob = cumsum(gsd$p.mob) / sum(gsd$p.mob)

  if(Map == T){
    par(mfcol = c(1,2))
    plot(gsd$size_mid, gsd$p,
         col = "darkgreen",
         type = "o",
         xlab = "grain size (mm)",
         ylab = "proportion in class",
         ylim = c(0, max(gsd$p)),
         log = "x")
    lines(gsd$size_mid, gsd$p.stab, lty = 2, col = "firebrick")
    lines(gsd$size_mid, gsd$p.mob, lty = 2, col = "orange")

    plot(gsd$size_class, gsd$cdf,
         col = "darkgreen",
         type = "o",
         ylim = c(0,1),
         xlab = "grain size (mm)",
         ylab = "proportion finer",
         log = "x")
    lines(gsd$size_class, gsd$cdf.stab, lty = 2, col = "firebrick")
    lines(gsd$size_class, gsd$cdf.mob, lty = 2, col = "orange")
    legend("topleft",
           cex = 0.9,
           bty = "n",
           col = c("orange", "darkgreen", "firebrick"),
           lty = c(2,1,2),
           legend = c("Mobile", "Surface", "Stable"))
  }
  return(gsd)
}


#' Solve for water depth
#'
#' \code{solve_Q} uses a midpoint solution approach to find water depth that
#' produces the user-specified stream discharge in a rectangular channel with
#' a known channel width.
#' @param Q discharge contained by the channel (m3/s)
#' @param W water surface width for the channel (m)
#' @param S reach average channel gradient (m/m)
#' @param D84 84th percentile of the bed surface grain size distribution (mm)
#' @export
solve_Q = function(Q, W, S, D84){
  tol = 0.0001
  bounds = c(0.01*Q^0.3, Q^0.3)
  test.d = mean(bounds)
  test.Q = W * test.d * velocity(test.d, S, D84)
  test.stat = (test.Q - Q) / Q

  while(abs(test.stat) > tol){
    if(test.stat < 0){
      bounds[1] = test.d
    }else{
      bounds[2] = test.d
    }
    test.d = mean(bounds)
    test.Q = W * test.d * velocity(test.d, S, D84)
    test.stat = (test.Q - Q) / Q
  }
  return(test.d)
}


#' Simulate a log-normal flood series
#'
#' \code{sim_flood} is a function to generate a log normal flood distribution
#' based on specifying the mean annual peak flow and an estimate of the maximum
#' probable flood
#'
#' @param Qmean the mean annual maximum flood (m3/s)
#' @param Qmax the maximum probably flood for the  system (m3/s)
#' @param n the length of the flood record to simulate
#' @param Map (optional) LOGICAL, when TRUE, the distributions are graphed
#' @export
sim_flood = function(Qmean, Qmax, n = 1000, Map = F){
  mean.log = log10(Qmean)
  max.log = log10(Qmax)
  sd = (max.log - mean.log) / 3 # FMI?
  # sd = (max.log - mean.log) * 0.8  # FMI set to
  Q = 10^rnorm(n, mean.log, sd)
  if(Map == TRUE){
    par(mfrow = c(1,2))
    ranked.Q = sort(Q, decreasing = T)
    RI = (n + 1) / 1:n
    plot(RI, ranked.Q,
         cex = 0.4,
         main = "flood frequency plot",
         xlab = "return period (yrs)",
         ylab = "discharge (cumecs)",
         col = "blue",
         log = "x")

    plot(density(Q),
         col = "blue",
         main = paste("sd = ", as.character(round(sd, digits = 2))))
  }
  return(Q)
}

#' Approximate relative bank strength
#'
#' \code{est_mu} uses an empirical function to estimate the relative bank
#' strength produced by vegetation with a known effective rooting depth (H)
#' (after Eaton, 2006) and a known average channel depth (d). The function only
#' applies when H < 0.9d.
#'
#' @param H effective rooting depth for riparian vegetation (m)
#' @param d mean hydraulic depth at channel-forming flow (m)
#' @export
est_mu = function(H,d){
  a = 0.85
  b = 0.87
  mu = 1 / (1 - a*(H/d))^b
  return(mu)
}

#' Run the STOCHASIM model
#'
#' \code{run_stochasim} executes the stochastic channel simulator code The
#' model predicts channel widening when a threshold shear stress is exceeded
#' and otherwise predicts channel narrowing due to vegetation colonization.
#' The length of the simulation is determined by the length of the input vector
#' containing the random sequence of annual peak flows (which can be created
#' using \code{sim_flood} or generated by the user).
#'
#' @param slope reach average channel slope (m/m)
#' @param gsd a data frame containing information on the bed surface GSD (see \code{sim_gsd})
#' @param Q a vector of floods to be used sequentially in the simulation (m3/s)
#' @param H (optional) estimated effective rooting depth, affecting bank stability (m)
#' @export
run_stochasim = function(gsd, Q, slope, H = 0, W = NA, D = NA, Qmaf = NA, veg.rate = 0.1){
  # gsd = gsd
  # Q = Qvals
  # slope = S
  # H = Hd
  # W = Wd
  # D = Dd

  #constants
  g = 9.81
  rho = 1000
  rho.sed = 2650 # good for granite (maybe adjust for limestone basins)
  a = 3.5  #hydraulic geometry coeff
  b = 0.5  #hydraulic geometry exponent
  # a = 4.05  # new ones for Victoria Cr - see get_hydraulicvars.R
  # b = 0.31  #
  # veg.rate = 0.1  #rate of annual gravel bar colonization (fraction)
  stab.prop = 0.06  #proportion of bed stable at threshold
  flood.dur = 12*60*60 #assume floods last for 1/2 day

  # sizes = approx(x = gsd$cdf, y = gsd$size_class, xout = c(0.5, 0.84, 0.99))
  # D50 = sizes$y[1]
  # D84 = sizes$y[2]
  # Dmax = sizes$y[3]

  # Replaces linear interp w/ empirical piecewise linear midpoint - MC
  getgrains = function(x){rep(gsd$size_mid[gsd$size_mid == x], gsd$p[gsd$size_mid == x]*100)}
  GrainSize = unlist(sapply(gsd$size_mid, getgrains))
  D50 = quantile(GrainSize, 0.50, type = 5)
  D84 = quantile(GrainSize, 0.84, type = 5)
  Dmax = quantile(GrainSize, 0.99, type = 5)
  # ----------------------------

  mean.Q = ifelse(is.na(Qmaf), 0.3 * mean(Q), Qmaf)  # Mean annual flow
  # width.min = a*(0.3 * mean.Q)^b # mean.Q - approx is good enough for this model
  width.min = a*(mean.Q)^b # bank width at MAF

  # find the critical threshold for gsd by converging
  tol = 0.001
  abound = c(0.5, 2.0)
  a.test = mean(abound)
  stress.threshold = a.test*crit_stress(Dmax, D50)  #max stress that can be maintained (no riparian veg)
  gsd = bed_stability(gsd, stress.threshold)
  test = sum(gsd$p.stab)
  converg = (test - stab.prop) / stab.prop
  counter = 0
  while(abs(converg) > tol){
    counter = counter + 1
    if(converg > 0){
      abound[1] = a.test #do this if depth.test > depth.target
    }else{
      abound[2] = a.test #do this if depth.target > depth.test
    }
    a.test = mean(abound)
    stress.threshold = a.test*crit_stress(Dmax, D50)  #max stress that can be maintained (no riparian veg)
    gsd = bed_stability(gsd, stress.threshold)
    test = sum(gsd$p.stab)
    converg = (test - stab.prop) / stab.prop
  }

  depth.threshold = stress.threshold / (g*rho*slope)  #max depth that can be maintained (no veg)

  if(H > 0){   #if mu > 0, adjust threshold and depth
    H2 = ifelse(H*0.9 > depth.threshold, depth.threshold, H) # modify initial conditions -- MC
    bounds = c(1, est_mu(H2,depth.threshold))*depth.threshold
    depth.test = mean(bounds)
    depth.target = stress.threshold * est_mu(H2,depth.test) / (g*rho*slope)
    converg = (depth.test - depth.target) / depth.target
    counter = 0
    while(abs(converg) > tol){
      counter = counter + 1
      if(converg > 0){
        bounds[2] = depth.test #do this if depth.test > depth.target
      }else{
        bounds[1] = depth.test #do this if depth.target > depth.test
      }
      depth.test = mean(bounds)
      depth.target = stress.threshold * est_mu(H2,depth.test) / (g*rho*slope)
      converg = (depth.test - depth.target) / depth.target
    }
    stress.threshold = a.test*crit_stress(Dmax, D50) * est_mu(H2,depth.test)
    depth.threshold = depth.test
  }

  vel.threshold = velocity(depth.threshold, slope, D84)  #velocity associated with threshold depth
  gsd = bed_stability(gsd, stress.threshold)
  print(paste("proportion of the bed that is stable at threshold stress is: ",as.character(sum(gsd$p.stab))))

  nsim = length(Q)
  t = seq(1,nsim,1) #time domain
  width = numeric(length = nsim)
  erosion = numeric(length = nsim)
  depth = numeric(length = nsim)
  transport = numeric(length = nsim)

  # Added override => if width provided, use it, otherwise calculate from relations MC
  if(is.na(W)){width[1] = a * mean.Q ^ b  } else{ width[1] = W  }  # doesn't really matter much...
  # width = W
  erosion[1] = NA
  depth[1] = NA
  transport[1] = NA

  for(i in seq(2, nsim,1)){
    #calculate the flood shear stress
    depth[i] = solve_Q(Q[i], width[i-1],slope, D84)
    tau = stress(depth[i], slope)
    vel = velocity(depth[i], slope, D84)

    if(tau > stress.threshold){
      #reduce the shear stress by widening the channel
      width[i] = Q[i] / (depth.threshold*vel.threshold)
      erosion[i] = width[i] - width[i-1]
      transport[i] = bedload(depth[i], D50, D84, slope) * width[i-1] * flood.dur/2 +
        bedload(depth.threshold, D50, D84, slope) * width[i] * flood.dur/2
    }else{
      erosion[i] = NA
      width[i] = width.min + (1 - veg.rate)*(width[i-1] - width.min)
      transport[i] = bedload(depth[i], D50, D84, slope) * width[i-1] * flood.dur
    }
  }
  results = data.frame(t,Q, width, depth, transport, erosion)
  results = results[order(results$Q, decreasing = T), ]
  stable = character(length = nsim)
  stable[is.na(results$erosion)] = "stable"
  stable[!is.na(results$erosion)] = "unstable"
  results$stable = as.factor(stable)
  rank = seq(1, nsim, 1)
  results$RI = (nsim + 1)/rank
  results = results[order(results$t, decreasing = F), ]

  return(results[-1,])  #remove the inital row with NAs
}


#' Plot erosion data from stochasim
#'
#' \code{plot_erosion} takes the results produced by running the STOCHASIM
#' model and produces a plot of the bank erosion magnitude and frequency
#'
#' @param results a data frame containing the output from \code{run_stochasim}
#' @param Est (optional) LOGICAL, TRUE returns RI of bank erosion
#' @export
plot_erosion = function(results, Est = F){
  par(mfrow = c(1,2), mar = c(4.5, 4, 1, 1), oma = c(1, 1, 1, 1))
  plot(results$Q, results$erosion,
       log = "x",
       col = "blue",
       xlab = expression(Q~(m^3/s)),
       ylab = "bank erosion (m)")
  abline(v = mean(results$Q), lty = 2, col = "red")

  erosion.ranked = sort(results$erosion, decreasing = T, na.last = T)
  nsim = length(results$Q)
  RI = (nsim + 1)/ seq(1:nsim)
  RI.erosion = RI[which(erosion.ranked == min(erosion.ranked, na.rm = T))]
  plot(RI, erosion.ranked,
       col = "blue",
       ylab = "bank erosion (m)",
       xlab = "return period (yrs)",
       log = "x")
  abline(v = 2.33, lty = 2, col = "red")
  if(Est == T) return(RI.erosion)
}


#' Estimate the effective discharge based on bedload transport
#'
#' \code{est_effective} takes results produced by running STOCHASIM and
#' estimates the effective discharge (i.e. the one that does the most geom
#' -orphic work. The function can also produce a plot of the proportion of
#' the total work associated with a given flow.
#'
#' @param results a data frame containing the output from \code{run_stochasim}
#' @param Map (optional) LOGICAL, when TRUE, the transp. frequency dist. is graphed
#' @export
est_effective = function(results, Map = F){
  # results = x
  results = results[order(results$Q),]
  sum.wk = cumsum(results$transport)/sum(results$transport)
  Xsum = seq(min(results$Q),
              max(results$Q),
              length.out = 40)
  sum = approx(y = sum.wk, x = results$Q, xout = Xsum)
  Xp = Xsum[-1] - 0.5*(Xsum[2] - Xsum[1])
  Yp = diff(sum$y)
  loess.fit = loess(Yp ~ Xp, span = 0.35, degree = 2)
  Ypred = predict( loess.fit, new = data.frame(Xp = results$Q) )
  RIp = approx(y = results$RI, x = results$Q, xout = Xp)
  #flow50 = approx(x = sum.wk, y = results$Q, xout = 0.5)
  EffPeak = results$Q[which(Ypred == max(Ypred, na.rm = T))]
  RIpeak = approx(y = results$RI, x = results$Q, xout = EffPeak)
  effective = data.frame(EffPeak, RIpeak$y)
  colnames(effective)= c("Q", "RI")
  if (Map == T){
    par(mfrow = c(1,1), mar = c(4, 4.5, 1, 1), oma = c(1, 1, 1, 1))
    # plot(results$Q, results$transport,
    #      col = "blue",
    #      log = "xy",
    #      xlab = expression(Q~(m^3/s)),
    #      ylab = expression(Q[b]~(m^3/event)))
    plot(RIp$y, Yp,
         log = "x",
         col = "blue",
         xlab = "return period (yrs)",
         ylab = "prop. of cumulative transp.")
    lines(results$RI, Ypred,
          col = "blue",
          lty = 2)
    # abline(v = 2.33, col = "red", lty = 2)
    abline(v = effective$RI, col = "red", lty = 2)
  }
 return(effective)
}


#' Plot time series of discharge and width
#'
#' \code{plot_timeseries} takes results produced by running STOCHASIM and
#' produces plots of the discharge series and of the channel width series
#' simulated by the model. The interquartile range for each distribution is
#' shown on the plots using dashed red lines.
#'
#' @param results a data frame containing the output from \code{run_stochasim}
#' @export
plot_timeseries = function(results){
  n = min(c(length(results$Q), 1000))
  k = floor(n/2)
  par(mfrow = c(2,1), mar = c(1, 4.5, 1, 1), oma = c(1, 0, 0, 1))
  plot(seq(k,n), results$Q[k:n],
       type = "l",
       log = "y",
       col = "blue",
       lwd = 0.5,
       ylab = expression(Q~(m^3/s)),
       xlab = "",
       ylim = c(min(results$Q), max(results$Q)),
       xaxt = "n")
  abline(h = quantile(results$Q, probs = c(0.25, 0.75)),
         col = "red",
         lty = 2,
         lwd = 0.5)
  plot(seq(k,n), results$width[k:n],
       col = "blue",
       ylab = expression(width~(m)),
       xlab = "",
       lwd = 0.5,
       ylim = c(min(results$width), max(results$width)),
       type = "l")
  abline(h = quantile(results$width, probs = c(0.25, 0.75)),
         col = "red",
         lty = 2,
         lwd = 0.5)
}

#' Plot flood frequency distribution
#'
#' \code{plot_floods} takes results produced by running STOCHASIM and
#' produces a typical flood frequency plot. The data are color coded by
#' whether or not the channel was unstable during the flood event. the limits
#' of the stable and unstable flood populations are indicated with a dashed
#' line.
#'
#' @param results a data frame containing the output from \code{run_stochasim}
#' @export
plot_floods = function(results){
  par(mfrow = c(1,1), mar = c(4, 4.5, 1, 1), oma = c(1, 1, 1, 1))
  index = results$stable == "stable"
  plot(results$RI[index], 0.98*results$Q[index],
       log = "x",
       pch = 20,
       cex = 0.5,
       ylim = c(0, max(results$Q)),
       xlim = c(1, max(results$RI)),
       xlab = "return period (yrs)",
       ylab = expression(Q~(m^3/s)),
       col = "blue")
  abline(v = max(results$RI[index]), col = "blue", lty = 2)
  index = results$stable == "unstable"
  points(results$RI[index], 1.02*results$Q[index],
         pch = 20,
         cex = 0.5,
         col = "firebrick")
  abline(v = min(results$RI[index]), col = "firebrick", lty = 2)
  legend("bottomright",
         cex = 0.9,
         bty = "n",
         col = c("firebrick","blue"),
         pch = c(20,20),
         legend = c("unstable", "stable"))
}

#' Estimate the formative discharge for a channel based on bank erosion
#'
#'  \code{est_formative} takes results produced by running STOCHASIM and
#' estimates the formative discharge, defined to be the discharge that does the
#' most bank erosion over time. Smaller flows can erode the banks, but these
#' ones are responsible for the majority of the bank erosion
#'
#' @param results a data frame containing the output from \code{run_stochasim}
#' @param Map (optional) LOGICAL, when TRUE, the bank erosion frequency dist. is graphed
#' @export
est_formative = function(results, Map = F){
  #results = results1
  results = results[order(results$Q),]
  index = results$stable == "unstable"
  X = results$Q[index]
  Y = results$erosion[index]
  RI = results$RI[index]
  sum = approx(y = cumsum(Y), x = X,
               xout = seq(min(X), max(X),length.out = 40))
  Xp = sum$x[-1]
  Yp = diff(sum$y)
  loess.fit = loess(Yp ~ Xp, span = 0.35, degree = 2)
  Ypred = predict( loess.fit, new = data.frame(Xp = X) )
  RIp = approx(y = results$RI, x = results$Q, xout = Xp)
  #flow50 = approx(x = sum.wk, y = results$Q, xout = 0.5)
  form.Q = X[which(Ypred == max(Ypred, na.rm = T))]
  form.RI = approx(y = RI, x = X, xout = form.Q)
  formative = data.frame(form.Q, form.RI$y)
  colnames(formative) = c("Q", "RI")

  if (Map == T){
    par(mfrow = c(1,1), mar = c(4, 4.5, 1, 1), oma = c(1, 1, 1, 1))
    plot(RIp$y, Yp,
         log = "x",
         col = "blue",
         xlim = c(1, max(RIp$y)),
         xlab = "return period (yrs)",
         ylab = "prop. of cumulative bank erosion")
    lines(RI, Ypred,
          col = "blue",
          lty = 2)

    abline(v = formative$RI, col = "red", lty = 2)
  }
  return(formative)
}
