#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Equilibrium-based size frequency estimator function (ebSizeFreqEstimator_function_CVt) ##
#
# Version100 by Nils Krueck, 19 April 2023
# Edited by Sarah Willington
#   Changes: function used for runs with fixed mortality and changing CVt
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

ebSizeFreqEstimator_function_sarah <- function(Linf,CV.Linf,t0,tmax,Lbins,L.classes,
                                               k,M, trunc.param, run.name) {
    
  #@@@@@@@@@@@@@@@@@
  # Input required
  #@@@@@@@@@@@@@@@@@

 # #key life history parameters
 #  Linf <- 49.5    # von Bertalanffy Linfinity (average length of fully grown individuals)
 #  t0 <- 0         # von Bertalanffy t0 (theoretical length at age 0)
 #  tmax <- 20      # maximum age
 #  CV.Linf <- seq(0.1,0.1 , length=tmax) #used for CV of length at age
 #  L.bins <- 1     # specify length bins in cm
 #  L.classes <- seq(0,Linf*1.5,L.bins)+0.5*L.bins # generates lengths per bin, i.e. length classes# growth rate
 #  k <- 0.2        # von Bertalanffy growth coefficient
 #  M <- 0.3       # instantaneous rate of matural mortality (per year)
 #  trunc.param <-  0.66 ## Sarah - include this parameter somewhere above, outside of the function. Make sure you run scenarios with 0, 0.25 and 0.66 truncation. 0.66*Linf truncation reflects life-history invariant of maturation size, so we then only look at adults.
 #  run.name <- ""

  # load relevant libraries
  library(tidyverse)
  library(nls2)
  library(scales)
  
  
  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  # Generate key modelling input and population #####
  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  
  R0 <- 100       # arbitrary number of recruits in unfished population
  ns <- length(L.classes) # number of length classes
  
  
  # calculate age at length
  age.at.length <- suppressWarnings(-log(-(L.classes-Linf)/Linf)/k+t0) # von Bertalanffy growth function (VBGF) to estimate age at length
  age.at.length[is.na(age.at.length)|!is.finite(age.at.length)] <- tmax # assign tmax to undefined values
  vbgf.formula <- length ~  Linf *(1 - exp(-k * (age - t0))) # generate VBGF formula for modelling
  vbgf.vals <- data.frame(Linf = Linf,k = k,t0 = t0) # specify key VBGF parameters from above
  vbgf.dat <- data.frame(age=age.at.length,length=L.classes) # generate input data as specified above (idealized)
  vbgf.model <- nls2(vbgf.formula,data = vbgf.dat,start = vbgf.vals, algorithm = "brute-force") # capture VGBF in model format
  ages <- seq(1,tmax,1) # specify ages for age-length predictions
  na <- length(ages)
  dages <- seq(0,tmax,0.001) # specify decimal ages for more precise age-length predictions
  length.at.age <- predict(vbgf.model,newdata=data.frame(age=ages-0.5),interval="prediction",alpha=c(0.95)) # length at age class
  length.at.dage <- predict(vbgf.model,newdata=data.frame(age=dages))  # length at decimal ages
  
  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  # Length frequencies by age matrix @
  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  
  ns <- length(L.classes)
  lf.age.mat <- matrix(rep(0,max(ages)*ns),nrow=ns,ncol=max(ages)) # generate age by length matrix (normalised version)
  lf.age.mat.raw <- matrix(rep(0,max(ages)*ns),nrow=ns,ncol=max(ages)) # generate age by length matrix
  for(ai in ages){ # loop over ages
    nprobs_length <- dnorm(L.classes,mean=length.at.age[ai],sd = length.at.age[ai]*CV.Linf[ai]) # generate length frequencies at age I EDITED THIS TO GET CV CHANGE WITH AGE
    lf.age.mat.raw[,ai] <- nprobs_length # save in matrix
    lf.age.mat[,ai] <- nprobs_length/sum(nprobs_length) # normalize to 1, then save in matrix 
  } # end of age loop
  
  
  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  # Prepare model and plotting #
  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  
  alpha.val <- 0.5 # color transparency for plotting
  ssa <- exp(-M)     # natural survival (converted from instantaneous M), assumed constant for all ages
  Mk <- M/k        # key life history invariant
  rMk <- round(Mk*10)/10 # rounded 
  Mk.ratio <- round(Mk, 3)
  print(paste0("Mk ratio = ",rMk))
  
  nprsa <- ssa^ages # relative numbers of survivors per age class
  
  if(identical(ssa,rev(ssa))=="FALSE"){ # if survival varies by age class, then we'll have to loop over to determine survivers
    nprsa <- 1; # numbers per recruit
    for(si in 2:length(ages)){
     # if(si < length(ssa)){
        nprsa[si] <- nprsa[si-1] * ssa[si-1]
     # }else if(si == length(ssa)){
    #    nprsa[si] <- nprsa[si-1] * ssa[si-1]/(1-ssa[si])}
    }
  }
  
  N0sa <- nprsa*R0 # virgin abundance 
  N0ss.mat <- lf.age.mat * t(matrix(rep(N0sa,ns),na,ns))
  N0ss <- rowSums(N0ss.mat)
  N0ss.scaled <- N0ss/max(N0ss)
  
  
  ## setting truncated distribution
  
  truncated.threshold.length <- trunc.param * Linf
  L.classes.truncated <- L.classes[L.classes>truncated.threshold.length]
  sdist.threshold <- length(L.classes)-length(L.classes.truncated)+1
  N0ss.truncated <- N0ss[sdist.threshold:length(N0ss)]
  N0ss.truncated.scaled <- N0ss.truncated/max(N0ss.truncated)
  #plot(N0ss.truncated.scaled, type = 'l')
  
  
  ## create lists to save outputs
  sim.pop <- list() # create list to save as outputs from calculations
  sim.pops <- list()
  save.names <- matrix()
  plot_list <- list() #save the size dist plots to this list
  df_list <- list() # save the df of each plot
    
    #@@@@@@@@@@@
    # Plotting
    #@@@@@@@@@@@
  
  #set threshold for filtering and scaling length class
  max.threshold <- 0.99 * sum(N0ss,na.rm = TRUE)
  Max.Length.p99 <- min(L.classes[which(cumsum(N0ss)>max.threshold)]) 
  
  # create a ggplot version to save - added by sarah 
  df <- data.frame(length=L.classes.truncated, number=N0ss.truncated.scaled, Mk.ratio=rMk, Mk.ratio.select=Mk.ratio, L.infinity=Linf, M=M, K=k, sens.run.name=sens.run.name)
  df <- df %>% 
        unite(group.col, Mk.ratio, L.infinity, sep = "_", remove = FALSE) %>% 
    filter(length<=Max.Length.p99) %>% 
    mutate(length.scaled = length/Max.Length.p99)
      
  bn <- ggplot(df, aes(x = (length), y=number)) +
        geom_line(linewidth=1) +
        labs(x = "Length class (cm)", y = "Relative number") +
        ggtitle(paste0("M/k ratio = ",rMk,"\nLinf = ",Linf, "\nM = ",M, "\nK = ",k))+
        scale_x_continuous(breaks = seq(0, length(L.classes.truncated), by = 5))+
        theme_classic() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
              plot.title = element_text(size = 11))
      
   # Save the plot and df to the list
  plot_list[[1]] <- bn  # Record the current plot
  df_list[[1]] <- df #record the df used to make the plot
      

    
  ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  
    #@@@@@@@@@@@@@
    # save output
    #@@@@@@@@@@@@@
    
    # set thresholds for max and median
    max.threshold.unfished <- 0.99 * sum(N0ss,na.rm = TRUE) # set to 1% of all individuals
    median.threshold.unfished <- 0.5 * sum(N0ss,na.rm = TRUE)
    
    sim.pop$Unfished.Mean.Length <- sum(L.classes.truncated*N0ss.truncated)/sum(N0ss.truncated)
    sim.pop$Unfished.SD.Length <- sqrt(abs(sum((L.classes.truncated-sim.pop$Unfished.Mean.Length)**2*N0ss.truncated)/(sum(N0ss.truncated)-1))) # added absolute value function here as some truncations were resulting in a negative value here which gave NaN 
    sim.pop$Unfished.Max.Length.p99 <- min(L.classes[which(cumsum(N0ss)>max.threshold.unfished)])
    sim.pop$Unfished.Median.Length.p50 <- min(L.classes[which(cumsum(N0ss)>median.threshold.unfished)])
    sim.pop$Length.classes <- L.classes.truncated
    sim.pop$Unfished.Number <- N0ss.truncated
    sim.pop$Unfished.Lengths.At.Age <- N0ss.mat
    sim.pop$Lengths.At.Age.Rel.Frequency <- lf.age.mat
    sim.pop$Mk.ratio <- Mk.ratio
    
    sim.pops[[1]] <- sim.pop
  
  # assign scenario names to list
  names(sim.pops) <- save.names
  
  # output from model
  return(list(sim.pops,plot_list, df_list))

} # end of function
