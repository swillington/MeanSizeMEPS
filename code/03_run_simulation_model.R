#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Script to run the equilibrium-based size frequency estimator (ebSizeFreqEstimator_function_MAIN) ##
# Version 100 by Nils Krueck, 19 April 2023 
# edited 8th april 2024 by Sarah Willington to update data loading and writing code
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

rm(list=ls()) # clear workspace

# run simulation loop ####

#libraries
library(tidyverse)
library(here)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#set file paths to load and write data
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

##NOTE for each run - update file paths, trunc.param (line 42), run.name (line 44) and check correct function (line 104)

project.dir <- here::here()
mk.inpath <- file.path(paste(project.dir, 
                                 "03_results",  
                                 "output1", 
                                 "01_baseline_scenarios_alfonsino.csv", sep = "/"))
mk.outpath <- file.path(paste(project.dir, 
                                  "03_results",
                                 "output1", 
                                 "08_sim_output_baseline_alfonsino_fulldist.csv", sep = "/"))
mk.sdist.outpath <- file.path(paste(project.dir, 
                                    "03_results", 
                                    "output1",
                                    "alfonsino_sdist_csv", sep = "/"))

mk.sdist.filename <- "_baseline_alfonsino_fulldist"   ### CHECK HERE

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# set size distribution threshold 
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# truncation parameter
trunc.param = 0
# add run name - change this based on sdist truncation e.g. sdist truncated 25%
run.name <- "full dist"

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# load data and create lists/vectors to store model output
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

#load data
mk <- read_csv(mk.inpath)
mk <- mk %>% 
  mutate(Linf = Linf/10) %>%  # Linf is in mm in data, need to convert to cm

# allocate vectors and lists

mean.length <- matrix() #save mean lengths
max.length <- matrix() #save max lengths
median.length <- matrix() #save median lengths
sd.length <- matrix() #save sd of length
Mk.ratio <- matrix() #save m/k ratio

mk.results <- list() #save function results
plot.list <- list() #save the size dist plots to this list
df.list <-  list() # save the df used to make each plot to this list 

nmk <- dim(mk)[1]

#@@@@@@@@@@@@@@@@@@
#run function loop
#@@@@@@@@@@@@@@@@@@

for(imk in 1:nmk){
    
  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  # Input required - key life history parameters
  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    
  
    Linf <- mk$Linf[imk]                           # von Bertalanffy Linfinity (average length of fully grown individuals)
    CV.Linf <- seq(0.1,0.1 , length=mk$tmax[imk])  #changing cv of length at age
   
    t0 <- mk$t0[imk]        # von Bertalanffy t0 (theoretical length at age 0)
    tmax <- mk$tmax[imk]    # maximum age
    
    # length bins 
    L.bins <- 1 # specify length bins in cm
    
    # generate length classes, i.e. mid lengths per bin, to analyse growth rate
    L.classes <- seq(0,Linf*1.5,L.bins)+0.5*L.bins 
    
    # growth rate
    k <- mk$K[imk]    # von Bertalanffy K (k = 0.2 is default when want m/k=1.5)
    
    # natural mortality
    #Hoenig M used as starting M value (initial value, from Prince et al. 2023) and increases/decreases are based off this 
    M <- mk$M[imk] # instantaneous rate of natural mortality (per year) (M = 0.3 is default when want m/k=1.5)
   
    
    #@@@@@@@@@@@@@@@@
    # load function @
    #@@@@@@@@@@@@@@@@
    
    source("01_code/functions/ebSizeFreqEstimator_function_MAIN.R") 
    
    #@@@@@@@@@@@@@@@
    # run function @
    #@@@@@@@@@@@@@@@
    
    sim.pops <- ebSizeFreqEstimator_function_sarah(Linf,CV.Linf,t0,tmax,Lbins,L.classes,
                                                   k,M, trunc.param, run.name)

    mean.length[imk] <- sim.pops[[1]][[1]]$Unfished.Mean.Length
    max.length[imk] <- sim.pops[[1]][[1]]$Unfished.Max.Length.p99
    sd.length[imk] <- sim.pops[[1]][[1]]$Unfished.SD.Length
    median.length[imk] <- sim.pops[[1]][[1]]$Unfished.Median.Length.p50
    Mk.ratio[imk] <- sim.pops[[1]][[1]]$Mk.ratio
    
    plot_list <- sim.pops[[2]]
    df_list <- sim.pops[[3]]
    mk.results[[imk]] <- sim.pops[[1]]
    plot.list[imk] <- plot_list
    df.list[imk] <- df_list
    
} # end mk data loop

new.data <- cbind(mk, mean.length, max.length, median.length, sd.length,  Mk.ratio) ## need to make sure mk and mean.length/max.length are the same length to bind

baseline.mean <- new.data$mean.length[new.data$scenario == "baseline"]
baseline.max <- new.data$max.length[new.data$scenario == "baseline"]

new.data <- new.data %>% 
  mutate(
    change.mean = mean.length - baseline.mean ,
    change.mean.percent = (change.mean/baseline.mean)*100,
    change.max =  max.length - baseline.max ,
    change.max.percent = (change.max/baseline.mean)*100,
    CV.mean = sd.length / mean.length,
    run.name=run.name)

#write_csv(new.data, mk.outpath)

# save the sdists for plotting

# save to csv
for (i in seq_along(df.list)) {
  file_name <- paste0(mk.sdist.outpath,"/", "08_sim_output_sdist", i, mk.sdist.filename, ".csv")
  write_csv(df.list[[i]], file_name)
}
