#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Script to run the equilibrium-based size frequency estimator (ebSizeFreqEstimator_function_SENS) ##
# Version 100 by Nils Krueck, 19 April 2023 
# edited 8th april 2024 to accommodate changes to the function ebSizeFreqEstimator_function_MAIN and update data loading and writing code
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

rm(list=ls()) # clear workspace

# run simulation loop ####

#libraries
library(tidyverse)
library(here)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#set file paths to load and write data
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

##NOTE for each run - update file paths, trunc.param (line 42), run.name (line 44) and check correct function (line 126)

project.dir <- here::here()
mk.inpath <- file.path(paste(project.dir, 
                                 "03_results", 
                                 "output2", 
                                 "01_Linf_scenarios_moonwrasse.csv", sep = "/"))
mk.outpath <- file.path(paste(project.dir, 
                                 "03_results", 
                                 "output3", 
                                 "08_sim_output_Linf_cM_fCV_moonwrasse_fulldist.csv", sep = "/")) ### CHECK HERE
mk.sdist.outpath <- file.path(paste(project.dir, 
                              "03_results",
                              "output3",
                              "sdist_dfs_moonwrasse", sep = "/"))

mk.sdist.filename <- "_cM_fCV_moonwrasse_fulldist"   ### CHECK HERE

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# set size distribution parameters
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# truncation parameter        ### CHECK HERE
trunc.param = 0
# add run name - change this based on sdist truncation e.g. sdist truncated 25%
run.name <- "full dist"

#sens run name
sens.run.name <- "hjuvM fCV"  ### CHECK HERE

#list of m/k ratios for plotting size distributions ### CHECK HERE
#alfonsino
#mk.list <- c(0.335,0.766,1.345,1.788,2.017,5.375,0.575,1.009,2.682,3.583,4.034,4.705)  

#swordfish
#mk.list <- c(0.281,0.422,0.482, 0.643, 0.985,1.286, 1.313, 1.969, 2.626, 3.001,3.939, 4.501)

#moonwrasse
mk.list <- c(0.337, 0.449, 0.769, 0.897, 1.177, 1.347, 1.569, 2.019, 2.691, 3.586, 4.037, 5.380)

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

sdist.list <- list() #list for saving sdists

nmk <- dim(mk)[1]

#@@@@@@@@@@@@@@@@@@
#run function loop
#@@@@@@@@@@@@@@@@@@

for(imk in 1:nmk){
    
  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  # Input required - key life history parameters
  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    
  
    Linf <- mk$Linf[imk]                           # von Bertalanffy Linfinity (average length of fully grown individuals)
    CV.Linf <- seq(0.1,0.1 , length=mk$tmax[imk]) ### CHECK HERE #changing cv of length at age            ### CHECK HERE
   
    t0 <- mk$t0[imk]        # von Bertalanffy t0 (theoretical length at age 0)
    tmax <- mk$tmax[imk]    # maximum age
    
    # length bins 
    L.bins <- 1 # specify length bins in cm
    
    # generate length classes, i.e. mid lengths per bin, to analyse growth rate
    L.classes <- seq(0,Linf*1.5,L.bins)+0.5*L.bins 
    
    # growth rate
    k <- mk$K[imk]    # von Bertalanffy K (k = 0.2 is default when want m/k=1.5)
    
    # natural mortality - constant M
    #Hoenig M used as starting M value (initial value, from Prince et al. 2023) and increases/decreases are based off this
    #M <- mk$M[imk] # instantaneous rate of natural mortality (per year) (M = 0.3 is default when want m/k=1.5)
    
    # set M for age class  ### CHECK HERE
    M.constant <- mk$M[imk]
    # Create the vector
    M <- rep(M.constant, tmax)
    M[1] <- 2 * M.constant
    #M[2] <- 1.2 * M.constant #[2:4] for alfonsino (A50@7.4), [2] for swordfish (A50@4.42), comment this line for moonwrasse
    
    
    #@@@@@@@@@@@@@@@@
    # load function @
    #@@@@@@@@@@@@@@@@
    
    source("01_code/functions/ebSizeFreqEstimator_function_SENS_mortality.R") ### CHECK HERE
    
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
    run.name=run.name, 
    sens.run.name=sens.run.name)

write_csv(new.data, mk.outpath)


# save the sdists for plotting
# loop over each m/k ratio
for (value in mk.list) {
  mk_filtered <- lapply(df.list, function(df) {
    filtered_df <- filter(df, Mk.ratio.select == value)
    if (nrow(filtered_df) > 0) filtered_df else NULL
  })

  # Remove NULL elements from the list and store the filtered data frames
  sdist.list[[as.character(value)]] <- mk_filtered[!sapply(mk_filtered, is.null)]

} #end loop

#unnest
sdist.list <- purrr::flatten(sdist.list)

# save to csv
for (i in seq_along(sdist.list)) {
  file_name <- paste0(mk.sdist.outpath,"/", "08_sim_output_sdist", i, mk.sdist.filename, ".csv")
  write_csv(sdist.list[[i]], file_name)
}



