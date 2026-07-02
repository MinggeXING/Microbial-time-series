# NO.1 AM_Cd0 START ------
# This 'EDM_test.R' file quickly go through the whole process of running rEDM on a test dataset, before using the real data. The test dataset is from Chang et al. 2021 Ecol Lett (Ricker model). Please refer to Chang et al. on how the data set is generated. The test dataset used is in 'data/result20191024_0_0_0_.csv'.  
# For real data, the analysis is conducted in “EDM_local.R” and "EDM_AM_Cd0" file

# Preparation ----
# check version of the installed rEDM package
packageVersion("rEDM") # note Chang et al. 2021 Ecol Lett use v1.2.3. Other versions may be incompatible. 

seed <- 49563
set.seed(seed)

# load intact functions from Chang et al. 
source('code/Demo_MDR_function.R')  

# Load original dataset ----
da.range <- 1:90 # Subsample for data analysis
out.sample <- T # T/F for out-of-sample forecast
if(out.sample){nout <- 2}else{nout <- 0}  # number of out-of-sample

# First look at the first 6 species over 90 time steps ----
df <- read.csv('data/filtered_for_EDM_network/ASV_tab.AM.Cd0.abs.EDM.csv',header=T,stringsAsFactors=F)
p1 <- ggplot(data.frame(df), aes(x = 1:nrow(df))) +
  geom_line(aes(y = ASV55, color = "ASV55"), linewidth = 0.7) +
  geom_line(aes(y = ASV153, color = "ASV153"), linewidth = 0.7) +
  geom_line(aes(y = ASV178, color = "ASV178"), linewidth = 0.7) +
  geom_line(aes(y = ASV256, color = "ASV256"), linewidth = 0.7) +
  geom_line(aes(y = ASV81, color = "ASV81"), linewidth = 0.7) +
  geom_line(aes(y = ASV4, color = "ASV4"), linewidth = 0.7) +
  labs(
    title = "AM_Cd0 (first 90*6)",
    x = "New Time",
    y = "number of DNA copies"
  ) +
  scale_color_manual(
    name = NULL,
    values = c(
      "ASV55" = "#BD210F", "ASV153" = "#F15D00", "ASV178" = "#DDB500",
      "ASV256" = "#20D6B5", "ASV81" = "#0092E0", "ASV4" = "#202023"
    )
  ) +
  theme(
    plot.title = element_text(
      size = 20,
      hjust = 0.5,
      margin = margin(b = 15)
    ),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    panel.border = element_rect(
      colour = "black", 
      fill = NA, 
      linewidth = 2),
    
    legend.position = c(0.995, 0.98),          
    legend.justification = c(1, 1),           
    legend.key.width = unit(1, "cm"),       
    legend.text = element_text(size = 14),
    legend.spacing.y = unit(0.2, "cm"),      
    legend.background = element_blank(),     
    legend.margin = margin(t = 0, r = 5, b = 0, l = 5) 
  ) +
  guides(linewidth = "none")
print(p1)



(da.name <- 'ASV_tab.AM.Cd0.abs.EDM')
do <- read.csv('data/filtered_for_EDM_network/ASV_tab.AM.Cd0.abs.EDM.csv',header=T,stringsAsFactors=F)
dot <- do[da.range,1] # get time column
do <- do[da.range,-1] # remove the first column (time) from the data frame
ndo <- nrow(do) # number of time steps
nin <- ndo-nout # library sample size

# take the 149 columns (total number of ASV),filted do data frame have 149 columns
do[, 2] #see column 2 "New time" 
ncol(do) #see total columns number
do_subset <- do[, 3:151]#only take all 'ASV' columns
do <- do_subset

# ~~~~ plot logistic map ----
# plot the normalized in-sample dataset at time t. Include V1 to V6.
p2 <- ggplot(data.frame(do), aes(x = 1:nrow(do))) +
     geom_line(aes(y = ASV55, color = "ASV55"), linewidth = 0.7) +
     geom_line(aes(y = ASV153, color = "ASV153"), linewidth = 0.7) +
     geom_line(aes(y = ASV178, color = "ASV178"), linewidth = 0.7) +
     geom_line(aes(y = ASV256, color = "ASV256"), linewidth = 0.7) +
     geom_line(aes(y = ASV81, color = "ASV81"), linewidth = 0.7) +
     geom_line(aes(y = ASV4, color = "ASV4"), linewidth = 0.7) +
     labs(
       title = "AM_Cd0 (first 90*6)",
       x = "New Time",
       y = "number of DNA copies"
  ) +
  scale_color_manual(
    name = NULL,
    values = c(
      "ASV55" = "#BD210F", "ASV153" = "#F15D00", "ASV178" = "#DDB500",
      "ASV256" = "#20D6B5", "ASV81" = "#0092E0", "ASV4" = "#202023"
    )
  ) +
  theme(
    plot.title = element_text(
      size = 20,
      hjust = 0.5,
      margin = margin(b = 15)
    ),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    panel.border = element_rect(
      colour = "black", 
      fill = NA, 
      linewidth = 2),
    
    legend.position = c(0.995, 0.98),          
    legend.justification = c(1, 1),          
    legend.key.width = unit(1, "cm"),       
    legend.text = element_text(size = 14),
    legend.spacing.y = unit(0.2, "cm"),      
    legend.background = element_blank(),     
    legend.margin = margin(t = 0, r = 5, b = 0, l = 5)  
  ) +
  guides(linewidth = "none")

print(p2)



# For real data, we may have 90 time steps and 100+ 'species' for each treatment with the same matrix structure. The 'species' can be ASV/OTU ID, taxon (at phylum, order, class, family, genus, or species level), depending on nature of data and research objectives. Also, we may focus on overlapping 'species' across treatments, if we want to compare the results across treatments.

dto <- apply(do,1,sum) 
# sum of each row (time step) in the data frame

dpo <- do*repmat(dto^-1,1,ncol(do)) 
# normalize the rows of the matrix 'do' so that each row sums to 1(?). 

#This is often done in data processing to scale the data and make comparisons across different rows more meaningful. Similar to using rarefied 16S rRNA gene sequencing data(?).
# repmat() is the own function by sourcing 'Demo_MDR_function.R'. 
# repmat() makes a matrix with repeated column or row.

# check and shall see all values in the rowsum column equal to 1. 
as.data.frame(dpo) %>% mutate(rowsum = rowSums(across(everything())))


# Exclusion of rare species ----
# Here we skip this step otherwise all species will be omitted.
# We shall consider the followings for real dataset.
# *Threshold is upon decision*.
# pcri <- 0;bcri <- 10^-3; # criteria for selecting species based on proportion of presnece (pcri) and mean relative abundance (bri) 
# doind2 <- (apply(dpo,2,mean,na.rm=T)>(bcri))&((apply(do>0,2,sum,na.rm=T)/nrow(do))>pcri) # index for selected species 
# exsp2 <- setdiff(1:ncol(do),which(doind2))   # index for rare species 
# do <- do[,-exsp2]                            # Dataset excluded rare species


(nsp <- ncol(do)) # number of species
ndo <- nrow(do) # number of time steps
nin <- ndo-nout # number of in-sample time steps

# Mean and SD of each node/species ----
do.mean <- apply(do[1:nin,],2,mean,na.rm=T) 
# mean of each column (species abundance) in the first 88 rows (in-sample)

do.sd <- apply(do[1:nin,],2,sd,na.rm=T) 
# sd of each column (species abundance) in the first 88 rows (in-sample)

# Construct a sd(i,j) matrix ----
# (dimension = nsp*nsp, i.e., 149*149 here).
# if sd(i,j)>1, meaning j varies more than i
# if sd(i,j)<1, meaning i varies more than j
# if sd(i,j)=1, meaning i and j vary the same
dosdM <- repmat(c(do.sd)^-1,1,nsp)*repmat(c(do.sd),nsp,1) 

# check the sd(i,j)
dosdM; dim(dosdM)

# In-sample ----
d <- do[1:(nin-1),] # In-sample dataset at time t (time 1-87)
d_tp1 <- do[2:(nin),] # In-sample dataset at time t+1 (time 2-88)

# ~~ normalization (z-score) ----
ds <- (d-repmat(do.mean,nrow(d),1))*repmat(do.sd,nrow(d),1)^-1 
# Normalized in-sample dataset at time t (i.e., *z-score = (x-mean)/sd)

ds_tp1 <- (d_tp1-repmat(do.mean,nrow(d_tp1),1))*repmat(do.sd,nrow(d_tp1),1)^-1 
# Normalized in-sample dataset at time t+1 (i.e., z-score)


# ~~~~ plot logistic map ----
# plot the normalized in-sample dataset at time t. Include V1 to V6.
p3 <- ggplot(data.frame(ds), aes(x = 1:nrow(ds))) +
  geom_line(aes(y = ASV55, color = "ASV55"), linewidth = 0.7) +
  geom_line(aes(y = ASV153, color = "ASV153"), linewidth = 0.7) +
  geom_line(aes(y = ASV178, color = "ASV178"), linewidth = 0.7) +
  geom_line(aes(y = ASV256, color = "ASV256"), linewidth = 0.7) +
  geom_line(aes(y = ASV81, color = "ASV81"), linewidth = 0.7) +
  geom_line(aes(y = ASV4, color = "ASV4"), linewidth = 0.7) +
  labs(
    title = "In-sample (T1~87*6) AM_Cd0",
    x = "Time t",
    y = "Normalized"
  ) +
  scale_color_manual(
    name = NULL,
    values = c(
      "ASV55" = "#BD210F", "ASV153" = "#F15D00", "ASV178" = "#DDB500",
      "ASV256" = "#20D6B5", "ASV81" = "#0092E0", "ASV4" = "#202023"
    )
  ) +
  theme(
    plot.title = element_text(
      size = 20,
      hjust = 0.5,
      margin = margin(b = 15)
    ),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    panel.border = element_rect(
      colour = "black", 
      fill = NA, 
      linewidth = 2),
  
    legend.position = c(0.995, 0.98),          
    legend.justification = c(1, 1),           
    legend.key.width = unit(1, "cm"),     
    legend.text = element_text(size = 14),
    legend.spacing.y = unit(0.2, "cm"),      
    legend.background = element_blank(),      
    legend.margin = margin(t = 0, r = 5, b = 0, l = 5)  
  ) +
  guides(linewidth = "none")

print(p3)

# plot the normalized in-sample dataset at time t+1. Include V1 to V6.
p4 <- ggplot(data.frame(ds_tp1), aes(x = 1:nrow(ds_tp1))) +
  geom_line(aes(y = ASV55, color = "ASV55"), linewidth = 0.7) +
  geom_line(aes(y = ASV153, color = "ASV153"), linewidth = 0.7) +
  geom_line(aes(y = ASV178, color = "ASV178"), linewidth = 0.7) +
  geom_line(aes(y = ASV256, color = "ASV256"), linewidth = 0.7) +
  geom_line(aes(y = ASV81, color = "ASV81"), linewidth = 0.7) +
  geom_line(aes(y = ASV4, color = "ASV4"), linewidth = 0.7) +
  labs(
    title = "In-sample (T2~88*6) AM_Cd0",
    x = "Time t+1",
    y = "Normalized"
  ) +
  scale_color_manual(
    name = NULL,
    values = c(
      "ASV55" = "#BD210F", "ASV153" = "#F15D00", "ASV178" = "#DDB500",
      "ASV256" = "#20D6B5", "ASV81" = "#0092E0", "ASV4" = "#202023"
    )
  ) +
  theme(
    plot.title = element_text(
      size = 20,
      hjust = 0.5,
      margin = margin(b = 15)
    ),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    panel.border = element_rect(
      colour = "black", 
      fill = NA, 
      linewidth = 2),
    
    legend.position = c(0.995, 0.98),          
    legend.justification = c(1, 1),           
    legend.key.width = unit(1, "cm"),      
    legend.text = element_text(size = 14),
    legend.spacing.y = unit(0.2, "cm"),     
    legend.background = element_blank(),      
    legend.margin = margin(t = 0, r = 5, b = 0, l = 5)  
  ) +
  guides(linewidth = "none")

print(p4)


# do lagged point plot for t and t+1, for V1
plot_data <- data.frame(
  t = ds[,1],        # x
  t_plus_1 = ds_tp1[,1]  # y
)


p5 <- ggplot(plot_data, aes(x = t, y = t_plus_1)) +
  geom_point() +  
  labs(
    x = "t",      
    y = "t+1",  
    title = "Lagged Point Plot for species 1 (AM_Cd0)"  
  ) +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5,margin = margin(b = 15)),
    geom_point(aes(y = V1), size = 1.2),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    panel.border = element_rect(colour = "black",  fill = NA, linewidth = 2))


print(p5)

# Out-sample ----
if(out.sample|nout!=0){
  d.test <- do[nin:(ndo-1),]                 # Out-of-sample dataset at time t 
  dt_tp1 <- do[(nin+1):ndo,]                 # Out-of-sample dataset at time t+1
  ds.test <- (d.test-repmat(do.mean,nrow(d.test),1))*repmat(do.sd,nrow(d.test),1)^-1 
  # Normalized out-of-sample dataset at time t
  
  dst_tp1 <- (dt_tp1-repmat(do.mean,nrow(dt_tp1),1))*repmat(do.sd,nrow(dt_tp1),1)^-1 
  # Normalized out-of-sample dataset at time t+1
}else{d.test <- dt_tp1 <- dst_tp1 <- ds.test <- NULL}

# Compiled data at time t -> '1-87' + '88-89'
ds.all <- rbind(ds,ds.test)
dim(ds.all) 
# 89 rows and 149 columns. Since we need to have lagged dataset, number of time steps of ds.all is 89, although ndo = 90.

# Finding optimal embedding dimension (Ed) and nonlinearity parameter ----
#############################################################
# Find the optimal embedding dimension & nonlinearity parameter for each variable 
# based on univariate simplex projection and S-map, respectively

# Univariate simplex projection
Emax <- 10
cri <- 'rmse' # model selection 
Ed <- NULL
forecast_skill_simplex <- NULL
for(i in 1:ncol(ds)){
  spx.i <- simplex(ds[,i],E=2:Emax)
  Ed <- c(Ed,spx.i[which.min(spx.i[,cri])[1],'E'])
  forecast_skill_simplex <- c(forecast_skill_simplex,spx.i[which.min(spx.i[,cri])[1],'rho'])
}
Ed # The optimal embedding dimension for each variable
forecast_skill_simplex # Forecast skills for each variable based on simplex projection

######################################################################
# Finding causal variables by CCM ----
# Find causal variables by CCM analysis for multiview embedding
# Warning: It is time consuming for calculating the causation for each node
# CCM causality test for all node pairs 
# do.CCM <- F 
if(do.CCM){ 
  ccm.out <- ccm.fast.demo(ds, Epair=T,cri=cri,Emax=Emax)
  ccm.sig <- ccm.out[['ccm.sig']]
  ccm.rho <- ccm.out[['ccm.rho']]
  if(save){
# To avoid overwrite the original files, we save them with different names, 'XXX_NEW'.
    write.csv(ccm.sig, file.path("out", paste("ccm_sig_", da.name, "_nin", nin, "_demo_NEW.csv", sep="")), row.names=FALSE)
    
    write.csv(ccm.rho, file.path('out', paste('ccm_rho_',da.name,'_nin',nin,'_demo_NEW.csv',sep='')),row.names=F)
  }
}

ccm.sig <- read.csv(file.path('out',paste('ccm_sig_',da.name,'_nin',nin,'_demo_NEW.csv',sep='')),header=T,stringsAsFactors = F)
ccm.rho <- read.csv(file.path('out',paste('ccm_rho_',da.name,'_nin',nin,'_demo_NEW.csv',sep='')),header=T,stringsAsFactors = F)

# ~~~~ meaning of the two tables? ----


######################################################################
# Multiview embedding ----
# Perform multiview embedding analysis for each node/species
# Warning: It is time consuming for running multview embedding for each node/species. 
# ** ~ 1 min for one node using intel CORE i7, ThinkPad P1 in balanced mode. 
# do.multiview <- F
if(do.multiview){
  esele_lag <- esim.lag.demo(ds,ccm.rho,ccm.sig,Ed,kmax=10000,kn=100,max_lag=3,Emax=Emax)
  # To avoid overwrite the original files, we save them with different names, 'XXX_NEW'.
  if(save){write.csv(esele_lag,file.path('out', paste('eseleLag_',da.name,'_nin',nin,'_demo_NEW.csv',sep='')),row.names=F)}
}

esele <- read.csv(file.path('out', paste('eseleLag_',da.name,'_nin',nin,'_demo_NEW.csv',sep='')),header=T)


# So far so good ----


####################################################
## The computation of multiview distance
dmatrix.mv <- mvdist.demo(ds,ds.all,esele)
dmatrix.train.mvx <- dmatrix.mv[['dmatrix.train.mvx']]
dmatrix.test.mvx <- dmatrix.mv[['dmatrix.test.mvx']]

save.image(file = "AM_Cd0_workspace.RData")

# Then use server
# NO.1 AM_Cd0  END -------

######################################################
######################################################

# NO.2 AM_Cd15 START ------
# This 'EDM_test.R' file quickly go through the whole process of running rEDM on a test dataset, before using the real data. The test dataset is from Chang et al. 2021 Ecol Lett (Ricker model). Please refer to Chang et al. on how the data set is generated. The test dataset used is in 'data/result20191024_0_0_0_.csv'.  
# For real data, the analysis is conducted in “EDM_local.R” and "EDM_AM_Cd15" file

# Preparation ----
# check version of the installed rEDM package
packageVersion("rEDM") # note Chang et al. 2021 Ecol Lett use v1.2.3. Other versions may be incompatible. 

seed <- 49563
set.seed(seed)

# load intact functions from Chang et al. 
source('code/Demo_MDR_function.R')  

# Load original dataset ----
da.range <- 1:90 # Subsample for data analysis
out.sample <- T # T/F for out-of-sample forecast
if(out.sample){nout <- 2}else{nout <- 0}  # number of out-of-sample

# First look at the first 6 species over 90 time steps ----
df <- read.csv('data/filtered_for_EDM_network/ASV_tab.AM.Cd15.abs.EDM.csv',header=T,stringsAsFactors=F)
p1 <- ggplot(data.frame(df), aes(x = 1:nrow(df))) +
  geom_line(aes(y = ASV55, color = "ASV55"), linewidth = 0.7) +
  geom_line(aes(y = ASV153, color = "ASV153"), linewidth = 0.7) +
  geom_line(aes(y = ASV178, color = "ASV178"), linewidth = 0.7) +
  geom_line(aes(y = ASV256, color = "ASV256"), linewidth = 0.7) +
  geom_line(aes(y = ASV81, color = "ASV81"), linewidth = 0.7) +
  geom_line(aes(y = ASV4, color = "ASV4"), linewidth = 0.7) +
  labs(
    title = "AM_Cd15 (first 90*6)",
    x = "New Time",
    y = "number of DNA copies"
  ) +
  scale_color_manual(
    name = NULL,
    values = c(
      "ASV55" = "#BD210F", "ASV153" = "#F15D00", "ASV178" = "#DDB500",
      "ASV256" = "#20D6B5", "ASV81" = "#0092E0", "ASV4" = "#202023"
    )
  ) +
  theme(
    plot.title = element_text(
      size = 20,
      hjust = 0.5,
      margin = margin(b = 15)
    ),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    panel.border = element_rect(
      colour = "black", 
      fill = NA, 
      linewidth = 2),
 
    legend.position = c(0.995, 0.98),         
    legend.justification = c(1, 1),           
    legend.key.width = unit(1, "cm"),       
    legend.text = element_text(size = 14),
    legend.spacing.y = unit(0.2, "cm"),       
    legend.background = element_blank(),     
    legend.margin = margin(t = 0, r = 5, b = 0, l = 5)  
  ) +
  guides(linewidth = "none")
print(p1)


(da.name <- 'ASV_tab.AM.Cd15.abs.EDM')
do <- read.csv('data/filtered_for_EDM_network/ASV_tab.AM.Cd15.abs.EDM.csv',header=T,stringsAsFactors=F)
dot <- do[da.range,2] # get time column
do <- do[da.range,-2] # remove the first 2 column (time) from the data frame
ndo <- nrow(do) # number of time steps
nin <- ndo-nout # library sample size

# take the 130 columns (total number of ASV),filted do data frame have 130 columns
do[, 2] #see column 2 "New time"
ncol(do) #see total columns 
do_subset <- do[, 3:132]#only take all 'ASV' columns
do <- do_subset

# ~~~~ plot logistic map ----
# plot the normalized in-sample dataset at time t. Include V1 to V6.
p2 <- ggplot(data.frame(do), aes(x = 1:nrow(do))) +
  geom_line(aes(y = ASV55, color = "ASV55"), linewidth = 0.7) +
  geom_line(aes(y = ASV153, color = "ASV153"), linewidth = 0.7) +
  geom_line(aes(y = ASV178, color = "ASV178"), linewidth = 0.7) +
  geom_line(aes(y = ASV256, color = "ASV256"), linewidth = 0.7) +
  geom_line(aes(y = ASV81, color = "ASV81"), linewidth = 0.7) +
  geom_line(aes(y = ASV4, color = "ASV4"), linewidth = 0.7) +
  labs(
    title = "AM_Cd15 (first 90*6)",
    x = "New Time",
    y = "number of DNA copies"
  ) +
  scale_color_manual(
    name = NULL,
    values = c(
      "ASV55" = "#BD210F", "ASV153" = "#F15D00", "ASV178" = "#DDB500",
      "ASV256" = "#20D6B5", "ASV81" = "#0092E0", "ASV4" = "#202023"
    )
  ) +
  theme(
    plot.title = element_text(
      size = 20,
      hjust = 0.5,
      margin = margin(b = 15)
    ),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    panel.border = element_rect(
      colour = "black", 
      fill = NA, 
      linewidth = 2),

    legend.position = c(0.995, 0.98),         
    legend.justification = c(1, 1),           
    legend.key.width = unit(1, "cm"),      
    legend.text = element_text(size = 14),
    legend.spacing.y = unit(0.2, "cm"),      
    legend.background = element_blank(),     
    legend.margin = margin(t = 0, r = 5, b = 0, l = 5)  
  ) +
  guides(linewidth = "none")

print(p2)



# For real data, we may have 90 time steps and 100+ 'species' for each treatment with the same matrix structure. The 'species' can be ASV/OTU ID, taxon (at phylum, order, class, family, genus, or species level), depending on nature of data and research objectives. Also, we may focus on overlapping 'species' across treatments, if we want to compare the results across treatments.

dto <- apply(do,1,sum) # sum of each row (time step) in the data frame
dpo <- do*repmat(dto^-1,1,ncol(do)) # normalize the rows of the matrix 'do' so that each row sums to 1(?). This is often done in data processing to scale the data and make comparisons across different rows more meaningful. Similar to using rarefied 16S rRNA gene sequencing data(?).
# repmat() is the own function by sourcing 'Demo_MDR_function.R'. 
# repmat() makes a matrix with repeated column or row.

# check and shall see all values in the rowsum column equal to 1. 
as.data.frame(dpo) %>% mutate(rowsum = rowSums(across(everything())))


# Exclusion of rare species ----
# Here we skip this step otherwise all species will be omitted.
# We shall consider the followings for real dataset.
# *Threshold is upon decision*.
# pcri <- 0;bcri <- 10^-3; # criteria for selecting species based on proportion of presnece (pcri) and mean relative abundance (bri) 
# doind2 <- (apply(dpo,2,mean,na.rm=T)>(bcri))&((apply(do>0,2,sum,na.rm=T)/nrow(do))>pcri) # index for selected species 
# exsp2 <- setdiff(1:ncol(do),which(doind2))   # index for rare species 
# do <- do[,-exsp2]                            # Dataset excluded rare species


(nsp <- ncol(do)) # number of species
ndo <- nrow(do) # number of time steps
nin <- ndo-nout # number of in-sample time steps

# Mean and SD of each node/species ----
do.mean <- apply(do[1:nin,],2,mean,na.rm=T) # mean of each column (species abundance) in the first 88 rows (in-sample)
do.sd <- apply(do[1:nin,],2,sd,na.rm=T) # sd of each column (species abundance) in the first 88 rows (in-sample)

# Construct a sd(i,j) matrix ----
# (dimension = nsp*nsp, i.e., 6*6 here).
# if sd(i,j)>1, meaning j varies more than i
# if sd(i,j)<1, meaning i varies more than j
# if sd(i,j)=1, meaning i and j vary the same
dosdM <- repmat(c(do.sd)^-1,1,nsp)*repmat(c(do.sd),nsp,1) 

# check the sd(i,j)
dosdM; dim(dosdM)

# In-sample ----
d <- do[1:(nin-1),] # In-sample dataset at time t (time 1-87)
d_tp1 <- do[2:(nin),] # In-sample dataset at time t+1 (time 2-88)

# ~~ normalization (z-score) ----
ds <- (d-repmat(do.mean,nrow(d),1))*repmat(do.sd,nrow(d),1)^-1 # Normalized in-sample dataset at time t (i.e., *z-score = (x-mean)/sd)
ds_tp1 <- (d_tp1-repmat(do.mean,nrow(d_tp1),1))*repmat(do.sd,nrow(d_tp1),1)^-1 # Normalized in-sample dataset at time t+1 (i.e., z-score)

# ~~~~ plot logistic map ----
# plot the normalized in-sample dataset at time t. Include V1 to V6.
p3 <- ggplot(data.frame(ds), aes(x = 1:nrow(ds))) +
  geom_line(aes(y = ASV55, color = "ASV55"), linewidth = 0.7) +
  geom_line(aes(y = ASV153, color = "ASV153"), linewidth = 0.7) +
  geom_line(aes(y = ASV178, color = "ASV178"), linewidth = 0.7) +
  geom_line(aes(y = ASV256, color = "ASV256"), linewidth = 0.7) +
  geom_line(aes(y = ASV81, color = "ASV81"), linewidth = 0.7) +
  geom_line(aes(y = ASV4, color = "ASV4"), linewidth = 0.7) +
  labs(
    title = "In-sample (T1~87*6) AM_Cd15",
    x = "Time t",
    y = "Normalized"
  ) +
  scale_color_manual(
    name = NULL,
    values = c(
      "ASV55" = "#BD210F", "ASV153" = "#F15D00", "ASV178" = "#DDB500",
      "ASV256" = "#20D6B5", "ASV81" = "#0092E0", "ASV4" = "#202023"
    )
  ) +
  theme(
    plot.title = element_text(
      size = 20,
      hjust = 0.5,
      margin = margin(b = 15)
    ),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    panel.border = element_rect(
      colour = "black", 
      fill = NA, 
      linewidth = 2),
 
    legend.position = c(0.995, 0.98),          
    legend.justification = c(1, 1),          
    legend.key.width = unit(1, "cm"),       
    legend.text = element_text(size = 14),
    legend.spacing.y = unit(0.2, "cm"),       
    legend.background = element_blank(),    
    legend.margin = margin(t = 0, r = 5, b = 0, l = 5)  
  ) +
  guides(linewidth = "none")

print(p3)


# we can do t+1 as well, but since the pattern is the same so we skip. 
# but we can do lagged point plot for t and t+1, for V1

plot_data <- data.frame(
  t = ds[,1],        
  t_plus_1 = ds_tp1[,1]  
)


p4 <- ggplot(plot_data, aes(x = t, y = t_plus_1)) +
  geom_point() +  
  labs(
    x = "t",     
    y = "t+1",   
    title = "Lagged Point Plot for species 1 (AM_Cd15)" 
  ) +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5,margin = margin(b = 15)),
    geom_point(aes(y = V1), size = 1.2),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    panel.border = element_rect(colour = "black",  fill = NA, linewidth = 2))


print(p4)


# Out-sample ----
if(out.sample|nout!=0){
  d.test <- do[nin:(ndo-1),]                 # Out-of-sample dataset at time t 
  dt_tp1 <- do[(nin+1):ndo,]                 # Out-of-sample dataset at time t+1
  ds.test <- (d.test-repmat(do.mean,nrow(d.test),1))*repmat(do.sd,nrow(d.test),1)^-1 # Normalized out-of-sample dataset at time t
  dst_tp1 <- (dt_tp1-repmat(do.mean,nrow(dt_tp1),1))*repmat(do.sd,nrow(dt_tp1),1)^-1 # Normalized out-of-sample dataset at time t+1
}else{d.test <- dt_tp1 <- dst_tp1 <- ds.test <- NULL}

# Compiled data at time t -> '1-87' + '88-89'
ds.all <- rbind(ds,ds.test)
dim(ds.all) # 89 rows and 130 columns. Since we need to have lagged dataset, number of time steps of ds.all is 89, although ndo = 30.

# Finding optimal embedding dimension (Ed) and nonlinearity parameter ----
#############################################################
# Find the optimal embedding dimension & nonlinearity parameter for each variable 
# based on univariate simplex projection and S-map, respectively

# Univariate simplex projection
Emax <- 10
cri <- 'rmse' # model selection 
Ed <- NULL
forecast_skill_simplex <- NULL
for(i in 1:ncol(ds)){
  spx.i <- simplex(ds[,i],E=2:Emax)
  Ed <- c(Ed,spx.i[which.min(spx.i[,cri])[1],'E'])
  forecast_skill_simplex <- c(forecast_skill_simplex,spx.i[which.min(spx.i[,cri])[1],'rho'])
}
Ed # The optimal embedding dimension for each variable
forecast_skill_simplex # Forecast skills for each variable based on simplex projection

######################################################################
# Finding causal variables by CCM ----
# Find causal variables by CCM analysis for multiview embedding
# Warning: It is time consuming for calculating the causation for each node
# CCM causality test for all node pairs 
# do.CCM <- F 
if(do.CCM){ 
  ccm.out <- ccm.fast.demo(ds, Epair=T,cri=cri,Emax=Emax)
  ccm.sig <- ccm.out[['ccm.sig']]
  ccm.rho <- ccm.out[['ccm.rho']]
  if(save){
    # To avoid overwrite the original files, we save them with different names, 'XXX_NEW'.
    write.csv(ccm.sig, file.path("out", paste("ccm_sig_", da.name, "_nin", nin, "_demo_NEW.csv", sep="")), row.names=FALSE)
    
    write.csv(ccm.rho, file.path('out', paste('ccm_rho_',da.name,'_nin',nin,'_demo_NEW.csv',sep='')),row.names=F)
  }
}

ccm.sig <- read.csv(file.path('out',paste('ccm_sig_',da.name,'_nin',nin,'_demo_NEW.csv',sep='')),header=T,stringsAsFactors = F)
ccm.rho <- read.csv(file.path('out',paste('ccm_rho_',da.name,'_nin',nin,'_demo_NEW.csv',sep='')),header=T,stringsAsFactors = F)

# ~~~~ meaning of the two tables? ----


######################################################################
# Multiview embedding ----
# Perform multiview embedding analysis for each node/species
# Warning: It is time consuming for running multview embedding for each node/species. 
# ** ~ 1 min for one node using intel CORE i7, ThinkPad P1 in balanced mode. 
# do.multiview <- F
if(do.multiview){
  esele_lag <- esim.lag.demo(ds,ccm.rho,ccm.sig,Ed,kmax=10000,kn=100,max_lag=3,Emax=Emax)
  # To avoid overwrite the original files, we save them with different names, 'XXX_NEW'.
  if(save){write.csv(esele_lag,file.path('out', paste('eseleLag_',da.name,'_nin',nin,'_demo_NEW.csv',sep='')),row.names=F)}
}

esele <- read.csv(file.path('out', paste('eseleLag_',da.name,'_nin',nin,'_demo_NEW.csv',sep='')),header=T)


# So far so good ----


####################################################
## The computation of multiview distance
dmatrix.mv <- mvdist.demo(ds,ds.all,esele)
dmatrix.train.mvx <- dmatrix.mv[['dmatrix.train.mvx']]
dmatrix.test.mvx <- dmatrix.mv[['dmatrix.test.mvx']]

save.image(file = "AM_Cd15_workspace.RData")

# Then use server
# NO.2 AM_Cd15  END -------

######################################################
######################################################

# NO.3 AM_Cd2 START ------
# This 'EDM_test.R' file quickly go through the whole process of running rEDM on a test dataset, before using the real data. The test dataset is from Chang et al. 2021 Ecol Lett (Ricker model). Please refer to Chang et al. on how the data set is generated. The test dataset used is in 'data/result20191024_0_0_0_.csv'.  
# For real data, the analysis is conducted “EDM_local.R” and in "EDM_AM_Cd2.R" file

# Preparation ----
# check version of the installed rEDM package
packageVersion("rEDM") # note Chang et al. 2021 Ecol Lett use v1.2.3. Other versions may be incompatible. 

seed <- 49563
set.seed(seed)

# load intact functions from Chang et al. 
source('code/Demo_MDR_function.R')  

# Load original dataset ----
da.range <- 1:90 # Subsample for data analysis
out.sample <- T # T/F for out-of-sample forecast
if(out.sample){nout <- 2}else{nout <- 0}  # number of out-of-sample

# First look at the first 6 species over 90 time steps ----
df <- read.csv('data/filtered_for_EDM_network/ASV_tab.AM.Cd2.abs.EDM.csv',header=T,stringsAsFactors=F)
p1 <- ggplot(data.frame(df), aes(x = 1:nrow(df))) +
  geom_line(aes(y = ASV55, color = "ASV55"), linewidth = 0.7) +
  geom_line(aes(y = ASV153, color = "ASV153"), linewidth = 0.7) +
  geom_line(aes(y = ASV178, color = "ASV178"), linewidth = 0.7) +
  geom_line(aes(y = ASV1211, color = "ASV1211"), linewidth = 0.7) +
  geom_line(aes(y = ASV358, color = "ASV358"), linewidth = 0.7) +
  geom_line(aes(y = ASV4, color = "ASV4"), linewidth = 0.7) +
  labs(
    title = "AM_Cd2 (first 90*6)",
    x = "New Time",
    y = "number of DNA copies"
  ) +
  scale_color_manual(
    name = NULL,
    values = c(
      "ASV55" = "#BD210F", "ASV153" = "#F15D00", "ASV178" = "#DDB500",
      "ASV1211" = "#20D6B5", "ASV358" = "#0092E0", "ASV4" = "#202023"
    )
  ) +
  theme(
    plot.title = element_text(
      size = 20,
      hjust = 0.5,
      margin = margin(b = 15)
    ),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    panel.border = element_rect(
      colour = "black", 
      fill = NA, 
      linewidth = 2),

    legend.position = c(0.995, 0.98),         
    legend.justification = c(1, 1),         
    legend.key.width = unit(1, "cm"),       
    legend.text = element_text(size = 14),
    legend.spacing.y = unit(0.2, "cm"),       
    legend.background = element_blank(),     
    legend.margin = margin(t = 0, r = 5, b = 0, l = 5)  
  ) +
  guides(linewidth = "none")
print(p1)



(da.name <- 'ASV_tab.AM.Cd2.abs.EDM')
do <- read.csv('data/filtered_for_EDM_network/ASV_tab.AM.Cd2.abs.EDM.csv',header=T,stringsAsFactors=F)
dot <- do[da.range,2] # get time column
do <- do[da.range,-2] # remove the first 2 column (time) from the data frame
ndo <- nrow(do) # number of time steps
nin <- ndo-nout # library sample size

# take the 141 columns (total number of ASV),filted do data frame have 141 columns
do[, 2] #see column 2 "New time"
ncol(do) #see total columns 
do_subset <- do[, 3:143]#only take all 'ASV' columns
do <- do_subset

# ~~~~ plot logistic map ----
# plot the normalized in-sample dataset at time t. Include V1 to V6.
p2 <- ggplot(data.frame(do), aes(x = 1:nrow(do))) +
  geom_line(aes(y = ASV55, color = "ASV55"), linewidth = 0.7) +
  geom_line(aes(y = ASV153, color = "ASV153"), linewidth = 0.7) +
  geom_line(aes(y = ASV178, color = "ASV178"), linewidth = 0.7) +
  geom_line(aes(y = ASV1211, color = "ASV1211"), linewidth = 0.7) +
  geom_line(aes(y = ASV358, color = "ASV358"), linewidth = 0.7) +
  geom_line(aes(y = ASV4, color = "ASV4"), linewidth = 0.7) +
  labs(
    title = "AM_Cd2 (first 90*6)",
    x = "New Time",
    y = "number of DNA copies"
  ) +
  scale_color_manual(
    name = NULL,
    values = c(
      "ASV55" = "#BD210F", "ASV153" = "#F15D00", "ASV178" = "#DDB500",
      "ASV1211" = "#20D6B5", "ASV358" = "#0092E0", "ASV4" = "#202023"
    )
  ) +
  theme(
    plot.title = element_text(
      size = 20,
      hjust = 0.5,
      margin = margin(b = 15)
    ),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    panel.border = element_rect(
      colour = "black", 
      fill = NA, 
      linewidth = 2),
  
    legend.position = c(0.995, 0.98),          
    legend.justification = c(1, 1),           
    legend.key.width = unit(1, "cm"),       
    legend.text = element_text(size = 14),
    legend.spacing.y = unit(0.2, "cm"),       
    legend.background = element_blank(),      
    legend.margin = margin(t = 0, r = 5, b = 0, l = 5)  
  ) +
  guides(linewidth = "none")

print(p2)



# For real data, we may have 90 time steps and 100+ 'species' for each treatment with the same matrix structure. The 'species' can be ASV/OTU ID, taxon (at phylum, order, class, family, genus, or species level), depending on nature of data and research objectives. Also, we may focus on overlapping 'species' across treatments, if we want to compare the results across treatments.

dto <- apply(do,1,sum) # sum of each row (time step) in the data frame
dpo <- do*repmat(dto^-1,1,ncol(do)) # normalize the rows of the matrix 'do' so that each row sums to 1(?). This is often done in data processing to scale the data and make comparisons across different rows more meaningful. Similar to using rarefied 16S rRNA gene sequencing data(?).
# repmat() is the own function by sourcing 'Demo_MDR_function.R'. 
# repmat() makes a matrix with repeated column or row.

# check and shall see all values in the rowsum column equal to 1. 
as.data.frame(dpo) %>% mutate(rowsum = rowSums(across(everything())))


# Exclusion of rare species ----
# Here we skip this step otherwise all species will be omitted.
# We shall consider the followings for real dataset.
# *Threshold is upon decision*.
# pcri <- 0;bcri <- 10^-3; # criteria for selecting species based on proportion of presnece (pcri) and mean relative abundance (bri) 
# doind2 <- (apply(dpo,2,mean,na.rm=T)>(bcri))&((apply(do>0,2,sum,na.rm=T)/nrow(do))>pcri) # index for selected species 
# exsp2 <- setdiff(1:ncol(do),which(doind2))   # index for rare species 
# do <- do[,-exsp2]                            # Dataset excluded rare species


(nsp <- ncol(do)) # number of species
ndo <- nrow(do) # number of time steps
nin <- ndo-nout # number of in-sample time steps

# Mean and SD of each node/species ----
do.mean <- apply(do[1:nin,],2,mean,na.rm=T) # mean of each column (species abundance) in the first 88 rows (in-sample)
do.sd <- apply(do[1:nin,],2,sd,na.rm=T) # sd of each column (species abundance) in the first 88 rows (in-sample)

# Construct a sd(i,j) matrix ----
# (dimension = nsp*nsp, i.e., 6*6 here).
# if sd(i,j)>1, meaning j varies more than i
# if sd(i,j)<1, meaning i varies more than j
# if sd(i,j)=1, meaning i and j vary the same
dosdM <- repmat(c(do.sd)^-1,1,nsp)*repmat(c(do.sd),nsp,1) 

# check the sd(i,j)
dosdM; dim(dosdM)

# In-sample ----
d <- do[1:(nin-1),] # In-sample dataset at time t (time 1-87)
d_tp1 <- do[2:(nin),] # In-sample dataset at time t+1 (time 2-88)

# ~~ normalization (z-score) ----
ds <- (d-repmat(do.mean,nrow(d),1))*repmat(do.sd,nrow(d),1)^-1 # Normalized in-sample dataset at time t (i.e., *z-score = (x-mean)/sd)
ds_tp1 <- (d_tp1-repmat(do.mean,nrow(d_tp1),1))*repmat(do.sd,nrow(d_tp1),1)^-1 # Normalized in-sample dataset at time t+1 (i.e., z-score)

# ~~~~ plot logistic map ----
# plot the normalized in-sample dataset at time t. Include V1 to V6.
p3 <- ggplot(data.frame(ds), aes(x = 1:nrow(ds))) +
  geom_line(aes(y = ASV55, color = "ASV55"), linewidth = 0.7) +
  geom_line(aes(y = ASV153, color = "ASV153"), linewidth = 0.7) +
  geom_line(aes(y = ASV178, color = "ASV178"), linewidth = 0.7) +
  geom_line(aes(y = ASV1211, color = "ASV1211"), linewidth = 0.7) +
  geom_line(aes(y = ASV358, color = "ASV358"), linewidth = 0.7) +
  geom_line(aes(y = ASV4, color = "ASV4"), linewidth = 0.7) +
  labs(
    title = "In-sample (T1~87*6) AM_Cd2",
    x = "Time t",
    y = "Normalized"
  ) +
  scale_color_manual(
    name = NULL,
    values = c(
      "ASV55" = "#BD210F", "ASV153" = "#F15D00", "ASV178" = "#DDB500",
      "ASV1211" = "#20D6B5", "ASV358" = "#0092E0", "ASV4" = "#202023"
    )
  ) +
  theme(
    plot.title = element_text(
      size = 20,
      hjust = 0.5,
      margin = margin(b = 15)
    ),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    panel.border = element_rect(
      colour = "black", 
      fill = NA, 
      linewidth = 2),

    legend.position = c(0.995, 0.98),         
    legend.justification = c(1, 1),           
    legend.key.width = unit(1, "cm"),      
    legend.text = element_text(size = 14),
    legend.spacing.y = unit(0.2, "cm"),       
    legend.background = element_blank(),      
    legend.margin = margin(t = 0, r = 5, b = 0, l = 5)  
  ) +
  guides(linewidth = "none")

print(p3)


# we can do t+1 as well, but since the pattern is the same so we skip. 
# but we can do lagged point plot for t and t+1, for V1

plot_data <- data.frame(
  t = ds[,1],        
  t_plus_1 = ds_tp1[,1]  
)


p4 <- ggplot(plot_data, aes(x = t, y = t_plus_1)) +
  geom_point() +  
  labs(
    x = "t",      
    y = "t+1",    
    title = "Lagged Point Plot for species 1 (AM_Cd2)"  
  ) +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5,margin = margin(b = 15)),
    geom_point(aes(y = V1), size = 1.2),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    panel.border = element_rect(colour = "black",  fill = NA, linewidth = 2))


print(p4)


# Out-sample ----
if(out.sample|nout!=0){
  d.test <- do[nin:(ndo-1),]                 # Out-of-sample dataset at time t 
  dt_tp1 <- do[(nin+1):ndo,]                 # Out-of-sample dataset at time t+1
  ds.test <- (d.test-repmat(do.mean,nrow(d.test),1))*repmat(do.sd,nrow(d.test),1)^-1 # Normalized out-of-sample dataset at time t
  dst_tp1 <- (dt_tp1-repmat(do.mean,nrow(dt_tp1),1))*repmat(do.sd,nrow(dt_tp1),1)^-1 # Normalized out-of-sample dataset at time t+1
}else{d.test <- dt_tp1 <- dst_tp1 <- ds.test <- NULL}

# Compiled data at time t -> '1-87' + '88-89'
ds.all <- rbind(ds,ds.test)
dim(ds.all) # 89 rows and 141 columns. Since we need to have lagged dataset, number of time steps of ds.all is 89, although ndo = 30.

# Finding optimal embedding dimension (Ed) and nonlinearity parameter ----
#############################################################
# Find the optimal embedding dimension & nonlinearity parameter for each variable 
# based on univariate simplex projection and S-map, respectively

# Univariate simplex projection
Emax <- 10
cri <- 'rmse' # model selection 
Ed <- NULL
forecast_skill_simplex <- NULL
for(i in 1:ncol(ds)){
  spx.i <- simplex(ds[,i],E=2:Emax)
  Ed <- c(Ed,spx.i[which.min(spx.i[,cri])[1],'E'])
  forecast_skill_simplex <- c(forecast_skill_simplex,spx.i[which.min(spx.i[,cri])[1],'rho'])
}
Ed # The optimal embedding dimension for each variable
forecast_skill_simplex # Forecast skills for each variable based on simplex projection

######################################################################
# Finding causal variables by CCM ----
# Find causal variables by CCM analysis for multiview embedding
# Warning: It is time consuming for calculating the causation for each node
# CCM causality test for all node pairs 
# do.CCM <- F 
if(do.CCM){ 
  ccm.out <- ccm.fast.demo(ds, Epair=T,cri=cri,Emax=Emax)
  ccm.sig <- ccm.out[['ccm.sig']]
  ccm.rho <- ccm.out[['ccm.rho']]
  if(save){
    # To avoid overwrite the original files, we save them with different names, 'XXX_NEW'.
    write.csv(ccm.sig, file.path("out", paste("ccm_sig_", da.name, "_nin", nin, "_demo_NEW.csv", sep="")), row.names=FALSE)
    
    write.csv(ccm.rho, file.path('out', paste('ccm_rho_',da.name,'_nin',nin,'_demo_NEW.csv',sep='')),row.names=F)
  }
}

ccm.sig <- read.csv(file.path('out',paste('ccm_sig_',da.name,'_nin',nin,'_demo_NEW.csv',sep='')),header=T,stringsAsFactors = F)
ccm.rho <- read.csv(file.path('out',paste('ccm_rho_',da.name,'_nin',nin,'_demo_NEW.csv',sep='')),header=T,stringsAsFactors = F)

# ~~~~ meaning of the two tables? ----


######################################################################
# Multiview embedding ----
# Perform multiview embedding analysis for each node/species
# Warning: It is time consuming for running multview embedding for each node/species. 
# ** ~ 1 min for one node using intel CORE i7, ThinkPad P1 in balanced mode. 
# do.multiview <- F
if(do.multiview){
  esele_lag <- esim.lag.demo(ds,ccm.rho,ccm.sig,Ed,kmax=10000,kn=100,max_lag=3,Emax=Emax)
  # To avoid overwrite the original files, we save them with different names, 'XXX_NEW'.
  if(save){write.csv(esele_lag,file.path('out', paste('eseleLag_',da.name,'_nin',nin,'_demo_NEW.csv',sep='')),row.names=F)}
}

esele <- read.csv(file.path('out', paste('eseleLag_',da.name,'_nin',nin,'_demo_NEW.csv',sep='')),header=T)


# So far so good ----


####################################################
## The computation of multiview distance
dmatrix.mv <- mvdist.demo(ds,ds.all,esele)
dmatrix.train.mvx <- dmatrix.mv[['dmatrix.train.mvx']]
dmatrix.test.mvx <- dmatrix.mv[['dmatrix.test.mvx']]

save.image(file = "AM_Cd2_workspace.RData")

# Then use server
# NO.3 AM_Cd2  END -------

######################################################
######################################################
# NO.4 AM_Cd5 START ------
# The 'EDM_test.R' file quickly go through the whole process of running rEDM on a test dataset, before using the real data. The test dataset is from Chang et al. 2021 Ecol Lett (Ricker model). Please refer to Chang et al. on how the data set is generated. The test dataset used is in 'data/result20191024_0_0_0_.csv'.  
# For real data, the analysis is conducted in “EDM_local.R” and "EDM_AM_Cd5.R" file

# Preparation ----
# check version of the installed rEDM package
packageVersion("rEDM") # note Chang et al. 2021 Ecol Lett use v1.2.3. Other versions may be incompatible. 

seed <- 49563
set.seed(seed)

# load intact functions from Chang et al. 
source('code/Demo_MDR_function.R')  

# Load original dataset ----
da.range <- 1:90 # Subsample for data analysis
out.sample <- T # T/F for out-of-sample forecast
if(out.sample){nout <- 2}else{nout <- 0}  # number of out-of-sample

# First look at the first 6 species over 90 time steps ----
df <- read.csv('data/filtered_for_EDM_network/ASV_tab.AM.Cd5.abs.EDM.csv',header=T,stringsAsFactors=F)
p1 <- ggplot(data.frame(df), aes(x = 1:nrow(df))) +
  geom_line(aes(y = ASV55, color = "ASV55"), linewidth = 0.7) +
  geom_line(aes(y = ASV153, color = "ASV153"), linewidth = 0.7) +
  geom_line(aes(y = ASV178, color = "ASV178"), linewidth = 0.7) +
  geom_line(aes(y = ASV277, color = "ASV277"), linewidth = 0.7) +
  geom_line(aes(y = ASV358, color = "ASV358"), linewidth = 0.7) +
  geom_line(aes(y = ASV4, color = "ASV4"), linewidth = 0.7) +
  labs(
    title = "AM_Cd5 (first 90*6)",
    x = "New Time",
    y = "number of DNA copies"
  ) +
  scale_color_manual(
    name = NULL,
    values = c(
      "ASV55" = "#BD210F", "ASV153" = "#F15D00", "ASV178" = "#DDB500",
      "ASV277" = "#20D6B5", "ASV358" = "#0092E0", "ASV4" = "#202023"
    )
  ) +
  theme(
    plot.title = element_text(
      size = 20,
      hjust = 0.5,
      margin = margin(b = 15)
    ),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    panel.border = element_rect(
      colour = "black", 
      fill = NA, 
      linewidth = 2),

    legend.position = c(0.995, 0.98),         
    legend.justification = c(1, 1),         
    legend.key.width = unit(1, "cm"),     
    legend.text = element_text(size = 14),
    legend.spacing.y = unit(0.2, "cm"),      
    legend.background = element_blank(),    
    legend.margin = margin(t = 0, r = 5, b = 0, l = 5) 
  ) +
  guides(linewidth = "none")
print(p1)

(da.name <- 'ASV_tab.AM.Cd5.abs.EDM')
do <- read.csv('data/filtered_for_EDM_network/ASV_tab.AM.Cd5.abs.EDM.csv',header=T,stringsAsFactors=F)
dot <- do[da.range,2] # get time column
do <- do[da.range,-2] # remove the first 2 column (time) from the data frame
ndo <- nrow(do) # number of time steps
nin <- ndo-nout # library sample size

# take the 135 columns (total number of ASV),filted do data frame have 135 columns
do[, 2] #see column 2 "New time"
ncol(do) #see total columns 
do_subset <- do[, 3:137]#only take all 'ASV' columns
do <- do_subset

# ~~~~ plot logistic map ----
# plot the normalized in-sample dataset at time t. Include V1 to V6.
p2 <- ggplot(data.frame(do), aes(x = 1:nrow(do))) +
  geom_line(aes(y = ASV55, color = "ASV55"), linewidth = 0.7) +
  geom_line(aes(y = ASV153, color = "ASV153"), linewidth = 0.7) +
  geom_line(aes(y = ASV178, color = "ASV178"), linewidth = 0.7) +
  geom_line(aes(y = ASV277, color = "ASV277"), linewidth = 0.7) +
  geom_line(aes(y = ASV358, color = "ASV358"), linewidth = 0.7) +
  geom_line(aes(y = ASV4, color = "ASV4"), linewidth = 0.7) +
  labs(
    title = "AM_Cd5 (first 90*6)",
    x = "New Time",
    y = "number of DNA copies"
  ) +
  scale_color_manual(
    name = NULL,
    values = c(
      "ASV55" = "#BD210F", "ASV153" = "#F15D00", "ASV178" = "#DDB500",
      "ASV277" = "#20D6B5", "ASV358" = "#0092E0", "ASV4" = "#202023"
    )
  ) +
  theme(
    plot.title = element_text(
      size = 20,
      hjust = 0.5,
      margin = margin(b = 15)
    ),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    panel.border = element_rect(
      colour = "black", 
      fill = NA, 
      linewidth = 2),
   
    legend.position = c(0.995, 0.98),         
    legend.justification = c(1, 1),          
    legend.key.width = unit(1, "cm"),       
    legend.text = element_text(size = 14),
    legend.spacing.y = unit(0.2, "cm"),     
    legend.background = element_blank(),     
    legend.margin = margin(t = 0, r = 5, b = 0, l = 5)  
  ) +
  guides(linewidth = "none")

print(p2)



# For real data, we may have 90 time steps and 100+ 'species' for each treatment with the same matrix structure. The 'species' can be ASV/OTU ID, taxon (at phylum, order, class, family, genus, or species level), depending on nature of data and research objectives. Also, we may focus on overlapping 'species' across treatments, if we want to compare the results across treatments.

dto <- apply(do,1,sum) # sum of each row (time step) in the data frame
dpo <- do*repmat(dto^-1,1,ncol(do)) # normalize the rows of the matrix 'do' so that each row sums to 1(?). This is often done in data processing to scale the data and make comparisons across different rows more meaningful. Similar to using rarefied 16S rRNA gene sequencing data(?).
# repmat() is the own function by sourcing 'Demo_MDR_function.R'. 
# repmat() makes a matrix with repeated column or row.

# check and shall see all values in the rowsum column equal to 1. 
as.data.frame(dpo) %>% mutate(rowsum = rowSums(across(everything())))


# Exclusion of rare species ----
# Here we skip this step otherwise all species will be omitted.
# We shall consider the followings for real dataset.
# *Threshold is upon decision*.
# pcri <- 0;bcri <- 10^-3; # criteria for selecting species based on proportion of presnece (pcri) and mean relative abundance (bri) 
# doind2 <- (apply(dpo,2,mean,na.rm=T)>(bcri))&((apply(do>0,2,sum,na.rm=T)/nrow(do))>pcri) # index for selected species 
# exsp2 <- setdiff(1:ncol(do),which(doind2))   # index for rare species 
# do <- do[,-exsp2]                            # Dataset excluded rare species


(nsp <- ncol(do)) # number of species
ndo <- nrow(do) # number of time steps
nin <- ndo-nout # number of in-sample time steps

# Mean and SD of each node/species ----
do.mean <- apply(do[1:nin,],2,mean,na.rm=T) # mean of each column (species abundance) in the first 88 rows (in-sample)
do.sd <- apply(do[1:nin,],2,sd,na.rm=T) # sd of each column (species abundance) in the first 88 rows (in-sample)

# Construct a sd(i,j) matrix ----
# (dimension = nsp*nsp, i.e., 6*6 here).
# if sd(i,j)>1, meaning j varies more than i
# if sd(i,j)<1, meaning i varies more than j
# if sd(i,j)=1, meaning i and j vary the same
dosdM <- repmat(c(do.sd)^-1,1,nsp)*repmat(c(do.sd),nsp,1) 

# check the sd(i,j)
dosdM; dim(dosdM)

# In-sample ----
d <- do[1:(nin-1),] # In-sample dataset at time t (time 1-87)
d_tp1 <- do[2:(nin),] # In-sample dataset at time t+1 (time 2-88)

# ~~ normalization (z-score) ----
ds <- (d-repmat(do.mean,nrow(d),1))*repmat(do.sd,nrow(d),1)^-1 # Normalized in-sample dataset at time t (i.e., *z-score = (x-mean)/sd)
ds_tp1 <- (d_tp1-repmat(do.mean,nrow(d_tp1),1))*repmat(do.sd,nrow(d_tp1),1)^-1 # Normalized in-sample dataset at time t+1 (i.e., z-score)

# ~~~~ plot logistic map ----
# plot the normalized in-sample dataset at time t. Include V1 to V6.
p3 <- ggplot(data.frame(ds), aes(x = 1:nrow(ds))) +
  geom_line(aes(y = ASV55, color = "ASV55"), linewidth = 0.7) +
  geom_line(aes(y = ASV153, color = "ASV153"), linewidth = 0.7) +
  geom_line(aes(y = ASV178, color = "ASV178"), linewidth = 0.7) +
  geom_line(aes(y = ASV277, color = "ASV277"), linewidth = 0.7) +
  geom_line(aes(y = ASV358, color = "ASV358"), linewidth = 0.7) +
  geom_line(aes(y = ASV4, color = "ASV4"), linewidth = 0.7) +
  labs(
    title = "AM_Cd5 (first 90*6)",
    x = "New Time",
    y = "number of DNA copies"
  ) +
  scale_color_manual(
    name = NULL,
    values = c(
      "ASV55" = "#BD210F", "ASV153" = "#F15D00", "ASV178" = "#DDB500",
      "ASV277" = "#20D6B5", "ASV358" = "#0092E0", "ASV4" = "#202023"
    )
  ) +
  theme(
    plot.title = element_text(
      size = 20,
      hjust = 0.5,
      margin = margin(b = 15)
    ),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    panel.border = element_rect(
      colour = "black", 
      fill = NA, 
      linewidth = 2),
    
    legend.position = c(0.995, 0.98),         
    legend.justification = c(1, 1),          
    legend.key.width = unit(1, "cm"),       
    legend.text = element_text(size = 14),
    legend.spacing.y = unit(0.2, "cm"),      
    legend.background = element_blank(),     
    legend.margin = margin(t = 0, r = 5, b = 0, l = 5)  
  ) +
  guides(linewidth = "none")

print(p3)


# we can do t+1 as well, but since the pattern is the same so we skip. 
# but we can do lagged point plot for t and t+1, for V1

plot_data <- data.frame(
  t = ds[,1],        
  t_plus_1 = ds_tp1[,1]  
)


p4 <- ggplot(plot_data, aes(x = t, y = t_plus_1)) +
  geom_point() +  
  labs(
    x = "t",      
    y = "t+1",    
    title = "Lagged Point Plot for species 1 (AM_Cd5)"  
  ) +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5,margin = margin(b = 15)),
    geom_point(aes(y = V1), size = 1.2),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    panel.border = element_rect(colour = "black",  fill = NA, linewidth = 2))


print(p4)


# Out-sample ----
if(out.sample|nout!=0){
  d.test <- do[nin:(ndo-1),]                 # Out-of-sample dataset at time t 
  dt_tp1 <- do[(nin+1):ndo,]                 # Out-of-sample dataset at time t+1
  ds.test <- (d.test-repmat(do.mean,nrow(d.test),1))*repmat(do.sd,nrow(d.test),1)^-1 # Normalized out-of-sample dataset at time t
  dst_tp1 <- (dt_tp1-repmat(do.mean,nrow(dt_tp1),1))*repmat(do.sd,nrow(dt_tp1),1)^-1 # Normalized out-of-sample dataset at time t+1
}else{d.test <- dt_tp1 <- dst_tp1 <- ds.test <- NULL}

# Compiled data at time t -> '1-87' + '88-89'
ds.all <- rbind(ds,ds.test)
dim(ds.all) # 89 rows and 141 columns. Since we need to have lagged dataset, number of time steps of ds.all is 89, although ndo = 30.

# Finding optimal embedding dimension (Ed) and nonlinearity parameter ----
#############################################################
# Find the optimal embedding dimension & nonlinearity parameter for each variable 
# based on univariate simplex projection and S-map, respectively

# Univariate simplex projection
Emax <- 10
cri <- 'rmse' # model selection 
Ed <- NULL
forecast_skill_simplex <- NULL
for(i in 1:ncol(ds)){
  spx.i <- simplex(ds[,i],E=2:Emax)
  Ed <- c(Ed,spx.i[which.min(spx.i[,cri])[1],'E'])
  forecast_skill_simplex <- c(forecast_skill_simplex,spx.i[which.min(spx.i[,cri])[1],'rho'])
}
Ed # The optimal embedding dimension for each variable
forecast_skill_simplex # Forecast skills for each variable based on simplex projection

######################################################################
# Finding causal variables by CCM ----
# Find causal variables by CCM analysis for multiview embedding
# Warning: It is time consuming for calculating the causation for each node
# CCM causality test for all node pairs 
# do.CCM <- F 
if(do.CCM){ 
  ccm.out <- ccm.fast.demo(ds, Epair=T,cri=cri,Emax=Emax)
  ccm.sig <- ccm.out[['ccm.sig']]
  ccm.rho <- ccm.out[['ccm.rho']]
  if(save){
    # To avoid overwrite the original files, we save them with different names, 'XXX_NEW'.
    write.csv(ccm.sig, file.path("out", paste("ccm_sig_", da.name, "_nin", nin, "_demo_NEW.csv", sep="")), row.names=FALSE)
    
    write.csv(ccm.rho, file.path('out', paste('ccm_rho_',da.name,'_nin',nin,'_demo_NEW.csv',sep='')),row.names=F)
  }
}

ccm.sig <- read.csv(file.path('out',paste('ccm_sig_',da.name,'_nin',nin,'_demo_NEW.csv',sep='')),header=T,stringsAsFactors = F)
ccm.rho <- read.csv(file.path('out',paste('ccm_rho_',da.name,'_nin',nin,'_demo_NEW.csv',sep='')),header=T,stringsAsFactors = F)

# ~~~~ meaning of the two tables? ----


######################################################################
# Multiview embedding ----
# Perform multiview embedding analysis for each node/species
# Warning: It is time consuming for running multview embedding for each node/species. 
# ** ~ 1 min for one node using intel CORE i7, ThinkPad P1 in balanced mode. 
# do.multiview <- F
if(do.multiview){
  esele_lag <- esim.lag.demo(ds,ccm.rho,ccm.sig,Ed,kmax=10000,kn=100,max_lag=3,Emax=Emax)
  # To avoid overwrite the original files, we save them with different names, 'XXX_NEW'.
  if(save){write.csv(esele_lag,file.path('out', paste('eseleLag_',da.name,'_nin',nin,'_demo_NEW.csv',sep='')),row.names=F)}
}

esele <- read.csv(file.path('out', paste('eseleLag_',da.name,'_nin',nin,'_demo_NEW.csv',sep='')),header=T)


# So far so good ----


####################################################
## The computation of multiview distance
dmatrix.mv <- mvdist.demo(ds,ds.all,esele)
dmatrix.train.mvx <- dmatrix.mv[['dmatrix.train.mvx']]
dmatrix.test.mvx <- dmatrix.mv[['dmatrix.test.mvx']]

save.image(file = "AM_Cd5_workspace.RData")

# Then use server
# NO.4 AM_Cd5  END -------

########################################################################
########################################################################

# NO.5 NM_Cd0 START ------
# This 'EDM_test.R' file quickly go through the whole process of running rEDM on a test dataset, before using the real data. The test dataset is from Chang et al. 2021 Ecol Lett (Ricker model). Please refer to Chang et al. on how the data set is generated. The test dataset used is in 'data/result20191024_0_0_0_.csv'.  
# For real data, the analysis is conducted in “EDM_local.R” and"EDM_NM_Cd0.R" file

# Preparation ----
# check version of the installed rEDM package
packageVersion("rEDM") # note Chang et al. 2021 Ecol Lett use v1.2.3. Other versions may be incompatible. 

seed <- 49563
set.seed(seed)

# load intact functions from Chang et al. 
source('code/Demo_MDR_function.R')  

# Load original dataset ----
da.range <- 1:45 # Subsample for data analysis
out.sample <- T # T/F for out-of-sample forecast
if(out.sample){nout <- 2}else{nout <- 0}  # number of out-of-sample

# First look at the first 6 species over 45 time steps ----
df <- read.csv('data/filtered_for_EDM_network/ASV_tab.NM.Cd0.abs.EDM.csv',header=T,stringsAsFactors=F)
p1 <- ggplot(data.frame(df), aes(x = 1:nrow(df))) +
  geom_line(aes(y = ASV55, color = "ASV55"), linewidth = 0.7) +
  geom_line(aes(y = ASV153, color = "ASV153"), linewidth = 0.7) +
  geom_line(aes(y = ASV178, color = "ASV178"), linewidth = 0.7) +
  geom_line(aes(y = ASV277, color = "ASV277"), linewidth = 0.7) +
  geom_line(aes(y = ASV2766, color = "ASV2766"), linewidth = 0.7) +
  geom_line(aes(y = ASV4, color = "ASV4"), linewidth = 0.7) +
  labs(
    title = "NM_Cd0 (first 45*6)",
    x = "New Time",
    y = "number of DNA copies"
  ) +
  scale_color_manual(
    name = NULL,
    values = c(
      "ASV55" = "#BD210F", "ASV153" = "#F15D00", "ASV178" = "#DDB500",
      "ASV277" = "#20D6B5", "ASV2766" = "#0092E0", "ASV4" = "#202023"
    )
  ) +
  theme(
    plot.title = element_text(
      size = 20,
      hjust = 0.5,
      margin = margin(b = 15)
    ),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    panel.border = element_rect(
      colour = "black", 
      fill = NA, 
      linewidth = 2),
    
    legend.position = c(0.995, 0.98),          
    legend.justification = c(1, 1),         
    legend.key.width = unit(1, "cm"),      
    legend.text = element_text(size = 14),
    legend.spacing.y = unit(0.2, "cm"),    
    legend.background = element_blank(),      
    legend.margin = margin(t = 0, r = 5, b = 0, l = 5)  
  ) +
  guides(linewidth = "none")
print(p1)
#202023


(da.name <- 'ASV_tab.NM.Cd0.abs.EDM')
do <- read.csv('data/filtered_for_EDM_network/ASV_tab.NM.Cd0.abs.EDM.csv',header=T,stringsAsFactors=F)
dot <- do[da.range,1] # get time column
do <- do[da.range,-1] # remove the first column (time) from the data frame
ndo <- nrow(do) # number of time steps
nin <- ndo-nout # library sample size

# take the 124 columns (total number of ASV),filted do data frame have 130 columns
do[, 2] #see column 2 "New time" 
ncol(do) #see total columns number
do_subset <- do[, 3:126]#only take all 'ASV' columns
do <- do_subset

# ~~~~ plot logistic map ----
# plot the normalized in-sample dataset at time t. Include V1 to V6.
p2 <- ggplot(data.frame(do), aes(x = 1:nrow(do))) +
  geom_line(aes(y = ASV55, color = "ASV55"), linewidth = 0.7) +
  geom_line(aes(y = ASV153, color = "ASV153"), linewidth = 0.7) +
  geom_line(aes(y = ASV178, color = "ASV178"), linewidth = 0.7) +
  geom_line(aes(y = ASV277, color = "ASV277"), linewidth = 0.7) +
  geom_line(aes(y = ASV2766, color = "ASV2766"), linewidth = 0.7) +
  geom_line(aes(y = ASV4, color = "ASV4"), linewidth = 0.7) +
  labs(
    title = "NM_Cd0 (first 45*6)",
    x = "New Time",
    y = "number of DNA copies"
  ) +
  scale_color_manual(
    name = NULL,
    values = c(
      "ASV55" = "#BD210F", "ASV153" = "#F15D00", "ASV178" = "#DDB500",
      "ASV277" = "#20D6B5", "ASV2766" = "#0092E0", "ASV4" = "#202023"
    )
  ) +
  theme(
    plot.title = element_text(
      size = 20,
      hjust = 0.5,
      margin = margin(b = 15)
    ),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    panel.border = element_rect(
      colour = "black", 
      fill = NA, 
      linewidth = 2),
    
    legend.position = c(0.995, 0.98),         
    legend.justification = c(1, 1),           
    legend.key.width = unit(1, "cm"),      
    legend.text = element_text(size = 14),
    legend.spacing.y = unit(0.2, "cm"),       
    legend.background = element_blank(),     
    legend.margin = margin(t = 0, r = 5, b = 0, l = 5)  
  ) +
  guides(linewidth = "none")

print(p2)



# For real data, we may have 45 time steps and 100+ 'species' for each treatment with the same matrix structure. The 'species' can be ASV/OTU ID, taxon (at phylum, order, class, family, genus, or species level), depending on nature of data and research objectives. Also, we may focus on overlapping 'species' across treatments, if we want to compare the results across treatments.

dto <- apply(do,1,sum) # sum of each row (time step) in the data frame
dpo <- do*repmat(dto^-1,1,ncol(do)) # normalize the rows of the matrix 'do' so that each row sums to 1(?). This is often done in data processing to scale the data and make comparisons across different rows more meaningful. Similar to using rarefied 16S rRNA gene sequencing data(?).
# repmat() is the own function by sourcing 'Demo_MDR_function.R'. 
# repmat() makes a matrix with repeated column or row.

# check and shall see all values in the rowsum column equal to 1. 
as.data.frame(dpo) %>% mutate(rowsum = rowSums(across(everything())))


# Exclusion of rare species ----
# Here we skip this step otherwise all species will be omitted.
# We shall consider the followings for real dataset.
# *Threshold is upon decision*.
# pcri <- 0;bcri <- 10^-3; # criteria for selecting species based on proportion of presnece (pcri) and mean relative abundance (bri) 
# doind2 <- (apply(dpo,2,mean,na.rm=T)>(bcri))&((apply(do>0,2,sum,na.rm=T)/nrow(do))>pcri) # index for selected species 
# exsp2 <- setdiff(1:ncol(do),which(doind2))   # index for rare species 
# do <- do[,-exsp2]                            # Dataset excluded rare species


(nsp <- ncol(do)) # number of species
ndo <- nrow(do) # number of time steps
nin <- ndo-nout # number of in-sample time steps

# Mean and SD of each node/species ----
do.mean <- apply(do[1:nin,],2,mean,na.rm=T) # mean of each column (species abundance) in the first 28 rows (in-sample)
do.sd <- apply(do[1:nin,],2,sd,na.rm=T) # sd of each column (species abundance) in the first 28 rows (in-sample)

# Construct a sd(i,j) matrix ----
# (dimension = nsp*nsp, i.e., 6*6 here).
# if sd(i,j)>1, meaning j varies more than i
# if sd(i,j)<1, meaning i varies more than j
# if sd(i,j)=1, meaning i and j vary the same
dosdM <- repmat(c(do.sd)^-1,1,nsp)*repmat(c(do.sd),nsp,1) 

# check the sd(i,j)
dosdM; dim(dosdM)

# In-sample ----
d <- do[1:(nin-1),] # In-sample dataset at time t (time 1-42)
d_tp1 <- do[2:(nin),] # In-sample dataset at time t+1 (time 2-43)

# ~~ normalization (z-score) ----
ds <- (d-repmat(do.mean,nrow(d),1))*repmat(do.sd,nrow(d),1)^-1 # Normalized in-sample dataset at time t (i.e., *z-score = (x-mean)/sd)
ds_tp1 <- (d_tp1-repmat(do.mean,nrow(d_tp1),1))*repmat(do.sd,nrow(d_tp1),1)^-1 # Normalized in-sample dataset at time t+1 (i.e., z-score)

# ~~~~ plot logistic map ----
# plot the normalized in-sample dataset at time t. Include V1 to V6.
p3 <- ggplot(data.frame(ds), aes(x = 1:nrow(ds))) +
  geom_line(aes(y = ASV55, color = "ASV55"), linewidth = 0.7) +
  geom_line(aes(y = ASV153, color = "ASV153"), linewidth = 0.7) +
  geom_line(aes(y = ASV178, color = "ASV178"), linewidth = 0.7) +
  geom_line(aes(y = ASV277, color = "ASV277"), linewidth = 0.7) +
  geom_line(aes(y = ASV2766, color = "ASV2766"), linewidth = 0.7) +
  geom_line(aes(y = ASV4, color = "ASV4"), linewidth = 0.7) +
  labs(
    title = "In-sample (T1~87*6) NM_Cd0",
    x = "Time t",
    y = "Normalized"
  ) +
  scale_color_manual(
    name = NULL,
    values = c(
      "ASV55" = "#BD210F", "ASV153" = "#F15D00", "ASV178" = "#DDB500",
      "ASV277" = "#20D6B5", "ASV2766" = "#0092E0", "ASV4" = "#202023"
    )
  ) +
  theme(
    plot.title = element_text(
      size = 20,
      hjust = 0.5,
      margin = margin(b = 15)
    ),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    panel.border = element_rect(
      colour = "black", 
      fill = NA, 
      linewidth = 2),
    
    legend.position = c(0.995, 0.98),          
    legend.justification = c(1, 1),        
    legend.key.width = unit(1, "cm"),       
    legend.text = element_text(size = 14),
    legend.spacing.y = unit(0.2, "cm"),      
    legend.background = element_blank(),     
    legend.margin = margin(t = 0, r = 5, b = 0, l = 5)  
  ) +
  guides(linewidth = "none")

print(p3)

# we can do t+1 as well, but since the pattern is the same so we skip. 
# but we can do lagged point plot for t and t+1, for V1

plot_data <- data.frame(
  t = ds[,1],     
  t_plus_1 = ds_tp1[,1]  
)


p4 <- ggplot(plot_data, aes(x = t, y = t_plus_1)) +
  geom_point() + 
  labs(
    x = "t",     
    y = "t+1",   
    title = "Lagged Point Plot for species 1 (NM_Cd0)"  
  ) +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5,margin = margin(b = 15)),
    geom_point(aes(y = V1), size = 1.2),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    panel.border = element_rect(colour = "black",  fill = NA, linewidth = 2))


print(p4)

# Out-sample ----
if(out.sample|nout!=0){
  d.test <- do[nin:(ndo-1),]                 # Out-of-sample dataset at time t 
  dt_tp1 <- do[(nin+1):ndo,]                 # Out-of-sample dataset at time t+1
  ds.test <- (d.test-repmat(do.mean,nrow(d.test),1))*repmat(do.sd,nrow(d.test),1)^-1 # Normalized out-of-sample dataset at time t
  dst_tp1 <- (dt_tp1-repmat(do.mean,nrow(dt_tp1),1))*repmat(do.sd,nrow(dt_tp1),1)^-1 # Normalized out-of-sample dataset at time t+1
}else{d.test <- dt_tp1 <- dst_tp1 <- ds.test <- NULL}

# Compiled data at time t -> '1-42' + '43-44'
ds.all <- rbind(ds,ds.test)
dim(ds.all) # 44 rows and 124 columns. Since we need to have lagged dataset, number of time steps of ds.all is 44, although ndo = 30.

# Finding optimal embedding dimension (Ed) and nonlinearity parameter ----
#############################################################
# Find the optimal embedding dimension & nonlinearity parameter for each variable 
# based on univariate simplex projection and S-map, respectively

# Univariate simplex projection
Emax <- 10
cri <- 'rmse' # model selection 
Ed <- NULL
forecast_skill_simplex <- NULL
for(i in 1:ncol(ds)){
  spx.i <- simplex(ds[,i],E=2:Emax)
  Ed <- c(Ed,spx.i[which.min(spx.i[,cri])[1],'E'])
  forecast_skill_simplex <- c(forecast_skill_simplex,spx.i[which.min(spx.i[,cri])[1],'rho'])
}
Ed # The optimal embedding dimension for each variable
forecast_skill_simplex # Forecast skills for each variable based on simplex projection

######################################################################
# Finding causal variables by CCM ----
# Find causal variables by CCM analysis for multiview embedding
# Warning: It is time consuming for calculating the causation for each node
# CCM causality test for all node pairs 
# do.CCM <- F 
if(do.CCM){ 
  ccm.out <- ccm.fast.demo(ds, Epair=T,cri=cri,Emax=Emax)
  ccm.sig <- ccm.out[['ccm.sig']]
  ccm.rho <- ccm.out[['ccm.rho']]
  if(save){
    # To avoid overwrite the original files, we save them with different names, 'XXX_NEW'.
    write.csv(ccm.sig, file.path("out", paste("ccm_sig_", da.name, "_nin", nin, "_demo_NEW.csv", sep="")), row.names=FALSE)
    
    write.csv(ccm.rho, file.path('out', paste('ccm_rho_',da.name,'_nin',nin,'_demo_NEW.csv',sep='')),row.names=F)
  }
}

ccm.sig <- read.csv(file.path('out',paste('ccm_sig_',da.name,'_nin',nin,'_demo_NEW.csv',sep='')),header=T,stringsAsFactors = F)
ccm.rho <- read.csv(file.path('out',paste('ccm_rho_',da.name,'_nin',nin,'_demo_NEW.csv',sep='')),header=T,stringsAsFactors = F)

# ~~~~ meaning of the two tables? ----


######################################################################
# Multiview embedding ----
# Perform multiview embedding analysis for each node/species
# Warning: It is time consuming for running multview embedding for each node/species. 
# ** ~ 1 min for one node using intel CORE i7, ThinkPad P1 in balanced mode. 
# do.multiview <- F
if(do.multiview){
  esele_lag <- esim.lag.demo(ds,ccm.rho,ccm.sig,Ed,kmax=10000,kn=100,max_lag=3,Emax=Emax)
  # To avoid overwrite the original files, we save them with different names, 'XXX_NEW'.
  if(save){write.csv(esele_lag,file.path('out', paste('eseleLag_',da.name,'_nin',nin,'_demo_NEW.csv',sep='')),row.names=F)}
}

esele <- read.csv(file.path('out', paste('eseleLag_',da.name,'_nin',nin,'_demo_NEW.csv',sep='')),header=T)


# So far so good ----


####################################################
## The computation of multiview distance
dmatrix.mv <- mvdist.demo(ds,ds.all,esele)
dmatrix.train.mvx <- dmatrix.mv[['dmatrix.train.mvx']]
dmatrix.test.mvx <- dmatrix.mv[['dmatrix.test.mvx']]

save.image(file = "NM_Cd0_workspace.RData")

# Then use server
# NO.5 NM_Cd0  END -------

######################################################
######################################################

# NO.6 NM_Cd2 START ------
# This 'EDM_test.R' file quickly go through the whole process of running rEDM on a test dataset, before using the real data. The test dataset is from Chang et al. 2021 Ecol Lett (Ricker model). Please refer to Chang et al. on how the data set is generated. The test dataset used is in 'data/result20191024_0_0_0_.csv'.  
# For real data, the analysis is conducted in “EDM_local.R” and "EDM_NM_Cd2.R" file

# Preparation ----
# check version of the installed rEDM package
packageVersion("rEDM") # note Chang et al. 2021 Ecol Lett use v1.2.3. Other versions may be incompatible. 

seed <- 49563
set.seed(seed)

# load intact functions from Chang et al. 
source('code/Demo_MDR_function.R')  

# Load original dataset ----
da.range <- 1:45 # Subsample for data analysis
out.sample <- T # T/F for out-of-sample forecast
if(out.sample){nout <- 2}else{nout <- 0}  # number of out-of-sample

# First look at the first 6 species over 45 time steps ----
df <- read.csv('data/filtered_for_EDM_network/ASV_tab.NM.Cd2.abs.EDM.csv',header=T,stringsAsFactors=F)
p1 <- ggplot(data.frame(df), aes(x = 1:nrow(df))) +
  geom_line(aes(y = ASV55, color = "ASV55"), linewidth = 0.7) +
  geom_line(aes(y = ASV153, color = "ASV153"), linewidth = 0.7) +
  geom_line(aes(y = ASV178, color = "ASV178"), linewidth = 0.7) +
  geom_line(aes(y = ASV277, color = "ASV277"), linewidth = 0.7) +
  geom_line(aes(y = ASV4, color = "ASV4"), linewidth = 0.7) +
  geom_line(aes(y = ASV4, color = "ASV4"), linewidth = 0.7) +
  labs(
    title = "NM_Cd2 (first 45*6)",
    x = "New Time",
    y = "number of DNA copies"
  ) +
  scale_color_manual(
    name = NULL,
    values = c(
      "ASV55" = "#BD210F", "ASV153" = "#F15D00", "ASV178" = "#DDB500",
      "ASV277" = "#20D6B5", "ASV4" = "#0092E0", "ASV4" = "#202023"
    )
  ) +
  theme(
    plot.title = element_text(
      size = 20,
      hjust = 0.5,
      margin = margin(b = 15)
    ),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    panel.border = element_rect(
      colour = "black", 
      fill = NA, 
      linewidth = 2),
    
    legend.position = c(0.995, 0.98),         
    legend.justification = c(1, 1),           
    legend.key.width = unit(1, "cm"),       
    legend.text = element_text(size = 14),
    legend.spacing.y = unit(0.2, "cm"),       
    legend.background = element_blank(),      
    legend.margin = margin(t = 0, r = 5, b = 0, l = 5)  
  ) +
  guides(linewidth = "none")
print(p1)


(da.name <- 'ASV_tab.NM.Cd2.abs.EDM')
do <- read.csv('data/filtered_for_EDM_network/ASV_tab.NM.Cd2.abs.EDM.csv',header=T,stringsAsFactors=F)
dot <- do[da.range,1] # get time column
do <- do[da.range,-1] # remove the first column (time) from the data frame
ndo <- nrow(do) # number of time steps
nin <- ndo-nout # library sample size

# take the 133 columns (total number of ASV),filted do data frame have 130 columns
do[, 2] #see column 2 "New time" 
ncol(do) #see total columns number
do_subset <- do[, 3:135]#only take all 'ASV' columns
do <- do_subset

# ~~~~ plot logistic map ----
# plot the normalized in-sample dataset at time t. Include V1 to V6.
p2 <- ggplot(data.frame(do), aes(x = 1:nrow(do))) +
  geom_line(aes(y = ASV55, color = "ASV55"), linewidth = 0.7) +
  geom_line(aes(y = ASV153, color = "ASV153"), linewidth = 0.7) +
  geom_line(aes(y = ASV178, color = "ASV178"), linewidth = 0.7) +
  geom_line(aes(y = ASV277, color = "ASV277"), linewidth = 0.7) +
  geom_line(aes(y = ASV4, color = "ASV4"), linewidth = 0.7) +
  geom_line(aes(y = ASV4, color = "ASV4"), linewidth = 0.7) +
  labs(
    title = "NM_Cd2 (first 45*6)",
    x = "New Time",
    y = "number of DNA copies"
  ) +
  scale_color_manual(
    name = NULL,
    values = c(
      "ASV55" = "#BD210F", "ASV153" = "#F15D00", "ASV178" = "#DDB500",
      "ASV277" = "#20D6B5", "ASV4" = "#0092E0", "ASV4" = "#202023"
    )
  ) +
  theme(
    plot.title = element_text(
      size = 20,
      hjust = 0.5,
      margin = margin(b = 15)
    ),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    panel.border = element_rect(
      colour = "black", 
      fill = NA, 
      linewidth = 2),
    
    legend.position = c(0.995, 0.98),          
    legend.justification = c(1, 1),          
    legend.key.width = unit(1, "cm"),      
    legend.text = element_text(size = 14),
    legend.spacing.y = unit(0.2, "cm"),       
    legend.background = element_blank(),     
    legend.margin = margin(t = 0, r = 5, b = 0, l = 5)  
  ) +
  guides(linewidth = "none")

print(p2)



# For real data, we may have 45 time steps and 100+ 'species' for each treatment with the same matrix structure. The 'species' can be ASV/OTU ID, taxon (at phylum, order, class, family, genus, or species level), depending on nature of data and research objectives. Also, we may focus on overlapping 'species' across treatments, if we want to compare the results across treatments.

dto <- apply(do,1,sum) # sum of each row (time step) in the data frame
dpo <- do*repmat(dto^-1,1,ncol(do)) # normalize the rows of the matrix 'do' so that each row sums to 1(?). This is often done in data processing to scale the data and make comparisons across different rows more meaningful. Similar to using rarefied 16S rRNA gene sequencing data(?).
# repmat() is the own function by sourcing 'Demo_MDR_function.R'. 
# repmat() makes a matrix with repeated column or row.

# check and shall see all values in the rowsum column equal to 1. 
as.data.frame(dpo) %>% mutate(rowsum = rowSums(across(everything())))


# Exclusion of rare species ----
# Here we skip this step otherwise all species will be omitted.
# We shall consider the followings for real dataset.
# *Threshold is upon decision*.
# pcri <- 0;bcri <- 10^-3; # criteria for selecting species based on proportion of presnece (pcri) and mean relative abundance (bri) 
# doind2 <- (apply(dpo,2,mean,na.rm=T)>(bcri))&((apply(do>0,2,sum,na.rm=T)/nrow(do))>pcri) # index for selected species 
# exsp2 <- setdiff(1:ncol(do),which(doind2))   # index for rare species 
# do <- do[,-exsp2]                            # Dataset excluded rare species


(nsp <- ncol(do)) # number of species
ndo <- nrow(do) # number of time steps
nin <- ndo-nout # number of in-sample time steps

# Mean and SD of each node/species ----
do.mean <- apply(do[1:nin,],2,mean,na.rm=T) # mean of each column (species abundance) in the first 28 rows (in-sample)
do.sd <- apply(do[1:nin,],2,sd,na.rm=T) # sd of each column (species abundance) in the first 28 rows (in-sample)

# Construct a sd(i,j) matrix ----
# (dimension = nsp*nsp, i.e., 6*6 here).
# if sd(i,j)>1, meaning j varies more than i
# if sd(i,j)<1, meaning i varies more than j
# if sd(i,j)=1, meaning i and j vary the same
dosdM <- repmat(c(do.sd)^-1,1,nsp)*repmat(c(do.sd),nsp,1) 

# check the sd(i,j)
dosdM; dim(dosdM)

# In-sample ----
d <- do[1:(nin-1),] # In-sample dataset at time t (time 1-42)
d_tp1 <- do[2:(nin),] # In-sample dataset at time t+1 (time 2-43)

# ~~ normalization (z-score) ----
ds <- (d-repmat(do.mean,nrow(d),1))*repmat(do.sd,nrow(d),1)^-1 # Normalized in-sample dataset at time t (i.e., *z-score = (x-mean)/sd)
ds_tp1 <- (d_tp1-repmat(do.mean,nrow(d_tp1),1))*repmat(do.sd,nrow(d_tp1),1)^-1 # Normalized in-sample dataset at time t+1 (i.e., z-score)

# ~~~~ plot logistic map ----
# plot the normalized in-sample dataset at time t. Include V1 to V6.
p3 <- ggplot(data.frame(ds), aes(x = 1:nrow(ds))) +
  geom_line(aes(y = ASV55, color = "ASV55"), linewidth = 0.7) +
  geom_line(aes(y = ASV153, color = "ASV153"), linewidth = 0.7) +
  geom_line(aes(y = ASV178, color = "ASV178"), linewidth = 0.7) +
  geom_line(aes(y = ASV277, color = "ASV277"), linewidth = 0.7) +
  geom_line(aes(y = ASV4, color = "ASV4"), linewidth = 0.7) +
  geom_line(aes(y = ASV4, color = "ASV4"), linewidth = 0.7) +
  labs(
    title = "In-sample (T1~87*6) NM_Cd2",
    x = "Time t",
    y = "Normalized"
  ) +
  scale_color_manual(
    name = NULL,
    values = c(
      "ASV55" = "#BD210F", "ASV153" = "#F15D00", "ASV178" = "#DDB500",
      "ASV277" = "#20D6B5", "ASV4" = "#0092E0", "ASV4" = "#202023"
    )
  ) +
  theme(
    plot.title = element_text(
      size = 20,
      hjust = 0.5,
      margin = margin(b = 15)
    ),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    panel.border = element_rect(
      colour = "black", 
      fill = NA, 
      linewidth = 2),
    
    legend.position = c(0.995, 0.98),         
    legend.justification = c(1, 1),          
    legend.key.width = unit(1, "cm"),       
    legend.text = element_text(size = 14),
    legend.spacing.y = unit(0.2, "cm"),       
    legend.background = element_blank(),      
    legend.margin = margin(t = 0, r = 5, b = 0, l = 5)  
  ) +
  guides(linewidth = "none")

print(p3)

# we can do t+1 as well, but since the pattern is the same so we skip. 
# but we can do lagged point plot for t and t+1, for V1

plot_data <- data.frame(
  t = ds[,1],       
  t_plus_1 = ds_tp1[,1] 
)

# ggplot2绘图
p4 <- ggplot(plot_data, aes(x = t, y = t_plus_1)) +
  geom_point() + 
  labs(
    x = "t",     
    y = "t+1",   
    title = "Lagged Point Plot for species 1 (NM_Cd2)"  
  ) +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5,margin = margin(b = 15)),
    geom_point(aes(y = V1), size = 1.2),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    panel.border = element_rect(colour = "black",  fill = NA, linewidth = 2))


print(p4)

# Out-sample ----
if(out.sample|nout!=0){
  d.test <- do[nin:(ndo-1),]                 # Out-of-sample dataset at time t 
  dt_tp1 <- do[(nin+1):ndo,]                 # Out-of-sample dataset at time t+1
  ds.test <- (d.test-repmat(do.mean,nrow(d.test),1))*repmat(do.sd,nrow(d.test),1)^-1 # Normalized out-of-sample dataset at time t
  dst_tp1 <- (dt_tp1-repmat(do.mean,nrow(dt_tp1),1))*repmat(do.sd,nrow(dt_tp1),1)^-1 # Normalized out-of-sample dataset at time t+1
}else{d.test <- dt_tp1 <- dst_tp1 <- ds.test <- NULL}

# Compiled data at time t -> '1-42' + '43-44'
ds.all <- rbind(ds,ds.test)
dim(ds.all) # 44 rows and 133 columns. Since we need to have lagged dataset, number of time steps of ds.all is 44, although ndo = 30.

# Finding optimal embedding dimension (Ed) and nonlinearity parameter ----
#############################################################
# Find the optimal embedding dimension & nonlinearity parameter for each variable 
# based on univariate simplex projection and S-map, respectively

# Univariate simplex projection
Emax <- 10
cri <- 'rmse' # model selection 
Ed <- NULL
forecast_skill_simplex <- NULL
for(i in 1:ncol(ds)){
  spx.i <- simplex(ds[,i],E=2:Emax)
  Ed <- c(Ed,spx.i[which.min(spx.i[,cri])[1],'E'])
  forecast_skill_simplex <- c(forecast_skill_simplex,spx.i[which.min(spx.i[,cri])[1],'rho'])
}
Ed # The optimal embedding dimension for each variable
forecast_skill_simplex # Forecast skills for each variable based on simplex projection

######################################################################
# Finding causal variables by CCM ----
# Find causal variables by CCM analysis for multiview embedding
# Warning: It is time consuming for calculating the causation for each node
# CCM causality test for all node pairs 
# do.CCM <- F 
if(do.CCM){ 
  ccm.out <- ccm.fast.demo(ds, Epair=T,cri=cri,Emax=Emax)
  ccm.sig <- ccm.out[['ccm.sig']]
  ccm.rho <- ccm.out[['ccm.rho']]
  if(save){
    # To avoid overwrite the original files, we save them with different names, 'XXX_NEW'.
    write.csv(ccm.sig, file.path("out", paste("ccm_sig_", da.name, "_nin", nin, "_demo_NEW.csv", sep="")), row.names=FALSE)
    
    write.csv(ccm.rho, file.path('out', paste('ccm_rho_',da.name,'_nin',nin,'_demo_NEW.csv',sep='')),row.names=F)
  }
}

ccm.sig <- read.csv(file.path('out',paste('ccm_sig_',da.name,'_nin',nin,'_demo_NEW.csv',sep='')),header=T,stringsAsFactors = F)
ccm.rho <- read.csv(file.path('out',paste('ccm_rho_',da.name,'_nin',nin,'_demo_NEW.csv',sep='')),header=T,stringsAsFactors = F)

# ~~~~ meaning of the two tables? ----


######################################################################
# Multiview embedding ----
# Perform multiview embedding analysis for each node/species
# Warning: It is time consuming for running multview embedding for each node/species. 
# ** ~ 1 min for one node using intel CORE i7, ThinkPad P1 in balanced mode. 
# do.multiview <- F
if(do.multiview){
  esele_lag <- esim.lag.demo(ds,ccm.rho,ccm.sig,Ed,kmax=10000,kn=100,max_lag=3,Emax=Emax)
  # To avoid overwrite the original files, we save them with different names, 'XXX_NEW'.
  if(save){write.csv(esele_lag,file.path('out', paste('eseleLag_',da.name,'_nin',nin,'_demo_NEW.csv',sep='')),row.names=F)}
}

esele <- read.csv(file.path('out', paste('eseleLag_',da.name,'_nin',nin,'_demo_NEW.csv',sep='')),header=T)


# So far so good ----


####################################################
## The computation of multiview distance
dmatrix.mv <- mvdist.demo(ds,ds.all,esele)
dmatrix.train.mvx <- dmatrix.mv[['dmatrix.train.mvx']]
dmatrix.test.mvx <- dmatrix.mv[['dmatrix.test.mvx']]

save.image(file = "NM_Cd2_workspace.RData")

# Then use server
# NO.6 NM_Cd2  END -------

######################################################
######################################################

# NO.7 NM_Cd5 START ------
# This 'EDM_test.R' file quickly go through the whole process of running rEDM on a test dataset, before using the real data. The test dataset is from Chang et al. 2021 Ecol Lett (Ricker model). Please refer to Chang et al. on how the data set is generated. The test dataset used is in 'data/result20191024_0_0_0_.csv'.  
# For real data, the analysis is conducted in “EDM_local.R” and "EDM_NM_Cd5.R" file

# Preparation ----
# check version of the installed rEDM package
packageVersion("rEDM") # note Chang et al. 2021 Ecol Lett use v1.2.3. Other versions may be incompatible. 

seed <- 49563
set.seed(seed)

# load intact functions from Chang et al. 
source('code/Demo_MDR_function.R')  

# Load original dataset ----
da.range <- 1:45 # Subsample for data analysis
out.sample <- T # T/F for out-of-sample forecast
if(out.sample){nout <- 2}else{nout <- 0}  # number of out-of-sample

# First look at the first 6 species over 45 time steps ----
df <- read.csv('data/filtered_for_EDM_network/ASV_tab.NM.Cd5.abs.EDM.csv',header=T,stringsAsFactors=F)
p1 <- ggplot(data.frame(df), aes(x = 1:nrow(df))) +
  geom_line(aes(y = ASV55, color = "ASV55"), linewidth = 0.7) +
  geom_line(aes(y = ASV153, color = "ASV153"), linewidth = 0.7) +
  geom_line(aes(y = ASV178, color = "ASV178"), linewidth = 0.7) +
  geom_line(aes(y = ASV277, color = "ASV277"), linewidth = 0.7) +
  geom_line(aes(y = ASV4, color = "ASV4"), linewidth = 0.7) +
  geom_line(aes(y = ASV4, color = "ASV4"), linewidth = 0.7) +
  labs(
    title = "NM_Cd5 (first 45*6)",
    x = "New Time",
    y = "number of DNA copies"
  ) +
  scale_color_manual(
    name = NULL,
    values = c(
      "ASV55" = "#BD210F", "ASV153" = "#F15D00", "ASV178" = "#DDB500",
      "ASV277" = "#20D6B5", "ASV4" = "#0092E0", "ASV4" = "#202023"
    )
  ) +
  theme(
    plot.title = element_text(
      size = 20,
      hjust = 0.5,
      margin = margin(b = 15)
    ),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    panel.border = element_rect(
      colour = "black", 
      fill = NA, 
      linewidth = 2),
    
    legend.position = c(0.995, 0.98),         
    legend.justification = c(1, 1),          
    legend.key.width = unit(1, "cm"),       
    legend.text = element_text(size = 14),
    legend.spacing.y = unit(0.2, "cm"),      
    legend.background = element_blank(),   
    legend.margin = margin(t = 0, r = 5, b = 0, l = 5)  
  ) +
  guides(linewidth = "none")
print(p1)



(da.name <- 'ASV_tab.NM.Cd5.abs.EDM')
do <- read.csv('data/filtered_for_EDM_network/ASV_tab.NM.Cd5.abs.EDM.csv',header=T,stringsAsFactors=F)
dot <- do[da.range,1] # get time column
do <- do[da.range,-1] # remove the first column (time) from the data frame
ndo <- nrow(do) # number of time steps
nin <- ndo-nout # library sample size

# take the 124 columns (total number of ASV),filted do data frame have 124 columns
do[, 2] #see column 2 "New time" 
ncol(do) #see total columns number
do_subset <- do[, 3:126]#only take all 'ASV' columns
do <- do_subset

# ~~~~ plot logistic map ----
# plot the normalized in-sample dataset at time t. Include V1 to V6.
p2 <- ggplot(data.frame(do), aes(x = 1:nrow(do))) +
  geom_line(aes(y = ASV55, color = "ASV55"), linewidth = 0.7) +
  geom_line(aes(y = ASV153, color = "ASV153"), linewidth = 0.7) +
  geom_line(aes(y = ASV178, color = "ASV178"), linewidth = 0.7) +
  geom_line(aes(y = ASV277, color = "ASV277"), linewidth = 0.7) +
  geom_line(aes(y = ASV4, color = "ASV4"), linewidth = 0.7) +
  geom_line(aes(y = ASV4, color = "ASV4"), linewidth = 0.7) +
  labs(
    title = "NM_Cd5 (first 45*6)",
    x = "New Time",
    y = "number of DNA copies"
  ) +
  scale_color_manual(
    name = NULL,
    values = c(
      "ASV55" = "#BD210F", "ASV153" = "#F15D00", "ASV178" = "#DDB500",
      "ASV277" = "#20D6B5", "ASV4" = "#0092E0", "ASV4" = "#202023"
    )
  ) +
  theme(
    plot.title = element_text(
      size = 20,
      hjust = 0.5,
      margin = margin(b = 15)
    ),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    panel.border = element_rect(
      colour = "black", 
      fill = NA, 
      linewidth = 2),
    
    legend.position = c(0.995, 0.98),          
    legend.justification = c(1, 1),           
    legend.key.width = unit(1, "cm"),      
    legend.text = element_text(size = 14),
    legend.spacing.y = unit(0.2, "cm"),       
    legend.background = element_blank(),     
    legend.margin = margin(t = 0, r = 5, b = 0, l = 5)  
  ) +
  guides(linewidth = "none")

print(p2)



# For real data, we may have 45 time steps and 100+ 'species' for each treatment with the same matrix structure. The 'species' can be ASV/OTU ID, taxon (at phylum, order, class, family, genus, or species level), depending on nature of data and research objectives. Also, we may focus on overlapping 'species' across treatments, if we want to compare the results across treatments.

dto <- apply(do,1,sum) # sum of each row (time step) in the data frame
dpo <- do*repmat(dto^-1,1,ncol(do)) # normalize the rows of the matrix 'do' so that each row sums to 1(?). This is often done in data processing to scale the data and make comparisons across different rows more meaningful. Similar to using rarefied 16S rRNA gene sequencing data(?).
# repmat() is the own function by sourcing 'Demo_MDR_function.R'. 
# repmat() makes a matrix with repeated column or row.

# check and shall see all values in the rowsum column equal to 1. 
as.data.frame(dpo) %>% mutate(rowsum = rowSums(across(everything())))


# Exclusion of rare species ----
# Here we skip this step otherwise all species will be omitted.
# We shall consider the followings for real dataset.
# *Threshold is upon decision*.
# pcri <- 0;bcri <- 10^-3; # criteria for selecting species based on proportion of presnece (pcri) and mean relative abundance (bri) 
# doind2 <- (apply(dpo,2,mean,na.rm=T)>(bcri))&((apply(do>0,2,sum,na.rm=T)/nrow(do))>pcri) # index for selected species 
# exsp2 <- setdiff(1:ncol(do),which(doind2))   # index for rare species 
# do <- do[,-exsp2]                            # Dataset excluded rare species


(nsp <- ncol(do)) # number of species
ndo <- nrow(do) # number of time steps
nin <- ndo-nout # number of in-sample time steps

# Mean and SD of each node/species ----
do.mean <- apply(do[1:nin,],2,mean,na.rm=T) # mean of each column (species abundance) in the first 28 rows (in-sample)
do.sd <- apply(do[1:nin,],2,sd,na.rm=T) # sd of each column (species abundance) in the first 28 rows (in-sample)

# Construct a sd(i,j) matrix ----
# (dimension = nsp*nsp, i.e., 6*6 here).
# if sd(i,j)>1, meaning j varies more than i
# if sd(i,j)<1, meaning i varies more than j
# if sd(i,j)=1, meaning i and j vary the same
dosdM <- repmat(c(do.sd)^-1,1,nsp)*repmat(c(do.sd),nsp,1) 

# check the sd(i,j)
dosdM; dim(dosdM)

# In-sample ----
d <- do[1:(nin-1),] # In-sample dataset at time t (time 1-42)
d_tp1 <- do[2:(nin),] # In-sample dataset at time t+1 (time 2-43)

# ~~ normalization (z-score) ----
ds <- (d-repmat(do.mean,nrow(d),1))*repmat(do.sd,nrow(d),1)^-1 # Normalized in-sample dataset at time t (i.e., *z-score = (x-mean)/sd)
ds_tp1 <- (d_tp1-repmat(do.mean,nrow(d_tp1),1))*repmat(do.sd,nrow(d_tp1),1)^-1 # Normalized in-sample dataset at time t+1 (i.e., z-score)

# ~~~~ plot logistic map ----
# plot the normalized in-sample dataset at time t. Include V1 to V6.
p3 <- ggplot(data.frame(ds), aes(x = 1:nrow(ds))) +
  geom_line(aes(y = ASV55, color = "ASV55"), linewidth = 0.7) +
  geom_line(aes(y = ASV153, color = "ASV153"), linewidth = 0.7) +
  geom_line(aes(y = ASV178, color = "ASV178"), linewidth = 0.7) +
  geom_line(aes(y = ASV277, color = "ASV277"), linewidth = 0.7) +
  geom_line(aes(y = ASV4, color = "ASV4"), linewidth = 0.7) +
  geom_line(aes(y = ASV4, color = "ASV4"), linewidth = 0.7) +
  labs(
    title = "In-sample (T1~87*6) NM_Cd5",
    x = "Time t",
    y = "Normalized"
  ) +
  scale_color_manual(
    name = NULL,
    values = c(
      "ASV55" = "#BD210F", "ASV153" = "#F15D00", "ASV178" = "#DDB500",
      "ASV277" = "#20D6B5", "ASV4" = "#0092E0", "ASV4" = "#202023"
    )
  ) +
  theme(
    plot.title = element_text(
      size = 20,
      hjust = 0.5,
      margin = margin(b = 15)
    ),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    panel.border = element_rect(
      colour = "black", 
      fill = NA, 
      linewidth = 2),
    
    legend.position = c(0.995, 0.98),          
    legend.justification = c(1, 1),         
    legend.key.width = unit(1, "cm"),       
    legend.text = element_text(size = 14),
    legend.spacing.y = unit(0.2, "cm"),     
    legend.background = element_blank(),     
    legend.margin = margin(t = 0, r = 5, b = 0, l = 5)  
  ) +
  guides(linewidth = "none")

print(p3)

# we can do t+1 as well, but since the pattern is the same so we skip. 
# but we can do lagged point plot for t and t+1, for V1

plot_data <- data.frame(
  t = ds[,1],        
  t_plus_1 = ds_tp1[,1]  
)


p4 <- ggplot(plot_data, aes(x = t, y = t_plus_1)) +
  geom_point() +  
  labs(
    x = "t",      
    y = "t+1",    
    title = "Lagged Point Plot for species 1 (NM_Cd5)"  
  ) +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5,margin = margin(b = 15)),
    geom_point(aes(y = V1), size = 1.2),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    panel.border = element_rect(colour = "black",  fill = NA, linewidth = 2))


print(p4)

# Out-sample ----
if(out.sample|nout!=0){
  d.test <- do[nin:(ndo-1),]                 # Out-of-sample dataset at time t 
  dt_tp1 <- do[(nin+1):ndo,]                 # Out-of-sample dataset at time t+1
  ds.test <- (d.test-repmat(do.mean,nrow(d.test),1))*repmat(do.sd,nrow(d.test),1)^-1 # Normalized out-of-sample dataset at time t
  dst_tp1 <- (dt_tp1-repmat(do.mean,nrow(dt_tp1),1))*repmat(do.sd,nrow(dt_tp1),1)^-1 # Normalized out-of-sample dataset at time t+1
}else{d.test <- dt_tp1 <- dst_tp1 <- ds.test <- NULL}

# Compiled data at time t -> '1-42' + '43-44'
ds.all <- rbind(ds,ds.test)
dim(ds.all) # 44 rows and 124 columns. Since we need to have lagged dataset, number of time steps of ds.all is 44, although ndo = 30.

# Finding optimal embedding dimension (Ed) and nonlinearity parameter ----
#############################################################
# Find the optimal embedding dimension & nonlinearity parameter for each variable 
# based on univariate simplex projection and S-map, respectively

# Univariate simplex projection
Emax <- 10
cri <- 'rmse' # model selection 
Ed <- NULL
forecast_skill_simplex <- NULL
for(i in 1:ncol(ds)){
  spx.i <- simplex(ds[,i],E=2:Emax)
  Ed <- c(Ed,spx.i[which.min(spx.i[,cri])[1],'E'])
  forecast_skill_simplex <- c(forecast_skill_simplex,spx.i[which.min(spx.i[,cri])[1],'rho'])
}
Ed # The optimal embedding dimension for each variable
forecast_skill_simplex # Forecast skills for each variable based on simplex projection

######################################################################
# Finding causal variables by CCM ----
# Find causal variables by CCM analysis for multiview embedding
# Warning: It is time consuming for calculating the causation for each node
# CCM causality test for all node pairs 
# do.CCM <- F 
if(do.CCM){ 
  ccm.out <- ccm.fast.demo(ds, Epair=T,cri=cri,Emax=Emax)
  ccm.sig <- ccm.out[['ccm.sig']]
  ccm.rho <- ccm.out[['ccm.rho']]
  if(save){
    # To avoid overwrite the original files, we save them with different names, 'XXX_NEW'.
    write.csv(ccm.sig, file.path("out", paste("ccm_sig_", da.name, "_nin", nin, "_demo_NEW.csv", sep="")), row.names=FALSE)
    
    write.csv(ccm.rho, file.path('out', paste('ccm_rho_',da.name,'_nin',nin,'_demo_NEW.csv',sep='')),row.names=F)
  }
}

ccm.sig <- read.csv(file.path('out',paste('ccm_sig_',da.name,'_nin',nin,'_demo_NEW.csv',sep='')),header=T,stringsAsFactors = F)
ccm.rho <- read.csv(file.path('out',paste('ccm_rho_',da.name,'_nin',nin,'_demo_NEW.csv',sep='')),header=T,stringsAsFactors = F)

# ~~~~ meaning of the two tables? ----


######################################################################
# Multiview embedding ----
# Perform multiview embedding analysis for each node/species
# Warning: It is time consuming for running multview embedding for each node/species. 
# ** ~ 1 min for one node using intel CORE i7, ThinkPad P1 in balanced mode. 
# do.multiview <- F
if(do.multiview){
  esele_lag <- esim.lag.demo(ds,ccm.rho,ccm.sig,Ed,kmax=10000,kn=100,max_lag=3,Emax=Emax)
  # To avoid overwrite the original files, we save them with different names, 'XXX_NEW'.
  if(save){write.csv(esele_lag,file.path('out', paste('eseleLag_',da.name,'_nin',nin,'_demo_NEW.csv',sep='')),row.names=F)}
}

esele <- read.csv(file.path('out', paste('eseleLag_',da.name,'_nin',nin,'_demo_NEW.csv',sep='')),header=T)


# So far so good ----


####################################################
## The computation of multiview distance
dmatrix.mv <- mvdist.demo(ds,ds.all,esele)
dmatrix.train.mvx <- dmatrix.mv[['dmatrix.train.mvx']]
dmatrix.test.mvx <- dmatrix.mv[['dmatrix.test.mvx']]

save.image(file = "NM_Cd5_workspace.RData")

# Then use server
# NO.7 NM_Cd5  END -------

######################################################
######################################################

# NO.8 NM_Cd15 START ------

# This 'EDM_test.R' file quickly go through the whole process of running rEDM on a test dataset, before using the real data. The test dataset is from Chang et al. 2021 Ecol Lett (Ricker model). Please refer to Chang et al. on how the data set is generated. The test dataset used is in 'data/result20191024_0_0_0_.csv'.  
# For real data, the analysis is conducted in “EDM_local.R” and "EDM_NM_Cd15.R" file

# Preparation ----
# check version of the installed rEDM package
packageVersion("rEDM") # note Chang et al. 2021 Ecol Lett use v1.2.3. Other versions may be incompatible. 

seed <- 49563
set.seed(seed)

# load intact functions from Chang et al. 
source('code/Demo_MDR_function.R')  

# Load original dataset ----
da.range <- 1:45 # Subsample for data analysis
out.sample <- T # T/F for out-of-sample forecast
if(out.sample){nout <- 2}else{nout <- 0}  # number of out-of-sample

# First look at the first 6 species over 45 time steps ----
df <- read.csv('data/filtered_for_EDM_network/ASV_tab.NM.Cd15.abs.EDM.csv',header=T,stringsAsFactors=F)
p1 <- ggplot(data.frame(df), aes(x = 1:nrow(df))) +
  geom_line(aes(y = ASV55, color = "ASV55"), linewidth = 0.7) +
  geom_line(aes(y = ASV153, color = "ASV153"), linewidth = 0.7) +
  geom_line(aes(y = ASV178, color = "ASV178"), linewidth = 0.7) +
  geom_line(aes(y = ASV277, color = "ASV277"), linewidth = 0.7) +
  geom_line(aes(y = ASV4, color = "ASV4"), linewidth = 0.7) +
  geom_line(aes(y = ASV4, color = "ASV4"), linewidth = 0.7) +
  labs(
    title = "NM_Cd15 (first 45*6)",
    x = "New Time",
    y = "number of DNA copies"
  ) +
  scale_color_manual(
    name = NULL,
    values = c(
      "ASV55" = "#BD210F", "ASV153" = "#F15D00", "ASV178" = "#DDB500",
      "ASV277" = "#20D6B5", "ASV4" = "#0092E0", "ASV4" = "#202023"
    )
  ) +
  theme(
    plot.title = element_text(
      size = 20,
      hjust = 0.5,
      margin = margin(b = 15)
    ),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    panel.border = element_rect(
      colour = "black", 
      fill = NA, 
      linewidth = 2),
    
    legend.position = c(0.995, 0.98),          
    legend.justification = c(1, 1),        
    legend.key.width = unit(1, "cm"),       
    legend.text = element_text(size = 14),
    legend.spacing.y = unit(0.2, "cm"),       
    legend.background = element_blank(),     
    legend.margin = margin(t = 0, r = 5, b = 0, l = 5)  
  ) +
  guides(linewidth = "none")
print(p1)



(da.name <- 'ASV_tab.NM.Cd15.abs.EDM')
do <- read.csv('data/filtered_for_EDM_network/ASV_tab.NM.Cd15.abs.EDM.csv',header=T,stringsAsFactors=F)
dot <- do[da.range,1] # get time column
do <- do[da.range,-1] # remove the first column (time) from the data frame
ndo <- nrow(do) # number of time steps
nin <- ndo-nout # library sample size

# take the 123 columns (total number of ASV),filted do data frame have 123 columns
do[, 2] #see column 2 "New time" 
ncol(do) #see total columns number
do_subset <- do[, 3:125]#only take all 'ASV' columns
do <- do_subset

# ~~~~ plot logistic map ----
# plot the normalized in-sample dataset at time t. Include V1 to V6.
p2 <- ggplot(data.frame(do), aes(x = 1:nrow(do))) +
  geom_line(aes(y = ASV55, color = "ASV55"), linewidth = 0.7) +
  geom_line(aes(y = ASV153, color = "ASV153"), linewidth = 0.7) +
  geom_line(aes(y = ASV178, color = "ASV178"), linewidth = 0.7) +
  geom_line(aes(y = ASV277, color = "ASV277"), linewidth = 0.7) +
  geom_line(aes(y = ASV4, color = "ASV4"), linewidth = 0.7) +
  geom_line(aes(y = ASV4, color = "ASV4"), linewidth = 0.7) +
  labs(
    title = "NM_Cd15 (first 45*6)",
    x = "New Time",
    y = "number of DNA copies"
  ) +
  scale_color_manual(
    name = NULL,
    values = c(
      "ASV55" = "#BD210F", "ASV153" = "#F15D00", "ASV178" = "#DDB500",
      "ASV277" = "#20D6B5", "ASV4" = "#0092E0", "ASV4" = "#202023"
    )
  ) +
  theme(
    plot.title = element_text(
      size = 20,
      hjust = 0.5,
      margin = margin(b = 15)
    ),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    panel.border = element_rect(
      colour = "black", 
      fill = NA, 
      linewidth = 2),
    
    legend.position = c(0.995, 0.98),          
    legend.justification = c(1, 1),           
    legend.key.width = unit(1, "cm"),       
    legend.text = element_text(size = 14),
    legend.spacing.y = unit(0.2, "cm"),      
    legend.background = element_blank(),      
    legend.margin = margin(t = 0, r = 5, b = 0, l = 5) 
  ) +
  guides(linewidth = "none")

print(p2)



# For real data, we may have 45 time steps and 100+ 'species' for each treatment with the same matrix structure. The 'species' can be ASV/OTU ID, taxon (at phylum, order, class, family, genus, or species level), depending on nature of data and research objectives. Also, we may focus on overlapping 'species' across treatments, if we want to compare the results across treatments.

dto <- apply(do,1,sum) # sum of each row (time step) in the data frame
dpo <- do*repmat(dto^-1,1,ncol(do)) # normalize the rows of the matrix 'do' so that each row sums to 1(?). This is often done in data processing to scale the data and make comparisons across different rows more meaningful. Similar to using rarefied 16S rRNA gene sequencing data(?).
# repmat() is the own function by sourcing 'Demo_MDR_function.R'. 
# repmat() makes a matrix with repeated column or row.

# check and shall see all values in the rowsum column equal to 1. 
as.data.frame(dpo) %>% mutate(rowsum = rowSums(across(everything())))


# Exclusion of rare species ----
# Here we skip this step otherwise all species will be omitted.
# We shall consider the followings for real dataset.
# *Threshold is upon decision*.
# pcri <- 0;bcri <- 10^-3; # criteria for selecting species based on proportion of presnece (pcri) and mean relative abundance (bri) 
# doind2 <- (apply(dpo,2,mean,na.rm=T)>(bcri))&((apply(do>0,2,sum,na.rm=T)/nrow(do))>pcri) # index for selected species 
# exsp2 <- setdiff(1:ncol(do),which(doind2))   # index for rare species 
# do <- do[,-exsp2]                            # Dataset excluded rare species


(nsp <- ncol(do)) # number of species
ndo <- nrow(do) # number of time steps
nin <- ndo-nout # number of in-sample time steps

# Mean and SD of each node/species ----
do.mean <- apply(do[1:nin,],2,mean,na.rm=T) # mean of each column (species abundance) in the first 28 rows (in-sample)
do.sd <- apply(do[1:nin,],2,sd,na.rm=T) # sd of each column (species abundance) in the first 28 rows (in-sample)

# Construct a sd(i,j) matrix ----
# (dimension = nsp*nsp, i.e., 6*6 here).
# if sd(i,j)>1, meaning j varies more than i
# if sd(i,j)<1, meaning i varies more than j
# if sd(i,j)=1, meaning i and j vary the same
dosdM <- repmat(c(do.sd)^-1,1,nsp)*repmat(c(do.sd),nsp,1) 

# check the sd(i,j)
dosdM; dim(dosdM)

# In-sample ----
d <- do[1:(nin-1),] # In-sample dataset at time t (time 1-42)
d_tp1 <- do[2:(nin),] # In-sample dataset at time t+1 (time 2-43)

# ~~ normalization (z-score) ----
ds <- (d-repmat(do.mean,nrow(d),1))*repmat(do.sd,nrow(d),1)^-1 # Normalized in-sample dataset at time t (i.e., *z-score = (x-mean)/sd)
ds_tp1 <- (d_tp1-repmat(do.mean,nrow(d_tp1),1))*repmat(do.sd,nrow(d_tp1),1)^-1 # Normalized in-sample dataset at time t+1 (i.e., z-score)

# ~~~~ plot logistic map ----
# plot the normalized in-sample dataset at time t. Include V1 to V6.
p3 <- ggplot(data.frame(ds), aes(x = 1:nrow(ds))) +
  geom_line(aes(y = ASV55, color = "ASV55"), linewidth = 0.7) +
  geom_line(aes(y = ASV153, color = "ASV153"), linewidth = 0.7) +
  geom_line(aes(y = ASV178, color = "ASV178"), linewidth = 0.7) +
  geom_line(aes(y = ASV277, color = "ASV277"), linewidth = 0.7) +
  geom_line(aes(y = ASV4, color = "ASV4"), linewidth = 0.7) +
  geom_line(aes(y = ASV4, color = "ASV4"), linewidth = 0.7) +
  labs(
    title = "In-sample (T1~87*6) NM_Cd15",
    x = "Time t",
    y = "Normalized"
  ) +
  scale_color_manual(
    name = NULL,
    values = c(
      "ASV55" = "#BD210F", "ASV153" = "#F15D00", "ASV178" = "#DDB500",
      "ASV277" = "#20D6B5", "ASV4" = "#0092E0", "ASV4" = "#202023"
    )
  ) +
  theme(
    plot.title = element_text(
      size = 20,
      hjust = 0.5,
      margin = margin(b = 15)
    ),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    panel.border = element_rect(
      colour = "black", 
      fill = NA, 
      linewidth = 2),
    
    legend.position = c(0.995, 0.98),          
    legend.justification = c(1, 1),          
    legend.key.width = unit(1, "cm"),       
    legend.text = element_text(size = 14),
    legend.spacing.y = unit(0.2, "cm"),      
    legend.background = element_blank(),      
    legend.margin = margin(t = 0, r = 5, b = 0, l = 5)  
  ) +
  guides(linewidth = "none")

print(p3)

# we can do t+1 as well, but since the pattern is the same so we skip. 
# but we can do lagged point plot for t and t+1, for V1

plot_data <- data.frame(
  t = ds[,1],        
  t_plus_1 = ds_tp1[,1] 
)

# ggplot2绘图
p4 <- ggplot(plot_data, aes(x = t, y = t_plus_1)) +
  geom_point() +  
  labs(
    x = "t",      
    y = "t+1",   
    title = "Lagged Point Plot for species 1 (NM_Cd15)"  
  ) +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5,margin = margin(b = 15)),
    geom_point(aes(y = V1), size = 1.2),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    panel.border = element_rect(colour = "black",  fill = NA, linewidth = 2))


print(p4)

# Out-sample ----
if(out.sample|nout!=0){
  d.test <- do[nin:(ndo-1),]                 # Out-of-sample dataset at time t 
  dt_tp1 <- do[(nin+1):ndo,]                 # Out-of-sample dataset at time t+1
  ds.test <- (d.test-repmat(do.mean,nrow(d.test),1))*repmat(do.sd,nrow(d.test),1)^-1 # Normalized out-of-sample dataset at time t
  dst_tp1 <- (dt_tp1-repmat(do.mean,nrow(dt_tp1),1))*repmat(do.sd,nrow(dt_tp1),1)^-1 # Normalized out-of-sample dataset at time t+1
}else{d.test <- dt_tp1 <- dst_tp1 <- ds.test <- NULL}

# Compiled data at time t -> '1-42' + '43-44'
ds.all <- rbind(ds,ds.test)
dim(ds.all) # 44 rows and 123 columns. Since we need to have lagged dataset, number of time steps of ds.all is 44, although ndo = 30.

# Finding optimal embedding dimension (Ed) and nonlinearity parameter ----
#############################################################
# Find the optimal embedding dimension & nonlinearity parameter for each variable 
# based on univariate simplex projection and S-map, respectively

# Univariate simplex projection
Emax <- 10
cri <- 'rmse' # model selection 
Ed <- NULL
forecast_skill_simplex <- NULL
for(i in 1:ncol(ds)){
  spx.i <- simplex(ds[,i],E=2:Emax)
  Ed <- c(Ed,spx.i[which.min(spx.i[,cri])[1],'E'])
  forecast_skill_simplex <- c(forecast_skill_simplex,spx.i[which.min(spx.i[,cri])[1],'rho'])
}
Ed # The optimal embedding dimension for each variable
forecast_skill_simplex # Forecast skills for each variable based on simplex projection

######################################################################
# Finding causal variables by CCM ----
# Find causal variables by CCM analysis for multiview embedding
# Warning: It is time consuming for calculating the causation for each node
# CCM causality test for all node pairs 
# do.CCM <- F 
if(do.CCM){ 
  ccm.out <- ccm.fast.demo(ds, Epair=T,cri=cri,Emax=Emax)
  ccm.sig <- ccm.out[['ccm.sig']]
  ccm.rho <- ccm.out[['ccm.rho']]
  if(save){
    # To avoid overwrite the original files, we save them with different names, 'XXX_NEW'.
    write.csv(ccm.sig, file.path("out", paste("ccm_sig_", da.name, "_nin", nin, "_demo_NEW.csv", sep="")), row.names=FALSE)
    
    write.csv(ccm.rho, file.path('out', paste('ccm_rho_',da.name,'_nin',nin,'_demo_NEW.csv',sep='')),row.names=F)
  }
}

ccm.sig <- read.csv(file.path('out',paste('ccm_sig_',da.name,'_nin',nin,'_demo_NEW.csv',sep='')),header=T,stringsAsFactors = F)
ccm.rho <- read.csv(file.path('out',paste('ccm_rho_',da.name,'_nin',nin,'_demo_NEW.csv',sep='')),header=T,stringsAsFactors = F)

# ~~~~ meaning of the two tables? ----


######################################################################
# Multiview embedding ----
# Perform multiview embedding analysis for each node/species
# Warning: It is time consuming for running multview embedding for each node/species. 
# ** ~ 1 min for one node using intel CORE i7, ThinkPad P1 in balanced mode. 
# do.multiview <- F
if(do.multiview){
  esele_lag <- esim.lag.demo(ds,ccm.rho,ccm.sig,Ed,kmax=10000,kn=100,max_lag=3,Emax=Emax)
  # To avoid overwrite the original files, we save them with different names, 'XXX_NEW'.
  if(save){write.csv(esele_lag,file.path('out', paste('eseleLag_',da.name,'_nin',nin,'_demo_NEW.csv',sep='')),row.names=F)}
}

esele <- read.csv(file.path('out', paste('eseleLag_',da.name,'_nin',nin,'_demo_NEW.csv',sep='')),header=T)


# So far so good ----


####################################################
## The computation of multiview distance
dmatrix.mv <- mvdist.demo(ds,ds.all,esele)
dmatrix.train.mvx <- dmatrix.mv[['dmatrix.train.mvx']]
dmatrix.test.mvx <- dmatrix.mv[['dmatrix.test.mvx']]

save.image(file = "NM_Cd15_workspace.RData")

# Then use server
# NO.8 NM_Cd15  END -------