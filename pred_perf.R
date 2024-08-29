
########### PREDICTIVE PERFORMANCE: REPETITIONS OF TREE, SCALE_FREE, HUB, CLUSTERED GRAPHS WITH 4TH SIM REGIME ###########
library(igraph)
library(vars)
library(forecast)
library(reshape2)
library(ggplot2)
library(dplyr)
#set.seed(12) used for 0.4 density
#seeds <- sample(seq(1,1000),50) used for 0.4 density
## Regime 4
set.seed(16)
seeds <- sample(1:500,20,replace=FALSE)
#SBM parameters
#probmat <- matrix(c(.7,.2,.1,.7),nrow = 2) for 0.4 density
probmat <- matrix(c(.2,.02,.02,.2),nrow = 2)
rmse_mat <- matrix(ncol = 20, nrow = 20)
colnames(rmse_mat) <- c("gnarnei_tr","gnar_tr","var_tr","ar_tr",
                        "gnarnei_sf","gnar_sf","var_sf","ar_sf",
                        "gnarnei_sbm","gnar_sbm","var_sbm","ar_sbm",
                        "gnarnei_hub","gnar_hub","var_hub","ar_hub",
                        "gnarnei_er","gnar_er","var_er","ar_er")
if (is.null(colnames(rmse_mat))) {
  stop("Column names of rmse_mat are not set correctly")
}
globalalpha <- TRUE
                        
for (i in 1:20){
  print(i)
  set.seed(seeds[i])
  
  # Tree network
  tr <- make_tree(n = 20, children = 3, mode = 'undirected')
  net_tr <- igraphtoGNAR(tr) 
  
  #Scale-free network
  sf <-  barabasi.game(n = 20, m = 4, directed = FALSE)
  net_sf <- igraphtoGNAR(sf) 
  
  # Sbm network
  sbm <- sample_sbm(sum(20), pref.matrix = probmat, block.sizes = c(10,10))
  net_sbm <- igraphtoGNAR(sbm)
  
  # Hub network
  hub <- make_star(20, mode = "undirected")
  net_hub <- igraphtoGNAR(hub)
  
  # Random network
  er <- erdos.renyi.game(20,p.or.m = 38,type = "gnm",directed = FALSE) 
  net_er <- igraphtoGNAR(er)
  
  # Simulate network data based on the GNAR model
  # Normalize the data
  
  ts_tr <- GNARsim(n = 200, net=net_tr, alphaParams = list(rep(0.2,20), rep(0.4,20), rep(-0.6,20)), betaParams = list(c(0.3,.1),c(0.1,.1),c(-0.2,.3)))
  #tsn_tr <- apply(ts_tr, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_tr <- ts_tr[1:199, ]
  
  ts_sf <- GNARsim(n = 200, net = net_sf, alphaParams = list(rep(0.2,20), rep(0.4,20), rep(-0.6,20)), betaParams = list(c(0.3,.1),c(0.1,.1),c(-0.2,.3)))
  #tsn_sf <- apply(ts_sf, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_sf <- ts_sf[1:199, ]
  
  ts_sbm <- GNARsim(n = 200, net = net_sbm, alphaParams = list(rep(0.2,20), rep(0.4,20), rep(-0.6,20)), betaParams = list(c(0.3,.1),c(0.1,.1),c(-0.2,.3)))
  #tsn_sbm <- apply(ts_sbm, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_sbm <- ts_sbm[1:199, ]
  
  ts_hub <- GNARsim(n = 200, net = net_hub, alphaParams = list(rep(0.2,20), rep(0.4,20), rep(-0.6,20)), betaParams = list(c(0.3,.1),c(0.1,.1),c(-0.2,.3)))
  #tsn_hub <- apply(ts_hub, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_hub <- ts_hub[1:199, ]
  
  ts_er <- GNARsim(n = 200, net = net_er, alphaParams = list(rep(0.2,20), rep(0.4,20), rep(-0.6,20)), betaParams = list(c(0.3,.1),c(0.1,.1),c(-0.2,.3)))
  #tsn_er <- apply(ts_er, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_er <- ts_er[1:199, ]
  
  for(j in c("tr","sf","sbm", "hub", "er")){
    datasim_train <- get(paste("datasim_train_",j,sep = ""))
    ts <- get(paste("ts_",j,sep = ""))
    net <- get(paste("net_",j,sep = ""))
    
    # Fit GNAR(3,[2,2,2])
    fit_pred <- predict(GNARfit(vts = datasim_train, net = net,
                                alphaOrder = 3, betaOrder = rep(2,3),
                                globalalpha = globalalpha))
    rmse_mat[i,paste("gnarnei_",j,sep = "")] <- sqrt(mean((fit_pred-ts[200,])^2))
    print(c("seed: ",i," gnar nei done"))
    
    # Fit GNAR(3,[0,0,0])
    simplevar_pred <- predict(GNARfit(vts = datasim_train, net = net,
                                      alphaOrder = 3, betaOrder = rep(0,3),
                                      globalalpha = globalalpha))
    rmse_mat[i,paste("gnar_",j,sep = "")] <- sqrt(mean((simplevar_pred-ts[200,])^2))
    print(c("seed: ",i," gnar nonei done"))
    
    # Fit VAR(1)
    #simtrainvar <- t(datasim_train)
    datasim_train[is.na(datasim_train)] <- 0
    varforecast <- predict(restrict(VAR(datasim_train,p=3,type = "none")),n.ahead=1)
    getfcst <- function(x){return(x[1])}
    varfor <- unlist(lapply(varforecast$fcst, getfcst))
    rmse_mat[i,paste("var_",j,sep = "")] <- sqrt(mean((varfor-ts[200,])^2))
    print(c("seed: ",i," var done"))
    
    # Fit simple AR max lag 3
    simple_ar <- apply(datasim_train, 2, function(x){forecast(auto.arima(x,d=0,D=0,max.p = 3,max.q = 0,max.P = 0,max.Q = 0,stationary = TRUE,seasonal = FALSE,
                                                                       ic="bic",allowmean = FALSE,allowdrift = FALSE,trace = FALSE),h=1)$mean})
    rmse_mat[i,paste("ar_",j,sep = "")] <- sqrt(mean((simple_ar-ts[200,])^2))
    print(c("seed: ",i," ar done"))
  }
}

library(patchwork)
## PLOTS RESULTS
# side by side for all models
rmse_df <- as.data.frame(rmse_mat)
colnames(rmse_df) <- rep(c("GNAR(3,[2,2,2])","GNAR(3,[0,0,0])","VAR","AR"),5)
rmse_melt <- melt(rmse_df)
#rmse_melt["model"] <- c(rep("GRG",200))#,rep("SBM",200))
rmse_melt$variable <- as.factor(rmse_melt$variable)
colnames(rmse_melt) <- c("m","rmse")
p1 <- ggplot(rmse_melt[1:80, ], aes(x=m, y=rmse)) +
  geom_boxplot(fill="red")+ theme(axis.text.x = element_text(angle = 60, hjust = 1))+xlab(" ")+ylim(0.4,1.8)+ggtitle("Tree")+theme(plot.title = element_text(hjust = 0.5))
p2 <- ggplot(rmse_melt[81:160, ], aes(x=m, y=rmse)) +
  geom_boxplot(fill="yellowgreen")+ theme(axis.text.x = element_text(angle = 60, hjust = 1))+xlab(" ")+ylim(0.4,1.8)+ggtitle("Scale-free")+theme(plot.title = element_text(hjust = 0.5))+ylab("")
p3 <- ggplot(rmse_melt[161:240, ], aes(x=m, y=rmse)) +
  geom_boxplot(fill="cornflowerblue")+ theme(axis.text.x = element_text(angle = 60, hjust = 1))+xlab(" ")+ylim(0.4,1.8)+ggtitle("SBM")+
  theme(plot.title = element_text(hjust = 0.5))+ylab("") 
p4 <- ggplot(rmse_melt[241:320, ], aes(x=m, y=rmse)) +
  geom_boxplot(fill="pink")+ theme(axis.text.x = element_text(angle = 60, hjust = 1))+xlab(" ")+ylim(0.4,1.8)+ggtitle("Hub")+
  theme(plot.title = element_text(hjust = 0.5))+ylab("")

p5 <- ggplot(rmse_melt[321:400, ], aes(x=m, y=rmse)) +
  geom_boxplot(fill="purple")+ theme(axis.text.x = element_text(angle = 60, hjust = 1))+xlab(" ")+ylim(0.4,1.8)+ggtitle("Random")+
  theme(plot.title = element_text(hjust = 0.5))+ylab("")


p1+p2+p3+p4+p5


############ Residual analysis
# Tree network
tr <- make_tree(n = 20, children = 3, mode = 'undirected')
net_tr <- igraphtoGNAR(tr) 
ts_tr <- GNARsim(n = 200, net=net_tr, alphaParams = list(rep(0.2,20), rep(0.4,20), rep(-0.6,20)), betaParams = list(c(0.3,.1),c(0.1,.1),c(-0.2,.3)))
#tsn_tr <- apply(ts_tr, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
datasim_train_tr <- ts_tr[1:199, ]
model_tr <- GNARfit(vts = datasim_train_tr, net = net_tr,
                            alphaOrder = 3, betaOrder = rep(2,3),
                            globalalpha = globalalpha)
tr_resi <- residuals(model_tr)[,1]
layout(matrix(c(1,2),2,1))
plot(ts(residuals(model_tr)[,1]), ylab='tree GNAR model residuals')
hist(residuals(model_tr)[,1], main = '', xlab= 'tree GNAR model residuals')


# scale-free
sf <-  barabasi.game(n = 20, m = 4, directed = FALSE)
net_sf <- igraphtoGNAR(sf) 
ts_sf <- GNARsim(n = 200, net=net_sf, alphaParams = list(rep(0.2,20), rep(0.4,20), rep(-0.6,20)), betaParams = list(c(0.3,.1),c(0.1,.1),c(-0.2,.3)))
#tsn_tr <- apply(ts_tr, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
datasim_train_sf <- ts_sf[1:199, ]
model_sf <- GNARfit(vts = datasim_train_sf, net = net_sf,
                    alphaOrder = 3, betaOrder = rep(2,3),
                    globalalpha = globalalpha)
sf_resi <- residuals(model_sf)[,1]
layout(matrix(c(1,2),2,1))
plot(ts(residuals(model_tr)[,1]), ylab='tree GNAR model residuals')
hist(residuals(model_tr)[,1], main = '', xlab= 'tree GNAR model residuals')




set.seed(16)
seeds <- sample(1:500,20,replace=FALSE)
#SBM parameters
#probmat <- matrix(c(.7,.2,.1,.7),nrow = 2) for 0.4 density
probmat <- matrix(c(.2,.02,.02,.2),nrow = 2)
resi_mat <- matrix(ncol = 5, nrow = 196)
colnames(resi_mat) <- c("resi_tr","resi_sf","resi_sbm","resi_hub", "resi_er")
if (is.null(colnames(rmse_mat))) {
  stop("Column names of rmse_mat are not set correctly")
}
globalalpha <- TRUE

# Tree network
tr <- make_tree(n = 20, children = 3, mode = 'undirected')
net_tr <- igraphtoGNAR(tr) 
#Scale-free network
sf <-  barabasi.game(n = 20, m = 4, directed = FALSE)
net_sf <- igraphtoGNAR(sf) 
# Sbm network
sbm <- sample_sbm(sum(20), pref.matrix = probmat, block.sizes = c(10,10))
net_sbm <- igraphtoGNAR(sbm)
# Hub network
hub <- make_star(20, mode = "undirected")
net_hub <- igraphtoGNAR(hub)
# Random network
er <- erdos.renyi.game(20,p.or.m = 38,type = "gnm",directed = FALSE) 
net_er <- igraphtoGNAR(er)
# Simulate network data based on the GNAR model
ts_tr <- GNARsim(n = 200, net=net_tr, alphaParams = list(rep(0.2,20), rep(0.4,20), rep(-0.6,20)), betaParams = list(c(0.3,.1),c(0.1,.1),c(-0.2,.3)))
#tsn_tr <- apply(ts_tr, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
datasim_train_tr <- ts_tr[1:199, ]
  
ts_sf <- GNARsim(n = 200, net = net_sf, alphaParams = list(rep(0.2,20), rep(0.4,20), rep(-0.6,20)), betaParams = list(c(0.3,.1),c(0.1,.1),c(-0.2,.3)))
#tsn_sf <- apply(ts_sf, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
datasim_train_sf <- ts_sf[1:199, ]
  
ts_sbm <- GNARsim(n = 200, net = net_sbm, alphaParams = list(rep(0.2,20), rep(0.4,20), rep(-0.6,20)), betaParams = list(c(0.3,.1),c(0.1,.1),c(-0.2,.3)))
#tsn_sbm <- apply(ts_sbm, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
datasim_train_sbm <- ts_sbm[1:199, ]
  
ts_hub <- GNARsim(n = 200, net = net_hub, alphaParams = list(rep(0.2,20), rep(0.4,20), rep(-0.6,20)), betaParams = list(c(0.3,.1),c(0.1,.1),c(-0.2,.3)))
#tsn_hub <- apply(ts_hub, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
datasim_train_hub <- ts_hub[1:199, ]
  
ts_er <- GNARsim(n = 200, net = net_er, alphaParams = list(rep(0.2,20), rep(0.4,20), rep(-0.6,20)), betaParams = list(c(0.3,.1),c(0.1,.1),c(-0.2,.3)))
#tsn_er <- apply(ts_er, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
datasim_train_er <- ts_er[1:199, ]
  
  for(j in c("tr","sf","sbm", "hub", "er")){
    datasim_train <- get(paste("datasim_train_",j,sep = ""))
    ts <- get(paste("ts_",j,sep = ""))
    net <- get(paste("net_",j,sep = ""))
    
    # Fit GNAR(3,[2,2,2])
    model <- GNARfit(vts = datasim_train, net = net,
                                alphaOrder = 3, betaOrder = rep(2,3),
                                globalalpha = globalalpha)
    resi_mat[, paste("resi_",j,sep = "")] <- residuals(model)[,1]
  }


# Convert resi_mat to a data frame for easier plotting
resi_df <- as.data.frame(resi_mat)
colnames(resi_df) <- rep(c("Tree","Scale-free","Clustered","Hub", "Random"))
# Melt the data frame for ggplot2
library(reshape2)
resi_melted <- melt(resi_df, variable.name = "Network_Type", value.name = "Residuals")

# Plotting the box plot
library(ggplot2)
ggplot(resi_melted, aes(x = Network_Type, y = Residuals, fill = Network_Type)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Box Plot of Residuals for Different Network Types", x = "Network Type", y = "Residuals")


# Adding row numbers as time points
resi_df$Time <- 1:nrow(resi_df)

# Melt the data frame for ggplot2
resi_melted_line <- melt(resi_df, id.vars = "Time", variable.name = "Network_Type", value.name = "Residuals")

# Plotting the line plot
ggplot(resi_melted_line, aes(x = Time, y = Residuals, color = Network_Type)) +
  geom_line() +
  theme_minimal() +
  labs(title = "Line Plot of Residuals Over Time for Different Network Types", x = "Time", y = "Residuals")

# Plotting histograms for each network type in the same plot
ggplot(resi_melted, aes(x = Residuals, fill = Network_Type)) +
  geom_histogram(alpha = 0.5, position = "identity", bins = 30) +
  theme_minimal() +
  labs(title = "Histogram of Residuals for Different Network Types", x = "Residuals", y = "Frequency") +
  facet_wrap(~ Network_Type, scales = "free") # Faceting to show individual histograms for each type

# Perform the Shapiro-Wilk test for each column in resi_mat
shapiro_results <- apply(resi_mat, 2, shapiro.test)

# Extract and display the p-values from the Shapiro-Wilk test
shapiro_pvalues <- sapply(shapiro_results, function(result) result$p.value)

# Display the p-values for each network type
print(shapiro_pvalues)


####### Regime 1:  GNAR(1,[1]), alpha = 0.2, beta = 0.3
set.seed(16)
seeds <- sample(1:500,20,replace=FALSE)
probmat <- matrix(c(.2,.02,.02,.2),nrow = 2)
rmse_mat <- matrix(ncol = 10, nrow = 20)
colnames(rmse_mat) <- c("gnarnei_tr","gnar_tr",
                        "gnarnei_sf","gnar_sf",
                        "gnarnei_sbm","gnar_sbm",
                        "gnarnei_hub","gnar_hub",
                        "gnarnei_er","gnar_er")
if (is.null(colnames(rmse_mat))) {
  stop("Column names of rmse_mat are not set correctly")
}
globalalpha <- TRUE

for (i in 1:20){
  print(i)
  set.seed(seeds[i])
  
  # Tree network
  tr <- make_tree(n = 20, children = 3, mode = 'undirected')
  net_tr <- igraphtoGNAR(tr) 
  
  #Scale-free network
  sf <-  barabasi.game(n = 20, m = 4, directed = FALSE)
  net_sf <- igraphtoGNAR(sf) 
  
  # Sbm network
  sbm <- sample_sbm(sum(20), pref.matrix = probmat, block.sizes = c(10,10))
  net_sbm <- igraphtoGNAR(sbm)
  
  # Hub network
  hub <- make_star(20, mode = "undirected")
  net_hub <- igraphtoGNAR(hub)
  
  # Random network
  er <- erdos.renyi.game(20,p.or.m = 38,type = "gnm",directed = FALSE) 
  net_er <- igraphtoGNAR(er)
  
  # Simulate network data based on the GNAR model
  # Normalize the data
  
  ts_tr <- GNARsim(n = 200, net=net_tr, alphaParams = list(rep(0.2,20)), betaParams = list(c(0.3)))
  #tsn_tr <- apply(ts_tr, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_tr <- tsn_tr[1:199, ]
  
  ts_sf <- GNARsim(n = 200, net = net_sf, alphaParams = list(rep(0.2,20)), betaParams = list(c(0.3)))
  #tsn_sf <- apply(ts_sf, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_sf <- ts_sf[1:199, ]
  
  ts_sbm <- GNARsim(n = 200, net = net_sbm, alphaParams = list(rep(0.2,20)), betaParams = list(c(0.3)))
  #tsn_sbm <- apply(ts_sbm, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_sbm <- ts_sbm[1:199, ]
  
  ts_hub <- GNARsim(n = 200, net = net_hub, alphaParams = list(rep(0.2,20)), betaParams = list(c(0.3)))
  #tsn_hub <- apply(ts_hub, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_hub <- ts_hub[1:199, ]
  
  ts_er <- GNARsim(n = 200, net = net_er, alphaParams = list(rep(0.2,20)), betaParams = list(c(0.3)))
  #tsn_er <- apply(ts_er, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_er <- ts_er[1:199, ]
  
  for(j in c("tr","sf","sbm", "hub", "er")){
    datasim_train <- get(paste("datasim_train_",j,sep = ""))
    ts <- get(paste("ts_",j,sep = ""))
    net <- get(paste("net_",j,sep = ""))
    
    # Fit GNAR(3,[2,2,2])
    fit_pred <- predict(GNARfit(vts = datasim_train, net = net,
                                alphaOrder = 1, betaOrder = c(1),
                                globalalpha = globalalpha))
    rmse_mat[i,paste("gnarnei_",j,sep = "")] <- sqrt(mean((fit_pred-ts[200,])^2))
    print(c("seed: ",i," gnar nei done"))
    
    # Fit GNAR(3,[0,0,0])
    simplevar_pred <- predict(GNARfit(vts = datasim_train, net = net,
                                      alphaOrder = 1, betaOrder = c(0),
                                      globalalpha = globalalpha))
    rmse_mat[i,paste("gnar_",j,sep = "")] <- sqrt(mean((simplevar_pred-ts[200,])^2))
    print(c("seed: ",i," gnar nonei done"))
    
    # Fit VAR(1)
    #simtrainvar <- t(datasim_train)
    #datasim_train[is.na(datasim_train)] <- 0
    #varforecast <- predict(restrict(VAR(datasim_train,p=3,type = "none")),n.ahead=1)
    #getfcst <- function(x){return(x[1])}
    #varfor <- unlist(lapply(varforecast$fcst, getfcst))
    #rmse_mat[i,paste("var_",j,sep = "")] <- sqrt(mean((varfor-tsn[200,])^2))
    #print(c("seed: ",i," var done"))
    
    # Fit simple AR max lag 3
    #simple_ar <- apply(datasim_train, 2, function(x){forecast(auto.arima(x,d=0,D=0,max.p = 3,max.q = 0,max.P = 0,max.Q = 0,stationary = TRUE,seasonal = FALSE,
      #                                                                   ic="bic",allowmean = FALSE,allowdrift = FALSE,trace = FALSE),h=1)$mean})
    #rmse_mat[i,paste("ar_",j,sep = "")] <- sqrt(mean((simple_ar-ts[200,])^2))
    #print(c("seed: ",i," ar done"))
  }
}

library(patchwork)
## PLOTS RESULTS
# side by side for all models
rmse_df <- as.data.frame(rmse_mat)
#colnames(rmse_df) <- rep(c("GNAR(3,[2,2,2])","GNAR(3,[0,0,0])","VAR","AR"),4)
rmse_melt <- melt(rmse_df)
#rmse_melt["model"] <- c(rep("GRG",200))#,rep("SBM",200))
rmse_melt$variable <- as.factor(rmse_melt$variable)
colnames(rmse_melt) <- c("m","rmse")
p1 <- ggplot(rmse_melt[1:40, ], aes(x=m, y=rmse)) +
  geom_boxplot(fill="red")+ theme(axis.text.x = element_text(angle = 60, hjust = 1))+xlab(" ")+ylim(0.25,1.5)+ggtitle("Tree")+theme(plot.title = element_text(hjust = 0.5))
p2 <- ggplot(rmse_melt[41:80, ], aes(x=m, y=rmse)) +
  geom_boxplot(fill="yellowgreen")+ theme(axis.text.x = element_text(angle = 60, hjust = 1))+xlab(" ")+ylim(0.25,1.5)+ggtitle("Scale-free")+theme(plot.title = element_text(hjust = 0.5))+ylab("")
p3 <- ggplot(rmse_melt[81:120, ], aes(x=m, y=rmse)) +
  geom_boxplot(fill="cornflowerblue")+ theme(axis.text.x = element_text(angle = 60, hjust = 1))+xlab(" ")+ylim(0.25,1.5)+ggtitle("SBM")+
  theme(plot.title = element_text(hjust = 0.5))+ylab("") 
p4 <- ggplot(rmse_melt[121:160, ], aes(x=m, y=rmse)) +
  geom_boxplot(fill="pink")+ theme(axis.text.x = element_text(angle = 60, hjust = 1))+xlab(" ")+ylim(0.25,1.5)+ggtitle("Hub")+
  theme(plot.title = element_text(hjust = 0.5))+ylab("")
p5 <- ggplot(rmse_melt[161:200, ], aes(x=m, y=rmse)) +
  geom_boxplot(fill="purple")+ theme(axis.text.x = element_text(angle = 60, hjust = 1))+xlab(" ")+ylim(0.25,1.5)+ggtitle("Random")+
  theme(plot.title = element_text(hjust = 0.5))+ylab("")


p1+p2+p3+p4+p5


####### Regime 2:  GNAR(1,[2]), alpha = 0.2, beta = (0.3,0.4)
set.seed(16)
seeds <- sample(1:500,20,replace=FALSE)
probmat <- matrix(c(.2,.02,.02,.2),nrow = 2)
rmse_mat <- matrix(ncol = 10, nrow = 20)
colnames(rmse_mat) <- c("gnarnei_tr","gnar_tr",
                        "gnarnei_sf","gnar_sf",
                        "gnarnei_sbm","gnar_sbm",
                        "gnarnei_hub","gnar_hub",
                        "gnarnei_er","gnar_er")
if (is.null(colnames(rmse_mat))) {
  stop("Column names of rmse_mat are not set correctly")
}
globalalpha <- TRUE

for (i in 1:20){
  print(i)
  set.seed(seeds[i])
  
  # Tree network
  tr <- make_tree(n = 20, children = 3, mode = 'undirected')
  net_tr <- igraphtoGNAR(tr) 
  
  #Scale-free network
  sf <-  barabasi.game(n = 20, m = 4, directed = FALSE)
  net_sf <- igraphtoGNAR(sf) 
  
  # Sbm network
  sbm <- sample_sbm(sum(20), pref.matrix = probmat, block.sizes = c(10,10))
  net_sbm <- igraphtoGNAR(sbm)
  
  # Hub network
  hub <- make_star(20, mode = "undirected")
  net_hub <- igraphtoGNAR(hub)
  
  # Random network
  er <- erdos.renyi.game(20,p.or.m = 38,type = "gnm",directed = FALSE) 
  net_er <- igraphtoGNAR(er)
  
  # Simulate network data based on the GNAR model
  # Normalize the data
  
  ts_tr <- GNARsim(n = 200, net=net_tr, alphaParams = list(rep(0.2,20)), betaParams = list(c(0.3, 0.4)))
  #tsn_tr <- apply(ts_tr, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_tr <- ts_tr[1:199, ]
  
  ts_sf <- GNARsim(n = 200, net = net_sf, alphaParams = list(rep(0.2,20)), betaParams = list(c(0.3, 0.4)))
  #tsn_sf <- apply(ts_sf, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_sf <- ts_sf[1:199, ]
  
  ts_sbm <- GNARsim(n = 200, net = net_sbm, alphaParams = list(rep(0.2,20)), betaParams = list(c(0.3, 0.4)))
  #tsn_sbm <- apply(ts_sbm, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_sbm <- ts_sbm[1:199, ]
  
  ts_hub <- GNARsim(n = 200, net = net_hub, alphaParams = list(rep(0.2,20)), betaParams = list(c(0.3, 0.4)))
  #tsn_hub <- apply(ts_hub, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_hub <- ts_hub[1:199, ]
  
  ts_er <- GNARsim(n = 200, net = net_er, alphaParams = list(rep(0.2,20)), betaParams = list(c(0.3, 0.4)))
  #tsn_er <- apply(ts_er, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_er <- ts_er[1:199, ]
  
  for(j in c("tr","sf","sbm", "hub", "er")){
    datasim_train <- get(paste("datasim_train_",j,sep = ""))
    ts <- get(paste("ts_",j,sep = ""))
    net <- get(paste("net_",j,sep = ""))
    
    # Fit GNAR(3,[2,2,2])
    fit_pred <- predict(GNARfit(vts = datasim_train, net = net,
                                alphaOrder = 1, betaOrder = c(2),
                                globalalpha = globalalpha))
    rmse_mat[i,paste("gnarnei_",j,sep = "")] <- sqrt(mean((fit_pred-ts[200,])^2))
    print(c("seed: ",i," gnar nei done"))
    
    # Fit GNAR(3,[0,0,0])
    simplevar_pred <- predict(GNARfit(vts = datasim_train, net = net,
                                      alphaOrder = 1, betaOrder = c(0),
                                      globalalpha = globalalpha))
    rmse_mat[i,paste("gnar_",j,sep = "")] <- sqrt(mean((simplevar_pred-ts[200,])^2))
    print(c("seed: ",i," gnar nonei done"))
    
    # Fit VAR(1)
    #simtrainvar <- t(datasim_train)
    #datasim_train[is.na(datasim_train)] <- 0
    #varforecast <- predict(restrict(VAR(datasim_train,p=3,type = "none")),n.ahead=1)
    #getfcst <- function(x){return(x[1])}
    #varfor <- unlist(lapply(varforecast$fcst, getfcst))
    #rmse_mat[i,paste("var_",j,sep = "")] <- sqrt(mean((varfor-tsn[200,])^2))
    #print(c("seed: ",i," var done"))
    
    # Fit simple AR max lag 3
    #simple_ar <- apply(datasim_train, 2, function(x){forecast(auto.arima(x,d=0,D=0,max.p = 3,max.q = 0,max.P = 0,max.Q = 0,stationary = TRUE,seasonal = FALSE,
      #                                                                   ic="bic",allowmean = FALSE,allowdrift = FALSE,trace = FALSE),h=1)$mean})
    #rmse_mat[i,paste("ar_",j,sep = "")] <- sqrt(mean((simple_ar-tsn[200,])^2))
   # print(c("seed: ",i," ar done"))
  }
}

library(patchwork)
## PLOTS RESULTS
# side by side for all models
rmse_df <- as.data.frame(rmse_mat)
#colnames(rmse_df) <- rep(c("GNAR(3,[2,2,2])","GNAR(3,[0,0,0])","VAR","AR"),4)
rmse_melt <- melt(rmse_df)
#rmse_melt["model"] <- c(rep("GRG",200))#,rep("SBM",200))
rmse_melt$variable <- as.factor(rmse_melt$variable)
colnames(rmse_melt) <- c("m","rmse")
p1 <- ggplot(rmse_melt[1:40, ], aes(x=m, y=rmse)) +
  geom_boxplot(fill="red")+ theme(axis.text.x = element_text(angle = 60, hjust = 1))+xlab(" ")+ylim(0.25,1.5)+ggtitle("Tree")+theme(plot.title = element_text(hjust = 0.5))
p2 <- ggplot(rmse_melt[41:80, ], aes(x=m, y=rmse)) +
  geom_boxplot(fill="yellowgreen")+ theme(axis.text.x = element_text(angle = 60, hjust = 1))+xlab(" ")+ylim(0.25,1.5)+ggtitle("Scale-free")+theme(plot.title = element_text(hjust = 0.5))+ylab("")
p3 <- ggplot(rmse_melt[81:120, ], aes(x=m, y=rmse)) +
  geom_boxplot(fill="cornflowerblue")+ theme(axis.text.x = element_text(angle = 60, hjust = 1))+xlab(" ")+ylim(0.25,1.5)+ggtitle("SBM")+
  theme(plot.title = element_text(hjust = 0.5))+ylab("") 
p4 <- ggplot(rmse_melt[121:160, ], aes(x=m, y=rmse)) +
  geom_boxplot(fill="pink")+ theme(axis.text.x = element_text(angle = 60, hjust = 1))+xlab(" ")+ylim(0.25,1.5)+ggtitle("Hub")+
  theme(plot.title = element_text(hjust = 0.5))+ylab("")
p5 <- ggplot(rmse_melt[161:200, ], aes(x=m, y=rmse)) +
  geom_boxplot(fill="purple")+ theme(axis.text.x = element_text(angle = 60, hjust = 1))+xlab(" ")+ylim(0.25,1.5)+ggtitle("Random")+
  theme(plot.title = element_text(hjust = 0.5))+ylab("")


p1+p2+p3+p4+p5

####### Regime 3:  GNAR(3,[1,1,1]), alpha = c(0.2, 0.4, -0.6), beta = (0.2,0.1,-0.2)
set.seed(16)
seeds <- sample(1:500,20,replace=FALSE)
probmat <- matrix(c(.2,.02,.02,.2),nrow = 2)
rmse_mat <- matrix(ncol = 10, nrow = 20)
colnames(rmse_mat) <- c("gnarnei_tr","gnar_tr",
                        "gnarnei_sf","gnar_sf",
                        "gnarnei_sbm","gnar_sbm",
                        "gnarnei_hub","gnar_hub",
                        "gnarnei_er","gnar_er")
if (is.null(colnames(rmse_mat))) {
  stop("Column names of rmse_mat are not set correctly")
}
globalalpha <- TRUE

for (i in 1:20){
  print(i)
  set.seed(seeds[i])
  
  # Tree network
  tr <- make_tree(n = 20, children = 3, mode = 'undirected')
  net_tr <- igraphtoGNAR(tr) 
  
  #Scale-free network
  sf <-  barabasi.game(n = 20, m = 4, directed = FALSE)
  net_sf <- igraphtoGNAR(sf) 
  
  # Sbm network
  sbm <- sample_sbm(sum(20), pref.matrix = probmat, block.sizes = c(10,10))
  net_sbm <- igraphtoGNAR(sbm)
  
  # Hub network
  hub <- make_star(20, mode = "undirected")
  net_hub <- igraphtoGNAR(hub)
  
  # Random network
  er <- erdos.renyi.game(20,p.or.m = 38,type = "gnm",directed = FALSE) 
  net_er <- igraphtoGNAR(er)
  
  # Simulate network data based on the GNAR model
  # Normalize the data
  
  ts_tr <- GNARsim(n = 200, net=net_tr, alphaParams = list(rep(0.2,20), rep(0.4,20), rep(-0.6,20)), betaParams =list(c(0.2),c(0.1),c(-.2)))
  #tsn_tr <- apply(ts_tr, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_tr <- ts_tr[1:199, ]
  
  ts_sf <- GNARsim(n = 200, net = net_sf, alphaParams = list(rep(0.2,20), rep(0.4,20), rep(-0.6,20)), betaParams = list(c(0.2),c(0.1),c(-.2)))
  #tsn_sf <- apply(ts_sf, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_sf <- ts_sf[1:199, ]
  
  ts_sbm <- GNARsim(n = 200, net = net_sbm, alphaParams = list(rep(0.2,20), rep(0.4,20), rep(-0.6,20)), betaParams = list(c(0.2),c(0.1),c(-.2)))
  #tsn_sbm <- apply(ts_sbm, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_sbm <- ts_sbm[1:199, ]
  
  ts_hub <- GNARsim(n = 200, net = net_hub, alphaParams = list(rep(0.2,20), rep(0.4,20), rep(-0.6,20)), betaParams = list(c(0.2),c(0.1),c(-.2)))
  #tsn_hub <- apply(ts_hub, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_hub <- ts_hub[1:199, ]
  
  ts_er <- GNARsim(n = 200, net = net_er, alphaParams = list(rep(0.2,20), rep(0.4,20), rep(-0.6,20)), betaParams = list(c(0.2),c(0.1),c(-.2)))
  #tsn_er <- apply(ts_er, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_er <- ts_er[1:199, ]
  
  for(j in c("tr","sf","sbm", "hub", "er")){
    datasim_train <- get(paste("datasim_train_",j,sep = ""))
    ts <- get(paste("ts_",j,sep = ""))
    net <- get(paste("net_",j,sep = ""))
    
    # Fit GNAR(3,[2,2,2])
    fit_pred <- predict(GNARfit(vts = datasim_train, net = net,
                                alphaOrder = 3, betaOrder = rep(1,3),
                                globalalpha = globalalpha))
    rmse_mat[i,paste("gnarnei_",j,sep = "")] <- sqrt(mean((fit_pred-ts[200,])^2))
    print(c("seed: ",i," gnar nei done"))
    
    # Fit GNAR(3,[0,0,0])
    simplevar_pred <- predict(GNARfit(vts = datasim_train, net = net,
                                      alphaOrder = 3, betaOrder = rep(0,3),
                                      globalalpha = globalalpha))
    rmse_mat[i,paste("gnar_",j,sep = "")] <- sqrt(mean((simplevar_pred-ts[200,])^2))
    print(c("seed: ",i," gnar nonei done"))
    
    # Fit VAR(1)
    #simtrainvar <- t(datasim_train)
   # datasim_train[is.na(datasim_train)] <- 0
    #varforecast <- predict(restrict(VAR(datasim_train,p=3,type = "none")),n.ahead=1)
    #getfcst <- function(x){return(x[1])}
    #varfor <- unlist(lapply(varforecast$fcst, getfcst))
    #rmse_mat[i,paste("var_",j,sep = "")] <- sqrt(mean((varfor-tsn[200,])^2))
   # print(c("seed: ",i," var done"))
    
    # Fit simple AR max lag 3
    #simple_ar <- apply(datasim_train, 2, function(x){forecast(auto.arima(x,d=0,D=0,max.p = 3,max.q = 0,max.P = 0,max.Q = 0,stationary = TRUE,seasonal = FALSE,
    #                                                                     ic="bic",allowmean = FALSE,allowdrift = FALSE,trace = FALSE),h=1)$mean})
    #rmse_mat[i,paste("ar_",j,sep = "")] <- sqrt(mean((simple_ar-tsn[200,])^2))
    #print(c("seed: ",i," ar done"))
  }
}

## PLOTS RESULTS
# side by side for all models
rmse_df <- as.data.frame(rmse_mat)
#colnames(rmse_df) <- rep(c("GNAR(3,[2,2,2])","GNAR(3,[0,0,0])","VAR","AR"),4)
rmse_melt <- melt(rmse_df)
#rmse_melt["model"] <- c(rep("GRG",200))#,rep("SBM",200))
rmse_melt$variable <- as.factor(rmse_melt$variable)
colnames(rmse_melt) <- c("m","rmse")
p1 <- ggplot(rmse_melt[1:40, ], aes(x=m, y=rmse)) +
  geom_boxplot(fill="red")+ theme(axis.text.x = element_text(angle = 60, hjust = 1))+xlab(" ")+ylim(0.25,1.5)+ggtitle("Tree")+theme(plot.title = element_text(hjust = 0.5))
p2 <- ggplot(rmse_melt[41:80, ], aes(x=m, y=rmse)) +
  geom_boxplot(fill="yellowgreen")+ theme(axis.text.x = element_text(angle = 60, hjust = 1))+xlab(" ")+ylim(0.25,1.5)+ggtitle("Scale-free")+theme(plot.title = element_text(hjust = 0.5))+ylab("")
p3 <- ggplot(rmse_melt[81:120, ], aes(x=m, y=rmse)) +
  geom_boxplot(fill="cornflowerblue")+ theme(axis.text.x = element_text(angle = 60, hjust = 1))+xlab(" ")+ylim(0.25,1.5)+ggtitle("SBM")+
  theme(plot.title = element_text(hjust = 0.5))+ylab("") 
p4 <- ggplot(rmse_melt[121:160, ], aes(x=m, y=rmse)) +
  geom_boxplot(fill="pink")+ theme(axis.text.x = element_text(angle = 60, hjust = 1))+xlab(" ")+ylim(0.25,1.5)+ggtitle("Hub")+
  theme(plot.title = element_text(hjust = 0.5))+ylab("")
p5 <- ggplot(rmse_melt[161:200, ], aes(x=m, y=rmse)) +
  geom_boxplot(fill="purple")+ theme(axis.text.x = element_text(angle = 60, hjust = 1))+xlab(" ")+ylim(0.25,1.5)+ggtitle("Random")+
  theme(plot.title = element_text(hjust = 0.5))+ylab("")



p1+p2+p3+p4+p5


############### Model comparison ###################################
## number of nodes is 90 for large network
# Regime 4 : GNAR(3, [3,3,3]
set.seed(16)
seeds <- sample(1:500,20,replace=FALSE)
#SBM parameters
#probmat <- matrix(c(.7,.2,.1,.7),nrow = 2) for 0.4 density
probmat <- matrix(c(.1,.02,.02,.1),nrow = 2)
rmse_mat <- matrix(ncol = 10, nrow = 20)
colnames(rmse_mat) <- c("gnarnei_tr","gnar_tr",
                        "gnarnei_sf","gnar_sf",
                        "gnarnei_sbm","gnar_sbm",
                        "gnarnei_hub","gnar_hub",
                        "gnarnei_er","gnar_er")
if (is.null(colnames(rmse_mat))) {
  stop("Column names of rmse_mat are not set correctly")
}
globalalpha <- TRUE

for (i in 1:20){
  print(i)
  set.seed(seeds[i])
  
  # Tree network
  tr <- make_tree(n = 90, children = 4, mode = 'undirected') #density 0.02
  net_tr <- igraphtoGNAR(tr) 
  
  #Scale-free network
  sf <-  barabasi.game(n = 90, m = 2, directed = FALSE) # 0.04419476
  net_sf <- igraphtoGNAR(sf) 
  
  # Sbm network
  sbm <- sample_sbm(sum(90), pref.matrix = probmat, block.sizes = c(45,45)) # 0.04868914
  net_sbm <- igraphtoGNAR(sbm)
  
  # Hub network
  hub <- create_hub_network(num_hubs = 15, nodes_per_hub = 5) # 0.04494382
  net_hub <- igraphtoGNAR(hub)
  
  # Random network
  er <- erdos.renyi.game(90,p.or.m = 161,type = "gnm",directed = FALSE)  # 0.04
  net_er <- igraphtoGNAR(er)
  
  # Simulate network data based on the GNAR model
  # Normalize the data
  
  ts_tr <- GNARsim(n = 200, net=net_tr, alphaParams = list(rep(0.2,90), rep(0.4,90), rep(-0.6,90)), betaParams = list(c(0.3,.1, 0.05),c(0.1,.1, 0.05),c(-0.2,.3, 0.05)))
  #tsn_tr <- apply(ts_tr, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_tr <- ts_tr[1:199, ]
  
  ts_sf <- GNARsim(n = 200, net = net_sf, alphaParams = list(rep(0.2,90), rep(0.4,90), rep(-0.6,90)), betaParams = list(c(0.3,.1, 0.05),c(0.1,.1, 0.05),c(-0.2,.3, 0.05)))
  #tsn_sf <- apply(ts_sf, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_sf <- ts_sf[1:199, ]
  
  ts_sbm <- GNARsim(n = 200, net = net_sbm, alphaParams = list(rep(0.2,90), rep(0.4,90), rep(-0.6,90)), betaParams = list(c(0.3,.1, 0.05),c(0.1,.1, 0.05),c(-0.2,.3, 0.05)))
  #tsn_sbm <- apply(ts_sbm, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_sbm <- ts_sbm[1:199, ]
  
  ts_hub <- GNARsim(n = 200, net = net_hub, alphaParams = list(rep(0.2,90), rep(0.4,90), rep(-0.6,90)), betaParams = list(c(0.3,.1, 0.05),c(0.1,.1, 0.05),c(-0.2,.3, 0.05)))
  #tsn_hub <- apply(ts_hub, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_hub <- ts_hub[1:199, ]
  
  ts_er <- GNARsim(n = 200, net = net_er, alphaParams = list(rep(0.2,90), rep(0.4,90), rep(-0.6,90)), betaParams = list(c(0.3,.1, 0.05),c(0.1,.1, 0.05),c(-0.2,.3, 0.05)))
  #tsn_er <- apply(ts_er, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_er <- ts_er[1:199, ]
  
  for(j in c("tr","sf","sbm", "hub", "er")){
    datasim_train <- get(paste("datasim_train_",j,sep = ""))
    ts <- get(paste("ts_",j,sep = ""))
    net <- get(paste("net_",j,sep = ""))
    
    # Fit GNAR(3,[2,2,2])
    fit_pred <- predict(GNARfit(vts = datasim_train, net = net,
                                alphaOrder = 3, betaOrder = rep(3,3),
                                globalalpha = globalalpha))
    rmse_mat[i,paste("gnarnei_",j,sep = "")] <- sqrt(mean((fit_pred-ts[200,])^2))
    print(c("seed: ",i," gnar nei done"))
    
    # Fit GNAR(3,[0,0,0])
    simplevar_pred <- predict(GNARfit(vts = datasim_train, net = net,
                                      alphaOrder = 3, betaOrder = rep(0,3),
                                      globalalpha = globalalpha))
    rmse_mat[i,paste("gnar_",j,sep = "")] <- sqrt(mean((simplevar_pred-ts[200,])^2))
    print(c("seed: ",i," gnar nonei done"))
    
    
  }
}

library(patchwork)
## PLOTS RESULTS
# side by side for all models
rmse_df <- as.data.frame(rmse_mat)
#colnames(rmse_df) <- rep(c("GNAR(3,[2,2,2])","GNAR(3,[0,0,0])","VAR","AR"),4)
rmse_melt <- melt(rmse_df)
#rmse_melt["model"] <- c(rep("GRG",200))#,rep("SBM",200))
rmse_melt$variable <- as.factor(rmse_melt$variable)
colnames(rmse_melt) <- c("m","rmse")
p1 <- ggplot(rmse_melt[1:40, ], aes(x=m, y=rmse)) +
  geom_boxplot(fill="red")+ theme(axis.text.x = element_text(angle = 60, hjust = 1))+xlab(" ")+ylim(0.8,1.25)+ggtitle("Tree")+theme(plot.title = element_text(hjust = 0.5))
p2 <- ggplot(rmse_melt[41:80, ], aes(x=m, y=rmse)) +
  geom_boxplot(fill="yellowgreen")+ theme(axis.text.x = element_text(angle = 60, hjust = 1))+xlab(" ")+ylim(0.8,1.25)+ggtitle("Scale-free")+theme(plot.title = element_text(hjust = 0.5))+ylab("")
p3 <- ggplot(rmse_melt[81:120, ], aes(x=m, y=rmse)) +
  geom_boxplot(fill="cornflowerblue")+ theme(axis.text.x = element_text(angle = 60, hjust = 1))+xlab(" ")+ylim(0.8,1.25)+ggtitle("SBM")+
  theme(plot.title = element_text(hjust = 0.5))+ylab("") 
p4 <- ggplot(rmse_melt[121:160, ], aes(x=m, y=rmse)) +
  geom_boxplot(fill="pink")+ theme(axis.text.x = element_text(angle = 60, hjust = 1))+xlab(" ")+ylim(0.8,1.25)+ggtitle("Hub")+
  theme(plot.title = element_text(hjust = 0.5))+ylab("")
p5 <- ggplot(rmse_melt[161:200, ], aes(x=m, y=rmse)) +
  geom_boxplot(fill="purple")+ theme(axis.text.x = element_text(angle = 60, hjust = 1))+xlab(" ")+ylim(0.8,1.25)+ggtitle("Random")+
  theme(plot.title = element_text(hjust = 0.5))+ylab("")


p1+p2+p3+p4+p5


######### Regime 5: GNAR(4, [1,1,1,1])

probmat <- matrix(c(.1,.02,.02,.1),nrow = 2)
rmse_mat <- matrix(ncol = 10, nrow = 20)
colnames(rmse_mat) <- c("gnarnei_tr","gnar_tr",
                        "gnarnei_sf","gnar_sf",
                        "gnarnei_sbm","gnar_sbm",
                        "gnarnei_hub","gnar_hub",
                        "gnarnei_er","gnar_er")
if (is.null(colnames(rmse_mat))) {
  stop("Column names of rmse_mat are not set correctly")
}
globalalpha <- TRUE

for (i in 1:20){
  print(i)
  set.seed(seeds[i])
  
  # Tree network
  tr <- make_tree(n = 90, children = 4, mode = 'undirected') #density 0.02
  net_tr <- igraphtoGNAR(tr) 
  
  #Scale-free network
  sf <-  barabasi.game(n = 90, m = 2, directed = FALSE) # 0.04419476
  net_sf <- igraphtoGNAR(sf) 
  
  # Sbm network
  sbm <- sample_sbm(sum(90), pref.matrix = probmat, block.sizes = c(45,45)) # 0.04868914
  net_sbm <- igraphtoGNAR(sbm)
  
  # Hub network
  hub <- create_hub_network(num_hubs = 15, nodes_per_hub = 5) # 0.04494382
  net_hub <- igraphtoGNAR(hub)
  
  # Random network
  er <- erdos.renyi.game(90,p.or.m = 161,type = "gnm",directed = FALSE)  # 0.04
  net_er <- igraphtoGNAR(er)
  
  # Simulate network data based on the GNAR model
  # Normalize the data
  
  ts_tr <- GNARsim(n = 200, net=net_tr, alphaParams = list(rep(-0.6,90), rep(-0.4,90), rep(-0.2,90), rep(-0.1,90)), betaParams = list(c(0.1),c(0.1),c(0.3),c(0.05)))
  #tsn_tr <- apply(ts_tr, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_tr <- ts_tr[1:199, ]
  
  ts_sf <- GNARsim(n = 200, net = net_sf, alphaParams = list(rep(-0.6,90), rep(-0.4,90), rep(-0.2,90), rep(-0.1,90)), betaParams = list(c(0.1),c(0.1),c(0.3),c(0.05)))
  #tsn_sf <- apply(ts_sf, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_sf <- ts_sf[1:199, ]
  
  ts_sbm <- GNARsim(n = 200, net = net_sbm, alphaParams = list(rep(-0.6,90), rep(-0.4,90), rep(-0.2,90), rep(-0.1,90)), betaParams = list(c(0.1),c(0.1),c(0.3),c(0.05)))
  #tsn_sbm <- apply(ts_sbm, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_sbm <- ts_sbm[1:199, ]
  
  ts_hub <- GNARsim(n = 200, net = net_hub, alphaParams = list(rep(-0.6,90), rep(-0.4,90), rep(-0.2,90), rep(-0.1,90)), betaParams = list(c(0.1),c(0.1),c(0.3),c(0.05)))
  #tsn_hub <- apply(ts_hub, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_hub <- ts_hub[1:199, ]
  
  ts_er <- GNARsim(n = 200, net = net_er, alphaParams = list(rep(-0.6,90), rep(-0.4,90), rep(-0.2,90), rep(-0.1,90)), betaParams = list(c(0.1),c(0.1),c(0.3),c(0.05)))
  #tsn_er <- apply(ts_er, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_er <- ts_er[1:199, ]
  
  for(j in c("tr","sf","sbm", "hub", "er")){
    datasim_train <- get(paste("datasim_train_",j,sep = ""))
    ts <- get(paste("ts_",j,sep = ""))
    net <- get(paste("net_",j,sep = ""))
    
    # Fit GNAR(3,[2,2,2])
    fit_pred <- predict(GNARfit(vts = datasim_train, net = net,
                                alphaOrder = 4, betaOrder = rep(1,4),
                                globalalpha = globalalpha))
    rmse_mat[i,paste("gnarnei_",j,sep = "")] <- sqrt(mean((fit_pred-ts[200,])^2))
    print(c("seed: ",i," gnar nei done"))
    
    # Fit GNAR(3,[0,0,0])
    simplevar_pred <- predict(GNARfit(vts = datasim_train, net = net,
                                      alphaOrder = 4, betaOrder = rep(0,4),
                                      globalalpha = globalalpha))
    rmse_mat[i,paste("gnar_",j,sep = "")] <- sqrt(mean((simplevar_pred-ts[200,])^2))
    print(c("seed: ",i," gnar nonei done"))
    
    
  }
}

library(patchwork)
## PLOTS RESULTS
# side by side for all models
rmse_df <- as.data.frame(rmse_mat)
#colnames(rmse_df) <- rep(c("GNAR(3,[2,2,2])","GNAR(3,[0,0,0])","VAR","AR"),4)
rmse_melt <- melt(rmse_df)
#rmse_melt["model"] <- c(rep("GRG",200))#,rep("SBM",200))
rmse_melt$variable <- as.factor(rmse_melt$variable)
colnames(rmse_melt) <- c("m","rmse")
p1 <- ggplot(rmse_melt[1:40, ], aes(x=m, y=rmse)) +
  geom_boxplot(fill="red")+ theme(axis.text.x = element_text(angle = 60, hjust = 1))+xlab(" ")+ylim(0.85,1.25)+ggtitle("Tree")+theme(plot.title = element_text(hjust = 0.5))
p2 <- ggplot(rmse_melt[41:80, ], aes(x=m, y=rmse)) +
  geom_boxplot(fill="yellowgreen")+ theme(axis.text.x = element_text(angle = 60, hjust = 1))+xlab(" ")+ylim(0.85,1.25)+ggtitle("Scale-free")+theme(plot.title = element_text(hjust = 0.5))+ylab("")
p3 <- ggplot(rmse_melt[81:120, ], aes(x=m, y=rmse)) +
  geom_boxplot(fill="cornflowerblue")+ theme(axis.text.x = element_text(angle = 60, hjust = 1))+xlab(" ")+ylim(0.85,1.25)+ggtitle("SBM")+
  theme(plot.title = element_text(hjust = 0.5))+ylab("") 
p4 <- ggplot(rmse_melt[121:160, ], aes(x=m, y=rmse)) +
  geom_boxplot(fill="pink")+ theme(axis.text.x = element_text(angle = 60, hjust = 1))+xlab(" ")+ylim(0.85,1.25)+ggtitle("Hub")+
  theme(plot.title = element_text(hjust = 0.5))+ylab("")
p5 <- ggplot(rmse_melt[161:200, ], aes(x=m, y=rmse)) +
  geom_boxplot(fill="purple")+ theme(axis.text.x = element_text(angle = 60, hjust = 1))+xlab(" ")+ylim(0.85,1.25)+ggtitle("Random")+
  theme(plot.title = element_text(hjust = 0.5))+ylab("")

dev.off()
p1+p2+p3+p4+p5

########## regime 6 : GNAR(4,[2,2,2,2])

probmat <- matrix(c(.1,.02,.02,.1),nrow = 2)
rmse_mat <- matrix(ncol = 10, nrow = 20)
colnames(rmse_mat) <- c("gnarnei_tr","gnar_tr",
                        "gnarnei_sf","gnar_sf",
                        "gnarnei_sbm","gnar_sbm",
                        "gnarnei_hub","gnar_hub",
                        "gnarnei_er","gnar_er")
if (is.null(colnames(rmse_mat))) {
  stop("Column names of rmse_mat are not set correctly")
}
globalalpha <- TRUE

for (i in 1:20){
  print(i)
  set.seed(seeds[i])
  
  # Tree network
  tr <- make_tree(n = 90, children = 4, mode = 'undirected') #density 0.02
  net_tr <- igraphtoGNAR(tr) 
  
  #Scale-free network
  sf <-  barabasi.game(n = 90, m = 2, directed = FALSE) # 0.04419476
  net_sf <- igraphtoGNAR(sf) 
  
  # Sbm network
  sbm <- sample_sbm(sum(90), pref.matrix = probmat, block.sizes = c(45,45)) # 0.04868914
  net_sbm <- igraphtoGNAR(sbm)
  
  # Hub network
  hub <- create_hub_network(num_hubs = 15, nodes_per_hub = 5) # 0.04494382
  net_hub <- igraphtoGNAR(hub)
  
  # Random network
  er <- erdos.renyi.game(90,p.or.m = 161,type = "gnm",directed = FALSE)  # 0.04
  net_er <- igraphtoGNAR(er)
  
  # Simulate network data based on the GNAR model
  # Normalize the data
  
  ts_tr <- GNARsim(n = 200, net=net_tr, alphaParams = list(rep(-0.6,90), rep(-0.4,90), rep(-0.2,90), rep(-0.1,90)), betaParams = list(c(0.4,0.4),c(0.3,-0.4),c(0.5, -0.3),c(0.05,-0.1)))
  #tsn_tr <- apply(ts_tr, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_tr <- ts_tr[1:199, ]
  
  ts_sf <- GNARsim(n = 200, net = net_sf, alphaParams = list(rep(-0.6,90), rep(-0.4,90), rep(-0.2,90), rep(-0.1,90)), betaParams = list(c(0.4,0.4),c(0.3,-0.4),c(0.5, -0.3),c(0.05,-0.1)))
  #tsn_sf <- apply(ts_sf, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_sf <- ts_sf[1:199, ]
  
  ts_sbm <- GNARsim(n = 200, net = net_sbm, alphaParams = list(rep(-0.6,90), rep(-0.4,90), rep(-0.2,90), rep(-0.1,90)), betaParams = list(c(0.4,0.4),c(0.3,-0.4),c(0.5, -0.3),c(0.05,-0.1)))
  #tsn_sbm <- apply(ts_sbm, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_sbm <- ts_sbm[1:199, ]
  
  ts_hub <- GNARsim(n = 200, net = net_hub, alphaParams = list(rep(-0.6,90), rep(-0.4,90), rep(-0.2,90), rep(-0.1,90)), betaParams = list(c(0.4,0.4),c(0.3,-0.4),c(0.5, -0.3),c(0.05,-0.1)))
  #tsn_hub <- apply(ts_hub, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_hub <- ts_hub[1:199, ]
  
  ts_er <- GNARsim(n = 200, net = net_er, alphaParams = list(rep(-0.6,90), rep(-0.4,90), rep(-0.2,90), rep(-0.1,90)), betaParams = list(c(0.4,0.4),c(0.3,-0.4),c(0.5, -0.3),c(0.05,-0.1)))
  #tsn_er <- apply(ts_er, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_er <- ts_er[1:199, ]
  
  for(j in c("tr","sf","sbm", "hub", "er")){
    datasim_train <- get(paste("datasim_train_",j,sep = ""))
    ts <- get(paste("ts_",j,sep = ""))
    net <- get(paste("net_",j,sep = ""))
    
    # Fit GNAR(3,[2,2,2])
    fit_pred <- predict(GNARfit(vts = datasim_train, net = net,
                                alphaOrder = 4, betaOrder = rep(2,4),
                                globalalpha = globalalpha))
    rmse_mat[i,paste("gnarnei_",j,sep = "")] <- sqrt(mean((fit_pred-ts[200,])^2))
    print(c("seed: ",i," gnar nei done"))
    
    # Fit GNAR(3,[0,0,0])
    simplevar_pred <- predict(GNARfit(vts = datasim_train, net = net,
                                      alphaOrder = 4, betaOrder = rep(0,4),
                                      globalalpha = globalalpha))
    rmse_mat[i,paste("gnar_",j,sep = "")] <- sqrt(mean((simplevar_pred-ts[200,])^2))
    print(c("seed: ",i," gnar nonei done"))
    
    
  }
}

library(patchwork)
## PLOTS RESULTS
# side by side for all models
rmse_df <- as.data.frame(rmse_mat)
#colnames(rmse_df) <- rep(c("GNAR(3,[2,2,2])","GNAR(3,[0,0,0])","VAR","AR"),4)
rmse_melt <- melt(rmse_df)
#rmse_melt["model"] <- c(rep("GRG",200))#,rep("SBM",200))
rmse_melt$variable <- as.factor(rmse_melt$variable)
colnames(rmse_melt) <- c("m","rmse")
p1 <- ggplot(rmse_melt[1:40, ], aes(x=m, y=rmse)) +
  geom_boxplot(fill="red")+ theme(axis.text.x = element_text(angle = 60, hjust = 1))+xlab(" ")+ylim(0.75,12)+ggtitle("Tree")+theme(plot.title = element_text(hjust = 0.5))
p2 <- ggplot(rmse_melt[41:80, ], aes(x=m, y=rmse)) +
  geom_boxplot(fill="yellowgreen")+ theme(axis.text.x = element_text(angle = 60, hjust = 1))+xlab(" ")+ylim(0.75,1.25)+ggtitle("Scale-free")+theme(plot.title = element_text(hjust = 0.5))+ylab("")
p3 <- ggplot(rmse_melt[81:120, ], aes(x=m, y=rmse)) +
  geom_boxplot(fill="cornflowerblue")+ theme(axis.text.x = element_text(angle = 60, hjust = 1))+xlab(" ")+ylim(0.75,1.25)+ggtitle("SBM")+
  theme(plot.title = element_text(hjust = 0.5))+ylab("") 
p4 <- ggplot(rmse_melt[121:160, ], aes(x=m, y=rmse)) +
  geom_boxplot(fill="pink")+ theme(axis.text.x = element_text(angle = 60, hjust = 1))+xlab(" ")+ylim(0.75,1.25)+ggtitle("Hub")+
  theme(plot.title = element_text(hjust = 0.5))+ylab("")
p5 <- ggplot(rmse_melt[161:200, ], aes(x=m, y=rmse)) +
  geom_boxplot(fill="purple")+ theme(axis.text.x = element_text(angle = 60, hjust = 1))+xlab(" ")+ylim(0.75,1.25)+ggtitle("Random")+
  theme(plot.title = element_text(hjust = 0.5))+ylab("")


p1+p2+p3+p4+p5






#####################################################################################
# consider a range of different densities and produce different seeded networks for
# each fixed density, the produce ts for their edges and see distribution of rmse
#####################################################################################


# Tree Network

set.seed(16)
seeds <- sample(1:500,50,replace=FALSE)
rmse_list_nei<- vector(mode = "list", length = 4)
rmse_list_nonei<- vector(mode = "list", length = 4)
ind <- 1
globalalpha = TRUE

for (c in c(2,3,4,5)){
  for (i in 1:50){
    print(c(c," ",i))
    set.seed(seeds[i])
    tr <- make_tree(n = 20, children = c, mode = 'undirected')
    net_tr <- igraphtoGNAR(tr) 

    # Regime 4: GNAR(3,[2,2,2])
    alpha_par_4 <- list(rep(.2,20),rep(.4,20),rep(-0.6,20))
    beta_par_4 <- list(c(0.3,.1),c(0.1,.1),c(-0.2,.3))
    simdata <- GNARsim(n=200,net=net_tr,alphaParams=alpha_par_4,betaParams=beta_par_4)
    #datasim_train <- apply(simdata, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  
    fit_pred <- predict(GNARfit(vts = simdata[1:199,], net=net_tr, alphaOrder = 3, betaOrder = rep(2,3), globalalpha = globalalpha))
    rmse_list_nei[[ind]][i] <- sqrt(mean((fit_pred-simdata[200,])^2))
    
    simplevar_pred <- predict(GNARfit(vts = simdata[1:199,], net=net_tr, alphaOrder = 3, betaOrder = rep(0,3), globalalpha = globalalpha))
    rmse_list_nonei[[ind]][i] <- sqrt(mean((simplevar_pred-simdata[200, ])^2))
    
  }
  ind <- ind+1
}

rmse_df_nei <- do.call(cbind, rmse_list_nei)
rmse_df_nonei <- do.call(cbind, rmse_list_nonei)

library(ggplot2)
# side by side
rmse_df_nonei <- as.data.frame(rmse_df_nonei)
rmse_df_nonei["model"] <- rep("no neighbour",50)
colnames(rmse_df_nonei) <- c("2","3","4","5","model")
rmse_df_nei <- as.data.frame(rmse_df_nei)
rmse_df_nei["model"] <- rep("neighbour",50)
colnames(rmse_df_nei) <- c("2","3","4","5","model")
rmseall <- rbind(rmse_df_nei,rmse_df_nonei)
rmseall_melt <- melt(rmseall)
colnames(rmseall_melt) <- c("model","children_nodes","rmse")
dev.off()
ggplot(rmseall_melt, aes(x=children_nodes, y=rmse, fill=model)) +
  geom_boxplot()+ggtitle("Tree")+theme(plot.title = element_text(hjust = 0.5))


# ER model

# nedges=38,76,114,152
set.seed(16)
seeds <- sample(1:500,50,replace=FALSE)
rmse_list_nei<- vector(mode = "list", length = 4)
rmse_list_nonei<- vector(mode = "list", length = 4)
ind <- 1
for (ne in c(38,76,114,152)){
  for (i in 1:50){
    print(c(ne," ",i))
    set.seed(seeds[i])
    er <- erdos.renyi.game(20,p.or.m = ne,type = "gnm",directed = FALSE) 
    net_er <- igraphtoGNAR(er)
    
    # Regime 4: GNAR(3,[2,2,2])
    alpha_par_4 <- list(rep(.2,20),rep(.4,20),rep(-0.6,20))
    beta_par_4 <- list(c(0.3,.1),c(0.1,.1),c(-0.2,.3))
    simdata <- GNARsim(n=200,net=net_er,alphaParams=alpha_par_4,betaParams=beta_par_4
                             )
    #datasim_train <- apply(simdata, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
    fit_pred <- predict(GNARfit(vts = simdata[1:199,], net=net_er, alphaOrder = 3, betaOrder = rep(2,3), globalalpha = globalalpha))
    rmse_list_nei[[ind]][i] <- sqrt(mean((fit_pred-simdata[200,])^2))
    
    simplevar_pred <- predict(GNARfit(vts = simdata[1:199,], net=net_er, alphaOrder = 3, betaOrder = rep(0,3), globalalpha = globalalpha))
    rmse_list_nonei[[ind]][i] <- sqrt(mean((simplevar_pred-simdata[200, ])^2))
    
  }
  ind <- ind+1
}

rmse_df_nei <- do.call(cbind, rmse_list_nei)
rmse_df_nonei <- do.call(cbind, rmse_list_nonei)

# side by side
rmse_df_nonei <- as.data.frame(rmse_df_nonei)
rmse_df_nonei["model"] <- rep("no neighbour",50)
colnames(rmse_df_nonei) <- c("0.1","0.2","0.3","0.4","model")
summary(rmse_df_nonei)

rmse_df_nei <- as.data.frame(rmse_df_nei)
rmse_df_nei["model"] <- rep("neighbour",50)
colnames(rmse_df_nei) <- c("0.1","0.2","0.3","0.4","model")
summary(rmse_df_nei)

rmseall <- rbind(rmse_df_nei,rmse_df_nonei)


rmseall_melt <- melt(rmseall)
colnames(rmseall_melt) <- c("model","density","rmse")
ggplot(rmseall_melt, aes(x=density, y=rmse, fill=model)) +
  geom_boxplot()+ggtitle("ER random network")+theme(plot.title = element_text(hjust = 0.5))



# Scaled-free Network

set.seed(16)
seeds <- sample(1:500,50,replace=FALSE)
rmse_list_nei<- vector(mode = "list", length = 4)
rmse_list_nonei<- vector(mode = "list", length = 4)
ind <- 1
for (m in c(1,2,3,4)){
  for (i in 1:50){
    print(c(m," ",i))
    set.seed(seeds[i])
    sf <- barabasi.game(n = 20, m = m, directed = FALSE)
    net_sf <- igraphtoGNAR(sf)
    
    # Regime 4: GNAR(3,[2,2,2])
    alpha_par_4 <- list(rep(.2,20),rep(.4,20),rep(-0.6,20))
    beta_par_4 <- list(c(0.3,.1),c(0.1,.1),c(-0.2,.3))
    simdata <- GNARsim(n=200,net=net_er,alphaParams=alpha_par_4,betaParams=beta_par_4
    )
    #datasim_train <- apply(simdata, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
    fit_pred <- predict(GNARfit(vts = simdata[1:199,], net=net_sf, alphaOrder = 3, betaOrder = rep(2,3), globalalpha = globalalpha))
    rmse_list_nei[[ind]][i] <- sqrt(mean((fit_pred-simdata[200,])^2))
    
    simplevar_pred <- predict(GNARfit(vts = simdata[1:199,], net=net_sf, alphaOrder = 3, betaOrder = rep(0,3), globalalpha = globalalpha))
    rmse_list_nonei[[ind]][i] <- sqrt(mean((simplevar_pred-simdata[200, ])^2))
    
  }
  ind <- ind+1
}

rmse_df_nei <- do.call(cbind, rmse_list_nei)
rmse_df_nonei <- do.call(cbind, rmse_list_nonei)

# side by side
rmse_df_nonei <- as.data.frame(rmse_df_nonei)
rmse_df_nonei["model"] <- rep("no neighbour",50)
colnames(rmse_df_nonei) <- c("0.1","0.2","0.3","0.4","model")
summary(rmse_df_nonei)

rmse_df_nei <- as.data.frame(rmse_df_nei)
rmse_df_nei["model"] <- rep("neighbour",50)
colnames(rmse_df_nei) <- c("0.1","0.2","0.3","0.4","model")
summary(rmse_df_nei)

rmseall <- rbind(rmse_df_nei,rmse_df_nonei)


rmseall_melt <- melt(rmseall)
colnames(rmseall_melt) <- c("model","density","rmse")
ggplot(rmseall_melt, aes(x=density, y=rmse, fill=model)) +
  geom_boxplot()+ggtitle("Scale-free network")+theme(plot.title = element_text(hjust = 0.5))



# Clustered Network

set.seed(16)
seeds <- sample(1:500,50,replace=FALSE)
probmat <- list(matrix(c(.2,.02,.02,.2),nrow = 2),matrix(c(.4,.05,.05,.4),nrow = 2),
                matrix(c(.5,.1,.1,.5),nrow = 2),matrix(c(.7,.2,.2,.7),nrow = 2))
rmse_list_nei_sbm<- vector(mode = "list", length = 4)
rmse_list_nonei_sbm<- vector(mode = "list", length = 4)

for (ne in 1:4){
  for (i in 1:50){
    print(c(ne," ",i))
    set.seed(seeds[i])
    sbm <- sample_sbm(20,probmat[[ne]],c(10,10),directed = FALSE)
    net_sbm <- igraphtoGNAR(sbm)
    
    # Regime 4: GNAR(3,[2,2,2])
    alpha_par_4 <- list(rep(.2,20),rep(.4,20),rep(-0.6,20))
    beta_par_4 <- list(c(0.3,.1),c(0.1,.1),c(-0.2,.3))
    simdata <- GNARsim(n=200,net=net_sbm,alphaParams=alpha_par_4,betaParams=beta_par_4
                             )
    #datasim_train <- apply(simdata, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
    fit_pred <- predict(GNARfit(vts = simdata[1:199,], net=net_sbm, alphaOrder = 3, betaOrder = rep(2,3), globalalpha = globalalpha))
    rmse_list_nei_sbm[[ne]][i] <- sqrt(mean((fit_pred-simdata[200,])^2))
    
    simplevar_pred <- predict(GNARfit(vts = simdata[1:199,], net=net_sbm, alphaOrder = 3, betaOrder = rep(0,3), globalalpha = globalalpha))
    rmse_list_nonei_sbm[[ne]][i] <- sqrt(mean((simplevar_pred-simdata[200, ])^2))
    
  }
}

# side by side
rmse_df_nei <- do.call(cbind, rmse_list_nei_sbm)
rmse_df_nonei <- do.call(cbind, rmse_list_nonei_sbm)
rmse_df_nonei <- as.data.frame(rmse_df_nonei)
rmse_df_nonei["model"] <- rep("no neighbour",50)
colnames(rmse_df_nonei) <- c("0.1","0.2","0.3","0.4","model")
summary(rmse_df_nonei)

rmse_df_nei <- as.data.frame(rmse_df_nei)
rmse_df_nei["model"] <- rep("neighbour",50)
colnames(rmse_df_nei) <- c("0.1","0.2","0.3","0.4","model")
summary(rmse_df_nei)

rmseall <- rbind(rmse_df_nei,rmse_df_nonei)


rmseall_melt <- melt(rmseall)
colnames(rmseall_melt) <- c("model","density","rmse")
ggplot(rmseall_melt, aes(x=density, y=rmse, fill=model)) +
  geom_boxplot()+ggtitle("Clustered network")+theme(plot.title = element_text(hjust = 0.5))



# Hub Network

set.seed(16)
seeds <- sample(1:500,50,replace=FALSE)
rmse_list_nei_hub<- vector(mode = "list", length = 4)
rmse_list_nonei_hub<- vector(mode = "list", length = 4)
commat <- list(c(2,9), c(4,4), c(5,3))
globalalpha = TRUE


for (ne in 1:4) {
  for (i in 1:50) {
    print(c(ne, " ", i))
    set.seed(seeds[i])
    
    if (ne == 1) {
      hub <- make_star(20, "undirected")
    } else {
      hub <- create_hub_network(num_hubs = commat[[ne-1]][1], nodes_per_hub = commat[[ne-1]][2])
    }
    
    net_hub <- igraphtoGNAR(hub)
    
    # Regime 4: GNAR(3,[2,2,2])
    alpha_par_4 <- list(rep(.2, 20), rep(.4, 20), rep(-0.6, 20))
    beta_par_4 <- list(c(0.3, .1), c(0.1, .1), c(-0.2, .3))
    simdata <- GNARsim(n=200, net=net_hub, alphaParams=alpha_par_4, betaParams=beta_par_4)
    #datasim_train <- apply(simdata, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
    fit_pred <- predict(GNARfit(vts = simdata[1:199, ], net=net_hub, alphaOrder = 3, betaOrder = rep(2, 3), globalalpha = globalalpha))
    rmse_list_nei_hub[[ne]][i] <- sqrt(mean((fit_pred - simdata[200, ])^2))
    
    simplevar_pred <- predict(GNARfit(vts = simdata[1:199, ], net=net_hub, alphaOrder = 3, betaOrder = rep(0, 3), globalalpha = globalalpha))
    rmse_list_nonei_hub[[ne]][i] <- sqrt(mean((simplevar_pred - simdata[200, ])^2))
  }
}


# side by side
rmse_df_nei <- do.call(cbind, rmse_list_nei_hub)
rmse_df_nonei <- do.call(cbind, rmse_list_nonei_hub)
rmse_df_nonei <- as.data.frame(rmse_df_nonei)
rmse_df_nonei["model"] <- rep("no neighbour",50)
colnames(rmse_df_nonei) <- c("1 hub+19 nodes/hub", "2 hub+9 nodes/hub","4 hubs+4 nodes/hub","5 hubs+3 nodes/hub","model")
summary(rmse_df_nonei)

rmse_df_nei <- as.data.frame(rmse_df_nei)
rmse_df_nei["model"] <- rep("neighbour",50)
colnames(rmse_df_nei) <- c("1 hub+19 nodes/hub", "2 hub+9 nodes/hub","4 hubs+4 nodes/hub","5 hubs+3 nodes/hub","model")
summary(rmse_df_nei)

rmseall <- rbind(rmse_df_nei,rmse_df_nonei)


rmseall_melt <- melt(rmseall)
colnames(rmseall_melt) <- c("model","Combination","rmse")
ggplot(rmseall_melt, aes(x=Combination, y=rmse, fill=model)) +
  geom_boxplot()+ggtitle("Hub")+theme(plot.title = element_text(hjust = 0.5))



rmse_list_nei_hub1<- vector(mode = "list", length = 1)
rmse_list_nonei_hub1<- vector(mode = "list", length = 1)
globalalpha <- TRUE
for (i in 1:50){
  print(c(i))
  set.seed(seeds[i])
  hub <-  make_star(20,"undirected")
  net_hub <- igraphtoGNAR(hub)
  
  # Regime 4: GNAR(3,[2,2,2])
  alpha_par_4 <- list(rep(.2,20),rep(.4,20),rep(-0.6,20))
  beta_par_4 <- list(c(0.3,.1),c(0.1,.1),c(-0.2,.3))
  simdata <- GNARsim(n=200,net=net_hub,alphaParams=alpha_par_4,betaParams=beta_par_4
  )
  datasim_train <- apply(simdata, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  fit_pred <- predict(GNARfit(vts = datasim_train[1:199,], net=net_hub, alphaOrder = 3, betaOrder = rep(2,3), globalalpha = globalalpha))
  rmse_list_nei_hub1[[1]][i] <- sqrt(mean((fit_pred-simdata[200,])^2))
  
  simplevar_pred <- predict(GNARfit(vts = datasim_train[1:199,], net=net_hub, alphaOrder = 3, betaOrder = rep(0,3), globalalpha = globalalpha))
  rmse_list_nonei_hub1[[1]][i] <- sqrt(mean((simplevar_pred-simdata[200, ])^2))
  
}

# side by side
rmse_df_nei1 <- do.call(cbind, rmse_list_nei_hub1)
rmse_df_nonei1 <- do.call(cbind, rmse_list_nonei_hub1)
rmse_df_nonei1 <- as.data.frame(rmse_df_nonei1)
rmse_df_nonei1["model"] <- rep("no neighbour",50)
colnames(rmse_df_nonei1) <- c("1 hub + 19 nodes per hub","model")
summary(rmse_df_nonei1)

rmse_df_nei1 <- as.data.frame(rmse_df_nei1)
rmse_df_nei1["model"] <- rep("neighbour",50)
colnames(rmse_df_nei1) <- c("1 hub + 19 nodes per hub","model")
summary(rmse_df_nei1)

rmseall <- rbind(rmse_df_nei1,rmse_df_nonei1, rmse_df_nei, rmse_df_nonei)

rmseall_melt <- melt(rmseall)
colnames(rmseall_melt) <- c("model","Combination","rmse")
ggplot(rmseall_melt, aes(x=Combination, y=rmse, fill=model)) +
  geom_boxplot()+ggtitle("Hub")+theme(plot.title = element_text(hjust = 0.5))

################### Hub Network with density changing ##################


set.seed(16)
seeds <- sample(1:500,50,replace=FALSE)
rmse_list_nei_hub<- vector(mode = "list", length = 6)
rmse_list_nonei_hub<- vector(mode = "list", length = 6)
commat <- list(c(4,20), c(4,10), c(4,5), c(4,4), c(4,3), c(4,2))
globalalpha = TRUE

for (ne in 1:6){
  for (i in 1:50){
    print(c(ne," ",i))
    set.seed(seeds[i])
    hub <- create_hub_network(num_hubs = commat[[ne]][1], nodes_per_hub = commat[[ne]][2])
    net_hub <- igraphtoGNAR(hub)
    
    # Regime 4: GNAR(3,[2,2,2])
    n <- commat[[ne]][1]*commat[[ne]][2] + commat[[ne]][1]
    alpha_par_4 <- list(rep(.2,n),rep(.4,n),rep(-0.6,n))
    beta_par_4 <- list(c(0.3,.1),c(0.1,.1),c(-0.2,.3))
    simdata <- GNARsim(n=200,net=net_hub,alphaParams=alpha_par_4,betaParams=beta_par_4
    )
    #datasim_train <- apply(simdata, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
    fit_pred <- predict(GNARfit(vts = simdata[1:199,], net=net_hub, alphaOrder = 3, betaOrder = rep(2,3), globalalpha = globalalpha))
    rmse_list_nei_hub[[ne]][i] <- sqrt(mean((fit_pred-simdata[200,])^2))
    
    simplevar_pred <- predict(GNARfit(vts = simdata[1:199,], net=net_hub, alphaOrder = 3, betaOrder = rep(0,3), globalalpha = globalalpha))
    rmse_list_nonei_hub[[ne]][i] <- sqrt(mean((simplevar_pred-simdata[200, ])^2))
    
  }
}


# side by side
rmse_df_nei <- do.call(cbind, rmse_list_nei_hub)
rmse_df_nonei <- do.call(cbind, rmse_list_nonei_hub)
rmse_df_nonei <- as.data.frame(rmse_df_nonei)
rmse_df_nonei["model"] <- rep("no neighbour",50)
colnames(rmse_df_nonei) <- c("0.02","0.05","0.09","0.12","0.15","0.21","model")
summary(rmse_df_nonei)

rmse_df_nei <- as.data.frame(rmse_df_nei)
rmse_df_nei["model"] <- rep("neighbour",50)
colnames(rmse_df_nei) <- c("0.02","0.05","0.09","0.12","0.15","0.21","model")
summary(rmse_df_nei)

rmseall <- rbind(rmse_df_nei,rmse_df_nonei)


rmseall_melt <- melt(rmseall)
colnames(rmseall_melt) <- c("model","density","rmse")
ggplot(rmseall_melt, aes(x=density, y=rmse, fill=model)) +
  geom_boxplot()+ggtitle("Hub")+theme(plot.title = element_text(hjust = 0.5))



##############################################################################
############################ Robustness Analysis #############################
##############################################################################


############################## Noise Injection ###############################


evaluate_noise_impact <- function(net, alphaParams, betaParams,  repetitions = 10) {
  rmse_results <- matrix(0, ncol = 2, nrow = repetitions)
  colnames(rmse_results) <- c("original_rmse", "noisy_rmse")
  
  for (rep in 1:repetitions) {
    set.seed(42 + rep) # Different seed for each repetition
    net_gnar <- igraphtoGNAR(net)
    ts_data <- GNARsim(n = 200, net = net_gnar, alphaParams = alphaParams, betaParams = betaParams)
    #tsn_data <- apply(ts_data, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
    train_data <- ts_data[1:199, ]
    
    # Add noise
    #noisy_data <- train_data + rnorm(n = length(train_data), mean = 0, sd = noise_level)
    noisy_data <- GNARsim_t_noise(n = 200, net = net_gnar, 
                                  alphaParams = list(rep(0.2, 20), rep(0.4, 20), rep(-0.6, 20)), 
                                  betaParams = list(c(0.3, .1), c(0.1, .1), c(-0.2, .3)))
    #tsn_noisy_data <- apply(noisy_data, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
    train_noisy_data <- noisy_data[1:199, ]
    
    # Fit models
    fit_gnar_original <- GNARfit(vts = train_data, net = net_gnar, alphaOrder = 3, betaOrder = rep(2, 3), globalalpha = TRUE)
    fit_gnar_noisy <- GNARfit(vts = train_noisy_data, net = net_gnar, alphaOrder = 3, betaOrder = rep(2, 3), globalalpha = TRUE)
    
    # Predict and calculate RMSE
    original_predictions <- predict(fit_gnar_original)
    noisy_predictions <- predict(fit_gnar_noisy)
    
    original_rmse <- calculate_rmse(ts_data[200, ], original_predictions)
    noisy_rmse <- calculate_rmse(ts_data[200, ], noisy_predictions)
    
    rmse_results[rep, "original_rmse"] <- original_rmse
    rmse_results[rep, "noisy_rmse"] <- noisy_rmse
  }
  
  return(rmse_results)
}

# Define parameters
alphaParams <- list(rep(0.2, 20), rep(0.4, 20), rep(-0.6, 20))
betaParams <- list(c(0.3, .1), c(0.1, .1), c(-0.2, .3))
#noise_level <- 0.5 # You can increase this to see the effect

# Evaluate noise impact on different network structures
results_tree <- evaluate_noise_impact(net_tree, alphaParams, betaParams,  repetitions = 10)
results_sf <- evaluate_noise_impact(net_sf, alphaParams, betaParams, repetitions = 10)
results_sbm <- evaluate_noise_impact(net_sbm, alphaParams, betaParams,  repetitions = 10)
results_hub <- evaluate_noise_impact(net_hub, alphaParams, betaParams, repetitions = 10)
results_er <- evaluate_noise_impact(net_er, alphaParams, betaParams, repetitions = 10)

# Calculate the mean RMSE for each network type
mean_results_tree <- colMeans(results_tree)
mean_results_sf <- colMeans(results_sf)
mean_results_sbm <- colMeans(results_sbm)
mean_results_hub <- colMeans(results_hub)
mean_results_er <- colMeans(results_er)

# Print results
print("Tree Network:")
print(mean_results_tree)
print("Scale-Free Network:")
print(mean_results_sf)
print("SBM Network:")
print(mean_results_sbm)
print("Hub Network:")
print(mean_results_hub)
print("Erdos-Renyi Network:")
print(mean_results_er)

# Combine the data into a data frame
rmse_data <- data.frame(
  RMSE = c(results_hub[,1], results_hub[,2]),
  Type = rep(c("Original", "t-Noise"), each = 10)
)

# Load ggplot2 library
library(ggplot2)

# Create the boxplot
ggplot(rmse_data, aes(x = Type, y = RMSE)) +
  geom_boxplot() +
  labs(title = "Comparison of RMSE: Original vs t-Noise Data under Hub network", 
       x = "Time series Data Type", 
       y = "RMSE") +
  theme_minimal()

######################### Missing Data #########################################


evaluate_missing_data_impact <- function(net, alphaParams, betaParams, missing_fraction = 0.1, repetitions = 10) {
  rmse_results <- matrix(0, ncol = 2, nrow = repetitions)
  colnames(rmse_results) <- c("original_rmse", "imputed_rmse")
  
  for (rep in 1:repetitions) {
    set.seed(42 + rep) # Different seed for each repetition
    net_gnar <- igraphtoGNAR(net)
    ts_data <- GNARsim(n = 200, net = net_gnar, alphaParams = alphaParams, betaParams = betaParams)
    tsn_data <- apply(ts_data, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
    train_data <- tsn_data[1:199, ]
    
    # Introduce missing data
    missing_indices <- sample(1:length(train_data), size = round(missing_fraction * length(train_data)), replace = FALSE)
    train_data_with_na <- train_data
    train_data_with_na[missing_indices] <- NA
    
    # Impute missing data
    train_data_imputed <- apply(train_data_with_na, 2, function(x) {
      x[is.na(x)] <- mean(x, na.rm = TRUE)
      return(x)
    })
    
    # Fit models
    fit_gnar_original <- GNARfit(vts = train_data, net = net_gnar, alphaOrder = 3, betaOrder = rep(2, 3), globalalpha = TRUE)
    fit_gnar_imputed <- GNARfit(vts = train_data_imputed, net = net_gnar, alphaOrder = 3, betaOrder = rep(2, 3), globalalpha = TRUE)
    
    # Predict and calculate RMSE
    original_predictions <- predict(fit_gnar_original)
    imputed_predictions <- predict(fit_gnar_imputed)
    
    original_rmse <- calculate_rmse(tsn_data[200, ], original_predictions)
    imputed_rmse <- calculate_rmse(tsn_data[200, ], imputed_predictions)
    
    rmse_results[rep, "original_rmse"] <- original_rmse
    rmse_results[rep, "imputed_rmse"] <- imputed_rmse
  }
  
  return(rmse_results)
}

# Define parameters
alphaParams <- list(rep(0.2, 20), rep(0.4, 20), rep(-0.6, 20))
betaParams <- list(c(0.3, .1), c(0.1, .1), c(-0.2, .3))
missing_fraction <- 0.1 # Fraction of data to be made missing
repetitions <- 10 # Number of repetitions

# Evaluate missing data impact on different network structures
results_tree <- evaluate_missing_data_impact(net_tree, alphaParams, betaParams, missing_fraction, repetitions)
results_sf <- evaluate_missing_data_impact(net_sf, alphaParams, betaParams, missing_fraction, repetitions)
results_sbm <- evaluate_missing_data_impact(net_sbm, alphaParams, betaParams, missing_fraction, repetitions)
results_hub <- evaluate_missing_data_impact(net_hub, alphaParams, betaParams, missing_fraction, repetitions)
results_er <- evaluate_missing_data_impact(net_er, alphaParams, betaParams, missing_fraction, repetitions)

# Calculate the mean RMSE for each network type
mean_results_tree <- colMeans(results_tree)
mean_results_sf <- colMeans(results_sf)
mean_results_sbm <- colMeans(results_sbm)
mean_results_hub <- colMeans(results_hub)
mean_results_er <- colMeans(results_er)

# Print results
print("Tree Network:")
print(mean_results_tree)
print("Scale-Free Network:")
print(mean_results_sf)
print("SBM Network:")
print(mean_results_sbm)
print("Hub Network:")
print(mean_results_hub)
print("Erdos-Renyi Network:")
print(mean_results_er)


##############################################################################
############################ Uncertainty Quantification ######################
##############################################################################

################ Prediction Intervals:
# Function to perform Moving Block Bootstrap
moving_block_bootstrap <- function(data, block_size) {
  n <- nrow(data)
  num_blocks <- ceiling(n / block_size)
  indices <- sample(1:(n - block_size + 1), size = num_blocks, replace = TRUE)
  resampled_data <- do.call(rbind, lapply(indices, function(start) data[start:(start + block_size - 1), ]))
  return(resampled_data[1:n, ])
}




######Define a function to perform bootstrapping and compute prediction intervals
compute_prediction_intervals_mbb <- function(fit_gnar, data, net, n_bootstrap = 1000, block_size = 5,confidence_level = 0.95) {
  predictions <- matrix(0, ncol = ncol(data), nrow = n_bootstrap)
  
  for (i in 1:n_bootstrap) {
    # Resample the data with replacement
    resampled_data <-  moving_block_bootstrap(data, block_size)
    
    # Fit the GNAR model to the resampled data
    fit_resampled <- GNARfit(vts = resampled_data, net = net, 
                             alphaOrder = 3, betaOrder = rep(2, 3), 
                             globalalpha = TRUE)
    
    # Generate predictions
    predictions[i, ] <- predict(fit_resampled)
  }
  
  # Compute prediction intervals
  lower_bound <- apply(predictions, 2, function(x) quantile(x, (1 - confidence_level) / 2))
  upper_bound <- apply(predictions, 2, function(x) quantile(x, 1 - (1 - confidence_level) / 2))
  
  return(list(lower_bound = lower_bound, upper_bound = upper_bound))
}


############ Function to evaluate prediction intervals

compute_prediction_intervals_mbb_repetitions <- function(data,actual , net, alphaOrder, betaOrder, n_bootstrap = 1000, block_size = 5, confidence_level = 0.95, repetitions = 10) {
  all_coverage_probabilities <- numeric(repetitions)
  
  for (rep in 1:repetitions) {
    set.seed(42 + rep) # Different seed for each repetition
    fit_gnar <- GNARfit(vts = data, net = net, alphaOrder = alphaOrder, betaOrder = betaOrder, globalalpha = TRUE)
    
    predictions <- matrix(0, ncol = ncol(data), nrow = n_bootstrap)
    for (i in 1:n_bootstrap) {
      resampled_data <- moving_block_bootstrap(data, block_size)
      fit_resampled <- GNARfit(vts = resampled_data, net = net, alphaOrder = alphaOrder, betaOrder = betaOrder, globalalpha = TRUE)
      predictions[i, ] <- predict(fit_resampled)
    }
    
    lower_bound <- apply(predictions, 2, function(x) quantile(x, (1 - confidence_level) / 2))
    upper_bound <- apply(predictions, 2, function(x) quantile(x, 1 - (1 - confidence_level) / 2))
    
    within_interval <- (actual  >= lower_bound) & (actual  <= upper_bound)
    coverage_probability <- mean(within_interval)
    all_coverage_probabilities[rep] <- coverage_probability
  }
  
  mean_coverage_probability <- mean(all_coverage_probabilities)
  return(list(mean_coverage_probability = mean_coverage_probability, all_coverage_probabilities = all_coverage_probabilities))
}


# Function to evaluate prediction intervals for a given network structure with repetitions
evaluate_network_structure_mbb_repetitions <- function(net, alphaParams, betaParams, n_bootstrap = 1000, block_size = 5, confidence_level = 0.95, repetitions = 10) {
  set.seed(42)
  net_gnar <- igraphtoGNAR(net)
  ts_data <- GNARsim(n = 200, net = net_gnar, alphaParams = alphaParams, betaParams = betaParams)
  tsn_data <- apply(ts_data, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  train_data <- tsn_data[1:199, ]
  actual <- tsn_data[200,]
  
  results <- compute_prediction_intervals_mbb_repetitions(train_data,actual , net_gnar, alphaOrder = 3, betaOrder = c(2, 2, 2), n_bootstrap, block_size, confidence_level, repetitions)
  
  return(results)
}


# Define parameters
alphaParams <- list(rep(0.2, 20), rep(0.4, 20), rep(-0.6, 20))
betaParams <- list(c(0.3, .1), c(0.1, .1), c(-0.2, .3))
n_bootstrap <- 1000
block_size <- 5
confidence_level <- 0.95
repetitions <- 10

# Evaluate prediction intervals for different network structures
results_tree <- evaluate_network_structure_mbb_repetitions(net_tree, alphaParams, betaParams, n_bootstrap, block_size, confidence_level, repetitions)
results_sf <- evaluate_network_structure_mbb_repetitions(net_sf, alphaParams, betaParams, n_bootstrap, block_size, confidence_level, repetitions)
results_sbm <- evaluate_network_structure_mbb_repetitions(net_sbm, alphaParams, betaParams, n_bootstrap, block_size, confidence_level, repetitions)
results_hub <- evaluate_network_structure_mbb_repetitions(net_hub, alphaParams, betaParams, n_bootstrap, block_size, confidence_level, repetitions)
results_er <- evaluate_network_structure_mbb_repetitions(net_er, alphaParams, betaParams, n_bootstrap, block_size, confidence_level, repetitions)

# Print results
print("Tree Network:")
print(results_tree)
print("Scale-Free Network:")
print(results_sf)
print("SBM Network:")
print(results_sbm)
print("Hub Network:")
print(results_hub)
print("Erdos-Renyi Network:")
print(results_er)



### Function to plot coverage probabilities across repetitions
plot_coverage_probabilities <- function(results, title) {
  coverage_probabilities <- results$all_coverage_probabilities
  hist(coverage_probabilities, breaks = 10, main = title, xlab = "Coverage Probability", ylab = "Frequency")
  abline(v = mean(coverage_probabilities), col = "red", lwd = 2)
  legend("topright", legend = paste("Mean:", round(mean(coverage_probabilities), 3)), col = "red", lwd = 2)
}

# Plot coverage probabilities for different network structures
plot_coverage_probabilities(results_tree, "Tree Network")
plot_coverage_probabilities(results_sf, "Scale-Free Network")
plot_coverage_probabilities(results_sbm, "SBM Network")
plot_coverage_probabilities(results_hub, "Hub Network")
plot_coverage_probabilities(results_er, "Erdos-Renyi Network")






#####################################################################################
# Model Misspecification 
#####################################################################################

###### frist analysis
probmat <- matrix(c(.7,.2,.2,.7),nrow = 2) #density 0.4
rmse_mat <- matrix(ncol = 2,nrow = 20)
colnames(rmse_mat) <- c("gnarnei_sbm&net_sf", "gnarnei_sf&net_sbm")
for (i in 1:20){
  print(i)
  set.seed(seeds[i])
  
  #Scale-free network
  sf <-  barabasi.game(n = 20, m = 4, directed = FALSE)
  net_sf <- igraphtoGNAR(sf) 
  
  # Sbm network
  sbm <- sample_sbm(sum(20), pref.matrix = probmat, block.sizes = c(10,10))
  net_sbm <- igraphtoGNAR(sbm)
  
  # Simulate network data based on the GNAR model
  # Normalize the data
  
  
  ts_sf <- GNARsim(n = 200, net = net_sf, alphaParams = list(rep(0.2,20), rep(0.4,20), rep(-0.6,20)), betaParams = list(c(0.3,.1),c(0.1,.1),c(-0.2,.3)))
  tsn_sf <- apply(ts_sf, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_sf <- tsn_sf[1:199, ]
  
  ts_sbm <- GNARsim(n = 200, net = net_sbm, alphaParams = list(rep(0.2,20), rep(0.4,20), rep(-0.6,20)), betaParams = list(c(0.3,.1),c(0.1,.1),c(-0.2,.3)))
  tsn_sbm <- apply(ts_sbm, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_sbm <- tsn_sbm[1:199, ]
  
  # Fit GNAR(3, [2,2,2]) model with sbm data and sf network
  fit_pred1 <- predict(GNARfit(vts = datasim_train_sbm, net = net_sf,
                               alphaOrder = 3, betaOrder = rep(2,3),
                               globalalpha = globalalpha))
  rmse_mat[i,paste("gnarnei_sbm&net_sf")] <- sqrt(mean((fit_pred1-tsn_sbm[200,])^2))
  print(c("seed: ",i," gnar nei1 done"))
  
  fit_pred2 <- predict(GNARfit(vts = datasim_train_sf, net = net_sbm,
                               alphaOrder = 3, betaOrder = rep(2,3),
                               globalalpha = globalalpha))
  rmse_mat[i,paste("gnarnei_sf&net_sbm")] <- sqrt(mean((fit_pred2-tsn_sf[200,])^2))
  print(c("seed: ",i," gnar nei2 done"))
  
}

# Convert the matrix to a data frame
rmse_df <- as.data.frame(rmse_mat)

# Reshape the data frame from wide to long format
rmse_melt <- melt(rmse_df, variable.name = "Network_Model", value.name = "RMSE")

# Create the box plot
ggplot(rmse_melt, aes(x = Network_Model, y = RMSE, fill = Network_Model)) +
  geom_boxplot() +
  ggtitle("RMSE Box Plot for Different Network Models") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Network and Model") +
  ylab("RMSE")

sf <-  barabasi.game(n = 20, m = 2, directed = FALSE)
edge_density(sf)
plot(sf)

################# Density 0.1 with scale-free Network Analysis #################

#probmat <- matrix(c(.7,.2,.2,.7),nrow = 2) #density 0.4
probmat <- matrix(c(.2,.02,.02,.2),nrow = 2)
rmse_mat <- matrix(ncol = 5,nrow = 20)
colnames(rmse_mat) <- c("tree&net_sf", "sf&net_sf", "sbm&net_sf", "hub&net_sf", "er&net_sf")

for (i in 1:20){
  print(i)
  set.seed(seeds[i])
  
  # Tree network
  tr <- make_tree(n = 20, children = 3, mode = 'undirected')
  net_tr <- igraphtoGNAR(tr) 
  
  #Scale-free network
  sf <-  barabasi.game(n = 20, m = 2, directed = FALSE)
  net_sf <- igraphtoGNAR(sf) 
  
  # Sbm network
  sbm <- sample_sbm(sum(20), pref.matrix = probmat, block.sizes = c(10,10))
  net_sbm <- igraphtoGNAR(sbm)
  
  # Hub network
  hub <- make_star(20, mode = "undirected")
  net_hub <- igraphtoGNAR(hub)
  
  # Random network
  er <- erdos.renyi.game(20,p.or.m = 38,type = "gnm",directed = FALSE) 
  net_er <- igraphtoGNAR(er)
  
  # Simulate network data based on the GNAR model
  # Normalize the data
  
  ts_tr <- GNARsim(n = 200, net=net_tr, alphaParams = list(rep(0.2,20), rep(0.4,20), rep(-0.6,20)), betaParams = list(c(0.3,.1),c(0.1,.1),c(-0.2,.3)))
  #tsn_tr <- apply(ts_tr, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_tr <- ts_tr[1:199, ]
  
  ts_sf <- GNARsim(n = 200, net = net_sf, alphaParams = list(rep(0.2,20), rep(0.4,20), rep(-0.6,20)), betaParams = list(c(0.3,.1),c(0.1,.1),c(-0.2,.3)))
  #tsn_sf <- apply(ts_sf, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_sf <- ts_sf[1:199, ]
  
  ts_sbm <- GNARsim(n = 200, net = net_sbm, alphaParams = list(rep(0.2,20), rep(0.4,20), rep(-0.6,20)), betaParams = list(c(0.3,.1),c(0.1,.1),c(-0.2,.3)))
  #tsn_sbm <- apply(ts_sbm, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_sbm <- ts_sbm[1:199, ]
  
  ts_hub <- GNARsim(n = 200, net = net_hub, alphaParams = list(rep(0.2,20), rep(0.4,20), rep(-0.6,20)), betaParams = list(c(0.3,.1),c(0.1,.1),c(-0.2,.3)))
  #tsn_hub <- apply(ts_hub, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_hub <- ts_hub[1:199, ]
  
  ts_er <- GNARsim(n = 200, net = net_er, alphaParams = list(rep(0.2,20), rep(0.4,20), rep(-0.6,20)), betaParams = list(c(0.3,.1),c(0.1,.1),c(-0.2,.3)))
  #tsn_er <- apply(ts_er, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_er <- ts_er[1:199, ]
  
  
  # Fit GNAR(3, [2,2,2]) model with tree data and sf network
  fit_pred1 <- predict(GNARfit(vts = datasim_train_tr, net = net_sf,
                               alphaOrder = 3, betaOrder = rep(2,3),
                               globalalpha = globalalpha))
  rmse_mat[i,paste("tree&net_sf")] <- sqrt(mean((fit_pred1-ts_tr[200,])^2))
  print(c("seed: ",i," gnar nei1 done"))
  
  
  # Fit GNAR(3, [2,2,2]) model with sf data and sf network
  fit_pred2 <- predict(GNARfit(vts = datasim_train_sf, net = net_sf,
                               alphaOrder = 3, betaOrder = rep(2,3),
                               globalalpha = globalalpha))
  rmse_mat[i,paste("sf&net_sf")] <- sqrt(mean((fit_pred2-ts_sf[200,])^2))
  print(c("seed: ",i," gnar nei2 done"))
  
  
  # Fit GNAR(3, [2,2,2]) model with sbm data and sf network
  fit_pred3 <- predict(GNARfit(vts = datasim_train_sbm, net = net_sf,
                               alphaOrder = 3, betaOrder = rep(2,3),
                               globalalpha = globalalpha))
  rmse_mat[i,paste("sbm&net_sf")] <- sqrt(mean((fit_pred3-ts_sbm[200,])^2))
  print(c("seed: ",i," gnar nei3 done"))
  
  
  # Fit GNAR(3, [2,2,2]) model with hub data and sf network
  fit_pred4 <- predict(GNARfit(vts = datasim_train_hub, net = net_sf,
                               alphaOrder = 3, betaOrder = rep(2,3),
                               globalalpha = globalalpha))
  rmse_mat[i,paste("hub&net_sf")] <- sqrt(mean((fit_pred4-ts_hub[200,])^2))
  print(c("seed: ",i," gnar nei4 done"))
  
  # Fit GNAR(3, [2,2,2]) model with random data and sf network
  fit_pred5 <- predict(GNARfit(vts = datasim_train_er, net = net_sf,
                               alphaOrder = 3, betaOrder = rep(2,3),
                               globalalpha = globalalpha))
  rmse_mat[i,paste("er&net_sf")] <- sqrt(mean((fit_pred5-ts_er[200,])^2))
  print(c("seed: ",i," gnar nei5 done"))
  
}


# Convert the matrix to a data frame
rmse_df <- as.data.frame(rmse_mat)
summary(rmse_df)

# Check the column names
colnames(rmse_df) <- c("Tree", "Sacle-free", "Clustered", "Hub", "Random")

# Reshape the data frame to long format
rmse_melt <- melt(rmse_df, variable.name = "Network_Time_seires", value.name = "RMSE")

dev.off()
# Create the box plot
ggplot(rmse_melt, aes(x = Network_Time_seires, y = RMSE, fill = Network_Time_seires)) +
  geom_boxplot() +
  ggtitle("RMSE Distribution for Scale-free Structures") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Network Time series") +
  ylab("RMSE") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


################# Density 0.1 with clustered Network Analysis #################

#probmat <- matrix(c(.7,.2,.2,.7),nrow = 2) #density 0.4
probmat <- matrix(c(.2,.02,.02,.2),nrow = 2)
rmse_mat <- matrix(ncol = 5,nrow = 20)
colnames(rmse_mat) <- c("tree&net_sbm", "sf&net_sbm", "sbm&net_sbm", "hub&net_sbm", "er&net_sbm")

for (i in 1:20){
  print(i)
  set.seed(seeds[i])
  
  # Tree network
  tr <- make_tree(n = 20, children = 3, mode = 'undirected')
  net_tr <- igraphtoGNAR(tr) 
  
  #Scale-free network
  sf <-  barabasi.game(n = 20, m = 2, directed = FALSE)
  net_sf <- igraphtoGNAR(sf) 
  
  # Sbm network
  sbm <- sample_sbm(sum(20), pref.matrix = probmat, block.sizes = c(10,10))
  net_sbm <- igraphtoGNAR(sbm)
  
  # Hub network
  hub <- make_star(20, mode = "undirected")
  net_hub <- igraphtoGNAR(hub)
  
  # Random network
  er <- erdos.renyi.game(20,p.or.m = 38,type = "gnm",directed = FALSE) 
  net_er <- igraphtoGNAR(er)
  
  # Simulate network data based on the GNAR model
  # Normalize the data
  
  ts_tr <- GNARsim(n = 200, net=net_tr, alphaParams = list(rep(0.2,20), rep(0.4,20), rep(-0.6,20)), betaParams = list(c(0.3,.1),c(0.1,.1),c(-0.2,.3)))
  #tsn_tr <- apply(ts_tr, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_tr <- ts_tr[1:199, ]
  
  ts_sf <- GNARsim(n = 200, net = net_sf, alphaParams = list(rep(0.2,20), rep(0.4,20), rep(-0.6,20)), betaParams = list(c(0.3,.1),c(0.1,.1),c(-0.2,.3)))
  #tsn_sf <- apply(ts_sf, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_sf <- ts_sf[1:199, ]
  
  ts_sbm <- GNARsim(n = 200, net = net_sbm, alphaParams = list(rep(0.2,20), rep(0.4,20), rep(-0.6,20)), betaParams = list(c(0.3,.1),c(0.1,.1),c(-0.2,.3)))
  #tsn_sbm <- apply(ts_sbm, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_sbm <- ts_sbm[1:199, ]
  
  ts_hub <- GNARsim(n = 200, net = net_hub, alphaParams = list(rep(0.2,20), rep(0.4,20), rep(-0.6,20)), betaParams = list(c(0.3,.1),c(0.1,.1),c(-0.2,.3)))
  #tsn_hub <- apply(ts_hub, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_hub <- ts_hub[1:199, ]
  
  ts_er <- GNARsim(n = 200, net = net_er, alphaParams = list(rep(0.2,20), rep(0.4,20), rep(-0.6,20)), betaParams = list(c(0.3,.1),c(0.1,.1),c(-0.2,.3)))
  #tsn_er <- apply(ts_er, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_er <- ts_er[1:199, ]
  
  
  # Fit GNAR(3, [2,2,2]) model with tree data and sbm network
  fit_pred1 <- predict(GNARfit(vts = datasim_train_tr, net = net_sbm,
                               alphaOrder = 3, betaOrder = rep(2,3),
                               globalalpha = globalalpha))
  rmse_mat[i,paste("tree&net_sbm")] <- sqrt(mean((fit_pred1-ts_tr[200,])^2))
  print(c("seed: ",i," gnar nei1 done"))
  
  
  # Fit GNAR(3, [2,2,2]) model with sf data and sbm network
  fit_pred2 <- predict(GNARfit(vts = datasim_train_sf, net = net_sbm,
                               alphaOrder = 3, betaOrder = rep(2,3),
                               globalalpha = globalalpha))
  rmse_mat[i,paste("sf&net_sbm")] <- sqrt(mean((fit_pred2-ts_sf[200,])^2))
  print(c("seed: ",i," gnar nei2 done"))
  
  
  # Fit GNAR(3, [2,2,2]) model with sbm data and sbm network
  fit_pred3 <- predict(GNARfit(vts = datasim_train_sbm, net = net_sbm,
                               alphaOrder = 3, betaOrder = rep(2,3),
                               globalalpha = globalalpha))
  rmse_mat[i,paste("sbm&net_sbm")] <- sqrt(mean((fit_pred3-ts_sbm[200,])^2))
  print(c("seed: ",i," gnar nei3 done"))
  
  
  # Fit GNAR(3, [2,2,2]) model with hub data and sbm network
  fit_pred4 <- predict(GNARfit(vts = datasim_train_hub, net = net_sbm,
                               alphaOrder = 3, betaOrder = rep(2,3),
                               globalalpha = globalalpha))
  rmse_mat[i,paste("hub&net_sbm")] <- sqrt(mean((fit_pred4-ts_hub[200,])^2))
  print(c("seed: ",i," gnar nei4 done"))
  
  # Fit GNAR(3, [2,2,2]) model with er data and tr network
  fit_pred5 <- predict(GNARfit(vts = datasim_train_er, net = net_sbm,
                               alphaOrder = 3, betaOrder = rep(2,3),
                               globalalpha = globalalpha))
  rmse_mat[i,paste("er&net_sbm")] <- sqrt(mean((fit_pred5-ts_er[200,])^2))
  print(c("seed: ",i," gnar nei5 done"))
  
}


# Convert the matrix to a data frame
rmse_df <- as.data.frame(rmse_mat)
summary(rmse_df)

# Check the column names
colnames(rmse_df) <- c("Tree", "Scale-free", "Clustered", "Hub", "Random")

# Reshape the data frame to long format
rmse_melt <- melt(rmse_df, variable.name = "Network_Model", value.name = "RMSE")

dev.off()
# Create the box plot
ggplot(rmse_melt, aes(x = Network_Model, y = RMSE, fill = Network_Model)) +
  geom_boxplot() +
  ggtitle("RMSE Distribution for Clustered Network Structures") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Network Time Series") +
  ylab("RMSE") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


################# Density 0.1 with tree Network Analysis #################

#probmat <- matrix(c(.7,.2,.2,.7),nrow = 2) #density 0.4
probmat <- matrix(c(.2,.02,.02,.2),nrow = 2)
rmse_mat <- matrix(ncol = 5,nrow = 20)
colnames(rmse_mat) <- c("tree&net_tr", "sf&net_tr", "sbm&net_tr", "hub&net_tr", "er&net_tr")

for (i in 1:20){
  print(i)
  set.seed(seeds[i])
  
  # Tree network
  tr <- make_tree(n = 20, children = 3, mode = 'undirected')
  net_tr <- igraphtoGNAR(tr) 
  
  #Scale-free network
  sf <-  barabasi.game(n = 20, m = 2, directed = FALSE)
  net_sf <- igraphtoGNAR(sf) 
  
  # Sbm network
  sbm <- sample_sbm(sum(20), pref.matrix = probmat, block.sizes = c(10,10))
  net_sbm <- igraphtoGNAR(sbm)
  
  # Hub network
  hub <- make_star(20, mode = "undirected")
  net_hub <- igraphtoGNAR(hub)
  
  # Random network
  er <- erdos.renyi.game(20,p.or.m = 38,type = "gnm",directed = FALSE) 
  net_er <- igraphtoGNAR(er)
  
  # Simulate network data based on the GNAR model
  # Normalize the data
  
  ts_tr <- GNARsim(n = 200, net=net_tr, alphaParams = list(rep(0.2,20), rep(0.4,20), rep(-0.6,20)), betaParams = list(c(0.3,.1),c(0.1,.1),c(-0.2,.3)))
  #tsn_tr <- apply(ts_tr, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_tr <- ts_tr[1:199, ]
  
  ts_sf <- GNARsim(n = 200, net = net_sf, alphaParams = list(rep(0.2,20), rep(0.4,20), rep(-0.6,20)), betaParams = list(c(0.3,.1),c(0.1,.1),c(-0.2,.3)))
  #tsn_sf <- apply(ts_sf, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_sf <- ts_sf[1:199, ]
  
  ts_sbm <- GNARsim(n = 200, net = net_sbm, alphaParams = list(rep(0.2,20), rep(0.4,20), rep(-0.6,20)), betaParams = list(c(0.3,.1),c(0.1,.1),c(-0.2,.3)))
  #tsn_sbm <- apply(ts_sbm, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_sbm <- ts_sbm[1:199, ]
  
  ts_hub <- GNARsim(n = 200, net = net_hub, alphaParams = list(rep(0.2,20), rep(0.4,20), rep(-0.6,20)), betaParams = list(c(0.3,.1),c(0.1,.1),c(-0.2,.3)))
  #tsn_hub <- apply(ts_hub, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_hub <- ts_hub[1:199, ]
  
  ts_er <- GNARsim(n = 200, net = net_er, alphaParams = list(rep(0.2,20), rep(0.4,20), rep(-0.6,20)), betaParams = list(c(0.3,.1),c(0.1,.1),c(-0.2,.3)))
  #tsn_er <- apply(ts_er, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_er <- ts_er[1:199, ]
  
  
  # Fit GNAR(3, [2,2,2]) model with tree data and tr network
  fit_pred1 <- predict(GNARfit(vts = datasim_train_tr, net = net_tr,
                               alphaOrder = 3, betaOrder = rep(2,3),
                               globalalpha = globalalpha))
  rmse_mat[i,paste("tree&net_tr")] <- sqrt(mean((fit_pred1-ts_tr[200,])^2))
  print(c("seed: ",i," gnar nei1 done"))
  
  
  # Fit GNAR(3, [2,2,2]) model with sf data and tr network
  fit_pred2 <- predict(GNARfit(vts = datasim_train_sf, net = net_tr,
                               alphaOrder = 3, betaOrder = rep(2,3),
                               globalalpha = globalalpha))
  rmse_mat[i,paste("sf&net_tr")] <- sqrt(mean((fit_pred2-ts_sf[200,])^2))
  print(c("seed: ",i," gnar nei2 done"))
  
  
  # Fit GNAR(3, [2,2,2]) model with sbm data and tr network
  fit_pred3 <- predict(GNARfit(vts = datasim_train_sbm, net = net_tr,
                               alphaOrder = 3, betaOrder = rep(2,3),
                               globalalpha = globalalpha))
  rmse_mat[i,paste("sbm&net_tr")] <- sqrt(mean((fit_pred3-ts_sbm[200,])^2))
  print(c("seed: ",i," gnar nei3 done"))
  
  
  # Fit GNAR(3, [2,2,2]) model with hub data and tr network
  fit_pred4 <- predict(GNARfit(vts = datasim_train_hub, net = net_tr,
                               alphaOrder = 3, betaOrder = rep(2,3),
                               globalalpha = globalalpha))
  rmse_mat[i,paste("hub&net_tr")] <- sqrt(mean((fit_pred4-ts_hub[200,])^2))
  print(c("seed: ",i," gnar nei4 done"))
  
  # Fit GNAR(3, [2,2,2]) model with er data and tr network
  fit_pred5 <- predict(GNARfit(vts = datasim_train_er, net = net_tr,
                               alphaOrder = 3, betaOrder = rep(2,3),
                               globalalpha = globalalpha))
  rmse_mat[i,paste("er&net_tr")] <- sqrt(mean((fit_pred5-ts_er[200,])^2))
  print(c("seed: ",i," gnar nei5 done"))
  
}


# Convert the matrix to a data frame
rmse_df <- as.data.frame(rmse_mat)
summary(rmse_df)

# Check the column names
colnames(rmse_df) <- c("tree", "sacle-free", "clustered", "hub", "random")

# Reshape the data frame to long format
rmse_melt <- melt(rmse_df, variable.name = "Network_Time_seires", value.name = "RMSE")

dev.off()
# Create the box plot
ggplot(rmse_melt, aes(x = Network_Time_seires, y = RMSE, fill = Network_Time_seires)) +
  geom_boxplot() +
  ggtitle("RMSE Distribution for Different Network Time series with Tree Structures") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Network Time series") +
  ylab("RMSE") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



################# Density 0.1 with hub Network Analysis #################

#probmat <- matrix(c(.7,.2,.2,.7),nrow = 2) #density 0.4
probmat <- matrix(c(.2,.02,.02,.2),nrow = 2)
rmse_mat <- matrix(ncol = 5,nrow = 20)
colnames(rmse_mat) <- c("tree&net_hub", "sf&net_hub", "sbm&net_hub", "hub&net_hub", "er&net_hub")

for (i in 1:20){
  print(i)
  set.seed(seeds[i])
  
  # Tree network
  tr <- make_tree(n = 20, children = 3, mode = 'undirected')
  net_tr <- igraphtoGNAR(tr) 
  
  #Scale-free network
  sf <-  barabasi.game(n = 20, m = 2, directed = FALSE)
  net_sf <- igraphtoGNAR(sf) 
  
  # Sbm network
  sbm <- sample_sbm(sum(20), pref.matrix = probmat, block.sizes = c(10,10))
  net_sbm <- igraphtoGNAR(sbm)
  
  # Hub network
  hub <- make_star(20, mode = "undirected")
  net_hub <- igraphtoGNAR(hub)
  
  # Random network density 0.1
  er <- erdos.renyi.game(20,p.or.m = 38,type = "gnm",directed = FALSE) 
  net_er <- igraphtoGNAR(er)

  
  # Simulate network data based on the GNAR model
  # Normalize the data
  
  ts_tr <- GNARsim(n = 200, net=net_tr, alphaParams = list(rep(0.2,20), rep(0.4,20), rep(-0.6,20)), betaParams = list(c(0.3,.1),c(0.1,.1),c(-0.2,.3)))
  #tsn_tr <- apply(ts_tr, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_tr <- ts_tr[1:199, ]
  
  ts_sf <- GNARsim(n = 200, net = net_sf, alphaParams = list(rep(0.2,20), rep(0.4,20), rep(-0.6,20)), betaParams = list(c(0.3,.1),c(0.1,.1),c(-0.2,.3)))
  #tsn_sf <- apply(ts_sf, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_sf <- ts_sf[1:199, ]
  
  ts_sbm <- GNARsim(n = 200, net = net_sbm, alphaParams = list(rep(0.2,20), rep(0.4,20), rep(-0.6,20)), betaParams = list(c(0.3,.1),c(0.1,.1),c(-0.2,.3)))
  #tsn_sbm <- apply(ts_sbm, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_sbm <- ts_sbm[1:199, ]
  
  ts_hub <- GNARsim(n = 200, net = net_hub, alphaParams = list(rep(0.2,20), rep(0.4,20), rep(-0.6,20)), betaParams = list(c(0.3,.1),c(0.1,.1),c(-0.2,.3)))
  #tsn_hub <- apply(ts_hub, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_hub <- ts_hub[1:199, ]
  
  ts_er <- GNARsim(n = 200, net = net_er, alphaParams = list(rep(0.2,20), rep(0.4,20), rep(-0.6,20)), betaParams = list(c(0.3,.1),c(0.1,.1),c(-0.2,.3)))
  #tsn_er <- apply(ts_er, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_er <- ts_er[1:199, ]
  
  
  
  # Fit GNAR(3, [2,2,2]) model with tree data and tr network
  fit_pred1 <- predict(GNARfit(vts = datasim_train_tr, net = net_hub,
                               alphaOrder = 3, betaOrder = rep(2,3),
                               globalalpha = globalalpha))
  rmse_mat[i,paste("tree&net_hub")] <- sqrt(mean((fit_pred1-ts_tr[200,])^2))
  print(c("seed: ",i," gnar nei1 done"))
  
  
  # Fit GNAR(3, [2,2,2]) model with sf data and tr network
  fit_pred2 <- predict(GNARfit(vts = datasim_train_sf, net = net_hub,
                               alphaOrder = 3, betaOrder = rep(2,3),
                               globalalpha = globalalpha))
  rmse_mat[i,paste("sf&net_hub")] <- sqrt(mean((fit_pred2-ts_sf[200,])^2))
  print(c("seed: ",i," gnar nei2 done"))
  
  
  # Fit GNAR(3, [2,2,2]) model with sbm data and tr network
  fit_pred3 <- predict(GNARfit(vts = datasim_train_sbm, net = net_hub,
                               alphaOrder = 3, betaOrder = rep(2,3),
                               globalalpha = globalalpha))
  rmse_mat[i,paste("sbm&net_hub")] <- sqrt(mean((fit_pred3-ts_sbm[200,])^2))
  print(c("seed: ",i," gnar nei3 done"))
  
  
  # Fit GNAR(3, [2,2,2]) model with hub data and tr network
  fit_pred4 <- predict(GNARfit(vts = datasim_train_hub, net = net_hub,
                               alphaOrder = 3, betaOrder = rep(2,3),
                               globalalpha = globalalpha))
  rmse_mat[i,paste("hub&net_hub")] <- sqrt(mean((fit_pred4-ts_hub[200,])^2))
  print(c("seed: ",i," gnar nei4 done"))
  
  # Fit GNAR(3, [2,2,2]) model with hub data and er network
  fit_pred5 <- predict(GNARfit(vts = datasim_train_er, net = net_hub,
                               alphaOrder = 3, betaOrder = rep(2,3),
                               globalalpha = globalalpha))
  rmse_mat[i,paste("er&net_hub")] <- sqrt(mean((fit_pred5-ts_er[200,])^2))
  print(c("seed: ",i," gnar nei5 done"))
  
}


# Convert the matrix to a data frame
rmse_df <- as.data.frame(rmse_mat)
summary(rmse_df)

# Check the column names
colnames(rmse_df) <- c("tree_net_hub", "sf_net_hub", "sbm_net_hub", "hub_net_hub", "er_net_hub")

# Reshape the data frame to long format
rmse_melt <- melt(rmse_df, variable.name = "Network_Model", value.name = "RMSE")

dev.off()
# Create the box plot
ggplot(rmse_melt, aes(x = Network_Model, y = RMSE, fill = Network_Model)) +
  geom_boxplot() +
  ggtitle("RMSE Distribution for Different Network Structures") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Network Structure") +
  ylab("RMSE") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

################# Density 0.1 with ramdom Network Analysis #################
# Regimes 4: GNAR(3,[2,2,2]), alpha = (0.2,0.4,-0.6), beta = ((0.3,0.1), (0.1,0.1), (-0.2,0.3))
#probmat <- matrix(c(.7,.2,.2,.7),nrow = 2) #density 0.4
probmat <- matrix(c(.2,.02,.02,.2),nrow = 2)
rmse_mat <- matrix(ncol = 5,nrow = 20)
colnames(rmse_mat) <- c("tree&net_er", "sf&net_er", "sbm&net_er", "hub&net_er", "er&net_er")

for (i in 1:20){
  print(i)
  set.seed(seeds[i])
  
  # Tree network
  tr <- make_tree(n = 20, children = 3, mode = 'undirected')
  net_tr <- igraphtoGNAR(tr) 
  
  #Scale-free network
  sf <-  barabasi.game(n = 20, m = 2, directed = FALSE)
  net_sf <- igraphtoGNAR(sf) 
  
  # Sbm network
  sbm <- sample_sbm(sum(20), pref.matrix = probmat, block.sizes = c(10,10))
  net_sbm <- igraphtoGNAR(sbm)
  
  # Hub network
  hub <- make_star(20, mode = "undirected")
  net_hub <- igraphtoGNAR(hub)
  
  # Random network density 0.1
  er <- erdos.renyi.game(20,p.or.m = 38,type = "gnm",directed = FALSE) 
  net_er <- igraphtoGNAR(er)
  
  # Simulate network data based on the GNAR model
  # Normalize the data
  
  ts_tr <- GNARsim(n = 200, net=net_tr, alphaParams = list(rep(0.2,20), rep(0.4,20), rep(-0.6,20)), betaParams = list(c(0.3,.1),c(0.1,.1),c(-0.2,.3)))
  #tsn_tr <- apply(ts_tr, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_tr <- ts_tr[1:199, ]
  
  ts_sf <- GNARsim(n = 200, net = net_sf, alphaParams = list(rep(0.2,20), rep(0.4,20), rep(-0.6,20)), betaParams = list(c(0.3,.1),c(0.1,.1),c(-0.2,.3)))
  #tsn_sf <- apply(ts_sf, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_sf <- ts_sf[1:199, ]
  
  ts_sbm <- GNARsim(n = 200, net = net_sbm, alphaParams = list(rep(0.2,20), rep(0.4,20), rep(-0.6,20)), betaParams = list(c(0.3,.1),c(0.1,.1),c(-0.2,.3)))
  #tsn_sbm <- apply(ts_sbm, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_sbm <- ts_sbm[1:199, ]
  
  ts_hub <- GNARsim(n = 200, net = net_hub, alphaParams = list(rep(0.2,20), rep(0.4,20), rep(-0.6,20)), betaParams = list(c(0.3,.1),c(0.1,.1),c(-0.2,.3)))
  #tsn_hub <- apply(ts_hub, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_hub <- ts_hub[1:199, ]
  
  ts_er <- GNARsim(n = 200, net = net_er, alphaParams = list(rep(0.2,20), rep(0.4,20), rep(-0.6,20)), betaParams = list(c(0.3,.1),c(0.1,.1),c(-0.2,.3)))
  #tsn_er <- apply(ts_er, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_er <- ts_er[1:199, ]
  
  
  # Fit GNAR(3, [2,2,2]) model with tree data and er network
  fit_pred1 <- predict(GNARfit(vts = datasim_train_tr, net = net_er,
                               alphaOrder = 3, betaOrder = rep(2,3),
                               globalalpha = globalalpha))
  rmse_mat[i,paste("tree&net_er")] <- sqrt(mean((fit_pred1-ts_tr[200,])^2))
  print(c("seed: ",i," gnar nei1 done"))
  
  
  # Fit GNAR(3, [2,2,2]) model with sf data and er network
  fit_pred2 <- predict(GNARfit(vts = datasim_train_sf, net = net_er,
                               alphaOrder = 3, betaOrder = rep(2,3),
                               globalalpha = globalalpha))
  rmse_mat[i,paste("sf&net_er")] <- sqrt(mean((fit_pred2-ts_sf[200,])^2))
  print(c("seed: ",i," gnar nei2 done"))
  
  
  # Fit GNAR(3, [2,2,2]) model with sbm data and er network
  fit_pred3 <- predict(GNARfit(vts = datasim_train_sbm, net = net_er,
                               alphaOrder = 3, betaOrder = rep(2,3),
                               globalalpha = globalalpha))
  rmse_mat[i,paste("sbm&net_er")] <- sqrt(mean((fit_pred3-ts_sbm[200,])^2))
  print(c("seed: ",i," gnar nei3 done"))
  
  
  # Fit GNAR(3, [2,2,2]) model with hub data and er network
  fit_pred4 <- predict(GNARfit(vts = datasim_train_hub, net = net_er,
                               alphaOrder = 3, betaOrder = rep(2,3),
                               globalalpha = globalalpha))
  rmse_mat[i,paste("hub&net_er")] <- sqrt(mean((fit_pred4-ts_hub[200,])^2))
  print(c("seed: ",i," gnar nei4 done"))
  
  # Fit GNAR(3, [2,2,2]) model with hub data and er network
  fit_pred5 <- predict(GNARfit(vts = datasim_train_er, net = net_er,
                               alphaOrder = 3, betaOrder = rep(2,3),
                               globalalpha = globalalpha))
  rmse_mat[i,paste("er&net_er")] <- sqrt(mean((fit_pred5-ts_er[200,])^2))
  print(c("seed: ",i," gnar nei5 done"))
  
}


# Convert the matrix to a data frame
rmse_df <- as.data.frame(rmse_mat)
summary(rmse_df)

# Check the column names
colnames(rmse_df) <- c("tree_net_er", "sf_net_er", "sbm_net_er", "hub_net_er", "er_net_er")

# Reshape the data frame to long format
rmse_melt <- melt(rmse_df, variable.name = "Network_Model", value.name = "RMSE")

dev.off()
# Create the box plot
ggplot(rmse_melt, aes(x = Network_Model, y = RMSE, fill = Network_Model)) +
  geom_boxplot() +
  ggtitle("RMSE Distribution for Different Network Structures") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Network Structure") +
  ylab("RMSE") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


################# Density 0.4 with Random Network Analysis #################

probmat <- matrix(c(.7,.2,.2,.7),nrow = 2) #density 0.4
#probmat <- matrix(c(.2,.02,.02,.2),nrow = 2)
rmse_mat <- matrix(ncol = 3, nrow = 20)
colnames(rmse_mat) <- c("sf&net_er", "sbm&net_er",  "er&net_er")

for (i in 1:20){
  print(i)
  set.seed(seeds[i])
  
  
  #Scale-free network
  sf <-  barabasi.game(n = 20, m = 4, directed = FALSE)
  net_sf <- igraphtoGNAR(sf) 
  
  # Sbm network
  sbm <- sample_sbm(sum(20), pref.matrix = probmat, block.sizes = c(10,10))
  net_sbm <- igraphtoGNAR(sbm)
  
  # Random network density 0.1
  er <- erdos.renyi.game(20,p.or.m = 76,type = "gnm",directed = FALSE) 
  net_er <- igraphtoGNAR(er)
  
  # Simulate network data based on the GNAR model
  # Normalize the data
  
  ts_sf <- GNARsim(n = 200, net = net_sf, alphaParams = list(rep(0.2,20), rep(0.4,20), rep(-0.6,20)), betaParams = list(c(0.3,.1),c(0.1,.1),c(-0.2,.3)))
  tsn_sf <- apply(ts_sf, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_sf <- tsn_sf[1:199, ]
  
  ts_sbm <- GNARsim(n = 200, net = net_sbm, alphaParams = list(rep(0.2,20), rep(0.4,20), rep(-0.6,20)), betaParams = list(c(0.3,.1),c(0.1,.1),c(-0.2,.3)))
  tsn_sbm <- apply(ts_sbm, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_sbm <- tsn_sbm[1:199, ]
  
  ts_er <- GNARsim(n = 200, net = net_er, alphaParams = list(rep(0.2,20), rep(0.4,20), rep(-0.6,20)), betaParams = list(c(0.3,.1),c(0.1,.1),c(-0.2,.3)))
  tsn_er <- apply(ts_er, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_er <- tsn_er[1:199, ]
  
  
  # Fit GNAR(3, [2,2,2]) model with sf data and er network
  fit_pred1 <- predict(GNARfit(vts = datasim_train_sf, net = net_er,
                               alphaOrder = 3, betaOrder = rep(2,3),
                               globalalpha = globalalpha))
  rmse_mat[i,paste("sf&net_er")] <- sqrt(mean((fit_pred1-tsn_sf[200,])^2))
  print(c("seed: ",i," gnar nei1 done"))
  
  
  # Fit GNAR(3, [2,2,2]) model with sbm data and er network
  fit_pred2 <- predict(GNARfit(vts = datasim_train_sbm, net = net_er,
                               alphaOrder = 3, betaOrder = rep(2,3),
                               globalalpha = globalalpha))
  rmse_mat[i,paste("sbm&net_er")] <- sqrt(mean((fit_pred2-tsn_sbm[200,])^2))
  print(c("seed: ",i," gnar nei2 done"))
  

  # Fit GNAR(3, [2,2,2]) model with hub data and er network
  fit_pred3 <- predict(GNARfit(vts = datasim_train_er, net = net_er,
                               alphaOrder = 3, betaOrder = rep(2,3),
                               globalalpha = globalalpha))
  rmse_mat[i,paste("er&net_er")] <- sqrt(mean((fit_pred3-tsn_er[200,])^2))
  print(c("seed: ",i," gnar nei3 done"))
  
}


# Convert the matrix to a data frame
rmse_df <- as.data.frame(rmse_mat)
summary(rmse_df)

# Check the column names
colnames(rmse_df) <- c( "sf_net_er", "sbm_net_er",  "er_net_er")

# Reshape the data frame to long format
rmse_melt <- melt(rmse_df, variable.name = "Network_Model", value.name = "RMSE")

dev.off()
# Create the box plot
ggplot(rmse_melt, aes(x = Network_Model, y = RMSE, fill = Network_Model)) +
  geom_boxplot() +
  ggtitle("RMSE Distribution for Different Network Structures") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Network Structure") +
  ylab("RMSE") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

################# Density 0.4 with Scale-free Network Analysis #################

probmat <- matrix(c(.7,.2,.2,.7),nrow = 2) #density 0.4
#probmat <- matrix(c(.2,.02,.02,.2),nrow = 2)
rmse_mat <- matrix(ncol = 3, nrow = 20)
colnames(rmse_mat) <- c("sf&net_sf", "sbm&net_sf",  "er&net_sf")

for (i in 1:20){
  print(i)
  set.seed(seeds[i])
  
  
  #Scale-free network
  sf <-  barabasi.game(n = 20, m = 4, directed = FALSE)
  net_sf <- igraphtoGNAR(sf) 
  
  # Sbm network
  sbm <- sample_sbm(sum(20), pref.matrix = probmat, block.sizes = c(10,10))
  net_sbm <- igraphtoGNAR(sbm)
  
  # Random network density 0.1
  er <- erdos.renyi.game(20,p.or.m = 76,type = "gnm",directed = FALSE) 
  net_er <- igraphtoGNAR(er)
  
  # Simulate network data based on the GNAR model
  # Normalize the data
  
  ts_sf <- GNARsim(n = 200, net = net_sf, alphaParams = list(rep(0.2,20), rep(0.4,20), rep(-0.6,20)), betaParams = list(c(0.3,.1),c(0.1,.1),c(-0.2,.3)))
  tsn_sf <- apply(ts_sf, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_sf <- tsn_sf[1:199, ]
  
  ts_sbm <- GNARsim(n = 200, net = net_sbm, alphaParams = list(rep(0.2,20), rep(0.4,20), rep(-0.6,20)), betaParams = list(c(0.3,.1),c(0.1,.1),c(-0.2,.3)))
  tsn_sbm <- apply(ts_sbm, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_sbm <- tsn_sbm[1:199, ]
  
  ts_er <- GNARsim(n = 200, net = net_er, alphaParams = list(rep(0.2,20), rep(0.4,20), rep(-0.6,20)), betaParams = list(c(0.3,.1),c(0.1,.1),c(-0.2,.3)))
  tsn_er <- apply(ts_er, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_er <- tsn_er[1:199, ]
  
  
  # Fit GNAR(3, [2,2,2]) model with sf data and sf network
  fit_pred1 <- predict(GNARfit(vts = datasim_train_sf, net = net_sf,
                               alphaOrder = 3, betaOrder = rep(2,3),
                               globalalpha = globalalpha))
  rmse_mat[i,paste("sf&net_sf")] <- sqrt(mean((fit_pred1-tsn_sf[200,])^2))
  print(c("seed: ",i," gnar nei1 done"))
  
  
  # Fit GNAR(3, [2,2,2]) model with sbm data and sf network
  fit_pred2 <- predict(GNARfit(vts = datasim_train_sbm, net = net_sf,
                               alphaOrder = 3, betaOrder = rep(2,3),
                               globalalpha = globalalpha))
  rmse_mat[i,paste("sbm&net_sf")] <- sqrt(mean((fit_pred2-tsn_sbm[200,])^2))
  print(c("seed: ",i," gnar nei2 done"))
  
  
  # Fit GNAR(3, [2,2,2]) model with hub data and sf network
  fit_pred3 <- predict(GNARfit(vts = datasim_train_er, net = net_sf,
                               alphaOrder = 3, betaOrder = rep(2,3),
                               globalalpha = globalalpha))
  rmse_mat[i,paste("er&net_sf")] <- sqrt(mean((fit_pred3-tsn_er[200,])^2))
  print(c("seed: ",i," gnar nei3 done"))
  
}


# Convert the matrix to a data frame
rmse_df <- as.data.frame(rmse_mat)
summary(rmse_df)

# Check the column names
colnames(rmse_df) <- c( "sf_net_sf", "sbm_net_sf",  "er_net_sf")

# Reshape the data frame to long format
rmse_melt <- melt(rmse_df, variable.name = "Network_Model", value.name = "RMSE")

dev.off()
# Create the box plot
ggplot(rmse_melt, aes(x = Network_Model, y = RMSE, fill = Network_Model)) +
  geom_boxplot() +
  ggtitle("RMSE Distribution for Different Network Structures") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Network Structure") +
  ylab("RMSE") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

################# Density 0.4 with Clustered Network Analysis #################
# Regimes 4: GNAR(3,[2,2,2]), alpha = (0.2,0.4,-0.6), beta = ((0.3,0.1), (0.1,0.1), (-0.2,0.3))
probmat <- matrix(c(.7,.2,.2,.7),nrow = 2) #density 0.4
#probmat <- matrix(c(.2,.02,.02,.2),nrow = 2)
rmse_mat <- matrix(ncol = 3, nrow = 20)
colnames(rmse_mat) <- c("sf&net_sbm", "sbm&net_sbm",  "er&net_sbm")

for (i in 1:20){
  print(i)
  set.seed(seeds[i])
  
  
  #Scale-free network
  sf <-  barabasi.game(n = 20, m = 4, directed = FALSE)
  net_sf <- igraphtoGNAR(sf) 
  
  # Sbm network
  sbm <- sample_sbm(sum(20), pref.matrix = probmat, block.sizes = c(10,10))
  net_sbm <- igraphtoGNAR(sbm)
  
  # Random network density 0.1
  er <- erdos.renyi.game(20,p.or.m = 76,type = "gnm",directed = FALSE) 
  net_er <- igraphtoGNAR(er)
  
  # Simulate network data based on the GNAR model
  # Normalize the data
  
  ts_sf <- GNARsim(n = 200, net = net_sf, alphaParams = list(rep(0.2,20), rep(0.4,20), rep(-0.6,20)), betaParams = list(c(0.3,.1),c(0.1,.1),c(-0.2,.3)))
  tsn_sf <- apply(ts_sf, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_sf <- tsn_sf[1:199, ]
  
  ts_sbm <- GNARsim(n = 200, net = net_sbm, alphaParams = list(rep(0.2,20), rep(0.4,20), rep(-0.6,20)), betaParams = list(c(0.3,.1),c(0.1,.1),c(-0.2,.3)))
  tsn_sbm <- apply(ts_sbm, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_sbm <- tsn_sbm[1:199, ]
  
  ts_er <- GNARsim(n = 200, net = net_er, alphaParams = list(rep(0.2,20), rep(0.4,20), rep(-0.6,20)), betaParams = list(c(0.3,.1),c(0.1,.1),c(-0.2,.3)))
  tsn_er <- apply(ts_er, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_er <- tsn_er[1:199, ]
  
  
  # Fit GNAR(3, [2,2,2]) model with sf data and sbm network
  fit_pred1 <- predict(GNARfit(vts = datasim_train_sf, net = net_sbm,
                               alphaOrder = 3, betaOrder = rep(2,3),
                               globalalpha = globalalpha))
  rmse_mat[i,paste("sf&net_sbm")] <- sqrt(mean((fit_pred1-tsn_sf[200,])^2))
  print(c("seed: ",i," gnar nei1 done"))
  
  
  # Fit GNAR(3, [2,2,2]) model with sbm data and sbm network
  fit_pred2 <- predict(GNARfit(vts = datasim_train_sbm, net = net_sbm,
                               alphaOrder = 3, betaOrder = rep(2,3),
                               globalalpha = globalalpha))
  rmse_mat[i,paste("sbm&net_sbm")] <- sqrt(mean((fit_pred2-tsn_sbm[200,])^2))
  print(c("seed: ",i," gnar nei2 done"))
  
  
  # Fit GNAR(3, [2,2,2]) model with hub data and sbm network
  fit_pred3 <- predict(GNARfit(vts = datasim_train_er, net = net_sbm,
                               alphaOrder = 3, betaOrder = rep(2,3),
                               globalalpha = globalalpha))
  rmse_mat[i,paste("er&net_sbm")] <- sqrt(mean((fit_pred3-tsn_er[200,])^2))
  print(c("seed: ",i," gnar nei3 done"))
  
}


# Convert the matrix to a data frame
rmse_df <- as.data.frame(rmse_mat)
summary(rmse_df)

# Check the column names
colnames(rmse_df) <- c( "sf_net_bm", "sbm_net_sbm",  "er_net_sbm")

# Reshape the data frame to long format
rmse_melt <- melt(rmse_df, variable.name = "Network_Model", value.name = "RMSE")

dev.off()
# Create the box plot
ggplot(rmse_melt, aes(x = Network_Model, y = RMSE, fill = Network_Model)) +
  geom_boxplot() +
  ggtitle("RMSE Distribution for Different Network Structures") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Network Structure") +
  ylab("RMSE") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


###############################################################################
################### Different simulation regimes analysis #####################
###############################################################################

################# Density 0.1 with ramdom Network Analysis #################
# Regime 1: GNAR(1,[1]), alpha = 0.2, beta = 0.3

#probmat <- matrix(c(.7,.2,.2,.7),nrow = 2) #density 0.4
probmat <- matrix(c(.2,.02,.02,.2),nrow = 2)
rmse_mat <- matrix(ncol = 5,nrow = 20)
colnames(rmse_mat) <- c("tree&net_er", "sf&net_er", "sbm&net_er", "hub&net_er", "er&net_er")

for (i in 1:20){
  print(i)
  set.seed(seeds[i])
  
  # Tree network
  tr <- make_tree(n = 20, children = 3, mode = 'undirected')
  net_tr <- igraphtoGNAR(tr) 
  
  #Scale-free network
  sf <-  barabasi.game(n = 20, m = 2, directed = FALSE)
  net_sf <- igraphtoGNAR(sf) 
  
  # Sbm network
  sbm <- sample_sbm(sum(20), pref.matrix = probmat, block.sizes = c(10,10))
  net_sbm <- igraphtoGNAR(sbm)
  
  # Hub network
  hub <- make_star(20, mode = "undirected")
  net_hub <- igraphtoGNAR(hub)
  
  # Random network density 0.1
  er <- erdos.renyi.game(20,p.or.m = 38,type = "gnm",directed = FALSE) 
  net_er <- igraphtoGNAR(er)
  
  # Simulate network data based on the GNAR model
  # Normalize the data
  
  ts_tr <- GNARsim(n = 200, net=net_tr, alphaParams = list(rep(0.2,20)), betaParams = list(c(0.3)))
  #tsn_tr <- apply(ts_tr, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_tr <- ts_tr[1:199, ]
  
  ts_sf <- GNARsim(n = 200, net = net_sf, alphaParams = list(rep(0.2,20)), betaParams =list(c(0.3)))
  #tsn_sf <- apply(ts_sf, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_sf <- ts_sf[1:199, ]
  
  ts_sbm <- GNARsim(n = 200, net = net_sbm, alphaParams =  list(rep(0.2,20)), betaParams = list(c(0.3)))
  #tsn_sbm <- apply(ts_sbm, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_sbm <- ts_sbm[1:199, ]
  
  ts_hub <- GNARsim(n = 200, net = net_hub, alphaParams =  list(rep(0.2,20)), betaParams =list(c(0.3)))
  #tsn_hub <- apply(ts_hub, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_hub <- ts_hub[1:199, ]
  
  ts_er <- GNARsim(n = 200, net = net_er, alphaParams =  list(rep(0.2,20)), betaParams = list(c(0.3)))
  #tsn_er <- apply(ts_er, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_er <- ts_er[1:199, ]
  
  
  # Fit GNAR(3, [2,2,2]) model with tree data and er network
  fit_pred1 <- predict(GNARfit(vts = datasim_train_tr, net = net_er,
                               alphaOrder = 1, betaOrder = c(1),
                               globalalpha = globalalpha))
  rmse_mat[i,paste("tree&net_er")] <- sqrt(mean((fit_pred1-ts_tr[200,])^2))
  print(c("seed: ",i," gnar nei1 done"))
  
  
  # Fit GNAR(3, [2,2,2]) model with sf data and er network
  fit_pred2 <- predict(GNARfit(vts = datasim_train_sf, net = net_er,
                               alphaOrder = 1, betaOrder = c(1),
                               globalalpha = globalalpha))
  rmse_mat[i,paste("sf&net_er")] <- sqrt(mean((fit_pred2-ts_sf[200,])^2))
  print(c("seed: ",i," gnar nei2 done"))
  
  
  # Fit GNAR(3, [2,2,2]) model with sbm data and er network
  fit_pred3 <- predict(GNARfit(vts = datasim_train_sbm, net = net_er,
                               alphaOrder = 1, betaOrder = c(1),
                               globalalpha = globalalpha))
  rmse_mat[i,paste("sbm&net_er")] <- sqrt(mean((fit_pred3-ts_sbm[200,])^2))
  print(c("seed: ",i," gnar nei3 done"))
  
  
  # Fit GNAR(3, [2,2,2]) model with hub data and er network
  fit_pred4 <- predict(GNARfit(vts = datasim_train_hub, net = net_er,
                               alphaOrder = 1, betaOrder = c(1),
                               globalalpha = globalalpha))
  rmse_mat[i,paste("hub&net_er")] <- sqrt(mean((fit_pred4-ts_hub[200,])^2))
  print(c("seed: ",i," gnar nei4 done"))
  
  # Fit GNAR(3, [2,2,2]) model with hub data and er network
  fit_pred5 <- predict(GNARfit(vts = datasim_train_er, net = net_er,
                               alphaOrder = 1, betaOrder = c(1),
                               globalalpha = globalalpha))
  rmse_mat[i,paste("er&net_er")] <- sqrt(mean((fit_pred5-ts_er[200,])^2))
  print(c("seed: ",i," gnar nei5 done"))
  
}


# Convert the matrix to a data frame
rmse_df <- as.data.frame(rmse_mat)
summary(rmse_df)

# Check the column names
colnames(rmse_df) <- c("tree_net_er", "sf_net_er", "sbm_net_er", "hub_net_er", "er_net_er")

# Reshape the data frame to long format
rmse_melt <- melt(rmse_df, variable.name = "Network_Model", value.name = "RMSE")

dev.off()
# Create the box plot
ggplot(rmse_melt, aes(x = Network_Model, y = RMSE, fill = Network_Model)) +
  geom_boxplot() +
  ggtitle("RMSE Distribution for Different Network Structures") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Network Structure") +
  ylab("RMSE") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

################## Regime 2: GNAR(1,[2]), alpha = 0.2, beta = (0.3, 0.4) ########

#probmat <- matrix(c(.7,.2,.2,.7),nrow = 2) #density 0.4
probmat <- matrix(c(.2,.02,.02,.2),nrow = 2)
rmse_mat <- matrix(ncol = 5,nrow = 20)
colnames(rmse_mat) <- c("tree&net_er", "sf&net_er", "sbm&net_er", "hub&net_er", "er&net_er")

for (i in 1:20){
  print(i)
  set.seed(seeds[i])
  
  # Tree network
  tr <- make_tree(n = 20, children = 3, mode = 'undirected')
  net_tr <- igraphtoGNAR(tr) 
  
  #Scale-free network
  sf <-  barabasi.game(n = 20, m = 2, directed = FALSE)
  net_sf <- igraphtoGNAR(sf) 
  
  # Sbm network
  sbm <- sample_sbm(sum(20), pref.matrix = probmat, block.sizes = c(10,10))
  net_sbm <- igraphtoGNAR(sbm)
  
  # Hub network
  hub <- make_star(20, mode = "undirected")
  net_hub <- igraphtoGNAR(hub)
  
  # Random network density 0.1
  er <- erdos.renyi.game(20,p.or.m = 38,type = "gnm",directed = FALSE) 
  net_er <- igraphtoGNAR(er)
  
  # Simulate network data based on the GNAR model
  # Normalize the data
  
  ts_tr <- GNARsim(n = 200, net=net_tr, alphaParams = list(rep(0.2,20)), betaParams = list(c(0.3,0.4)))
  tsn_tr <- apply(ts_tr, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_tr <- tsn_tr[1:199, ]
  
  ts_sf <- GNARsim(n = 200, net = net_sf, alphaParams = list(rep(0.2,20)), betaParams =list(c(0.3,0.4)))
  tsn_sf <- apply(ts_sf, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_sf <- tsn_sf[1:199, ]
  
  ts_sbm <- GNARsim(n = 200, net = net_sbm, alphaParams =  list(rep(0.2,20)), betaParams = list(c(0.3,0.4)))
  tsn_sbm <- apply(ts_sbm, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_sbm <- tsn_sbm[1:199, ]
  
  ts_hub <- GNARsim(n = 200, net = net_hub, alphaParams =  list(rep(0.2,20)), betaParams =list(c(0.3, 0.4)))
  tsn_hub <- apply(ts_hub, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_hub <- tsn_hub[1:199, ]
  
  ts_er <- GNARsim(n = 200, net = net_er, alphaParams =  list(rep(0.2,20)), betaParams = list(c(0.3,0.4)))
  tsn_er <- apply(ts_er, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_er <- tsn_er[1:199, ]
  
  
  # Fit GNAR(3, [2,2,2]) model with tree data and er network
  fit_pred1 <- predict(GNARfit(vts = datasim_train_tr, net = net_er,
                               alphaOrder = 1, betaOrder = c(2),
                               globalalpha = globalalpha))
  rmse_mat[i,paste("tree&net_er")] <- sqrt(mean((fit_pred1-tsn_tr[200,])^2))
  print(c("seed: ",i," gnar nei1 done"))
  
  
  # Fit GNAR(3, [2,2,2]) model with sf data and er network
  fit_pred2 <- predict(GNARfit(vts = datasim_train_sf, net = net_er,
                               alphaOrder = 1, betaOrder = c(2),
                               globalalpha = globalalpha))
  rmse_mat[i,paste("sf&net_er")] <- sqrt(mean((fit_pred2-tsn_sf[200,])^2))
  print(c("seed: ",i," gnar nei2 done"))
  
  
  # Fit GNAR(3, [2,2,2]) model with sbm data and er network
  fit_pred3 <- predict(GNARfit(vts = datasim_train_sbm, net = net_er,
                               alphaOrder = 1, betaOrder = c(2),
                               globalalpha = globalalpha))
  rmse_mat[i,paste("sbm&net_er")] <- sqrt(mean((fit_pred3-tsn_sbm[200,])^2))
  print(c("seed: ",i," gnar nei3 done"))
  
  
  # Fit GNAR(3, [2,2,2]) model with hub data and er network
  fit_pred4 <- predict(GNARfit(vts = datasim_train_hub, net = net_er,
                               alphaOrder = 1, betaOrder = c(2),
                               globalalpha = globalalpha))
  rmse_mat[i,paste("hub&net_er")] <- sqrt(mean((fit_pred4-tsn_hub[200,])^2))
  print(c("seed: ",i," gnar nei4 done"))
  
  # Fit GNAR(3, [2,2,2]) model with hub data and er network
  fit_pred5 <- predict(GNARfit(vts = datasim_train_er, net = net_er,
                               alphaOrder = 1, betaOrder = c(2),
                               globalalpha = globalalpha))
  rmse_mat[i,paste("er&net_er")] <- sqrt(mean((fit_pred5-tsn_er[200,])^2))
  print(c("seed: ",i," gnar nei5 done"))
  
}


# Convert the matrix to a data frame
rmse_df <- as.data.frame(rmse_mat)
summary(rmse_df)

# Check the column names
colnames(rmse_df) <- c("tree_net_er", "sf_net_er", "sbm_net_er", "hub_net_er", "er_net_er")

# Reshape the data frame to long format
rmse_melt <- melt(rmse_df, variable.name = "Network_Model", value.name = "RMSE")

dev.off()
# Create the box plot
ggplot(rmse_melt, aes(x = Network_Model, y = RMSE, fill = Network_Model)) +
  geom_boxplot() +
  ggtitle("RMSE Distribution for Different Network Structures") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Network Structure") +
  ylab("RMSE") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

################## Regime 3: GNAR(3,[1,1,1]), alpha = (0.2, 0.4,-0.6), beta = (0.2, 0.1, -0.2) ########

#probmat <- matrix(c(.7,.2,.2,.7),nrow = 2) #density 0.4
probmat <- matrix(c(.2,.02,.02,.2),nrow = 2)
rmse_mat <- matrix(ncol = 5,nrow = 20)
colnames(rmse_mat) <- c("tree&net_er", "sf&net_er", "sbm&net_er", "hub&net_er", "er&net_er")

for (i in 1:20){
  print(i)
  set.seed(seeds[i])
  
  # Tree network
  tr <- make_tree(n = 20, children = 3, mode = 'undirected')
  net_tr <- igraphtoGNAR(tr) 
  
  #Scale-free network
  sf <-  barabasi.game(n = 20, m = 2, directed = FALSE)
  net_sf <- igraphtoGNAR(sf) 
  
  # Sbm network
  sbm <- sample_sbm(sum(20), pref.matrix = probmat, block.sizes = c(10,10))
  net_sbm <- igraphtoGNAR(sbm)
  
  # Hub network
  hub <- make_star(20, mode = "undirected")
  net_hub <- igraphtoGNAR(hub)
  
  # Random network density 0.1
  er <- erdos.renyi.game(20,p.or.m = 38,type = "gnm",directed = FALSE) 
  net_er <- igraphtoGNAR(er)
  
  # Simulate network data based on the GNAR model
  # Normalize the data
  
  ts_tr <- GNARsim(n = 200, net=net_tr, alphaParams = list(rep(0.2,20), rep(0.4,20), rep(-0.6,20)), betaParams = list(c(0.2),c(0.1),c(-.2)))
  ts_tr <- apply(ts_tr, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_tr <- ts_tr[1:199, ]
  
  ts_sf <- GNARsim(n = 200, net = net_sf, alphaParams = list(rep(0.2,20), rep(0.4,20), rep(-0.6,20)), betaParams =list(c(0.2),c(0.1),c(-.2)))
  #tsn_sf <- apply(ts_sf, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_sf <- ts_sf[1:199, ]
  
  ts_sbm <- GNARsim(n = 200, net = net_sbm, alphaParams =  list(rep(0.2,20), rep(0.4,20), rep(-0.6,20)), betaParams = list(c(0.2),c(0.1),c(-.2)))
  #tsn_sbm <- apply(ts_sbm, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_sbm <- ts_sbm[1:199, ]
  
  ts_hub <- GNARsim(n = 200, net = net_hub, alphaParams =  list(rep(0.2,20), rep(0.4,20), rep(-0.6,20)), betaParams =list(c(0.2),c(0.1),c(-.2)))
  #tsn_hub <- apply(ts_hub, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_hub <- ts_hub[1:199, ]
  
  ts_er <- GNARsim(n = 200, net = net_er, alphaParams =  list(rep(0.2,20), rep(0.4,20), rep(-0.6,20)), betaParams = list(c(0.2),c(0.1),c(-.2)))
  #tsn_er <- apply(ts_er, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_er <- ts_er[1:199, ]
  
  
  # Fit GNAR(3, [2,2,2]) model with tree data and er network
  fit_pred1 <- predict(GNARfit(vts = datasim_train_tr, net = net_er,
                               alphaOrder = 3, betaOrder = rep(1,3),
                               globalalpha = globalalpha))
  rmse_mat[i,paste("tree&net_er")] <- sqrt(mean((fit_pred1-ts_tr[200,])^2))
  print(c("seed: ",i," gnar nei1 done"))
  
  
  # Fit GNAR(3, [2,2,2]) model with sf data and er network
  fit_pred2 <- predict(GNARfit(vts = datasim_train_sf, net = net_er,
                               alphaOrder = 3, betaOrder = rep(1,3),
                               globalalpha = globalalpha))
  rmse_mat[i,paste("sf&net_er")] <- sqrt(mean((fit_pred2-ts_sf[200,])^2))
  print(c("seed: ",i," gnar nei2 done"))
  
  
  # Fit GNAR(3, [2,2,2]) model with sbm data and er network
  fit_pred3 <- predict(GNARfit(vts = datasim_train_sbm, net = net_er,
                               alphaOrder = 3, betaOrder = rep(1,3),
                               globalalpha = globalalpha))
  rmse_mat[i,paste("sbm&net_er")] <- sqrt(mean((fit_pred3-ts_sbm[200,])^2))
  print(c("seed: ",i," gnar nei3 done"))
  
  
  # Fit GNAR(3, [2,2,2]) model with hub data and er network
  fit_pred4 <- predict(GNARfit(vts = datasim_train_hub, net = net_er,
                               alphaOrder = 3, betaOrder = rep(1,3),
                               globalalpha = globalalpha))
  rmse_mat[i,paste("hub&net_er")] <- sqrt(mean((fit_pred4-ts_hub[200,])^2))
  print(c("seed: ",i," gnar nei4 done"))
  
  # Fit GNAR(3, [2,2,2]) model with hub data and er network
  fit_pred5 <- predict(GNARfit(vts = datasim_train_er, net = net_er,
                               alphaOrder = 3, betaOrder = rep(1,3),
                               globalalpha = globalalpha))
  rmse_mat[i,paste("er&net_er")] <- sqrt(mean((fit_pred5-ts_er[200,])^2))
  print(c("seed: ",i," gnar nei5 done"))
  
}


# Convert the matrix to a data frame
rmse_df <- as.data.frame(rmse_mat)
summary(rmse_df)

# Check the column names
colnames(rmse_df) <- c("tree_net_er", "sf_net_er", "sbm_net_er", "hub_net_er", "er_net_er")

# Reshape the data frame to long format
rmse_melt <- melt(rmse_df, variable.name = "Network_Model", value.name = "RMSE")

dev.off()
# Create the box plot
ggplot(rmse_melt, aes(x = Network_Model, y = RMSE, fill = Network_Model)) +
  geom_boxplot() +
  ggtitle("RMSE Distribution for Different Network Structures") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Network Structure") +
  ylab("RMSE") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


################## Regime 5: GNAR(3,[2,0,0]), alpha = (0.2, 0.4,-0.6), beta = ((0.3, 0.4), 0,0) ########
# random network
#probmat <- matrix(c(.7,.2,.2,.7),nrow = 2) #density 0.4
probmat <- matrix(c(.2,.02,.02,.2),nrow = 2)
rmse_mat <- matrix(ncol = 5,nrow = 20)
colnames(rmse_mat) <- c("tree&net_er", "sf&net_er", "sbm&net_er", "hub&net_er", "er&net_er")

for (i in 1:20){
  print(i)
  set.seed(seeds[i])
  
  # Tree network
  tr <- make_tree(n = 20, children = 3, mode = 'undirected')
  net_tr <- igraphtoGNAR(tr) 
  
  #Scale-free network
  sf <-  barabasi.game(n = 20, m = 2, directed = FALSE)
  net_sf <- igraphtoGNAR(sf) 
  
  # Sbm network
  sbm <- sample_sbm(sum(20), pref.matrix = probmat, block.sizes = c(10,10))
  net_sbm <- igraphtoGNAR(sbm)
  
  # Hub network
  hub <- make_star(20, mode = "undirected")
  net_hub <- igraphtoGNAR(hub)
  
  # Random network density 0.1
  er <- erdos.renyi.game(20,p.or.m = 38,type = "gnm",directed = FALSE) 
  net_er <- igraphtoGNAR(er)
  
  # Simulate network data based on the GNAR model
  # Normalize the data
  
  ts_tr <- GNARsim(n = 200, net=net_tr, alphaParams = list(rep(0.2,20), rep(0.4,20), rep(-0.6,20)), betaParams = list(c(0.3,.4),c(0),c(0)))
  tsn_tr <- apply(ts_tr, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_tr <- tsn_tr[1:199, ]
  
  ts_sf <- GNARsim(n = 200, net = net_sf, alphaParams = list(rep(0.2,20), rep(0.4,20), rep(-0.6,20)), betaParams =list(c(0.3,.4),c(0),c(0)))
  tsn_sf <- apply(ts_sf, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_sf <- tsn_sf[1:199, ]
  
  ts_sbm <- GNARsim(n = 200, net = net_sbm, alphaParams =  list(rep(0.2,20), rep(0.4,20), rep(-0.6,20)), betaParams = list(c(0.3,.4),c(0),c(0)))
  tsn_sbm <- apply(ts_sbm, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_sbm <- tsn_sbm[1:199, ]
  
  ts_hub <- GNARsim(n = 200, net = net_hub, alphaParams =  list(rep(0.2,20), rep(0.4,20), rep(-0.6,20)), betaParams =list(c(0.3,.4),c(0),c(0)))
  tsn_hub <- apply(ts_hub, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_hub <- tsn_hub[1:199, ]
  
  ts_er <- GNARsim(n = 200, net = net_er, alphaParams =  list(rep(0.2,20), rep(0.4,20), rep(-0.6,20)), betaParams = list(c(0.3,.4),c(0),c(0)))
  tsn_er <- apply(ts_er, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_er <- tsn_er[1:199, ]
  
  
  # Fit GNAR(3, [2,2,2]) model with tree data and er network
  fit_pred1 <- predict(GNARfit(vts = datasim_train_tr, net = net_er,
                               alphaOrder = 3, betaOrder = c(2,0,0),
                               globalalpha = globalalpha))
  rmse_mat[i,paste("tree&net_er")] <- sqrt(mean((fit_pred1-tsn_tr[200,])^2))
  print(c("seed: ",i," gnar nei1 done"))
  
  
  # Fit GNAR(3, [2,2,2]) model with sf data and er network
  fit_pred2 <- predict(GNARfit(vts = datasim_train_sf, net = net_er,
                               alphaOrder = 3, betaOrder = c(2,0,0),
                               globalalpha = globalalpha))
  rmse_mat[i,paste("sf&net_er")] <- sqrt(mean((fit_pred2-tsn_sf[200,])^2))
  print(c("seed: ",i," gnar nei2 done"))
  
  
  # Fit GNAR(3, [2,2,2]) model with sbm data and er network
  fit_pred3 <- predict(GNARfit(vts = datasim_train_sbm, net = net_er,
                               alphaOrder = 3, betaOrder = c(2,0,0),
                               globalalpha = globalalpha))
  rmse_mat[i,paste("sbm&net_er")] <- sqrt(mean((fit_pred3-tsn_sbm[200,])^2))
  print(c("seed: ",i," gnar nei3 done"))
  
  
  # Fit GNAR(3, [2,2,2]) model with hub data and er network
  fit_pred4 <- predict(GNARfit(vts = datasim_train_hub, net = net_er,
                               alphaOrder = 3, betaOrder = c(2,0,0),
                               globalalpha = globalalpha))
  rmse_mat[i,paste("hub&net_er")] <- sqrt(mean((fit_pred4-tsn_hub[200,])^2))
  print(c("seed: ",i," gnar nei4 done"))
  
  # Fit GNAR(3, [2,2,2]) model with hub data and er network
  fit_pred5 <- predict(GNARfit(vts = datasim_train_er, net = net_er,
                               alphaOrder = 3, betaOrder = c(2,0,0),
                               globalalpha = globalalpha))
  rmse_mat[i,paste("er&net_er")] <- sqrt(mean((fit_pred5-tsn_er[200,])^2))
  print(c("seed: ",i," gnar nei5 done"))
  
}


# Convert the matrix to a data frame
rmse_df <- as.data.frame(rmse_mat)
summary(rmse_df)

# Check the column names
colnames(rmse_df) <- c("tree_net_er", "sf_net_er", "sbm_net_er", "hub_net_er", "er_net_er")

# Reshape the data frame to long format
rmse_melt <- melt(rmse_df, variable.name = "Network_Model", value.name = "RMSE")

dev.off()
# Create the box plot
ggplot(rmse_melt, aes(x = Network_Model, y = RMSE, fill = Network_Model)) +
  geom_boxplot() +
  ggtitle("RMSE Distribution for Different Network Structures") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Network Structure") +
  ylab("RMSE") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

############################### Regime 5 with clustered network #################
probmat <- matrix(c(.2,.02,.02,.2),nrow = 2)
rmse_mat <- matrix(ncol = 5,nrow = 20)
colnames(rmse_mat) <- c("tree&net_sbm", "sf&net_sbm", "sbm&net_sbm", "hub&net_sbm", "er&net_sbm")

for (i in 1:20){
  print(i)
  set.seed(seeds[i])
  
  # Tree network
  tr <- make_tree(n = 20, children = 3, mode = 'undirected')
  net_tr <- igraphtoGNAR(tr) 
  
  #Scale-free network
  sf <-  barabasi.game(n = 20, m = 2, directed = FALSE)
  net_sf <- igraphtoGNAR(sf) 
  
  # Sbm network
  sbm <- sample_sbm(sum(20), pref.matrix = probmat, block.sizes = c(10,10))
  net_sbm <- igraphtoGNAR(sbm)
  
  # Hub network
  hub <- make_star(20, mode = "undirected")
  net_hub <- igraphtoGNAR(hub)
  
  # Random network density 0.1
  er <- erdos.renyi.game(20,p.or.m = 38,type = "gnm",directed = FALSE) 
  net_er <- igraphtoGNAR(er)
  
  # Simulate network data based on the GNAR model
  # Normalize the data
  
  ts_tr <- GNARsim(n = 200, net=net_tr, alphaParams = list(rep(0.2,20), rep(0.4,20), rep(-0.6,20)), betaParams = list(c(0.3,.4),c(0),c(0)))
  #tsn_tr <- apply(ts_tr, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_tr <- ts_tr[1:199, ]
  
  ts_sf <- GNARsim(n = 200, net = net_sf, alphaParams = list(rep(0.2,20), rep(0.4,20), rep(-0.6,20)), betaParams =list(c(0.3,.4),c(0),c(0)))
  #tsn_sf <- apply(ts_sf, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_sf <- ts_sf[1:199, ]
  
  ts_sbm <- GNARsim(n = 200, net = net_sbm, alphaParams =  list(rep(0.2,20), rep(0.4,20), rep(-0.6,20)), betaParams = list(c(0.3,.4),c(0),c(0)))
  #tsn_sbm <- apply(ts_sbm, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_sbm <- ts_sbm[1:199, ]
  
  ts_hub <- GNARsim(n = 200, net = net_hub, alphaParams =  list(rep(0.2,20), rep(0.4,20), rep(-0.6,20)), betaParams =list(c(0.3,.4),c(0),c(0)))
  #tsn_hub <- apply(ts_hub, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_hub <- ts_hub[1:199, ]
  
  ts_er <- GNARsim(n = 200, net = net_er, alphaParams =  list(rep(0.2,20), rep(0.4,20), rep(-0.6,20)), betaParams = list(c(0.3,.4),c(0),c(0)))
  #tsn_er <- apply(ts_er, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_er <- ts_er[1:199, ]
  
  
  # Fit GNAR(3, [2,2,2]) model with tree data and sbm network
  fit_pred1 <- predict(GNARfit(vts = datasim_train_tr, net = net_sbm,
                               alphaOrder = 3, betaOrder = c(2,0,0),
                               globalalpha = globalalpha))
  rmse_mat[i,paste("tree&net_sbm")] <- sqrt(mean((fit_pred1-ts_tr[200,])^2))
  print(c("seed: ",i," gnar nei1 done"))
  
  
  # Fit GNAR(3, [2,2,2]) model with sf data and sbm network
  fit_pred2 <- predict(GNARfit(vts = datasim_train_sf, net = net_sbm,
                               alphaOrder = 3, betaOrder = c(2,0,0),
                               globalalpha = globalalpha))
  rmse_mat[i,paste("sf&net_sbm")] <- sqrt(mean((fit_pred2-ts_sf[200,])^2))
  print(c("seed: ",i," gnar nei2 done"))
  
  
  # Fit GNAR(3, [2,2,2]) model with sbm data and sbm network
  fit_pred3 <- predict(GNARfit(vts = datasim_train_sbm, net = net_sbm,
                               alphaOrder = 3, betaOrder = c(2,0,0),
                               globalalpha = globalalpha))
  rmse_mat[i,paste("sbm&net_sbm")] <- sqrt(mean((fit_pred3-ts_sbm[200,])^2))
  print(c("seed: ",i," gnar nei3 done"))
  
  
  # Fit GNAR(3, [2,2,2]) model with hub data and sbm network
  fit_pred4 <- predict(GNARfit(vts = datasim_train_hub, net = net_sbm,
                               alphaOrder = 3, betaOrder = c(2,0,0),
                               globalalpha = globalalpha))
  rmse_mat[i,paste("hub&net_sbm")] <- sqrt(mean((fit_pred4-ts_hub[200,])^2))
  print(c("seed: ",i," gnar nei4 done"))
  
  # Fit GNAR(3, [2,2,2]) model with hub data and sbm network
  fit_pred5 <- predict(GNARfit(vts = datasim_train_er, net = net_sbm,
                               alphaOrder = 3, betaOrder = c(2,0,0),
                               globalalpha = globalalpha))
  rmse_mat[i,paste("er&net_sbm")] <- sqrt(mean((fit_pred5-ts_er[200,])^2))
  print(c("seed: ",i," gnar nei5 done"))
  
}


# Convert the matrix to a data frame
rmse_df <- as.data.frame(rmse_mat)
summary(rmse_df)

# Check the column names
colnames(rmse_df) <- c("tree_net_sbm", "sf_net_sbm", "sbm_net_sbm", "hub_net_sbm", "er_net_sbm")

# Reshape the data frame to long format
rmse_melt <- melt(rmse_df, variable.name = "Network_Model", value.name = "RMSE")

dev.off()
# Create the box plot
ggplot(rmse_melt, aes(x = Network_Model, y = RMSE, fill = Network_Model)) +
  geom_boxplot() +
  ggtitle("RMSE Distribution for Different Network Structures") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Network Structure") +
  ylab("RMSE") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


############################### Regime 5 with scale-free network #################
probmat <- matrix(c(.2,.02,.02,.2),nrow = 2)
rmse_mat <- matrix(ncol = 5,nrow = 20)
colnames(rmse_mat) <- c("tree&net_sf", "sf&net_sf", "sbm&net_sf", "hub&net_sf", "er&net_sf")

for (i in 1:20){
  print(i)
  set.seed(seeds[i])
  
  # Tree network
  tr <- make_tree(n = 20, children = 3, mode = 'undirected')
  net_tr <- igraphtoGNAR(tr) 
  
  #Scale-free network
  sf <-  barabasi.game(n = 20, m = 2, directed = FALSE)
  net_sf <- igraphtoGNAR(sf) 
  
  # Sbm network
  sbm <- sample_sbm(sum(20), pref.matrix = probmat, block.sizes = c(10,10))
  net_sbm <- igraphtoGNAR(sbm)
  
  # Hub network
  hub <- make_star(20, mode = "undirected")
  net_hub <- igraphtoGNAR(hub)
  
  # Random network density 0.1
  er <- erdos.renyi.game(20,p.or.m = 38,type = "gnm",directed = FALSE) 
  net_er <- igraphtoGNAR(er)
  
  # Simulate network data based on the GNAR model
  # Normalize the data
  
  ts_tr <- GNARsim(n = 200, net=net_tr, alphaParams = list(rep(0.2,20), rep(0.4,20), rep(-0.6,20)), betaParams = list(c(0.3,.4),c(0),c(0)))
  #tsn_tr <- apply(ts_tr, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_tr <- ts_tr[1:199, ]
  
  ts_sf <- GNARsim(n = 200, net = net_sf, alphaParams = list(rep(0.2,20), rep(0.4,20), rep(-0.6,20)), betaParams =list(c(0.3,.4),c(0),c(0)))
  #tsn_sf <- apply(ts_sf, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_sf <- ts_sf[1:199, ]
  
  ts_sbm <- GNARsim(n = 200, net = net_sbm, alphaParams =  list(rep(0.2,20), rep(0.4,20), rep(-0.6,20)), betaParams = list(c(0.3,.4),c(0),c(0)))
  #tsn_sbm <- apply(ts_sbm, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_sbm <- ts_sbm[1:199, ]
  
  ts_hub <- GNARsim(n = 200, net = net_hub, alphaParams =  list(rep(0.2,20), rep(0.4,20), rep(-0.6,20)), betaParams =list(c(0.3,.4),c(0),c(0)))
  #tsn_hub <- apply(ts_hub, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_hub <- ts_hub[1:199, ]
  
  ts_er <- GNARsim(n = 200, net = net_er, alphaParams =  list(rep(0.2,20), rep(0.4,20), rep(-0.6,20)), betaParams = list(c(0.3,.4),c(0),c(0)))
  #tsn_er <- apply(ts_er, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_er <- ts_er[1:199, ]
  
  
  # Fit GNAR(3, [2,2,2]) model with tree data and sf network
  fit_pred1 <- predict(GNARfit(vts = datasim_train_tr, net = net_sf,
                               alphaOrder = 3, betaOrder = c(2,0,0),
                               globalalpha = globalalpha))
  rmse_mat[i,paste("tree&net_sf")] <- sqrt(mean((fit_pred1-ts_tr[200,])^2))
  print(c("seed: ",i," gnar nei1 done"))
  
  
  # Fit GNAR(3, [2,2,2]) model with sf data and sf network
  fit_pred2 <- predict(GNARfit(vts = datasim_train_sf, net = net_sf,
                               alphaOrder = 3, betaOrder = c(2,0,0),
                               globalalpha = globalalpha))
  rmse_mat[i,paste("sf&net_sf")] <- sqrt(mean((fit_pred2-ts_sf[200,])^2))
  print(c("seed: ",i," gnar nei2 done"))
  
  
  # Fit GNAR(3, [2,2,2]) model with sbm data and sf network
  fit_pred3 <- predict(GNARfit(vts = datasim_train_sbm, net = net_sf,
                               alphaOrder = 3, betaOrder = c(2,0,0),
                               globalalpha = globalalpha))
  rmse_mat[i,paste("sbm&net_sf")] <- sqrt(mean((fit_pred3-ts_sbm[200,])^2))
  print(c("seed: ",i," gnar nei3 done"))
  
  
  # Fit GNAR(3, [2,2,2]) model with hub data and sf network
  fit_pred4 <- predict(GNARfit(vts = datasim_train_hub, net = net_sf,
                               alphaOrder = 3, betaOrder = c(2,0,0),
                               globalalpha = globalalpha))
  rmse_mat[i,paste("hub&net_sf")] <- sqrt(mean((fit_pred4-ts_hub[200,])^2))
  print(c("seed: ",i," gnar nei4 done"))
  
  # Fit GNAR(3, [2,2,2]) model with hub data and sf network
  fit_pred5 <- predict(GNARfit(vts = datasim_train_er, net = net_sf,
                               alphaOrder = 3, betaOrder = c(2,0,0),
                               globalalpha = globalalpha))
  rmse_mat[i,paste("er&net_sf")] <- sqrt(mean((fit_pred5-ts_er[200,])^2))
  print(c("seed: ",i," gnar nei5 done"))
  
}


# Convert the matrix to a data frame
rmse_df <- as.data.frame(rmse_mat)
summary(rmse_df)

# Check the column names
colnames(rmse_df) <- c("tree_net_sf", "sf_net_sf", "sbm_net_sf", "hub_net_sf", "er_net_sf")

# Reshape the data frame to long format
rmse_melt <- melt(rmse_df, variable.name = "Network_Model", value.name = "RMSE")

dev.off()
# Create the box plot
ggplot(rmse_melt, aes(x = Network_Model, y = RMSE, fill = Network_Model)) +
  geom_boxplot() +
  ggtitle("RMSE Distribution for Different Network Structures") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Network Structure") +
  ylab("RMSE") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

############################### Regime 5 with tree network #################
probmat <- matrix(c(.2,.02,.02,.2),nrow = 2)
rmse_mat <- matrix(ncol = 5,nrow = 20)
colnames(rmse_mat) <- c("tree&net_tr", "sf&net_tr", "sbm&net_tr", "hub&net_tr", "er&net_tr")

for (i in 1:20){
  print(i)
  set.seed(seeds[i])
  
  # Tree network
  tr <- make_tree(n = 20, children = 3, mode = 'undirected')
  net_tr <- igraphtoGNAR(tr) 
  
  #Scale-free network
  sf <-  barabasi.game(n = 20, m = 2, directed = FALSE)
  net_sf <- igraphtoGNAR(sf) 
  
  # Sbm network
  sbm <- sample_sbm(sum(20), pref.matrix = probmat, block.sizes = c(10,10))
  net_sbm <- igraphtoGNAR(sbm)
  
  # Hub network
  hub <- make_star(20, mode = "undirected")
  net_hub <- igraphtoGNAR(hub)
  
  # Random network density 0.1
  er <- erdos.renyi.game(20,p.or.m = 38,type = "gnm",directed = FALSE) 
  net_er <- igraphtoGNAR(er)
  
  # Simulate network data based on the GNAR model
  # Normalize the data
  
  ts_tr <- GNARsim(n = 200, net=net_tr, alphaParams = list(rep(0.2,20), rep(0.4,20), rep(-0.6,20)), betaParams = list(c(0.3,.4),c(0),c(0)))
  tsn_tr <- apply(ts_tr, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_tr <- tsn_tr[1:199, ]
  
  ts_sf <- GNARsim(n = 200, net = net_sf, alphaParams = list(rep(0.2,20), rep(0.4,20), rep(-0.6,20)), betaParams =list(c(0.3,.4),c(0),c(0)))
  tsn_sf <- apply(ts_sf, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_sf <- tsn_sf[1:199, ]
  
  ts_sbm <- GNARsim(n = 200, net = net_sbm, alphaParams =  list(rep(0.2,20), rep(0.4,20), rep(-0.6,20)), betaParams = list(c(0.3,.4),c(0),c(0)))
  tsn_sbm <- apply(ts_sbm, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_sbm <- tsn_sbm[1:199, ]
  
  ts_hub <- GNARsim(n = 200, net = net_hub, alphaParams =  list(rep(0.2,20), rep(0.4,20), rep(-0.6,20)), betaParams =list(c(0.3,.4),c(0),c(0)))
  tsn_hub <- apply(ts_hub, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_hub <- tsn_hub[1:199, ]
  
  ts_er <- GNARsim(n = 200, net = net_er, alphaParams =  list(rep(0.2,20), rep(0.4,20), rep(-0.6,20)), betaParams = list(c(0.3,.4),c(0),c(0)))
  tsn_er <- apply(ts_er, 2, function(x) { x / sd(x[1:199], na.rm = TRUE) })
  datasim_train_er <- tsn_er[1:199, ]
  
  
  # Fit GNAR(3, [2,2,2]) model with tree data and tr network
  fit_pred1 <- predict(GNARfit(vts = datasim_train_tr, net = net_tr,
                               alphaOrder = 3, betaOrder = c(2,0,0),
                               globalalpha = globalalpha))
  rmse_mat[i,paste("tree&net_tr")] <- sqrt(mean((fit_pred1-tsn_tr[200,])^2))
  print(c("seed: ",i," gnar nei1 done"))
  
  
  # Fit GNAR(3, [2,2,2]) model with sf data and tr network
  fit_pred2 <- predict(GNARfit(vts = datasim_train_sf, net = net_tr,
                               alphaOrder = 3, betaOrder = c(2,0,0),
                               globalalpha = globalalpha))
  rmse_mat[i,paste("sf&net_tr")] <- sqrt(mean((fit_pred2-tsn_sf[200,])^2))
  print(c("seed: ",i," gnar nei2 done"))
  
  
  # Fit GNAR(3, [2,2,2]) model with sbm data and tr network
  fit_pred3 <- predict(GNARfit(vts = datasim_train_sbm, net = net_tr,
                               alphaOrder = 3, betaOrder = c(2,0,0),
                               globalalpha = globalalpha))
  rmse_mat[i,paste("sbm&net_tr")] <- sqrt(mean((fit_pred3-tsn_sbm[200,])^2))
  print(c("seed: ",i," gnar nei3 done"))
  
  
  # Fit GNAR(3, [2,2,2]) model with hub data and tr network
  fit_pred4 <- predict(GNARfit(vts = datasim_train_hub, net = net_tr,
                               alphaOrder = 3, betaOrder = c(2,0,0),
                               globalalpha = globalalpha))
  rmse_mat[i,paste("hub&net_tr")] <- sqrt(mean((fit_pred4-tsn_hub[200,])^2))
  print(c("seed: ",i," gnar nei4 done"))
  
  # Fit GNAR(3, [2,2,2]) model with hub data and tr network
  fit_pred5 <- predict(GNARfit(vts = datasim_train_er, net = net_tr,
                               alphaOrder = 3, betaOrder = c(2,0,0),
                               globalalpha = globalalpha))
  rmse_mat[i,paste("er&net_tr")] <- sqrt(mean((fit_pred5-tsn_er[200,])^2))
  print(c("seed: ",i," gnar nei5 done"))
  
}


# Convert the matrix to a data frame
rmse_df <- as.data.frame(rmse_mat)
summary(rmse_df)

# Check the column names
colnames(rmse_df) <- c("tree_net_tr", "sf_net_tr", "sbm_net_tr", "hub_net_tr", "er_net_tr")

# Reshape the data frame to long format
rmse_melt <- melt(rmse_df, variable.name = "Network_Model", value.name = "RMSE")

dev.off()
# Create the box plot
ggplot(rmse_melt, aes(x = Network_Model, y = RMSE, fill = Network_Model)) +
  geom_boxplot() +
  ggtitle("RMSE Distribution for Different Network Structures") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Network Structure") +
  ylab("RMSE") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


###############################################################################
############## (Large) Network with Regime 6, 7 and 8 ###########################
###############################################################################
# Function to create a hub network with multiple hubs
create_hub_network <- function(num_hubs, nodes_per_hub) {
  g <- make_empty_graph(directed = FALSE)
  
  # Track the current highest vertex ID to avoid overlaps
  current_vertex_id <- 0
  
  # Add hubs and their peripheral nodes
  for (i in 1:num_hubs) {
    # Create a star for each hub
    star <- make_star(nodes_per_hub + 1, mode = "undirected")
    
    # Relabel the vertices to avoid overlap
    V(star)$name <- as.character(current_vertex_id + seq_len(vcount(star)))
    current_vertex_id <- current_vertex_id + vcount(star)
    
    # Union the star graph with the main graph
    g <- g + star
  }
  # Connect the hubs with each other (fully connected hubs)
  hub_nodes <- seq(1, current_vertex_id, by = (nodes_per_hub + 1))
  for (i in 1:(length(hub_nodes) - 1)) {
    for (j in (i + 1):length(hub_nodes)) {
      g <- add_edges(g, c(hub_nodes[i], hub_nodes[j]))
    }
  }
  
  return(g)
}


