#==============================================================================#
# Mrode: Chapter 9
#==============================================================================#

# Remove all
  rm(list=ls())

#------------------------------------------------------------------------------#
# Set inputs
#------------------------------------------------------------------------------#

# read in functions for M and Phi matrix
  source(paste0("~/Documents/Programming/R/Animal_Breeding/",
                "Gota_Morota/RandomRegression/legendre.R"))
  source(paste0("~/Documents/Programming/R/Animal_Breeding/",
                "Gota_Morota/RandomRegression/stdtime.R"))
  source(paste0("~/Documents/Programming/R/Animal_Breeding/",
                "Gota_Morota/Pedigrees/createA.R"))
  
#==============================================================================#
# Mrode: Section 9.2
#==============================================================================#

#------------------------------------------------------------------------------#
# Data
#------------------------------------------------------------------------------#
  
# variances
  sigma2.u  <- 5.521
  sigma2.pe <- 8.470
  sigma2.e  <- 3.710
  alpha1 <- sigma2.e / sigma2.u
  alpha2 <- sigma2.e / sigma2.pe

# read in dataset
  y <- c(17.0, 18.6, 24.0, 20.0, 20.0, 15.6, 16.0, 13.0,  8.2,  8.0,
  		 23.0, 21.0, 18.0, 17.0, 16.2, 14.0, 14.2, 13.4, 11.8, 11.4,
  		                               10.4, 12.3, 13.2, 11.6,  8.4,
  		                   22.8, 22.4, 21.4, 18.8, 18.3, 16.2, 15.0,
  		 22.2, 20.0, 21.0, 23.0, 16.8, 11.0, 13.0, 17.0, 13.0, 12.6)
  DIM <- c(4, 38, 72, 106, 140, 174, 208, 242, 276, 310,
           4, 38, 72, 106, 140, 174, 208, 242, 276, 310,
                                174, 208, 242, 276, 310,
                      106, 140, 174, 208, 242, 276, 310,
           4, 38, 72, 106, 140, 174, 208, 242, 276, 310)
  HTD <- c(1,2,3,4,5,6,7,8,9,10,
  		   1,2,3,4,5,6,7,8,9,10,
  		   			 6,7,8,9,10,
  		         4,5,6,7,8,9,10,
  		   1,2,3,4,5,6,7,8,9,10)
  Animal <- c(rep(4, 10), rep(5, 10), rep(6, 5), rep(7, 7), rep(8, 10))
  
# put data together
  data <- data.frame(Animal, DIM, HTD, y)

# remove vectors
  rm(y); rm(DIM); rm(HTD); rm(Animal)
  
# Pedigree
  Animal <- 1:8
  Sire   <- c(0, 0, 0, 1, 3, 1, 3, 1)
  Dam    <- c(0, 0, 0, 2, 2, 5, 4, 7)

# Pedigree
  Ped <- data.frame(Animal, Sire, Dam)
  
# remove vectors
  rm(Animal); rm(Sire); rm(Dam)
  
# write data out for ASReml
  write.table(data, file=paste0("/Users/austinputz/Documents/Programming/R/", 
                                "Animal_Breeding/Mrode/ASReml/ch9.txt"),
        quote=FALSE, row.names=FALSE, col.names=TRUE)
  write.table(Ped, file=paste0("/Users/austinputz/Documents/Programming/R/", 
                                    "Animal_Breeding/Mrode/ASReml/ch9.ped"),
        quote=FALSE, row.names=FALSE, col.names=FALSE)
  
#------------------------------------------------------------------------------#
# Get matrices for Legendre polynomials
#------------------------------------------------------------------------------#

# get unique sorted DIM values
  DIM.uniq <- unique(data$DIM)

# get M matrix
  M <- stdtime(DIM.uniq, 4)
  rownames(M) <- DIM.uniq
  
# get Lambda matrix
  Lambda <- legendre(4, gengler=FALSE)

# get Phi matrix
  Phi <- as.data.frame(M %*% Lambda)
  
# add rownames as column to merge with
  Phi$DIM <- rownames(Phi)

# merge to dataset
  data.leg <- merge(data, Phi, by="DIM")

# sort it by Animal and DIM
  data.leg <- data.leg[order(data.leg$Animal, data.leg$DIM), 
                      c("Animal","DIM","HTD","y","V1","V2","V3","V4","V5")]
  
#------------------------------------------------------------------------------#
# Get matrices for MME
#------------------------------------------------------------------------------#
  
# change DIM to factor
  data.leg$DIM <- as.factor(data.leg$DIM)
  
# MME
  X  <- model.matrix(y ~ -1 + DIM + V1+V2+V3+V4+V5, data=data.leg)
  X  <- X[, c(1:9, 11:15)]
  
# split X into 2
  X1 <- model.matrix(y ~ -1 + DIM, data=data.leg)
  X1 <- X1[, c(1:9)]
  X2 <- model.matrix(y ~ -1 + V1+V2+V3+V4+V5, data=data.leg)
  
# Enter his X2 matrix
  X2pX2.book <- matrix(c(20.9996, -4.4261, 4.0568, -0.8441, 8.7149,
                      -4.4261, 24.6271, -4.7012, 11.1628, -3.0641,
                      4.0568, -4.7012, 31.0621, -6.6603, 19.0867,
                      -0.8441, 11.1628, -6.6603, 38.6470, -8.8550,
                      8.7149, -3.0641, 19.0867, -8.8550, 48.2930), 
                              byrow=TRUE, ncol=5, nrow=5)
  
  X.new <- matrix(0, nrow=14, ncol=14)
  X.new[1:9, 1:9] <- t(X1) %*% X1
  X.new[10:14, 10:14] <- X2pX2.book
  
# Create Q
  Q1 <- matrix(0, 
            ncol=length(unique(Ped$Animal))-length(unique(data.leg$Animal)), 
            nrow=nrow(data.leg))
  Q2 <- model.matrix(y ~ -1 + as.factor(Animal), data=data.leg)
  Q <- cbind(Q1, Q2)

# Create Z
  Z <- Q2

# create A (his A-1 in chapter 4 was incorrect)
  A    <- createA(Ped)
  Ainv <- solve(A)
  # Ainv <- matrix(c(2.50, 0.50, 0.00, -1.00, 0.50, -1.00, 0.50, -1.00,
  #                  0.50, 1.50, 0.00, -1.00, 0.00, 0.00, 0.00, 0.00,
  #                  0.0, 0.00, 1.83, 0.50, -0.67, 0.00, -1.00, 0.00,
  #                  -1.00, -1.00, 0.50, 2.50, 0.00, 0.00, -1.00, 0.00,
  #                  0.50, 0.00, -0.67, 0.00, 1.83, -1.00, 0.00, 0.00,
  #                  -1.00, 0.00, 0.00, 0.00, -1.00, 2.00, 0.00, 0.00,
  #                  0.50, 0.00, -1.00, -1.00, 0.00, 0.00, 2.50, -1.00,
  #                  -1.00, 0.00, 0.00, 0.00, 0.00, 0.00, -1.00, 2.00), 
  #     nrow=8, ncol=8) # This is from page 64 in chapter 4
  
#------------------------------------------------------------------------------#
# Solve MME function
#------------------------------------------------------------------------------#

# create function to solve MME
  fixedRegMME <- function(X, Q, Z, Ainv, y, alpha1, alpha2) {
    
    #----------------------------------------#
    # LHS
    #----------------------------------------#
    
    # multiply row 1
    XpX <- t(X) %*% X
    XpQ <- t(X) %*% Q
    XpZ <- t(X) %*% Z
    
    # multiply row 2
    QpX <- t(Q) %*% X
    QpQ <- (t(Q) %*% Q) + Ainv*alpha1
    QpZ <- t(Q) %*% Z
    
    # multiply row 3
    ZpX <- t(Z) %*% X
    ZpQ <- t(Z) %*% Q
    ZpZ <- (t(Z) %*% Z) + alpha2
    
    # concatenate top row
    top <- cbind(XpX, XpQ, XpZ)
    mid <- cbind(QpX, QpQ, QpZ)
    bot <- cbind(ZpX, ZpQ, ZpZ)
    
    # rbind all for LHS
    LHS <- rbind(top, mid, bot)
    
    #----------------------------------------#
    # RHS
    #----------------------------------------#
    
    # get RHS elements
    rhs.top <- t(X) %*% y
    rhs.mid <- t(Q) %*% y
    rhs.bot <- t(Z) %*% y
    
    # RHS
    RHS <- rbind(rhs.top, rhs.mid, rhs.bot)
    
    # solve MME
    solutions <- solve(LHS) %*% RHS
    
    # return them
    return(list(LHS=LHS, RHS=RHS, sols=solutions))
     
  }
  

# Calculate the output
  output <- fixedRegMME(X=X, Q=Q, Z=Z, Ainv=Ainv, y=data.leg$y, alpha1, alpha2)
  
# Get the LHS output
  output$LHS - output2$LHS
  
# Get the RHS output
  output$RHS - output2$RHS

# get solutions rounded to 4 digits like he has
  round(sols, 4)
  
#------------------------------------------------------------------------------#
# Plot lactation curve
#------------------------------------------------------------------------------#

  library(ggplot2)
  
# plot the lactation curve
  b2_hat <- matrix(c(16.3082, -0.5227, -0.1245, 0.5355, -0.4195))
  data.plot <- as.data.frame(as.matrix(Phi[, 1:5]) %*% b2_hat)
  data.plot$DIM <- rownames(data.plot)
  data.plot$DIM <- as.numeric(data.plot$DIM)
  
# plot lactation curve points
  ggplot(data=data.plot, aes(x=DIM, y=V1, group=1)) +
    geom_line(color="steelblue") +
    geom_point(color="steelblue", size=3) +
    ggtitle("Section 9.2: Fixed Regression Model") +
    xlab("Days in Milk (DIM)") +
    ylab("Estimates values for Leg Polynomials")
  
  
#==============================================================================#
# Mrode: Section 9.3 Random Regression Model
#==============================================================================#

# Remove all
  rm(list=ls())
  
  library(dplyr)
  library(tidyr)

# read in functions for M and Phi matrix
  source(paste0("~/Documents/Programming/R/Animal_Breeding/",
                "Gota_Morota/RandomRegression/legendre.R"))
  source(paste0("~/Documents/Programming/R/Animal_Breeding/",
                "Gota_Morota/RandomRegression/stdtime.R"))
  source(paste0("~/Documents/Programming/R/Animal_Breeding/",
                "Gota_Morota/Pedigrees/createA.R"))
  
#------------------------------------------------------------------------------#
# Data
#------------------------------------------------------------------------------#
  
# variances
  sigma2.u  <- 5.521
  sigma2.pe <- 8.470
  sigma2.e  <- 3.710
  alpha1 <- sigma2.e / sigma2.u
  alpha2 <- sigma2.e / sigma2.pe

# read in dataset
  y <- c(17.0, 18.6, 24.0, 20.0, 20.0, 15.6, 16.0, 13.0,  8.2,  8.0,
  		 23.0, 21.0, 18.0, 17.0, 16.2, 14.0, 14.2, 13.4, 11.8, 11.4,
  		                               10.4, 12.3, 13.2, 11.6,  8.4,
  		                   22.8, 22.4, 21.4, 18.8, 18.3, 16.2, 15.0,
  		 22.2, 20.0, 21.0, 23.0, 16.8, 11.0, 13.0, 17.0, 13.0, 12.6)
  DIM <- c(4, 38, 72, 106, 140, 174, 208, 242, 276, 310,
           4, 38, 72, 106, 140, 174, 208, 242, 276, 310,
                                174, 208, 242, 276, 310,
                      106, 140, 174, 208, 242, 276, 310,
           4, 38, 72, 106, 140, 174, 208, 242, 276, 310)
  HTD <- c(1,2,3,4,5,6,7,8,9,10,
  		   1,2,3,4,5,6,7,8,9,10,
  		   			 6,7,8,9,10,
  		         4,5,6,7,8,9,10,
  		   1,2,3,4,5,6,7,8,9,10)
  Animal <- c(rep(4, 10), rep(5, 10), rep(6, 5), rep(7, 7), rep(8, 10))
  
# put data together
  data <- data.frame(Animal, DIM, HTD, y)

# remove vectors
  rm(y); rm(DIM); rm(HTD); rm(Animal)
  
# Pedigree
  Animal <- 1:8
  Sire   <- c(0, 0, 0, 1, 3, 1, 3, 1)
  Dam    <- c(0, 0, 0, 2, 2, 5, 4, 7)

# Pedigree
  Ped <- data.frame(Animal, Sire, Dam)
  
# remove vectors
  rm(Animal); rm(Sire); rm(Dam)
  
#------------------------------------------------------------------------------#
# Get matrices for Legendre polynomials
#------------------------------------------------------------------------------#

# get unique sorted DIM values
  DIM.uniq <- unique(data$DIM)

# get M matrix
  M <- stdtime(DIM.uniq, 4)
  rownames(M) <- DIM.uniq
  
# get Lambda matrix
  Delta <- legendre(4)

# get Phi matrix
  Phi <- M %*% Delta
  
# Phi data frame
  Phi.df <- data.frame(Phi)
  
# add rownames as column to merge with
  Phi.df$DIM <- rownames(Phi)

# merge to dataset
  data.leg <- merge(data, Phi.df, by="DIM")

# sort it
  data.leg <- data.leg[order(data.leg$Animal, data.leg$DIM), 
                      c("Animal","DIM","HTD","y","X1","X2","X3","X4","X5")]
  
#------------------------------------------------------------------------------#
# Get matrices for MME
#------------------------------------------------------------------------------#
  
# create A and invert it
  A <- createA(Ped)
  Ainv <- solve(A)
  
# G matrix for breeding values
  G = matrix(c(3.297,  0.594, -1.381,
               0.594,  0.921, -0.289,
              -1.381, -0.289,  1.005), byrow=TRUE, nrow=3)
  
# P matrix for permanent env effects
  P = matrix(c(6.872, -0.254, -1.101,
              -0.254,  3.171,  0.167,
              -1.101,  0.167,  2.457), byrow=TRUE, nrow=3)
  
# Get R matrix for residuals
  R <- diag(rep(sigma2.e, nrow(data)))
  
# Invert the R matrix for MME
  Rinv <- solve(R)
  
# set HTD as factor for X matrix
  data.leg$HTD <- as.factor(data.leg$HTD)
  
# MME
  X  <- model.matrix(y ~ -1 + HTD + X1+X2+X3+X4+X5, data=data.leg)
  X <- X[, c(1:9, 11:15)] # had to cut down because of singularities... 
  
# Create Q
  Q1 <- matrix(0, 
          ncol=(length(unique(Ped$Animal))-length(unique(data.leg$Animal)))*ncol(Phi[, c(1:3)]), 
          nrow=nrow(data.leg))
  Q2 <- matrix(0,
          ncol=length(unique(data.leg$Animal))*ncol(Phi[, c(1:3)]),
          nrow=nrow(data.leg))
  Q2[1:10, 1:3]    <- as.matrix(data.leg[1:10, c("X1", "X2", "X3")])
  Q2[11:20, 4:6]   <- as.matrix(data.leg[11:20, c("X1", "X2", "X3")])
  Q2[21:25, 7:9]   <- as.matrix(data.leg[21:25, c("X1", "X2", "X3")])
  Q2[26:32, 10:12] <- as.matrix(data.leg[26:32, c("X1", "X2", "X3")])
  Q2[33:42, 13:15] <- as.matrix(data.leg[33:42, c("X1", "X2", "X3")])
  
# combine them together (side by side)
  Q <- cbind(Q1, Q2)
  
# set row and column names for Q
  colnames(Q) <- c(rep(1:8, each=3))
  rownames(Q) <- data.leg$Animal
  
# Create Z
  Z <- Q2
  colnames(Z) <- c(rep(4:8, each=3))
  rownames(Z) <- data.leg$Animal
  
# set I to kronecker by P
  I <- diag(length(unique(data.leg$Animal)))
  
#------------------------------------------------------------------------------#
# Solve MME function
#------------------------------------------------------------------------------#
  
# create function to solve MME
  RandomRegMME <- function(X, Q, Z, I, Rinv, Ainv, y, G, P) {
    
    #----------------------------------------#
    # LHS
    #----------------------------------------#
    
    # multiply row 1
    XpX <- t(X) %*% Rinv %*% X
    XpQ <- t(X) %*% Rinv %*% Q
    XpZ <- t(X) %*% Rinv %*% Z
    
    # multiply row 2
    QpX <-  t(Q) %*% Rinv %*% X
    QpQ <- (t(Q) %*% Rinv %*% Q) + (Ainv %x% G)
    QpZ <-  t(Q) %*% Rinv %*% Z
    
    # multiply row 3
    ZpX <-  t(Z) %*% Rinv %*% X
    ZpQ <-  t(Z) %*% Rinv %*% Q
    ZpZ <- (t(Z) %*% Rinv %*% Z) + (I %x% P)
    
    # concatenate top row
    top <- cbind(XpX, XpQ, XpZ)
    mid <- cbind(QpX, QpQ, QpZ)
    bot <- cbind(ZpX, ZpQ, ZpZ)
    
    # rbind all for LHS
    LHS <- rbind(top, mid, bot)
    
    #----------------------------------------#
    # RHS
    #----------------------------------------#
    
    # get RHS elements
    rhs.top <- t(X) %*% Rinv %*% y
    rhs.mid <- t(Q) %*% Rinv %*% y
    rhs.bot <- t(Z) %*% Rinv %*% y
    
    # RHS
    RHS <- rbind(rhs.top, rhs.mid, rhs.bot)
    
    # solve MME
    solutions <- solve(LHS) %*% RHS
    
    # return them
    return(list(LHS=LHS, RHS=RHS, sols=solutions))
     
  }
  
# get output from RR model
  output <- RandomRegMME(X=X, Q=Q, Z=Z, I=I, 
                                Rinv=Rinv, Ainv=Ainv, 
                                y=data.leg$y, G=G, P=P)
  
# Get LHS
  round(output$LHS, 2)
  round(output$RHS, 2)
  
# solutions
  round(output$sols, 2)
  
#------------------------------------------------------------------------------#
# Enter solutions since I can't solve for them
#------------------------------------------------------------------------------#
  
# add solutions
  b1 <- matrix(c(10.0862, 7.5908, 8.5601, 8.2430, 6.3161, 3.0101,
                 3.1085, 3.1718, 0.5044, 0.0000), ncol=1)
  b2 <- matrix(c(16.6384, -0.6253, -0.1346, 0.3479, -0.4218), ncol=1)
  an.sol <- matrix(c(-0.0583, 0.0552, -0.0442, 
                     -0.0728, -0.0305, -0.0244,
                     0.1311, -0.0247, 0.0686,
                     0.3445, 0.0063, -0.3164,
                     -0.4537, -0.0520, 0.2798,
                    -0.5485, 0.0730, 0.1946,
                    0.8518, -0.0095, -0.3131,
                    0.2209, 0.0127, -0.0174), byrow=TRUE, ncol=3)
  pe.sol <- matrix(c(-0.6487, -0.3601, -1.4718,
                     -0.7761, 0.1370, 0.9688,
                     -1.9927, 0.9851, -0.0693,
                      3.5188, -1.0510, -0.4048,
                     -0.1013, 0.2889, 0.9771), byrow=TRUE, ncol=3)
  
#------------------------------------------------------------------------------#
# Get genetic and PE variances and covariances
#------------------------------------------------------------------------------#
  
# Get genetic and pe variances and covariances
  covar.gen <- Phi[, 1:3] %*% G %*% t(Phi[, 1:3])
  covar.pe  <- Phi[, 1:3] %*% P %*% t(Phi[, 1:3])
  
# create df
  data.gen.pe <- data.frame(x=c(4,38,72,106,140,174,208,242,276,310),
                        gen.var=diag(covar.gen),
                        pe.var =diag(covar.pe))
  
  data.gen.pe.melt <- data.gen.pe %>% gather(Type, Value, gen.var:pe.var)
  
# plot bv variances (on diagonal)
  ggplot(data=data.gen.pe.melt, aes(x=x, y=Value, color=Type, group=Type)) +
    geom_line() +
    geom_point()
  

  
#----------------------------------------#
# new Phi matrix
#----------------------------------------#
  
# get M matrix
  M <- stdtime(6:310, 2)
  rownames(M) <- 6:310
  
# get Lambda matrix
  Delta <- legendre(2)
  
# get Phi matrix
  Phi.305 <- M %*% Delta
  
# sum columns
  apply(Phi.305, 2, sum)
  
#----------------------------------------#
# Get daily BV's
#----------------------------------------#
  
# 305 d yields
  bvs.305d <- apply(Phi.305 %*% t(an.sol), 2, sum)
  
  matrix(c(215.6655, 2.4414, -1.5561), byrow=TRUE, nrow=1) %*% 
    matrix(c(0.3445, 0.0063, -0.3164), ncol=1)
  
# daily BVs for each animal
  daily.bvs <- Phi.305 %*% t(an.sol)
  daily.bvs.df <- as.data.frame(daily.bvs)
  
# rename columns to animal numbers
  names(daily.bvs.df) <- 1:8
  
# melt into one column
  daily.bvs.melt <- daily.bvs.df %>% gather(Animal, Value)
  
# add DIM
  daily.bvs.melt$DIM <- rep(6:310, 8)
  
# plot individual animals
  plot(x=seq(6, 310, 1), y=daily.bvs[, 8])
  
# plot each individuals BV for each day
  ggplot(data=daily.bvs.melt, aes(x=DIM, y=Value, color=Animal, linetype=Animal, group=Animal)) +
    geom_line() +
    scale_shape_manual(values=1:8) +
    ggtitle("Daily EBVs for each Animal") +
    xlab("Days in Milk (DIM)") +
    ylab("EBV")
  
#==============================================================================#
# Plot lactation curve
#==============================================================================#
  
  source("~/Documents/Programming/R/Functions/plotCurve.R")

# Plot lactation curve from fixed effects, example from chapter 9
  x     <- seq(-1, 1, 0.01)  
  coefs <- c(16.6384, -0.6353, -0.1346, 0.3479, -0.4218)
  
# plot cure
  plotCurve(x=x, coefs=coefs)
  
  
  
#==============================================================================#
# For my pig study
#==============================================================================#
 
# 10th order
  coefs <- c(0.000, -0.4723, 0.8492, 0.8110, -0.3259, -0.7682, 
            -0.5576E-01, 0.3954, -0.8676E-01, -0.1890, -0.9973E-01)
  y <- (coefs[1] + x*coefs[2] + (x)^2*coefs[3] + (x)^3*coefs[4] + (x)^4*coefs[5] +
       (x)^5*coefs[6] + (x)^6*coefs[7] + (x)^7*coefs[8] + (x)^8*coefs[9] + 
        (x)^9*coefs[10]) 
  plot(x,y, type="l")
  
# 15th order
  coefs <- c(0.000, -0.4659, 0.8455, 0.8281, -0.3101, -0.7693, -0.5448E-01, 0.4085,
              -0.8680E-01, -0.1830, -0.9150E-01, -0.7782E-01, 
              -0.1060, -0.1392, -0.4406E-01, 0.1673)
  y <- (coefs[1] + x*coefs[2] + (x)^2*coefs[3] + (x)^3*coefs[4] + (x)^4*coefs[5] +
       (x)^5*coefs[6] + (x)^6*coefs[7] + (x)^7*coefs[8] + (x)^8*coefs[9] + 
        (x)^9*coefs[10] + (x)^10*coefs[11] + (x)^11*coefs[12] + (x)^12*coefs[13] +
       (x)^13*coefs[14] + (x)^14*coefs[15] + (x)^15*coefs[16])
  plot(x,y, type="l")
  
  
  
  
  
  
  


