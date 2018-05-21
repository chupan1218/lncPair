
setwd("C:\\Users\\panchu\\Desktop\\lncPair")

#read BRCA-specific lncRNA expression data
lncexpression <- read.csv('BRCAlncexpression.csv', sep = ",", header = T, row.names = 1, stringsAsFactors = F, check.names = F)

#divide the tumors into multiple cohorts according to the cancer subtype information
lncexpression1 <- lncexpression[, c(1:140)]    #Basal-like
lncexpression2 <- lncexpression[, c(141:205)]  #HER2-enriched
lncexpression3 <- lncexpression[, c(206:625)]  #Luminal A
lncexpression4 <- lncexpression[, c(626:799)]  #Luminal B

#calculate Pearson correlation coefficient matrix for each case
pcc1 <- cor(t(lncexpression1), t(lncexpression1))

pcc2 <- cor(t(lncexpression2), t(lncexpression2))

pcc3 <- cor(t(lncexpression3), t(lncexpression3))

pcc4 <- cor(t(lncexpression4), t(lncexpression4))


if(isSymmetric(pcc1)) {
    warning(paste("lncRNA Pearson correlation coefficient matrix interaction is symmetric!",
                  "but, only the Upper triangle matrix is used in this situtation.", sep="\n"))
    pcc1[lower.tri(pcc1)] <- 0		
    
    diag(pcc1) <- 0 #the diagonal elements are set to zero
  }


if(isSymmetric(pcc2)) {
  warning(paste("lncRNA Pearson correlation coefficient matrix interaction is symmetric!",
                "but, only the Upper triangle matrix is used in this situtation.", sep="\n"))
  
  pcc2[lower.tri(pcc2)] <- 0		
  
  diag(pcc2) <- 0 
}


if(isSymmetric(pcc3)) {
  warning(paste("lncRNA Pearson correlation coefficient matrix interaction is symmetric!",
                "but, only the Upper triangle matrix is used in this situtation.", sep="\n"))
  
  pcc3[lower.tri(pcc3)] <- 0		
  
  diag(pcc3) <- 0 
}


if(isSymmetric(pcc4)) {
  warning(paste("lncRNA Pearson correlation coefficient matrix interaction is symmetric!",
                "but, only the Upper triangle matrix is used in this situtation.", sep="\n"))
  
  pcc4[lower.tri(pcc4)] <- 0		
  
  diag(pcc4) <- 0 
}


#tranfer the upper triangle matrix into a vector for each case
p1<- c(); p2 <- c(); p3 <- c(); p4 <- c()
pairname <- c(); pairindex <- c()

for(i in 1:(nrow(lncexpression)-1)){
  name1 <- rownames(lncexpression)[i]
  index1 <- i
  
  for(j in (i+1):nrow(lncexpression)){
    p1 <- c(p1, pcc1[i,j])
    p2 <- c(p2, pcc2[i,j])
    p3 <- c(p3, pcc3[i,j])
    p4 <- c(p4, pcc4[i,j])
    
    name2 <- colnames(pcc1)[j]
    index2 <- j
    
    pairname <- c(pairname, paste(name1, name2, sep = ","))
    pairindex <- c(pairindex, paste(index1, index2, sep = ","))
  }
}



subtypes = c("pcc1","pcc2","pcc3","pcc4")

#the absolute value of the Pearson correlation coefficient
zmat <- matrix(0, nrow = length(p1), ncol = length(subtypes))
colnames(zmat) <- subtypes

#the Pearson correlation coefficient
mat <- matrix(0, nrow = length(p1), ncol = length(subtypes))
colnames(mat) <- subtypes


mat[,1] <- p1 
mat[,2] <- p2
mat[,3] <- p3
mat[,4] <- p4


zmat[,1] <- abs(p1)
zmat[,2] <- abs(p2)
zmat[,3] <- abs(p3)
zmat[,4] <- abs(p4)



#local fdr calculation   
library(locfdr)
library(splines)
library(MASS)
N_subtypes <- length(subtypes)
N <- nrow(zmat)
lfdr = matrix(0,N,N_subtypes)  #Local fdr
p0 = rep(0,N_subtypes)         #Null prior probability   
f0 = matrix(0,N,N_subtypes)    #Null prior density
p1f1 = matrix(0,N,N_subtypes)  #Non-null prior density * p1

#local fdr
for(i in 1:N_subtypes) {
  Mfdr = locfdr(zmat[,i]);
  lfdr[,i] = Mfdr$fdr;
  p0[i] = Mfdr$fp0[3,3]
  f0[,i] = predict(interpSpline(Mfdr$mat[,1],Mfdr$mat[,6]),zmat[,i])$y
  p1f1[,i] = predict(interpSpline(Mfdr$mat[,1],Mfdr$mat[,11]),zmat[,i])$y
}

density_0 <- matrix(0, nrow = N, ncol = N_subtypes)
for(i in 1:N_subtypes){
  density_0[,i] <- lfdr[,i]/p0[i]*(p0[i]*f0[,i]+p1f1[,i]) 
}

density_1 <- matrix(0, nrow = N, ncol = N_subtypes)
for(i in 1:N_subtypes){
  density_1[,i] <- (1-lfdr[,i])/(1-p0[i])*(p0[i]*f0[,i]+p1f1[,i])
}

#Calculate likelihood of zscores
J <- 2^N_subtypes
pivec <- rep(1/J,J)
likelihood <- matrix(0,nrow=N, ncol=J)
status = data.frame(c(rep(0,J/2),rep(1,J/2)))
for(coln in 2:N_subtypes){
  status <- cbind(status, rep(c(rep(0,J/(2^coln)),rep(1,J/(2^coln))),2^(coln-1)))
}
colnames(status) <- c(1:N_subtypes)

pioutput <- array(c(0, rep(1/J,J)))

for(i in 1:nrow(likelihood)) { 
  output <- matrix(0, nrow=J, ncol=N_subtypes)
  for(td in 1:N_subtypes){
    output[status[,td]==0,td] <- density_0[i,td]
    output[status[,td]==1,td] <- density_1[i,td]
  }
  likelihood[i,] <- apply(output,1, prod)
}


#iteration functions
myiteration <- function(pivec, tol=.0001){
  likelihoodnew <- likelihood*matrix(rep(pivec,N),byrow = TRUE,nrow = N)
  conditionalprobsnew <- likelihoodnew/apply(likelihoodnew,1,sum)
  
  
  marginalprobsnew <- matrix(0, nrow = N, ncol = 2)
  colnames(marginalprobsnew) <- c("positive probability", "negative probability")
                             #dimnames = list( pairindex,c("positive probability", "negative probability"))) # two probabilities, (1,1) and (0,0)
  
  for(i in 1:N){
    
    marginalprobsnew[i,1] <- conditionalprobsnew[i,which(rowSums(status)==N_subtypes)] #positive probability
    
    marginalprobsnew[i,2] <- conditionalprobsnew[i,which(rowSums(status)==0)] #negative probability
  }

  
  pivecnew <- (1/N)*apply(conditionalprobsnew,2,sum)
  stopyn = "n"
  if(max(abs(pivec - pivecnew))<tol) stopyn = "y"
  return(list(conditionalprobs= conditionalprobsnew, marginalprobs=marginalprobsnew, pivec = pivecnew, stopyn = stopyn))
}


#iteration
stopyn <- FALSE
itnumber <- 1
limit <- 100

while(!stopyn & itnumber <=  limit){
  result <- myiteration(pivec)
  pivec <- result$pivec
  if(result$stopyn == "y") stopyn <- TRUE
  flush.console()
  print(itnumber)
  print(result$pivec)
  pioutput <- rbind(pioutput, c(itnumber, result$pivec))
  
  itnumber <- itnumber + 1
}

#the posterior marginal probability of similar and dissimilar of lncRNA pairs
mar <- result$marginalprobs 


MAR <- mar #the posterior marginal probability matrix
MAT <- mat #the Pearson correlation coefficient matrix
ZMAT <- zmat #the absolute Pearson correlation coefficient matrix

MAR <- MAR[order(mar[,1],decreasing = T),]
MAT <- MAT[order(mar[,1],decreasing = T),]
ZMAT <- ZMAT[order(mar[,1],decreasing = T),]

pairindex1 <- as.matrix(pairindex[order(mar[,1],decreasing = T)])
pairname1 <- as.matrix(pairname[order(mar[,1],decreasing = T)])

OUTPUT <- cbind(pairindex1,pairname1,MAR,MAT)


write.csv(OUTPUT, 'OUTPUT.csv')























































































































# local fdr calculation   
library(locfdr)
library(splines)
library(MASS)
localfdr <- function(zmat, algorithms){ 
  N <- nrow(zmat)
  N_alg <-length(algorithms)
  
  lfdr=matrix(0,N,N_alg)  #Local fdr
  p0=rep(0,N_alg)         #Null prior probability   
  f0=matrix(0,N,N_alg)    #Null prior density
  p1f1=matrix(0,N,N_alg)  #Non-null prior density * p1
  
  # local fdr
  for(i in 1:N_alg) {
    Mfdr = locfdr(zmat[,i]);
    lfdr[,i] = Mfdr$fdr;
    p0[i] = Mfdr$fp0[3,3]
    f0[,i] = predict(interpSpline(Mfdr$mat[,1], Mfdr$mat[,6]), zmat[,i])$y
    p1f1[,i] = predict(interpSpline(Mfdr$mat[,1], Mfdr$mat[,11]), zmat[,i])$y
  }
  
  density_0 <- matrix(0, nrow = N, ncol = N_alg)
  for(i in 1:N_alg){
    density_0[,i] <- lfdr[,i]/p0[i]*(p0[i]*f0[,i]+p1f1[,i]) #realize p(zdjk|tdjk=0)
  }
  
  density_1 <- matrix(0, nrow = N, ncol = N_alg)
  for(i in 1:N_alg){
    density_1[,i] <- (1-lfdr[,i])/(1-p0[i])*(p0[i]*f0[,i]+p1f1[,i]) #p(zdjk|tdjk=1)
  }
  
  # Calculate likelihood of zscores     
  J <- 2^N_alg
  pivec <- rep(1/J,J) #initial probability for four states: (0.25, 0.25, 0.25, 0.25)
  likelihood <- matrix(0,nrow=N, ncol=J)
  status = data.frame(c(rep(0,J/2),rep(1,J/2))) #(0, 0, 1, 1)
  
  for(coln in 2:N_alg){
    status <- cbind(status, rep(c(rep(0,J/(2^coln)),rep(1,J/(2^coln))),2^(coln-1)))
  }
  colnames(status) <- c(1:N_alg)
  
  pioutput <- array(c(0, rep(1/J,J))) #iteration probability
  
  for(i in 1:nrow(likelihood)) { 
    output <- matrix(0, nrow=J, ncol=N_alg) # 4 by 2 matrix 
    for(td in 1:N_alg){
      output[status[,td]==0,td] <- density_0[i,td] #first column status==0 <- density_0[i,1], second column status==0 <- density_0[i,2]
      output[status[,td]==1,td] <- density_1[i,td] #first column status==1 <- density_1[i,1], second column status==1 <- density_0[i,2]
    }
    likelihood[i,] <- apply(output,1, prod) #prod is a function of product
  }
  
  return (list(status = status ,likelihood = likelihood)) #p(z1jk,...,zDjk|t1jk,...,tDjk)=p(z1jk|t1jk)...p(zDjk|tDjk)
}


data <- localfdr(zmat, algorithms)


#main functions 
iteration <- function(zmat, data, algorithms, tol=.0001){
  N <- nrow(zmat)
  N_alg <- length(algorithms)
  J <- 2^N_alg
  pivec <- rep(1/J,J) #initial probability: 0.25, 0.25, 0.25, 0.25
  
  
  likelihoodnew <- data$likelihood * matrix(rep(pivec,N),byrow=TRUE,nrow=N) # likelihood * [ 0.25 0.25 0.25 0.25]
  conditionalprobsnew <- likelihoodnew/apply(likelihoodnew,1,sum) # 
  
  
  marginalprobsnew <- matrix(0, nrow=N, ncol=2) # two probabilities, (1,1) and (0,0)
  
  for(i in 1:N){

    marginalprobsnew[i,1] <- conditionalprobsnew[i,which(rowSums(data$status)==N_alg)] #positive probability
    
    marginalprobsnew[i,2] <- conditionalprobsnew[i,which(rowSums(data$status)==0)] #negative probability
  }
  
  pivecnew <- (1/N)*apply(conditionalprobsnew,2,sum) #update pivec
  stopyn = "n"
  
  if(max(abs(pivec - pivecnew))<tol) stopyn = "y" 
  return(list(conditionalprobs= conditionalprobsnew, marginalprobs=marginalprobsnew, pivec = pivecnew, stopyn = stopyn))
}

### iteration ###
stopyn <- FALSE
itnumber <- 1
limit <- 100

while(!stopyn & itnumber <=  limit){
  result <- iteration(zmat, data, algorithms, tol=.0001)
  pivec <- result$pivec
  if(result$stopyn == "y") stopyn <- TRUE
  flush.console()
  print(itnumber)
  print(result$pivec)
  pioutput <- rbind(pioutput, c(itnumber, result$pivec))
  
  itnumber <- itnumber + 1
}


mar <- result$marginalprobs



