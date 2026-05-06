#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("Missing parameters.n", call.=FALSE)
}

outPath = args[1]
print("This is the visualisation-tool.")
print(paste("Outpath is",outPath,sep = " "))

estimators = rep(NA,length(args)-1)
if (length(args)>1) {
  for(i in 2:length(args)) {
    # save the names of the estimators
    estimators[i-1] = args[i]
  }
}

scalar1 <- function(x) {x / sum(x)}

#dieChoice history
dieChoiceHistory <-function(path2files, it)
{
  iterations = it
  dieHistory = rep(NA,iterations)
  
  for(i in 1:iterations){
    adv = read.csv2(paste(path2files,"Toy",i,".csv" ,sep=""))
    #dieHistory[i] = as.character(adv$best[1])
    dieHistory[i] = levels(adv$best[1])
  }
  #plot(c(1:iterations),as.factor(dieHistory), xlab="Iteration", ylab = "dice", type = "p", main = "Die Choice")
  plot(c(1:iterations),as.factor(dieHistory), xlab="Iteration", ylab = "dice", type = "p")
  mtext("Iteration", side=1, line=2,outer = FALSE,cex = 0.5)
}


# we need this to determine the plot ranges
findMaxScore <- function(path,estimators){
  ylimMax = -100000000000000
  ylimMin = 1000000000000000
  for(i in 1: length(estimators)){
    files <- list.files(path= paste(path, estimators[i], sep = ""), pattern="Toy.*.csv", full.names=T, recursive=FALSE)
    
    it = length(files)
    for(j in 1: it){
      adv = read.csv2(paste(path, estimators[i],"/Toy",j,".csv" ,sep=""))
      adv$totalProfit <- as.numeric(as.character(adv$totalProfit))
      
      pot <- adv$totalProfit[1]
      
      
      #adv$totalScore <- as.numeric(as.character(adv$totalScore))        # <----------------
      
      #pot <- adv$totalScore[1]        # <----------------
      
      if(ylimMax < as.numeric(pot)) {
        ylimMax = pot
      }
      if(ylimMin > as.numeric(pot)) {
        ylimMin = pot
      }
    }
  }
  print(ylimMax)
  re = c(ylimMin,ylimMax)
  return(re)
}

totalScoreHistory <-function(path2files, it, yl)
{
  iterations = it
  scoreHistory = rep(NA,iterations)
  profitHistory <- numeric(iterations)
  
  for(i in 1:iterations){
    adv = read.csv2(paste(path2files,"Toy",i,".csv" ,sep=""))
    

    adv$totalProfit <- as.numeric(as.character(adv$totalProfit))

    profitHistory[i] <- adv$totalProfit[1]
    
    
    #adv$totalScore <- as.numeric(as.character(adv$totalScore))
    
   # profitHistory[i] <- adv$totalScore[1]

  }
#  plot(c(1:iterations),scoreHistory, xlab="Iteration", ylab = "totalScore", type = "l", 
#       main = paste("Total Score: ", scoreHistory[iterations], sep=""))
  
  #print(paste("ylim ", ylimMax, sep = ""))
  plot(c(1:iterations),profitHistory, xlab="Iteration", ylab = "Profit", type = "l", 
       main = paste("Final Profit: ", profitHistory[iterations]," max: ",max(profitHistory), sep=""), ylim = c(yl[1],yl[2]*1.1))
}



totalScoreHistoryAllTogether <-function(path,yl){
  maxX = 0;
  for(j in 1:length(estimators)){
    files <- list.files(path= paste(outPath, estimators[j], sep = ""), pattern="Toy.*.csv", full.names=T, recursive=FALSE)
    if(maxX < length(files)){
      maxX = length(files)
    }
  }
  
  colvec = c()
  for(j in 1:length(estimators)){
    files <- list.files(path= paste(outPath, estimators[j], sep = ""), pattern="Toy.*.csv", full.names=T, recursive=FALSE)
    iterations = length(files)
    
    print(iterations)
    
    path=paste(outPath, estimators[j], "/", sep = "") 
    
    scoreHistory = rep(NA,iterations)
    profitHistory <- numeric(iterations)
    
    for(i in 1:iterations){
      adv = read.csv2(paste(path,"Toy",i,".csv" ,sep=""))
      
      
      adv$totalProfit <- as.numeric(as.character(adv$totalProfit))
      
      profitHistory[i] <- adv$totalProfit[1]
      
      #adv$totalScore <- as.numeric(as.character(adv$totalScore))        # <----------------
      
      #profitHistory[i] <- adv$totalScore[1]        # <----------------
    }
    
    
    if(1 == j){
      #plot(c(1:iterations),profitHistory, xlab="Iteration", ylab = "Profit", type = "l", 
      #   main = paste("Final Profit: ", profitHistory[iterations]," max: ",max(profitHistory), sep=""),xlim = c(0,maxX),ylim = c(yl[1],yl[2]*1.1), col = j+1)
      
      plot(c(1:iterations),profitHistory, xlab="Iteration", ylab = "Score", type = "l", 
           xlim = c(0,maxX),ylim = c(yl[1],yl[2]*1.1), col = j+1)
    } else {
      par(new = TRUE)
      plot(c(1:iterations), profitHistory, axes = FALSE, xlim = c(0,maxX),ylim = c(yl[1],yl[2]*1.1), xlab = "", ylab = "", type = "l", col = j+1)
    }
    colvec = c(colvec,j+1)
  }
  
  legend('topright', estimators, lty=c(1,1), lwd=c(2.5,2.5),col=colvec)        # <----------------------
  
  mtext("Profit", side=2, line=2,outer = FALSE)
  #mtext("Profit", side=2, line=2,outer = FALSE)
  mtext("Iteration", side=1, line=2,outer = FALSE, cex = 0.5)
}

dieFrequency <-function(path2files, it) {
  dieHistory = rep(NA,it)
  for(i in 1:it){
    adv = read.csv2(paste(path2files,"Toy",i,".csv" ,sep=""))
    dieHistory[i] = levels(adv$best[1])
    #dieHistory[i] = adv$best[1]
  }
  t = table(dieHistory)
  
  barplot(t)
}

totalObservations <-function(path2files, it) {
  sim = read.csv2(paste(path2files,"Toy",it,".csv" ,sep=""))
  n2 = scalar1(sim$totalObservations)
  barplot(n2,ylim=c(0,1.0))
}

plotAll <-function(path2files, estimator,ylimMax) {
  files <- list.files(path= paste(outPath, estimator, sep = ""), pattern="Toy.*.csv", full.names=T, recursive=FALSE)
 # starts with 0, so the last file will have index length(files)-1
  #it = length(files)-1
  
  it = length(files)
  
  path=paste(outPath, estimator, "/", sep = "") 
  
  dieChoiceHistory(path,it)
  totalScoreHistory(path,it,ylimMax)
  dieFrequency(path,it)
  totalObservations(path,it)
}

plotPHidden <-function(outPath) {
  jpeg(paste(outPath,"pHidden.jpg", sep = ""))
  files <- list.files(path= outPath, pattern="die.*.csv", full.names=T, recursive=FALSE)
  # starts with 0, so the last file will have index length(files)-1
  it = length(files)
  
  par(mfcol=c(2,it/2), mar=c(0.5,0.5,0.5,0.5) + c(4,3,4,3))
  for(i in 1:it){
    mydata = read.csv2(paste(outPath,"die", i, ".csv", sep=""), header = FALSE)
    mydata
    mydata$V2
    #  r = as.matrix(mydata)
    r = scalar1(mydata$V2)
    barplot(r, cex.main = 3, main = paste("die",i,sep=""), ylim = c(0,1))
    #mtext(paste("die",i,sep=""), side=4, line=1, las=1) 
    #mtext(text=mydata$best[1], side=4, line=1, las=1) 
  }
  dev.off()
}

#outPath = "/home/willy/Bachelorarbeit/Manager/test32/"

#plotPHidden(outPath)

#jpeg(paste(outPath,"Comparison.jpg",sep=""))
#pdf(paste(outPath,"Comparison.pdf",sep=""))

#print(length(estimators))
#par(mfcol=c(4,length(estimators)),mar = c(2,2,2,2), oma = c(1,1,3,1))

# looks in all estimators and finds the max score of all
#ylimMax = findMaxScore(outPath,estimators);

#title =""
#for(i in 1:length(estimators)) {
#  print(paste("plotting estimator ",estimators[i], sep=""))
  
#  plotAll(outPath, estimators[i],ylimMax)
#  title = paste(title,estimators[i],sep="              ")
#}
#mtext(title, side=3, line=1, las=1,outer = TRUE)

#dev.off()


plotAll2 <-function(path2files, estimator,ylimMax,i) {
  files <- list.files(path= paste(outPath, estimator, sep = ""), pattern="Toy.*.csv", full.names=T, recursive=FALSE)
  # starts with 0, so the last file will have index length(files)-1
  #it = length(files)-1
  
  it = length(files)
  
  path=paste(outPath, estimator, "/", sep = "") 
  
  dieChoiceHistory(path,it)
  
  if(1 == i) mtext("Die Choice", side=2, line=2,outer = FALSE)
  
  mtext(estimator, side=3, line=2,outer = FALSE)

  #totalScoreHistory(path,it,ylimMax)
  dieFrequency(path,it)
  if(1 == i) mtext("Die Frequency", side=2, line=2,outer = FALSE)
  
  totalObservations(path,it)
  if(1 == i) mtext("Total Observations", side=2, line=2,outer = FALSE)
}

pdf(paste(outPath,"ComparisonNEW.pdf",sep=""))



print(length(estimators))
par(mar = c(2,2,2,2), oma = c(1,0.1,3,1))
#layout(matrix(c(1,2,3,4,5,6,7,8,9,10,10,10), 4, 3, byrow = TRUE))

if(1 == length(estimators)){
  layout(matrix(c(1, 2, 3, 4, 4), 5, 1, byrow = TRUE))
}else if(2 == length(estimators)){
  layout(matrix(c(1,4, 2,5, 3,6, 7,7, 7,7), 5, 2, byrow = TRUE))
}else if(3 == length(estimators)){
  layout(matrix(c(1,4,7, 2,5,8, 3,6,9, 10,10,10, 10,10,10), 5, 3, byrow = TRUE))
} else if(4 == length(estimators)) {
  layout(matrix(c(1,4,7,10, 2,5,8,11, 3,6,9,12, 13,13,13,13, 13,13,13,13), 5, 4, byrow = TRUE))
}
par(mar = c(2,3,2,2), oma = c(1,1,3,1))
  

# looks in all estimators and finds the max score of all
ylimMax = findMaxScore(outPath,estimators);

for(i in 1:length(estimators)) {
  print(paste("plotting estimator ",estimators[i], sep=""))
  
  plotAll2(outPath, estimators[i],ylimMax,i)
  #title = paste(title,estimators[i],sep="              ")
  

}

totalScoreHistoryAllTogether(outPath,ylimMax)

#mtext(title, side=3, line=1, las=1,outer = TRUE)


dev.off()
