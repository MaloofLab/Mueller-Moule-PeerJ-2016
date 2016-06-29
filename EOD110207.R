get.excel2 <- function(filename)
  #function to convert a column based spreadsheet of measurements into a dataframe.
  #the first column should consist of grouping variable names, corresponding to the contents
  #of each row.  The first row of data should have a row name of "Data"
  {
     vars <- as.character(read.table(file=filename,sep=",",fill=TRUE,blank.lines.skip=FALSE)[,1])
     data.start <- match("Data", vars)
     data.table <- read.table(file=filename,sep=",",skip=data.start-1,fill=TRUE)[,-1]
     last <- ifelse(match("",vars)>data.start,data.start-1,match("",vars)-1) #where is the last grouping variable? 
     header <- read.table(file=filename, sep=",",fill=TRUE,as.is=TRUE)[1:last,-1]
     time.frame <- as.data.frame(data.table[,header[1,]=="time"] )#extract the time information
     colnames(time.frame) <- header[5,header[1,]=="time"]
     data.table <- data.table[,header[1,]!="time"] #remove time information from data table, leaving hyp info
     header <- header[,header[1,]!="time"] #remove time information from header
     data.na <- is.na(data.table)
     time.table <- time.frame[,match(header[5,],colnames(time.frame))] # expand time.frame so that there is one column per plant
     tmp.frame <- data.frame(line=col(data.na)[!data.na],
     						 time=time.table[!data.na],
     						 y=data.table[!data.na])

      for (i in (1:last))
       tmp.frame[,vars[i]] <- factor(unlist(header[i,])[tmp.frame$line])
     return(tmp.frame)
   }
   
data <- get.excel2("EODFR110207.csv")

data <- data[data$time>32,] #first time point is high

summary(data)

library(lattice)

xyplot(y~time|treat_gt,data=data,type="l",lwd=2,groups=Plant)

dr5.mean <- tapply(data$y,list(data$time,data$treat,data$gt),mean)

dr5.sem <- tapply(data$y,list(data$time,data$treat,data$gt),function(x) sd(x)/sqrt(length(x)))

pdf("col and yuccaQ DR5.pdf")#,width=10.5,height=8)

#par(mfrow=c(1,2))

op <- par(mar=par()$mar+.5)

ylim=c(min(dr5.mean-dr5.sem)*.8,max(dr5.mean+dr5.sem)*1.2)

matplot(x=as.numeric(rownames(dr5.mean)),y=cbind(dr5.mean[,,"Col"]+dr5.sem[,,"Col"],dr5.mean[,,"Col"]-dr5.sem[,,"Col"]),type="l",lwd=2,lty=3,col=c("red","black"),
		xlab="time after FR (min)",ylab="DR5::LUC Bioluminescence",main="Col",cex.lab=1.5,cex.main=1.5,xlim=c(0,1000),xaxp=c(0,1000,10),ylim=ylim)

matlines(x=as.numeric(rownames(dr5.mean[,,"Col"])),y=dr5.mean[,,"Col"],lwd=2,lty=2:1,col=c("red","black"))

legend("topright",legend=c("control","EOD-FR"),lwd=2,lty=1:2,col=c("black","red"))

matplot(x=as.numeric(rownames(dr5.mean)),y=cbind(dr5.mean[,,"yuccaQ"]+dr5.sem[,,"yuccaQ"],dr5.mean[,,"yuccaQ"]-dr5.sem[,,"yuccaQ"]),type="l",lwd=2,lty=3,col=c("red","black"),
		xlab="time after FR (min)",ylab="DR5::LUC Bioluminescence",main="yuccaQ",cex.lab=1.5,cex.main=1.5,xlim=c(0,1000),xaxp=c(0,1000,10),ylim=ylim)

matlines(x=as.numeric(rownames(dr5.mean[,,"yuccaQ"])),y=dr5.mean[,,"yuccaQ"],lwd=2,lty=2:1,col=c("red","black"))

legend("topright",legend=c("control","EOD-FR"),lwd=2,lty=1:2,col=c("black","red"))


dev.off()
