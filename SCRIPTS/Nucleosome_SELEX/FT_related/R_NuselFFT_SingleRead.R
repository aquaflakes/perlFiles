#library(sfsmisc)

setwd("./")
allfiles<- list.files(path="./FFT_SR/", pattern="*.txt", full.names=FALSE, recursive=FALSE)#read in all txt and fastq files
# files in different subfolder of the same name

filecnt<-length(allfiles)
#colour of lines
col1<-c("#000000","#333333","#666666","#888888","#AAAAAA","#0000FF","#2222FF","#5555FF","#8888FF","#BBBBFF","#FF0000","#FF2222","#FF5555","#FF8888","#FFBBBB","#00FF00","#22FF22","#55FF55","#88FF88","#BBFFBB")
pdf(file="./FFT_SR/FFT_plots.pdf", title="A vs position", paper="a4",  width = 0, height = 0) #width and height is necessary for positioning
newpg_cnt<-0
for (i in 1:(length(allfiles))) 
{
  if (newpg_cnt %% 12 == 0) {par(bg=NA, mfrow = c(4,3), mar=c(3.2,3,1,1))} #place 3 figure on 1 page
  newpg_cnt<-newpg_cnt+1 
  
  FFTresult1<- as.matrix(read.table(paste("./FFT_SR/",allfiles[i],sep=""), header=TRUE, sep = "\t", as.is=TRUE))
  FFTresult<- FFTresult1[1:length(FFTresult1[,1])/2,]  # only take half of the result set due to symmetry
  keys<- colnames(FFTresult)[2:length(colnames(FFTresult))]
  plot(FFTresult[,1][1:length(FFTresult[,1])/2],FFTresult[,2][1:length(FFTresult[,2])/2],xlim = c(0,0.55),
       ylim=c(2.5,9.5),
       type="n",
       xlab="spacial frequency /bp-1",mgp=c(2,1,0),
       ylab="" ,yaxt='n',bty="n",
       cex.lab=1.2, cex.axis=1, cex.main=0.9, cex.sub=1.1)
  title(ylab = "Power Spectra Density", mgp = c(0, 1, 0), cex.lab=1.2) # do it seperately so only affect Y axis)
  shift<-1
  for (key in keys){ 
    lines(FFTresult[,"keys"],(FFTresult[,key]-mean(FFTresult[,key]))/mean(FFTresult[,key])*6+3+(shift-1)*0.3, # normalize counts by dividing mean
          col=col1[shift]) #scale by 3 to view easier
    points(0.5,(2.6+(shift-1)*0.3), col=col1[shift], pch=key)	
    points(0.52,(2.6+(shift-1)*0.3), col=col1[shift], pch=substr(key,2,3))		
    shift<-shift+1
  }
  

  title(allfiles[i],cex.main=1);

}

dev.off()
