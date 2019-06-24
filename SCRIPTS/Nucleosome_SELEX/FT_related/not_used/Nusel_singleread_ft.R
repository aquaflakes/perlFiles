library(sfsmisc)

#setwd("/.")
allfiles<- list.files(path="./", pattern="*.txt", full.names=FALSE, recursive=FALSE)#read in all txt and fastq files
# files in different subfolder of the same name

filecnt<-length(allfiles)
#colour of lines
col1<-c("#000000","#333333","#666666","#888888","#AAAAAA","#0000FF","#2222FF","#5555FF","#8888FF","#BBBBFF","#FF0000","#FF2222","#FF5555","#FF8888","#FFBBBB","#00FF00","#22FF22","#55FF55","#88FF88","#BBFFBB")
pdf(file="./ntFreq_plots.pdf", title="A vs position", paper="a4",  width = 0, height = 0) #width and height is necessary for positioning
newpg_cnt<-0
for (i in 1:(length(allfiles))) 
{
	data <- read.table(paste("./",allfiles[i],sep=""), header=FALSE, sep="\t", row.names=1) #read in data from txt files
	dat <- as.list(as.data.frame(t(data)))# chage data frame to list
	keys <- names(dat) 

	if (newpg_cnt %% 4 == 0) {par(mfrow = c(4,1), mar=c(4,4,4,5))} #place 3 figure on 1 page
	newpg_cnt<-newpg_cnt+1
	peakheight <- as.matrix(read.table(paste("./",allfiles[i],sep=""), header=TRUE, sep = "\t", row.names = 1, as.is=TRUE))
	#peakArea <- as.matrix(read.table(paste("./Area/",allfiles[i],sep=""), header=TRUE, sep = "\t", row.names = 1, as.is=TRUE))
	barplot(peakheight[,1],col=col1,ylim = c(-0.02,0.13)) 	
	title(allfiles[i])

}

dev.off()
