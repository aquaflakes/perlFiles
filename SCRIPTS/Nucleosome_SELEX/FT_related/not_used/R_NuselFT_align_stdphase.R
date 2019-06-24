library(sfsmisc)

setwd("./")
allfiles<- list.files(path="./ADM", pattern="*.txt", full.names=FALSE, recursive=FALSE)#read in all txt and fastq files
# files in different subfolder of the same name

filecnt<-length(allfiles)
#colour of lines
col1<-c("#000000","#333333","#666666","#888888","#AAAAAA","#0000FF","#2222FF","#5555FF","#8888FF","#BBBBFF","#FF0000","#FF2222","#FF5555","#FF8888","#FFBBBB","#00FF00","#22FF22","#55FF55","#88FF88","#BBFFBB","#000000")
pdf(file="./ntFreq_plots.pdf", title="A vs position", paper="a4",  width = 0, height = 0) #width and height is necessary for positioning
newpg_cnt<-0
for (i in 1:(length(allfiles))) 
{
	data <- read.table(paste("./ADM/",allfiles[i],sep=""), header=FALSE, sep="\t", row.names=1) #read in data from txt files
	dat <- as.list(as.data.frame(t(data)))# chage data frame to list
	keys <- names(dat) 
	xkey<-keys[1]
	keys<-keys[-1] # remove the key of the first caption line


	#options(scipen=1) 
	if (newpg_cnt %% 4 == 0) {par(bg=NA, mfrow = c(4,3), mar=c(3.2,3,1,1))} #place 3 figure on 1 page
	newpg_cnt<-newpg_cnt+1
	par(xpd=T) #display legend out of plot area
	# plot(unlist(dat),type="n",xlim=c(0,max(sapply(dat,length))),xlab="position /bp",ylab="nt count")
	#title(allfiles[i])
	#ymax<-par('usr') #get max of y axis for the next plot
	#mapply(lines,dat,lty=1,col=col1) #col=seq_along(dat)
	#legend(par()$usr[2], par()$usr[4],keys,lty=1,col=col1)#c("#000000","#333333","#666666","#888888","#AAAAAA","#0000FF","#2222FF","#5555FF","#8888FF","#BBBBFF","#FF0000","#FF2222","#FF5555","#FF8888","#FFBBBB","#00FF00","#22FF22","#55FF55","#88FF88","#BBFFBB")) #col=seq_along(dat)
	#enlarge
	#par(xpd=F) #limit the plot out of plot area
	# plot(unlist(dat),type="n",xlim=c(1,max(sapply(dat,length))),ylim=c(0,0.33*ymax[4]),xlab="position /bp",ylab="nt count",font.main = 1, main = "(enlarge of dinucleotide counts)")
	# mapply(lines,dat,col=col1,lty=1)

	
	plot(dat[[xkey]],dat$TA,type="n",col=col1[1],ylim=c(0,0.3*21), xlim=c(0,100),xlab="position /bp",ylab="",
		yaxt='n',cex.lab=1.2, cex.axis=1, cex.main=1, cex.sub=1.1, bty="n",mgp=c(2,1,0))
	title(ylab = "Normalized Nucleotide Freq", mgp = c(0, 1, 0), cex.lab=1.2) # do it seperately so only affect Y axis

	colCnt<-1
	shift<-1
	for (key in keys){ 
		lines(dat[[xkey]],dat[[key]]*3+shift*0.3,col=col1[colCnt]) #scale by 3 to view easier
		points(100,(0.27+(shift-1)*0.3), col=col1[colCnt], pch=key)	
		points(104,(0.27+(shift-1)*0.3), col=col1[colCnt], pch=substr(key,2,3))		
		colCnt<-colCnt+1
		shift<-shift+1
	}
	title(allfiles[i],cex.main=1.4)
	#legend(par()$usr[2], par()$usr[4],keys,lty=1,col=col1)



	# also read in 1 file of FT result
	data <- read.table(paste("./FT/",allfiles[i],sep=""), header=FALSE, sep="\t", row.names=1) #read in data from txt files
	dat <- as.list(as.data.frame(t(data)))# chage data frame to list
	keys <- names(dat) 
	xkey<-keys[1]  # data of the caption line to x keys
	keys<-keys[-1] # remove the caption line, now keys contains only the data keys

	plot(dat[[xkey]],dat$TA,type="n",col=col1[1],ylim=c(0,0.08*21),
		main="after FT",xlab="spacial frequency /bp-1",mgp=c(2,1,0),
		ylab="" ,yaxt='n',bty="n",
		cex.lab=1.2, cex.axis=1, cex.main=1.1, cex.sub=1.1)
	title(ylab = "Power Spectra Density", mgp = c(0, 1, 0), cex.lab=1.2) # do it seperately so only affect Y axis
	colCnt<-1
	for (key in keys){ 
		lines(dat[[xkey]],dat[[key]]*3+colCnt*0.08,col=col1[colCnt]) #scale by 3 to view easier
		points(0.31,(0.085+(colCnt-1)*0.08), col=col1[colCnt], pch=key)	
		points(0.322,(0.085+(colCnt-1)*0.08), col=col1[colCnt], pch=substr(key,2,3))	
		colCnt<-colCnt+1
	}
	

	library(plotrix)
	phase <- as.matrix(read.table(paste("./phaseAng/",allfiles[i],sep=""), header=TRUE, sep = "\t", row.names = 1, as.is=TRUE))
	peakArea <- as.matrix(read.table(paste("./Area/",allfiles[i],sep=""), header=TRUE, sep = "\t", row.names = 1, as.is=TRUE))
	 	
	polar.plot(peakArea[,1]-0.05*max(peakArea[,1]), phase[,1]*180, # original input is -1 to 1
		radial.labels=NA, rp.type="s",
				radial.lim=c(0,max(peakArea[,1])),lwd=2,line.col=col1,point.col=col1, 
				point.symbols=substr(keys,2,3),#substring the second char of key
				 cex=1
				)
	par(new=TRUE)
	polar.plot(peakArea[,1]+0.05*max(peakArea[,1]), phase[,1]*180, # original input is -1 to 1
		main="amplitude and phase (10.2bp signal)",rp.type="s",radial.labels=NA,
				radial.lim=c(0,max(peakArea[,1])),lwd=2,line.col=col1,point.col=col1, 
				point.symbols=keys, #default the first char of key
				cex=1.7
				)	
	# x scales
	text(0,-par()$usr[2]*0.1,"0", cex=1.5, col="#008B9B")
	text(par()$usr[2],-par()$usr[2]*0.1,round(max(peakArea[,1]),2), cex=1.5, col="#008B9B")

}

dev.off()
