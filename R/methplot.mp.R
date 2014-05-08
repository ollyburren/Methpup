#!/usr/bin/RScript

library(ggplot2)
library(reshape)

## EXPECTS
#methplotfile
#mp.file
#where to write plot
#out.dir


args<-commandArgs(TRUE)
if(length(args) < 2){
  cat("Error incorrect number of args","\n",sep="")
  q()
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}


## subroutines

mode <- function(x){
	table_data<-table(x) 
	as.integer(names(subset(table_data, table_data==max(table_data))))
}


## this code is adapted from Xin Yang's methplot R package
## http://github.com/XinYang6699/Methplot
strReverse <- function(x){
	sapply(lapply(strsplit(x, NULL), rev), paste,collapse="")
}


## this code is adapted from Xin Yang's methplot R package
## http://github.com/XinYang6699/Methplot
plotdata<- function(x,x.title, legendpos="null") {
  x$nT <- nchar(gsub("C","",x$seq))
  x$prev <- strReverse(x$seq)
  x <- x[order(x$nT,x$prev),]
  n <- mode(nchar(x$seq))
  ## remove all those that don't have sequences that conform to the mode
  x<-x[nchar(x$seq)==n,]
  ss<-strsplit(as.character(x$seq),"")
  ss <- t(as.data.frame(ss))
  colnames(ss) <- sprintf("s%s",1:n)
  rownames(ss) <- NULL
  x <- cbind(x,ss)
  x$y2 <- cumsum(x$read_count)
  x$y2 <- 100*x$y2/max(x$y2)
  x$y1 <- c(0,x$y2[-nrow(x)])
  df <- lapply(1:n, function(i) {
    data.frame(x1=i-1, x2=i, y1=x$y1, y2=x$y2, base=x[,sprintf("s%s",i)])
  })
  df <- do.call("rbind",df)
  x2 <- melt(x[,c("nT","y1","y2")], id="nT")
   ggplot(df) +
    geom_rect(aes(xmin=x1,xmax=x2,ymin=y1,ymax=y2,fill=base)) +
    geom_text(aes(label=x2,x=x1+0.5,y=-2),size=4,angle=90,alpha=0.2,fontface=1) +
    scale_fill_brewer(palette="Greens") + xlab("Position") +
    ylab(paste("Proportion (total reads = ", sum(x$read_count),")", sep="")) +
    labs(title=x.title) + theme_set(theme_bw()) +
    scale_x_discrete(breaks=1:n, labels=sprintf("s%s", 1:n), expand=c(0.05, -0.4))+
    theme(legend.position=legendpos, axis.text.x=element_blank())

}

## read in file

mp.data<-read.table(file=mp.file,sep="\t",stringsAsFactors=FALSE)
names(mp.data)<-c('gene','seq','read_count','well')
## split by well (just in case)
lapply(split(mp.data,mp.data$well),function(w){
	lapply(split(w,w$gene),function(g){
		fname=paste('methplot',unique(g$well),unique(g$gene),'pdf',sep=".")
		ofile=paste(out.dir,fname,sep="/")
		p.title<-paste("Methylation Plot",paste(unique(g$well),unique(g$gene),sep=":"))
		pdf(ofile)
		gg <- plotdata(g, x.title=p.title,legendpos="right")
		print(gg)
		dev.off()
	})
})

print("Success")


