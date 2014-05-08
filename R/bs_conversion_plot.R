#!/usr/bin/RScript

library(reshape)
library(ggplot2)

## EXPECTS

in.file=
out.dir=

args<-commandArgs(TRUE)
if(length(args) < 2){
  cat("Error incorrect number of args","\n",sep="")
  q()
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}

#THRESH<-0.95
tab<-read.table(file=in.file,sep="\t",stringsAsFactors=FALSE)
names(tab)<-c('gene','seq','read_count','well')
tab$nT<-nchar(gsub("C","",tab$seq)) * tab$read_count
tab$nBase<-nchar(tab$seq) * tab$read_count

##for each gene and well combination calculate the conversion efficiency and store total reads

summ.tab<-do.call("rbind",lapply(split(tab,tab$well),function(w){
	do.call("rbind",lapply(split(w,w$gene),function(g){
		tmp<-list()
		tmp['well']<-unique(g$well)
		tmp['gene']<-unique(g$gene)
		tmp['total_reads']<-sum(g$read_count)
		tmp['conversion']<-sum(g$nT)/sum(g$nBase)
		#tmp['status']<-ifelse(tmp['conversion'] > THRESH,'PASS','FAIL')
		as.data.frame(tmp)
		}))
}))

summ.tab$plate<-as.integer(sub("^([^_]+)\\_.*$","\\1",summ.tab$well))
summ.tab$row<-sub("[^_]+\\_([A-Z]+)[0-9]+$","\\1",summ.tab$well)
summ.tab$col<-as.integer(sub("[^_]+\\_[A-Z]+([0-9]+)$","\\1",summ.tab$well))
summ.tab$col <- factor(summ.tab$col, levels = sort(unique(summ.tab$col)))
summ.tab<-summ.tab[order(summ.tab$plate,summ.tab$col,summ.tab$row),]
summ.tab$well<-factor(summ.tab$well,level=unique(summ.tab$well))

summ.tab<-summ.tab[summ.tab$total_reads>10,]

plot.df<-split(summ.tab,summ.tab$plate)
ofile=ofile=paste(out.dir,"bs_conversion.pdf",sep="/")
pdf(ofile,onefile=TRUE,height=8,width=14)
theme_set(theme_bw())
sapply(seq_along(plot.df),function(i){
	p<-plot.df[[i]]
	gg<-ggplot(p,aes(x=gene,y=conversion)) + geom_bar(position="dodge",stat="identity") 
	gg<-gg + geom_text(position=position_dodge(width=0.9),color="black",size=3,aes(label=total_reads,y=0.95),angle=-90) 
	gg<-gg + xlab("Percentage conversion") + guides(fill=guide_legend(title="95% Conversion")) + coord_cartesian(xlim=NULL,ylim=c(0.9,1))
	gg<-gg + facet_grid(row ~ col)
	gg<-gg + labs(title=paste("Plate",names(plot.df)[i])) + theme(axis.text.x=element_text(angle = -90, hjust = 0))
	print(gg)
	})
dev.off()
print("Success")
