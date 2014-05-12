#!/usr/bin/RScript

library(reshape)
library(ggplot2)

## GLOBAL VARIABLES

READ_COUNT_THRESH=300
BOWTIE_STEP_NAME='Bowtie'
PDF_FILE_NAME="gene_dropout.pdf"

## EXPECTS

## in.file=
## out.dir=

args<-commandArgs(TRUE)
if(length(args) < 2){
  cat("Error incorrect number of args","\n",sep="")
  q()
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}

## MAIN

tab<-read.table(file=in.file,sep="\t",stringsAsFactors=FALSE)
names(tab)<-c('gene','step','read_count','well')
tab$step<-factor(tab$step,levels=unique(tab$step))
tab<-tab[tab$well!="NOT_CALLED",]

## which well gene combos 

summ.tab<-do.call("rbind",lapply(split(tab,tab$well),function(w){
	do.call("rbind",lapply(split(w,w$gene),function(g){
		## assume demux is reference count
		ref_count<-g[g$step==levels(g$step)[1],]$read_count
		bt_count<-g[g$step==BOWTIE_STEP_NAME,]$read_count
		if(length(bt_count)==0)
			bt_count=0
		if(length(ref_count)==0)
			ref_count=0
		if(ref_count > READ_COUNT_THRESH){
			p.df<-do.call("rbind",lapply(levels(g$step)[-1],function(step){
				prop_drop<- 1-(g[g$step==step,]$read_count/ref_count)
				if(length(prop_drop)==0)
						prop_drop=1/length(levels(g$step)[-1])
				if(prop_drop<0)
						prop_drop=0
				data.frame(step=step,pro.dropped=prop_drop)
			}))
			p.df$gene<-unique(g$gene)
			p.df$well<-unique(g$well)
			p.df$total_reads<-character(length=nrow(p.df))
			p.df[1,]$total_reads<-bt_count
			rownames(p.df)<-NULL
			p.df
		}
	}))
}))
		

## for ease 

summ.tab$step<-factor(summ.tab$step,levels=unique(summ.tab$step))

summ.tab$plate<-as.integer(sub("^([^_]+)\\_.*$","\\1",summ.tab$well))
summ.tab$row<-sub("[^_]+\\_([A-Z]+)[0-9]+$","\\1",summ.tab$well)
summ.tab$col<-as.integer(sub("[^_]+\\_[A-Z]+([0-9]+)$","\\1",summ.tab$well))
summ.tab$col <- factor(summ.tab$col, levels = sort(unique(summ.tab$col)))
summ.tab<-summ.tab[order(summ.tab$plate,summ.tab$col,summ.tab$row),]
summ.tab$well<-factor(summ.tab$well,level=unique(summ.tab$well))

plot.df<-split(summ.tab,summ.tab$plate)
## there is only one of these per experiment so 
## name can be hardcoded
ofile=paste(out.dir,PDF_FILE_NAME,sep="/")
pdf(ofile,onefile=TRUE,height=8,width=14)
theme_set(theme_bw())
sapply(seq_along(plot.df),function(i){
	p<-plot.df[[i]]
	gg<-ggplot(p,aes(x=gene,y=pro.dropped,fill=step)) + geom_bar(stat="identity") 
	gg<-gg + geom_text(hjust=2,position="stack",color="black",size=3,aes(label=total_reads),angle=-90)
	gg<-gg + coord_cartesian(ylim=c(0,0.6))
	gg<-gg + facet_grid(row ~ col)
	gg<-gg + labs(title=paste("Plate",names(plot.df)[i])) + theme(axis.text.x=element_text(angle = -90, hjust = 0))
	print(gg)
	})
dev.off()

