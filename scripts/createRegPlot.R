require(ggplot2)

setwd("~/Documents/cours/Labo/ValPipelinePosfai")

###----------------###
#LOADING FUNCTIONS####
###----------------###

RequiredfunctionsDir<-"requiredFiles/"
source(paste0(RequiredfunctionsDir,"/SmoothRegressionPipeline.R"))


###-----------###
#LOADING DATA####
###-----------###

sampleAnnot<-lire("data/sampleAnnot.tsv")
expr<-lire("data/exprDat.norm.tsv")
TreeModel<-lire("data/SegmentsPosfai.tsv")
n<-100 # Cell Per segment
span<-0.75
method<-median #median, mean... usual function to fitted common segments expression
errorType<-T #F : sd / T : standard error (se)
  
###---------------###
#Script launching####
###---------------###

SmoothPipeline(expr,sampleAnnot,TreeModel,n,span,method,errorType) 
## Results are in results directory
# One directory by fate model
# One directory with a unique model based on common segments fusion


sampleFinal<-lire("results/SegmentFusion/SampleAnnotation.tsv")
exprFinal<-lire("results/SegmentFusion/ExpressionSample.tsv")
exprFinal[exprFinal<0]<-0
exprsdFinal<-lire("results/SegmentFusion/ExpressionSdSample.tsv")
sdmin<-exprFinal[,]-exprsdFinal[,]
sdmax<-exprFinal[,]+exprsdFinal[,]

genes<-row.names(exprFinal)
markers<-c("Nanog","Sox2","Pou5f1","Klf2","Gata6","Sox17","Pgdfra","Gata4","Cdx2","Gata2","Gata3","Klf6")

for(i in markers){
	sd<-as.data.frame(t(rbind(unlist(sdmin[i,]),unlist(sdmax[i,]))))
	colnames(sd)<-c("Min","Max")
	sd$Min[sd$Min<0]<-0
	# svg(filename = paste0("figs/",i,".svg"),width = 5,height = 3,bg = "transparent")
	print(ggplot(as.data.frame(t(exprFinal)), 
							 aes(x=as.numeric(as.character(sampleFinal$Pseudotime)), 
							 		y=as.numeric(unlist(exprFinal[i,])),color=sampleFinal$Lineage)) +  
					geom_ribbon(data = as.data.frame(sd),
											aes(ymin=unlist(sd$Min),ymax=sd$Max,x=sampleFinal$Pseudotime,fill=sampleFinal$Lineage), 
											alpha = 0.3,inherit.aes = F,color=NA)+
					geom_line(size=1.5)+
					scale_color_manual(values=c("red","green","blue"))+
					scale_fill_manual(values=c("red","green","blue"))+
					#scale_x_continuous(limits=xlim, expand = c(0,0))+
					#scale_y_continuous(limits=c(0,NA), expand = c(0,0))+
					guides(colour=FALSE,fill=FALSE)+
					theme(panel.background = element_rect(fill = "#EEEEEE",colour="black"),
								panel.grid.major = element_line(colour = NA),
								panel.grid.minor = element_line(colour = NA),
								# axis.title.x=element_text(),
								# axis.title.y=element_text(),
								plot.background = element_rect(fill = "transparent",colour=NA)
					)+
					xlab("Pseudotime")+
					ylab(paste0(i, " fitted expression"))
	)
	# dev.off()
}
