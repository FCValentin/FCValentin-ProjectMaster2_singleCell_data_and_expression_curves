###-------------###
#Home functions####
###-------------###

## Create a new directory
createDir<-function(file){ 
  if(!file.exists(file)){
    dir.create(file)
  }
}

# read a file
lire<-function(x, character=FALSE){
  if(character){
    d<-read.table(file = x,sep = "\t",header=T,row.names = 1,colClasses = "character",quote="")
  }else{
    d<-read.table(file = x,sep = "\t",header=T,row.names = 1,quote="")
  }
  return(d)
}

# write a file
ecrire<-function(x,file="default.tsv",headRow="Name",row.names=TRUE,col.names=TRUE){
  options(warn=-1) #Supress unecessary warning about append 
  if(row.names && col.names){
    write.table(x = paste0(headRow,"\t"),file = file,sep = "\t",eol="",quote=F,row.names=F,col.names=F)
    write.table(x=x,file=file,sep="\t", row.names = T, col.names = T, quote = FALSE,append=T)
  }else{
    write.table(x=x,file=file,sep="\t", row.names = row.names, col.names = col.names, quote = FALSE)
  }
  options(warn=0)
}


controlData<-function(sampleAnnot){
  #Lineage must be factor and organised by a continiuous factor levels (used number or letter to classify factor segments in time)
  #Pseudotime must be numeric
  if(class(sampleAnnot$Lineage)!="factor") sampleAnnot$Lineage<-as.factor(sampleAnnot$Lineage)
  if(class(sampleAnnot$Pseudotime)!="numeric") sampleAnnot$Pseudotime<-as.numeric(as.character(sampleAnnot$Pseudotime))
  return(sampleAnnot)
}

###-----------------###
#Automatised script####
###-----------------###


SmoothPipeline<-function(expr,sampleAnnot,TreeModel,n=100,span=0.75,method=median,errorType=FALSE){
  
  sampleAnnot<-controlData(sampleAnnot)
  #Fate model formatter
  Fate<-strsplit(as.character(TreeModel$Lineages),",")
  names(Fate)<-row.names(TreeModel)
  for(i in 1:nrow(TreeModel)){
    Fate[[row.names(TreeModel)[i]]]<-as.factor(Fate[[row.names(TreeModel)[i]]])
  }
  
  
  ### Sample selection by fate
  sample<-list()
  for(i in names(Fate)){
    sample[[i]]<-sampleAnnot[which(sampleAnnot$Lineage%in%Fate[[i]]),]
  }
  
  
  ## NewSamplesPseudotime
  Newsample<-list()
  for(s in names(sample)){
    NewSegments<-list()
    begin<-T
    Lineage<-c()
    for(i in levels(as.factor(as.character(sample[[s]]$Lineage)))){
      if(begin){
        NewSegments[[i]]<-round(seq(min(sampleAnnot$Pseudotime[which(sampleAnnot$Lineage==i)]),max(sampleAnnot$Pseudotime[which(sampleAnnot$Lineage==i)]),length.out=n),2)
        begin<-F
      }else{
        NewSegments[[i]]<-round(seq(max(sampleAnnot$Pseudotime[which(sampleAnnot$Lineage==last)]+0.01),(max(sampleAnnot$Pseudotime[which(sampleAnnot$Lineage==i)])-0.05),length.out=n),2)
      } 
      last<-i
      Lineage<-c(Lineage,rep(i,n))
    }
    Newsample[[s]]<-data.frame(unlist(NewSegments),Lineage)
    colnames(Newsample[[s]])[1]<-c("Pseudotime")
  }
  
  ## NewExpression
  NewExpression<-list()
  NewExpressionSd<-list()
  for(s in names(sample)){
    dataSegment<-as.matrix(expr[,row.names(sample[[s]])])
    NewExpression[[s]]<-as.data.frame(t(apply(dataSegment[,], 1, 
                                              function(x){ model<-loess(x~sample[[s]]$Pseudotime,span=span,degree=2)
                                              return(predict(model,newdata=as.numeric(as.character(Newsample[[s]]$Pseudotime))))           
                                              })))
    NewExpressionSd[[s]]<-as.data.frame(t(apply(dataSegment[,], 1, 
                                                function(x){ model<-loess(x~sample[[s]]$Pseudotime,span=span,degree=2)
                                                return(predict(model,newdata=as.numeric(as.character(Newsample[[s]]$Pseudotime)),se=T)$se)           
                                                })))
    #if(errorType) NewExpressionSd[[s]]<-NewExpressionSd[[s]]/sqrt(n)
    colnames(NewExpression[[s]])<-colnames(NewExpressionSd[[s]])<-row.names(Newsample[[s]])
    row.names(NewExpression[[s]])<-row.names(NewExpressionSd[[s]])<-row.names(expr)
  }
  
  
  ###------------------###
  #Export data by fate####
  ###------------------###
  
  createDir("results")
  for(s in names(sample)){
    createDir(paste("results/",s,sep=""))
    ecrire(NewExpression[[s]],paste(paste("results/",s,sep=""),"/ExpressionSample.tsv",sep=""))
    ecrire(NewExpressionSd[[s]],paste(paste("results/",s,sep=""),"/ExpressionSdSample.tsv",sep=""))
    ecrire(Newsample[[s]],paste(paste("results/",s,sep=""),"/SampleAnnotation.tsv",sep=""))
  }
  
  ###-----------------###
  #Aggregate segments####
  ###-----------------###
  
  #Keep samples specific to an unique segment
  UniqueSample<-t(data.frame(row.names=c("Pseudotime","Lineage")))
  UniqueExpr<-data.frame(row.names=row.names(expr))
  UniqueExprSd<-data.frame(row.names=row.names(expr))
  for(s in names(Newsample)){
    UniqueSample<-rbind(UniqueSample,Newsample[[s]][which(Newsample[[s]]$Lineage%in%TreeModel$leef),])
    UniqueExpr<-cbind(UniqueExpr,NewExpression[[s]][row.names(Newsample[[s]][which(Newsample[[s]]$Lineage%in%TreeModel$leef),])])
    UniqueExprSd<-cbind(UniqueExprSd,NewExpressionSd[[s]][row.names(Newsample[[s]][which(Newsample[[s]]$Lineage%in%TreeModel$leef),])])
  }
  
  # Find redundant samples between fate
  redundantSample<-Newsample
  redundantExpression<-NewExpression
  redundantExpressionSd<-NewExpressionSd
  for(s in names(redundantExpression)){
    redundantSample[[s]]<-redundantSample[[s]][which(!row.names(redundantSample[[s]])%in%row.names(UniqueSample)),]
    redundantExpression[[s]]<-redundantExpression[[s]][,which(!colnames(redundantExpression[[s]])%in%row.names(UniqueSample))]
    redundantExpressionSd[[s]]<-redundantExpressionSd[[s]][,which(!colnames(redundantExpressionSd[[s]])%in%row.names(UniqueSample))]
  }
  
  # Calculate mean and sd expr of redundant samples
  CommonSample<-t(data.frame(row.names=c("Pseudotime","Lineage")))
  CommonExpr<-data.frame(row.names=row.names(expr))
  CommonExprSd<-data.frame(row.names=row.names(expr))
  CommonFittedSd<-data.frame(row.names=row.names(expr))
  j<-1
  for(l in levels(sampleAnnot$Lineage)[!levels(sampleAnnot$Lineage)%in%TreeModel$leef]){
    SegmentExpr<-list()
    SegmentExprSd<-list()
    segments<-c()
    for(s in names(Fate)){
      if(l%in%Fate[[s]])  segments<-c(segments,s)
    }
    for(s in segments){
      SegmentExpr[[s]]<-data.frame(row.names=row.names(expr))
      SegmentExprSd[[s]]<-data.frame(row.names=row.names(expr))
      SegmentExpr[[s]]<-redundantExpression[[s]][,row.names(redundantSample[[s]][which(redundantSample[[s]]$Lineage%in%l),])]
      SegmentExprSd[[s]]<-redundantExpressionSd[[s]][,row.names(redundantSample[[s]][which(redundantSample[[s]]$Lineage%in%l),])]
    }
    for(i in 1:n){
      CommonExpr<-cbind(CommonExpr,apply(sapply(SegmentExpr,"[[",i),1,method))
      CommonExprSd<-cbind(CommonExprSd,apply(sapply(SegmentExprSd,"[[",i),1,method))
      CommonFittedSd<-cbind(CommonFittedSd,apply(sapply(SegmentExpr,"[[",i),1,sd))
    }
    CommonSample<-rbind(CommonSample,redundantSample[[s]][which(redundantSample[[s]]$Lineage%in%l),])
    colnames(CommonExpr)[j:(j+(n-1))]<-colnames(CommonExprSd)[j:(j+(n-1))]<-colnames(SegmentExpr[[s]])
    colnames(CommonFittedSd)[j:(j+(n-1))]<-colnames(SegmentExpr[[s]])
    j<-j+n
  }
  
  ## export sd and expr of new samples
  UniqueFittedSd<-UniqueExpr
  UniqueFittedSd[,]<-0
  FinalExpr<-cbind(CommonExpr,UniqueExpr)
  FinalExprSd<-cbind(CommonExprSd,UniqueExprSd)
  FinalFittedSd<-cbind(CommonFittedSd,UniqueFittedSd)
  FinalSample<-rbind(CommonSample,UniqueSample)
  createDir(paste("results/","SegmentFusion",sep=""))
  ecrire(FinalExpr,paste(paste("results/","SegmentFusion",sep=""),"/ExpressionSample.tsv",sep=""))
  ecrire(FinalExprSd,paste(paste("results/","SegmentFusion",sep=""),"/ExpressionSdSample.tsv",sep=""))
  ecrire(FinalFittedSd,paste(paste("results/","SegmentFusion",sep=""),"/FittedSdSample.tsv",sep=""))
  ecrire(FinalSample,paste(paste("results/","SegmentFusion",sep=""),"/SampleAnnotation.tsv",sep=""))
}

