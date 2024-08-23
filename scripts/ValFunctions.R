####-----------------###
### Complex HeatMap ####
####-----------------###

ComplexHeatMapLineage<-function(dataHeatMap,clusteringRowsH2,sampleAnnotBis,colHt,KinomeGroup){
  require(ComplexHeatmap)
  EchantMarqueursMorulaB1<-sampleAnnotBis[which(sampleAnnotBis$Branches == "3.Morula" | sampleAnnotBis$Branches == "4.B1"),]
  EchantMarqueursICM<-sampleAnnotBis[which(sampleAnnotBis$Branches == "5.Unspecified ICM"),]
  EchantMarqueursSpéPE<-sampleAnnotBis[which(sampleAnnotBis$Branches == "7.Primitive endorderm"),]
  EchantMarqueursSpéEPI<-sampleAnnotBis[which(sampleAnnotBis$Branches == "6.Epiblast"),]
  EchantMarqueursSpéTE<-sampleAnnotBis[which(sampleAnnotBis$Branches == "8.early TE" | sampleAnnotBis$Branches == "9.Trophectoderm"),]
  haMorulaB1<-HeatmapAnnotation(EchantMarqueursMorulaB1[,c(9,11)],col = list(State=c("25"="red","24"="blue","23"="black","4"="green"),
                                                                             Branches=c("3.Morula"="black","4.B1"="grey")))
  HeatMapMorulaB1<-Heatmap(dataHeatMap[,which(colnames(dataHeatMap)%in%rownames(EchantMarqueursMorulaB1))],name="MorulaB1",
                           row_title = "Gènes", column_title = "MorulaB1",row_names_gp = gpar(cex=0.1),
                           column_names_gp = gpar(cex=0.5),show_heatmap_legend=T,show_row_names=T,
                           show_column_names = F,col=colHt,top_annotation = haMorulaB1,column_dend_reorder=F,cluster_columns=F,
                           cluster_rows = clusteringRowsH2)
  
  haICM<-HeatmapAnnotation(EchantMarqueursICM[,c(9,11)],col = list(State=c("5"="cyan"),
                                                                   Branches=c("5.Unspecified ICM"="orange")))
  HeatMapICM<-Heatmap(dataHeatMap[,which(colnames(dataHeatMap)%in%rownames(EchantMarqueursICM))],name="ICM",
                      row_title = "Gènes", column_title = "ICM",row_names_gp = gpar(cex=0.1),
                      column_names_gp = gpar(cex=0.5),show_heatmap_legend=T,show_row_names=T,
                      show_column_names = F,col=colHt,top_annotation = haICM,column_dend_reorder=F,cluster_columns=F,
                      cluster_rows = clusteringRowsH2,width = unit(1, "cm"))
  
  haEPI<-HeatmapAnnotation(EchantMarqueursSpéEPI[,c(9,11)],col = list(State=c("6"="magenta","7"="pink","8"="orange"),
                                                                      Branches=c("6.Epiblast"="red")))
  HeatMapEPI<-Heatmap(dataHeatMap[,which(colnames(dataHeatMap)%in%rownames(EchantMarqueursSpéEPI))],name="EPI",
                      row_title = "Gènes", column_title = "EPI",row_names_gp = gpar(cex=0.1),
                      column_names_gp = gpar(cex=0.5),show_heatmap_legend=T,show_row_names=T,
                      show_column_names = F,col=colHt,top_annotation = haEPI,column_dend_reorder=F,cluster_columns=F,
                      cluster_rows = clusteringRowsH2)
  
  haPE<-HeatmapAnnotation(EchantMarqueursSpéPE[,c(9,11)],col = list(State=c("9"="magenta"),
                                                                    Branches=c("7.Primitive endorderm"="green")))
  HeatMapPE<-Heatmap(dataHeatMap[,which(colnames(dataHeatMap)%in%rownames(EchantMarqueursSpéPE))],name="PE",
                     row_title = "Gènes", column_title = "PE",row_names_gp = gpar(cex=0.1),
                     column_names_gp = gpar(cex=0.5),show_heatmap_legend=T,show_row_names=T,
                     show_column_names = F,col=colHt,top_annotation = haPE,column_dend_reorder=F,cluster_columns=F,
                     cluster_rows = clusteringRowsH2,width = unit(1, "cm"))
  
  haTE<-HeatmapAnnotation(EchantMarqueursSpéTE[,c(9,11)],col = list(State=c("22"="green","21"="cyan","20"="magenta","19"="pink","18"="orange","17"="yellow","16"="darkblue","15"="grey","13"="purple","12"="gold","11"="darkorange","10"="tomato"),
                                                                    Branches=c("8.early TE"="cyan","9.Trophectoderm"="blue")))
  HeatMapTE<-Heatmap(dataHeatMap[,which(colnames(dataHeatMap)%in%rownames(EchantMarqueursSpéTE))],name="TE",
                     row_title = "Gènes", column_title = "TE",row_names_gp = gpar(cex=0.1),
                     column_names_gp = gpar(cex=0.5),show_heatmap_legend=T,show_row_names=T,
                     show_column_names = F,col=colHt,top_annotation = haTE,column_dend_reorder=F,cluster_columns=F,
                     cluster_rows = clusteringRowsH2)
  
  HeatMap<-HeatMapMorulaB1+HeatMapICM+HeatMapEPI+HeatMapPE+HeatMapTE
  pdf(paste(paste("fig/HeatMapComplex",KinomeGroup,sep="_"),".pdf",sep=""),width=10,height=20)
  draw(HeatMap, row_title = paste("Genes Clustering",KinomeGroup,sep=" "), row_title_gp = gpar(col = "red"),
       column_title = "Sample Pseudotime", column_title_side = "bottom",gap = unit(0.2, "cm"))
  dev.off()
  ecrire(dataHeatMap,paste(paste("results/exprHeatMapComplex",KinomeGroup,sep="_"),".tsv",sep=""))
  ecrire(cutree(clusteringRowsH2,best.cutree(clusteringRowsH2)),paste(paste("results/BestCutreeClusterMapComplex",KinomeGroup,sep="_"),".tsv",sep=""))
}

ComplexHeatMapCluster<-function(dataHeatMap,sampleAnnotBis,colHt,KinomeGroup,split,...){
  require(ComplexHeatmap) # Package
  # Sample Selection
  EchantMarqueursMorulaB1<-sampleAnnotBis[which(sampleAnnotBis$Branches == "3.Morula" | sampleAnnotBis$Branches == "4.B1"),]
  EchantMarqueursICM<-sampleAnnotBis[which(sampleAnnotBis$Branches == "5.Unspecified ICM"),]
  EchantMarqueursSpéPE<-sampleAnnotBis[which(sampleAnnotBis$Branches == "7.Primitive endorderm"),]
  EchantMarqueursSpéEPI<-sampleAnnotBis[which(sampleAnnotBis$Branches == "6.Epiblast"),]
  EchantMarqueursSpéTE<-sampleAnnotBis[which(sampleAnnotBis$Branches == "8.early TE" | sampleAnnotBis$Branches == "9.Trophectoderm"),]
  
  #Annot and HeatMap
  haMorulaB1<-HeatmapAnnotation(EchantMarqueursMorulaB1[,c(9,11)],col = list(State=c("25"="red","24"="blue","23"="black","4"="green"),
                                                                             Branches=c("3.Morula"="black","4.B1"="grey")))
  dataHeatMap[,which(colnames(dataHeatMap)%in%rownames(EchantMarqueursMorulaB1))]
  HeatMapMorulaB1<-Heatmap(dataHeatMap[which(rn(dataHeatMap)%in%rn(dataHeatMap)[which(apply(dataHeatMap,1,sd)!=0)]),which(colnames(dataHeatMap)%in%rownames(EchantMarqueursMorulaB1))],name="MorulaB1",
                           row_title = "Gènes", column_title = "MorulaB1",row_names_gp = gpar(cex=0.1),
                           column_names_gp = gpar(cex=0.5),show_heatmap_legend=T,show_row_names=T,
                           show_column_names = F,col=colHt,top_annotation = haMorulaB1,column_dend_reorder=F,cluster_columns=F,
                           clustering_distance_rows = "pearson",clustering_method_rows="ward.D2",split=split, gap = unit(5, "mm"))
  
  haICM<-HeatmapAnnotation(EchantMarqueursICM[,c(9,11)],col = list(State=c("5"="cyan"),
                                                                   Branches=c("5.Unspecified ICM"="orange")))
  dataHeatMap[,which(colnames(dataHeatMap)%in%rownames(EchantMarqueursICM))]
  HeatMapICM<-Heatmap(dataHeatMap[which(rn(dataHeatMap)%in%rn(dataHeatMap)[which(apply(dataHeatMap,1,sd)!=0)]),which(colnames(dataHeatMap)%in%rownames(EchantMarqueursICM))],name="ICM",
                      row_title = "Gènes", column_title = "ICM",row_names_gp = gpar(cex=0.1),
                      column_names_gp = gpar(cex=0.5),show_heatmap_legend=T,show_row_names=T,
                      show_column_names = F,col=colHt,top_annotation = haICM,column_dend_reorder=F,cluster_columns=F,
                      clustering_distance_rows = "pearson",clustering_method_rows="ward.D2",width = unit(1, "cm"),split=split, gap = unit(5, "mm"),...)
  
  haEPI<-HeatmapAnnotation(EchantMarqueursSpéEPI[,c(9,11)],col = list(State=c("6"="magenta","7"="pink","8"="orange"),
                                                                      Branches=c("6.Epiblast"="red")))
  dataHeatMap[,which(colnames(dataHeatMap)%in%rownames(EchantMarqueursSpéEPI))]
  HeatMapEPI<-Heatmap(dataHeatMap[which(rn(dataHeatMap)%in%rn(dataHeatMap)[which(apply(dataHeatMap,1,sd)!=0)]),which(colnames(dataHeatMap)%in%rownames(EchantMarqueursSpéEPI))],name="EPI",
                      row_title = "Gènes", column_title = "EPI",row_names_gp = gpar(cex=0.1),
                      column_names_gp = gpar(cex=0.5),show_heatmap_legend=T,show_row_names=T,
                      show_column_names = F,col=colHt,top_annotation = haEPI,column_dend_reorder=F,cluster_columns=F,
                      clustering_distance_rows = "pearson",clustering_method_rows="ward.D2",split=split, gap = unit(5, "mm"))
  
  haPE<-HeatmapAnnotation(EchantMarqueursSpéPE[,c(9,11)],col = list(State=c("9"="magenta"),
                                                                    Branches=c("7.Primitive endorderm"="green")))
  dataHeatMap[,which(colnames(dataHeatMap)%in%rownames(EchantMarqueursSpéPE))]
  HeatMapPE<-Heatmap(dataHeatMap[which(rn(dataHeatMap)%in%rn(dataHeatMap)[which(apply(dataHeatMap,1,sd)!=0)]),which(colnames(dataHeatMap)%in%rownames(EchantMarqueursSpéPE))],name="PE",
                     row_title = "Gènes", column_title = "PE",row_names_gp = gpar(cex=0.1),
                     column_names_gp = gpar(cex=0.5),show_heatmap_legend=T,show_row_names=T,
                     show_column_names = F,col=colHt,top_annotation = haPE,column_dend_reorder=F,cluster_columns=F,
                     clustering_distance_rows = "pearson",clustering_method_rows="ward.D2",width = unit(1, "cm"),split=split, gap = unit(5, "mm"))
  
  haTE<-HeatmapAnnotation(EchantMarqueursSpéTE[,c(9,11)],col = list(State=c("22"="green","21"="cyan","20"="magenta","19"="pink","18"="orange","17"="yellow","16"="darkblue","15"="grey","13"="purple","12"="gold","11"="darkorange","10"="tomato"),
                                                                    Branches=c("8.early TE"="cyan","9.Trophectoderm"="blue")))
  dataHeatMap[,which(colnames(dataHeatMap)%in%rownames(EchantMarqueursSpéTE))]
  HeatMapTE<-Heatmap(dataHeatMap[which(rn(dataHeatMap)%in%rn(dataHeatMap)[which(apply(dataHeatMap,1,sd)!=0)]),which(colnames(dataHeatMap)%in%rownames(EchantMarqueursSpéTE))],name="TE",
                     row_title = "Gènes", column_title = "TE",row_names_gp = gpar(cex=0.1),
                     column_names_gp = gpar(cex=0.5),show_heatmap_legend=T,show_row_names=T,
                     show_column_names = F,col=colHt,top_annotation = haTE,column_dend_reorder=F,cluster_columns=F,
                     clustering_distance_rows = "pearson",clustering_method_rows="ward.D2",split=split, gap = unit(5, "mm"))
  
  HeatMap<-HeatMapMorulaB1+HeatMapICM+HeatMapEPI+HeatMapPE+HeatMapTE
  pdf(paste(paste("fig/HeatMapComplex",KinomeGroup,sep="_"),".pdf",sep=""))
  draw(HeatMap, row_title = paste("Genes Clustering",KinomeGroup,sep=" "), row_title_gp = gpar(col = "red"),
       column_title = "Sample Pseudotime", column_title_side = "bottom",gap = unit(0.2, "cm"))
  dev.off()
  ecrire(dataHeatMap,paste(paste("results/exprHeatMapComplex",KinomeGroup,sep="_"),".tsv",sep=""))
}


SmoothComplexHeatMap<-function(mat,sampleAnnotSmooth,colHt,GeneGroup,clusteringRowsH2=TRUE,split=NULL){
  ## mat : expression matrix
  ## sampleAnnotSmooth : sampleAnnot with Branches
  ## colHt : color for HeatMap
  ## GeneGroup : Group of gene
  ## split : How to split genes
  require(circlize)
  require(ComplexHeatmap)
  colHt<-colorRamp2(c(-2,0,2),c("blue","white","red"))
  clustering_method_rows<-NULL
  clustering_distance_rows<-NULL
  if(!is.null(split)){
    if(nrow(mat)!=length(split)){
      stop("Error in length of split : They are more or less genes in split than in expression matrix")
    }
    # We can't provided a single manual clustering with split
    clusteringRowsH2<-TRUE
    clustering_method_rows<-"ward.D2"
    clustering_distance_rows<-"minkowski"
  }
  require(ComplexHeatmap) # Package
  
  # Sample Selection
  EchantMarqueursMorulaB1<-sampleAnnotSmooth[which(sampleAnnotSmooth$Lineage == "3.Morula4.B1"),]
  EchantMarqueursICM<-sampleAnnotSmooth[which(sampleAnnotSmooth$Lineage == "5. Unspecified ICM"),]
  EchantMarqueursSpéPE<-sampleAnnotSmooth[which(sampleAnnotSmooth$Lineage == "7.Primitive endoderm"),]
  EchantMarqueursSpéEPI<-sampleAnnotSmooth[which(sampleAnnotSmooth$Lineage == "6.Epiblast"),]
  EchantMarqueursSpéTE<-sampleAnnotSmooth[which(sampleAnnotSmooth$Lineage == "9.Trophectoderm"),]
  
  #Annot and HeatMap
  annot<-as.data.frame(EchantMarqueursMorulaB1[,2])
  colnames(annot)<-"Lineage"
  haMorulaB1<-HeatmapAnnotation(annot,col = list(Lineage=c("3.Morula4.B1"="black")))
  HeatMapMorulaB1<-Heatmap(mat[,which(colnames(mat)%in%rownames(EchantMarqueursMorulaB1))],name="MorulaB1",
                           row_title = "Gènes", column_title = "Morula",row_names_gp = gpar(cex=0.1),
                           column_names_gp = gpar(cex=0.2),column_title_gp=gpar(cex=0.5),show_heatmap_legend=T,show_row_names=T,
                           show_column_names = F,col=colHt,top_annotation = haMorulaB1,column_dend_reorder=F,cluster_columns=F,
                           cluster_rows=clusteringRowsH2,clustering_distance_rows=clustering_distance_rows,clustering_method_rows=clustering_method_rows,split=split, gap = unit(5, "mm"))
  
  annot<-as.data.frame(EchantMarqueursICM[,2])
  colnames(annot)<-"Lineage"
  haICM<-HeatmapAnnotation(annot,col = list(Lineage=c("5. Unspecified ICM"="yellow")))
  HeatMapMorulaICM<-Heatmap(mat[,which(colnames(mat)%in%rownames(EchantMarqueursICM))],name="ICM",
                           row_title = "Gènes", column_title = "ICM",row_names_gp = gpar(cex=0.1),
                           column_names_gp = gpar(cex=0.2),column_title_gp=gpar(cex=0.5),show_heatmap_legend=T,show_row_names=T,
                           show_column_names = F,col=colHt,top_annotation = haICM,column_dend_reorder=F,cluster_columns=F,
                           cluster_rows=clusteringRowsH2,clustering_distance_rows=clustering_distance_rows,clustering_method_rows=clustering_method_rows,split=split, gap = unit(5, "mm"))
  
  
  
  annot<-as.data.frame(EchantMarqueursSpéEPI[,2])
  colnames(annot)<-"Lineage"
  haEPI<-HeatmapAnnotation(annot,col = list(Lineage=c("6.Epiblast"="red")))
  HeatMapEPI<-Heatmap(mat[,which(colnames(mat)%in%rownames(EchantMarqueursSpéEPI))],name="EPI",
                      row_title = "Gènes", column_title = "EPI",row_names_gp = gpar(cex=0.1),
                      column_names_gp = gpar(cex=0.2),column_title_gp=gpar(cex=0.5),show_heatmap_legend=T,show_row_names=T,
                      show_column_names = F,col=colHt,top_annotation = haEPI,column_dend_reorder=F,cluster_columns=F,
                      cluster_rows=clusteringRowsH2,clustering_distance_rows=clustering_distance_rows,clustering_method_rows=clustering_method_rows,split=split, gap = unit(5, "mm"))
  
  annot<-as.data.frame(EchantMarqueursSpéPE[,2])
  colnames(annot)<-"Lineage"
  haPE<-HeatmapAnnotation(annot,col = list(Lineage=c("7.Primitive endoderm"="green")))
  HeatMapPE<-Heatmap(mat[,which(colnames(mat)%in%rownames(EchantMarqueursSpéPE))],name="PE",
                     row_title = "Gènes", column_title = "PE",row_names_gp = gpar(cex=0.1),
                     column_names_gp = gpar(cex=0.2),column_title_gp=gpar(cex=0.5),show_heatmap_legend=T,show_row_names=T,
                     show_column_names = F,col=colHt,top_annotation = haPE,column_dend_reorder=F,cluster_columns=F,
                     cluster_rows=clusteringRowsH2,clustering_distance_rows=clustering_distance_rows,clustering_method_rows=clustering_method_rows,width = unit(1, "cm"),split=split, gap = unit(5, "mm"))
  
  annot<-as.data.frame(EchantMarqueursSpéTE[,2])
  colnames(annot)<-"Lineage"
  haTE<-HeatmapAnnotation(annot,col = list(Lineage=c("9.Trophectoderm"="blue")))
  HeatMapTE<-Heatmap(mat[,which(colnames(mat)%in%rownames(EchantMarqueursSpéTE))],name="TE",
                     row_title = "Gènes", column_title = "TE",column_title_gp=gpar(cex=0.5),row_names_gp = gpar(cex=0.1),
                     column_names_gp = gpar(cex=0.2),show_heatmap_legend=T,show_row_names=T,
                     show_column_names = F,col=colHt,top_annotation = haTE,column_dend_reorder=F,cluster_columns=F,
                     cluster_rows=clusteringRowsH2,clustering_distance_rows=clustering_distance_rows,clustering_method_rows=clustering_method_rows,split=split, gap = unit(5, "mm"))
  
  HeatMap<-HeatMapMorulaB1+HeatMapMorulaICM+HeatMapEPI+HeatMapPE+HeatMapTE
  pdf(paste(paste("fig/HeatMapComplex",GeneGroup,sep="_"),".pdf",sep=""),height = 10)
  draw(HeatMap, row_title = paste("Genes Clustering",GeneGroup,sep=" "), row_title_gp = gpar(col = "red"),
       column_title = "Sample Pseudotime", column_title_side = "bottom",gap = unit(0.2, "cm"))
  dev.off()
  # ecrire(mat,paste(paste("results/exprHeatMapComplex",GeneGroup,sep="_"),".tsv",sep=""))
}


SmoothComplexHeatMap2<-function(mat,sampleAnnotSmooth,colHt,GeneGroup,clusteringRowsH2=TRUE,split=NULL){
  ## mat : expression matrix
  ## sampleAnnotSmooth : sampleAnnot with Branches
  ## colHt : color for HeatMap
  ## GeneGroup : Group of gene
  ## split : How to split genes
  require(circlize)
  require(ComplexHeatmap)
  colHt<-colorRamp2(c(-2,0,2),c("blue","white","red"))
  clustering_method_rows<-NULL
  clustering_distance_rows<-NULL
  if(!is.null(split)){
    if(nrow(mat)!=length(split)){
      stop("Error in length of split : They are more or less genes in split than in expression matrix")
    }
    # We can't provided a single manual clustering with split
    clusteringRowsH2<-TRUE
    clustering_method_rows<-"ward.D2"
    clustering_distance_rows<-corrDist
  }
  require(ComplexHeatmap) # Package
  
  # Sample Selection
  EchantMarqueursSpéPE<-sampleAnnotSmooth[which(sampleAnnotSmooth$Lineage == "3.Early blastocyst"|sampleAnnotSmooth$Lineage == "1.Pre-morula"|sampleAnnotSmooth$Lineage == "2.Morula"|sampleAnnotSmooth$Lineage == "4.Inner cell mass"|sampleAnnotSmooth$Lineage == "7.Primitive endoderm"),]
  EchantMarqueursSpéEPI<-sampleAnnotSmooth[which(sampleAnnotSmooth$Lineage == "3.Early blastocyst"|sampleAnnotSmooth$Lineage == "1.Pre-morula"|sampleAnnotSmooth$Lineage == "2.Morula"|sampleAnnotSmooth$Lineage == "4.Inner cell mass"|sampleAnnotSmooth$Lineage == "6.Epiblast"),]
  EchantMarqueursSpéTE<-sampleAnnotSmooth[which(sampleAnnotSmooth$Lineage == "3.Early blastocyst"|sampleAnnotSmooth$Lineage == "1.Pre-morula"|sampleAnnotSmooth$Lineage == "2.Morula"|sampleAnnotSmooth$Lineage == "8.Trophectoderm"|sampleAnnotSmooth$Lineage == "5.Early trophectoderm"),]
  
  #Annot and HeatMap
  
  annot<-as.data.frame(EchantMarqueursSpéEPI[,"Lineage"])
  colnames(annot)<-"Lineage"
  haEPI<-HeatmapAnnotation(annot,col = list(Lineage=c( "1.Pre-morula"="white","3.Early blastocyst"="#f3a0b5","2.Morula"= "#808080","4.Inner cell mass"="#f3e646","6.Epiblast"="#e40522")))
  HeatMapEPI<-Heatmap(mat[,which(colnames(mat)%in%rownames(EchantMarqueursSpéEPI))],name="EPI",
                      row_title = "Gènes", column_title = "EPI",row_names_gp = gpar(cex=0.3),row_title_rot=0,
                      column_names_gp = gpar(cex=0.2),column_title_gp=gpar(cex=0.5),show_heatmap_legend=T,show_row_names=F,
                      show_column_names = F,col=colHt,top_annotation = haEPI,column_dend_reorder=F,cluster_columns=F,
                      cluster_rows=clusteringRowsH2,clustering_distance_rows=clustering_distance_rows,clustering_method_rows=clustering_method_rows,split=split, gap = unit(5, "mm"))
  
  annot<-as.data.frame(EchantMarqueursSpéPE[,"Lineage"])
  colnames(annot)<-"Lineage"
  haPE<-HeatmapAnnotation(annot,col = list(Lineage=c("1.Pre-morula"="white","3.Early blastocyst"="#f3a0b5","2.Morula"= "#808080","4.Inner cell mass"="#f3e646","7.Primitive endoderm"="#66b32e")))
  HeatMapPE<-Heatmap(mat[,which(colnames(mat)%in%rownames(EchantMarqueursSpéPE))],name="PE",
                     row_title = "Gènes", column_title = "PE",row_names_gp = gpar(cex=0.3),
                     column_names_gp = gpar(cex=0.2),column_title_gp=gpar(cex=0.5),show_heatmap_legend=T,show_row_names=F,
                     show_column_names = F,col=colHt,top_annotation = haPE,column_dend_reorder=F,cluster_columns=F,
                     cluster_rows=clusteringRowsH2,clustering_distance_rows=clustering_distance_rows,clustering_method_rows=clustering_method_rows,width = unit(1, "cm"),split=split, gap = unit(5, "mm"))

  annot<-as.data.frame(EchantMarqueursSpéTE[,"Lineage"])
  colnames(annot)<-"Lineage"
  haTE<-HeatmapAnnotation(annot,col = list(Lineage=c("1.Pre-morula"="white","3.Early blastocyst"="#f3a0b5","2.Morula"= "#808080","5.Early trophectoderm"="#3fb4e8","8.Trophectoderm"="#354999")))
  HeatMapTE<-Heatmap(mat[,which(colnames(mat)%in%rownames(EchantMarqueursSpéTE))],name="TE",
                     row_title = "Gènes", column_title = "TE",column_title_gp=gpar(cex=0.5),row_names_gp = gpar(cex=0.3),
                     column_names_gp = gpar(cex=0.2),show_heatmap_legend=T,show_row_names=T,
                     show_column_names = F,col=colHt,top_annotation = haTE,column_dend_reorder=F,cluster_columns=F,
                     cluster_rows=clusteringRowsH2,clustering_distance_rows=clustering_distance_rows,clustering_method_rows=clustering_method_rows,split=split, gap = unit(5, "mm"))
  
  HeatMap<-HeatMapEPI+HeatMapPE+HeatMapTE
  pdf(paste(paste("fig/HeatMapComplex",GeneGroup,sep="_"),".pdf",sep=""),height = 15)
  draw(HeatMap, row_title = "", row_title_gp = gpar(col = "red",fontsize(5)),
       column_title = "Sample Pseudotime", column_title_side = "bottom",gap = unit(0.2, "cm"))
  dev.off()
  # ecrire(mat,paste(paste("results/exprHeatMapComplex",GeneGroup,sep="_"),".tsv",sep=""))
}


SmoothComplexHeatMapSpeLineage<-function(mat,sampleAnnotSmooth,colHt,GeneGroup,clusteringRowsH2=TRUE,split=NULL){
  ## mat : expression matrix
  ## sampleAnnotSmooth : sampleAnnot with Branches, without Morula!
  ## colHt : color for HeatMap
  ## GeneGroup : Group of gene
  ## split : How to split genes
  require(circlize)
  require(ComplexHeatmap)
  colHt<-colorRamp2(c(-2,0,2),c("blue","white","red"))
  clustering_method_rows<-NULL
  clustering_distance_rows<-NULL
  if(!is.null(split)){
    if(nrow(mat)!=length(split)){
      stop("Error in length of split : They are more or less genes in split than in expression matrix")
    }
    # We can't provided a single manual clustering with split
    clusteringRowsH2<-TRUE
    clustering_method_rows<-"ward.D2"
    clustering_distance_rows<-"pearson"
  }
  require(ComplexHeatmap) # Package
  
  # Sample Selection
  EchantMarqueursSpéPE<-sampleAnnotSmooth[which(sampleAnnotSmooth$Lineage == "4.Inner cell mass"|sampleAnnotSmooth$Lineage == "7.Primitive endoderm"),]
  EchantMarqueursSpéEPI<-sampleAnnotSmooth[which(sampleAnnotSmooth$Lineage == "4.Inner cell mass"|sampleAnnotSmooth$Lineage == "6.Epiblast"),]
  EchantMarqueursSpéTE<-sampleAnnotSmooth[which(sampleAnnotSmooth$Lineage == "8.Trophectoderm"|sampleAnnotSmooth$Lineage == "5.Early trophectoderm"),]
  
  #Annot and HeatMap
  
  annot<-as.data.frame(EchantMarqueursSpéEPI[,2])
  colnames(annot)<-"Lineage"
  haEPI<-HeatmapAnnotation(annot,col = list(Lineage=c("4.Inner cell mass"="#f3e646","6.Epiblast"="#e40522")))
  HeatMapEPI<-Heatmap(mat[,which(colnames(mat)%in%rownames(EchantMarqueursSpéEPI))],name="EPI",
                      row_title = "Gènes", column_title = "EPI",row_names_gp = gpar(cex=0.1),row_title_rot=0,
                      column_names_gp = gpar(cex=0.2),column_title_gp=gpar(cex=0.5),show_heatmap_legend=T,show_row_names=T,
                      show_column_names = F,col=colHt,top_annotation = haEPI,column_dend_reorder=F,cluster_columns=F,
                      cluster_rows=clusteringRowsH2,clustering_distance_rows=clustering_distance_rows,clustering_method_rows=clustering_method_rows,split=split, gap = unit(5, "mm"))
  
  annot<-as.data.frame(EchantMarqueursSpéPE[,2])
  colnames(annot)<-"Lineage"
  haPE<-HeatmapAnnotation(annot,col = list(Lineage=c("4.Inner cell mass"="#f3e646","7.Primitive endoderm"="#66b32e")))
  HeatMapPE<-Heatmap(mat[,which(colnames(mat)%in%rownames(EchantMarqueursSpéPE))],name="PE",
                     row_title = "Gènes", column_title = "PE",row_names_gp = gpar(cex=0.1),
                     column_names_gp = gpar(cex=0.2),column_title_gp=gpar(cex=0.5),show_heatmap_legend=T,show_row_names=T,
                     show_column_names = F,col=colHt,top_annotation = haPE,column_dend_reorder=F,cluster_columns=F,
                     cluster_rows=clusteringRowsH2,clustering_distance_rows=clustering_distance_rows,clustering_method_rows=clustering_method_rows,width = unit(1, "cm"),split=split, gap = unit(5, "mm"))
  
  annot<-as.data.frame(EchantMarqueursSpéTE[,2])
  colnames(annot)<-"Lineage"
  haTE<-HeatmapAnnotation(annot,col = list(Lineage=c("5.Early trophectoderm"="#3fb4e8","8.Trophectoderm"="#354999")))
  HeatMapTE<-Heatmap(mat[,which(colnames(mat)%in%rownames(EchantMarqueursSpéTE))],name="TE",
                     row_title = "Gènes", column_title = "TE",column_title_gp=gpar(cex=0.5),row_names_gp = gpar(cex=0.1),
                     column_names_gp = gpar(cex=0.2),show_heatmap_legend=T,show_row_names=T,
                     show_column_names = F,col=colHt,top_annotation = haTE,column_dend_reorder=F,cluster_columns=F,
                     cluster_rows=clusteringRowsH2,clustering_distance_rows=clustering_distance_rows,clustering_method_rows=clustering_method_rows,split=split, gap = unit(5, "mm"))
  
  HeatMap<-HeatMapEPI+HeatMapPE+HeatMapTE
  pdf(paste(paste("fig/HeatMapComplex",GeneGroup,sep="_"),".pdf",sep=""),height = 15)
  draw(HeatMap, row_title = paste("Genes Clustering",GeneGroup,sep=" "), row_title_gp = gpar(col = "red",fontsize(5)),
       column_title = "Sample Pseudotime", column_title_side = "bottom",gap = unit(0.2, "cm"))
  dev.off()
  # ecrire(mat,paste(paste("results/exprHeatMapComplex",GeneGroup,sep="_"),".tsv",sep=""))
}


###----------------------------------------###
####  MOUSE Gene to HUMAN Gene converter  ####
###----------------------------------------###

convertMouseGeneList <- function(x){
  
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}

# 
# dataGenes<-lire("results/ComplexHeatMap/GenesLineage/exprHeatMapComplex.tsv")
# plot(y=exprDat.SC.Norm.log["NANOG",],x=sampleAnnotBis$Pseudotime,col=sampleAnnotBis$Color,xlab="Pseudotime",ylab="Expression")
# legend("bottomright", legend=levels(as.factor(as.character(sampleAnnotBis$Branches))),col=c("black","orange","red","green","cyan","blue"), lty=1:2, cex=0.8, text.font=4, bg='lightblue')

###----------------------------------------###
#### PLOT EXPRESSION GENE HEATMAP LINEAGE ####
###----------------------------------------###

PlotExprGene<-function(exprDat.SC.Norm.log,gene,EchantMarqueursEPI,EchantMarqueursPE,EchantMarqueursTE){
  require(ggplot2)
  require(cowplot)
  EPI<-ggplot(as.data.frame(exprDat.SC.Norm.log[gene,which(colnames(exprDat.SC.Norm.log)%in%rownames(EchantMarqueursEPI))]), aes(x=EchantMarqueursEPI$Pseudotime, y=exprDat.SC.Norm.log[gene,which(colnames(exprDat.SC.Norm.log)%in%rownames(EchantMarqueursEPI))],color=as.factor(as.character(EchantMarqueursEPI$BranchesVal)))) + 
    geom_point(shape=18)+ 
    ylim(0,20)+
    labs(title = "EPI\n", x = "Pseudotime", y = paste("Log2 Expr Gene",gene,sep=" "), color = "Lineage\n")+
    scale_color_manual(values=c("black","orange","red"))+
    geom_smooth(method=loess, se=T, linetype="solid",color="red",fullrange=T)+
    theme(legend.position = "none",axis.title.x=element_blank())
  PE<-ggplot(as.data.frame(exprDat.SC.Norm.log[gene,which(colnames(exprDat.SC.Norm.log)%in%rownames(EchantMarqueursPE))]), aes(x=EchantMarqueursPE$Pseudotime, y=exprDat.SC.Norm.log[gene,which(colnames(exprDat.SC.Norm.log)%in%rownames(EchantMarqueursPE))],color=as.factor(as.character(EchantMarqueursPE$BranchesVal)))) + 
    geom_point(shape=18)+ 
    ylim(0,20)+
    labs(title = "PE\n", x = "Pseudotime", y = "Log2 Expr Gene", color = "Lineage\n")+
    scale_color_manual(values=c("black","orange","green"))+
    geom_smooth(method=loess, se=T, linetype="solid",color="green",fullrange=T)+
    theme(legend.position = "none",axis.title.y=element_blank(),axis.title.x=element_blank())
  TE<-ggplot(as.data.frame(exprDat.SC.Norm.log[gene,which(colnames(exprDat.SC.Norm.log)%in%rownames(EchantMarqueursTE))]), aes(x=EchantMarqueursTE$Pseudotime, y=exprDat.SC.Norm.log[gene,which(colnames(exprDat.SC.Norm.log)%in%rownames(EchantMarqueursTE))],color=as.factor(as.character(EchantMarqueursTE$BranchesVal)))) + 
    geom_point(shape=18)+ 
    ylim(0,20)+
    labs(title = "TE\n", x = "Pseudotime", y = "Log2 Expr Gene", color = "Lineage\n")+
    scale_color_manual(values=c("black","cyan","blue"))+
    geom_smooth(method=loess, se=T, linetype="solid",color="blue",fullrange=T)+
    theme(legend.position = "none",axis.title.y=element_blank(),axis.title.x=element_blank())
  p<-ggdraw() +
    draw_plot(EPI, 0, 0, 0.33, 1) +
    draw_plot(PE, 0.33, 0, .33, 1) +
    draw_plot(TE, .66, 0, .33, 1)
  print(p)
}

PlotExprGeneCluster<-function(ClusterMeanExpr,EchantMarqueursEPI,EchantMarqueursPE,EchantMarqueursTE){
  require(ggplot2)
  require(cowplot)
  EPI<-ggplot(as.data.frame(ClusterMeanExpr[,which(colnames(ClusterMeanExpr)%in%rownames(EchantMarqueursEPI))]), aes(x=EchantMarqueursEPI$Pseudotime, y=ClusterMeanExpr[,which(colnames(ClusterMeanExpr)%in%rownames(EchantMarqueursEPI))],color=as.factor(as.character(EchantMarqueursEPI$BranchesVal)))) + 
    geom_point(shape=18)+ 
    ylim(0,20)+
    labs(title = "EPI\n", x = "Pseudotime", y = "Expr Genes Cluster", color = "Lineage\n")+
    scale_color_manual(values=c("black","orange","red"))+
    geom_smooth(method=loess, se=T, linetype="solid",color="red",fullrange=T)+
    theme(legend.position = "none",axis.title.x=element_blank())
  PE<-ggplot(as.data.frame(ClusterMeanExpr[,which(colnames(ClusterMeanExpr)%in%rownames(EchantMarqueursPE))]), aes(x=EchantMarqueursPE$Pseudotime, y=ClusterMeanExpr[,which(colnames(ClusterMeanExpr)%in%rownames(EchantMarqueursPE))],color=as.factor(as.character(EchantMarqueursPE$BranchesVal)))) + 
    geom_point(shape=18)+ 
    ylim(0,20)+
    labs(title = "PE\n", x = "Pseudotime", y = "Expr Gene", color = "Lineage\n")+
    scale_color_manual(values=c("black","orange","green"))+
    geom_smooth(method=loess, se=T, linetype="solid",color="green",fullrange=T)+
    theme(legend.position = "none",axis.title.y=element_blank(),axis.title.x=element_blank())
  TE<-ggplot(as.data.frame(ClusterMeanExpr[,which(colnames(ClusterMeanExpr)%in%rownames(EchantMarqueursTE))]), aes(x=EchantMarqueursTE$Pseudotime, y=ClusterMeanExpr[,which(colnames(ClusterMeanExpr)%in%rownames(EchantMarqueursTE))],color=as.factor(as.character(EchantMarqueursTE$BranchesVal)))) + 
    geom_point(shape=18)+ 
    ylim(0,20)+
    labs(title = "TE\n", x = "Pseudotime", y = "Expr Gene", color = "Lineage\n")+
    scale_color_manual(values=c("black","cyan","blue"))+
    geom_smooth(method=loess, se=T, linetype="solid",color="blue",fullrange=T)+
    theme(legend.position = "none",axis.title.y=element_blank(),axis.title.x=element_blank())
  p<-ggdraw() +
    draw_plot(EPI, 0, 0, 0.33, 1) +
    draw_plot(PE, 0.33, 0, .33, 1) +
    draw_plot(TE, .66, 0, .33, 1)
  print(p)
}


SmoothPlotExprGene<-function(ClusterMeanExpr,sampleAnnotSmooth){
  ## dataHeatMap : expression matrix
  ## sampleAnnotSmooth : sampleAnnot with Branches
  require(ggplot2)
  require(cowplot)
  
  # Sample Selection
  EchantMarqueursMorulaB1<-sampleAnnotSmooth[which(sampleAnnotSmooth$Lineage == "3.Morula4.B1"),]
  EchantMarqueursICM<-sampleAnnotSmooth[which(sampleAnnotSmooth$Lineage == "5. Unspecified ICM"),]
  EchantMarqueursSpéPE<-sampleAnnotSmooth[which(sampleAnnotSmooth$Lineage == "7.Primitive endoderm"),]
  EchantMarqueursSpéEPI<-sampleAnnotSmooth[which(sampleAnnotSmooth$Lineage == "6.Epiblast"),]
  EchantMarqueursSpéTE<-sampleAnnotSmooth[which(sampleAnnotSmooth$Lineage == "9.Trophectoderm"),]
  
  EchantMarqueursEPI<-rbind(rbind(EchantMarqueursMorulaB1,EchantMarqueursICM),EchantMarqueursSpéEPI)
  EPI<-ggplot(as.data.frame(ClusterMeanExpr[,which(colnames(ClusterMeanExpr)%in%rownames(EchantMarqueursEPI))]), aes(x=EchantMarqueursEPI$Pseudotime, y=ClusterMeanExpr[,which(colnames(ClusterMeanExpr)%in%rownames(EchantMarqueursEPI))],color=as.factor(as.character(EchantMarqueursEPI$Lineage)))) + 
    geom_point(shape=18)+ 
    ylim(0,20)+
    labs(title = "EPI\n", x = "Pseudotime", y = "Expr Genes Cluster", color = "Lineage\n")+
    scale_color_manual(values=c("black","yellow","red"))+
    theme(legend.position = "none",axis.title.x=element_blank())
  
  EchantMarqueursPE<-rbind(rbind(EchantMarqueursMorulaB1,EchantMarqueursICM),EchantMarqueursSpéPE)
  PE<-ggplot(as.data.frame(ClusterMeanExpr[,which(colnames(ClusterMeanExpr)%in%rownames(EchantMarqueursPE))]), aes(x=EchantMarqueursPE$Pseudotime, y=ClusterMeanExpr[,which(colnames(ClusterMeanExpr)%in%rownames(EchantMarqueursPE))],color=as.factor(as.character(EchantMarqueursPE$Lineage)))) + 
    geom_point(shape=18)+ 
    ylim(0,20)+
    labs(title = "PE\n", x = "Pseudotime", y = "Expr Gene", color = "Lineage\n")+
    scale_color_manual(values=c("black","yellow","green"))+
    theme(legend.position = "none",axis.title.y=element_blank(),axis.title.x=element_blank())
  
  EchantMarqueursTE<-rbind(EchantMarqueursMorulaB1,EchantMarqueursSpéTE)
  TE<-ggplot(as.data.frame(ClusterMeanExpr[,which(colnames(ClusterMeanExpr)%in%rownames(EchantMarqueursTE))]), aes(x=EchantMarqueursTE$Pseudotime, y=ClusterMeanExpr[,which(colnames(ClusterMeanExpr)%in%rownames(EchantMarqueursTE))],color=as.factor(as.character(EchantMarqueursTE$Lineage)))) + 
    geom_point(shape=18)+ 
    ylim(0,20)+
    labs(title = "TE\n", x = "Pseudotime", y = "Expr Gene", color = "Lineage\n")+
    scale_color_manual(values=c("black","blue"))+
    theme(legend.position = "none",axis.title.y=element_blank(),axis.title.x=element_blank())
  p<-ggdraw() +
    draw_plot(EPI, 0, 0, 0.33, 1) +
    draw_plot(PE, 0.33, 0, .33, 1) +
    draw_plot(TE, .66, 0, .33, 1)
  print(p)
}



###---------------------------------------------###
### New Expression Lineage-Pseudotime Specific ####
###---------------------------------------------###




SmoothCurveExpressionLineage<-function(scds,first_lineage_cell,last_lineage_cell,vectorPseudotime,DataMethodTransformation="N"){
  # scds : a CellDataSet Object obtain with Monocle
  # first_lineage_cell and last_lineage cell : two cells/sample that delimitated the lineage of interest. They are rownames of scds expression matrix "exprs(scds)"
  # DataMethodTransformation : Transformation of data, log2, vst, log or nothing by default.
  # vectorPseudotime : Vector of values for the pseudotime. (seq(0,100,length.out=100) : 100 values betweens 0 and 100 with same delta between each values)
  require(igraph) 
  require(scater)
  require(monocle)
  require(Matrix)
  if(!(DataMethodTransformation %in%c("N","vst","log10","log2"))){
    stop('DataMethodTransformation must be "vst", "log10","log2 or "N"')
  }
  if((!(first_lineage_cell %in%colnames(exprs(scds))))|(!(last_lineage_cell %in%colnames(exprs(scds))))){
    stop("At least one Cell don't exist")
  }
  if((!(is.numeric(vectorPseudotime)))|(length(vectorPseudotime)<10)){ # No numeric values or null object
    stop("vectorPseudotime must contains at least 10 numeric continuous values")
  }
  pr_graph_cell_proj_mst <- minSpanningTree(scds) # Tree caracterization
  path_to_ancestor <- shortest_paths(pr_graph_cell_proj_mst, first_lineage_cell, last_lineage_cell)
  path_to_ancestor <- names(unlist(path_to_ancestor$vpath)) # Path cells that describes the lineage
  scds2 <- scds[, row.names(pData(scds[,path_to_ancestor]))] #or just union(ancestor_cells, branch_cells)
  Pseudotime <- pData(scds2)$Pseudotime 
  pData <- pData(scds2)
  max_pseudotime_on_path <- max(pData[path_to_ancestor,]$Pseudotime) # Pseudotemps max branch
  pData$Pseudotime <- 100 * pData$Pseudotime / max_pseudotime_on_path # Adjusted Pseudotime
  pData$original_cell_id <- row.names(pData)
  exprs_data<-exprs(scds2) # Copy gene expression
  pData$Branch <- paste(first_lineage_cell,last_lineage_cell,sep="_to_") # Lineage 
  pData$State <- factor(pData$State)
  Size_Factor <- pData$Size_Factor
  fData <- fData(scds2)
  colnames(exprs_data) <- row.names(pData) # Create a pData annotation object from a cellDataSet Object
  options(warn=-1) ## Don't show warning (miss column gene_short_name in featureData that is necessary for others functions)
  new_scds <- newCellDataSet(as.matrix(exprs_data),
                             phenoData = new("AnnotatedDataFrame", data = pData),
                             featureData = new("AnnotatedDataFrame", data = fData),
                             expressionFamily=scds@expressionFamily,
                             lowerDetectionLimit=scds@lowerDetectionLimit)
  options(warn=0) 
  pData(new_scds)$State <- as.factor(pData(new_scds)$State)
  pData(new_scds)$Size_Factor <- Size_Factor
  new_scds@dispFitInfo <- scds@dispFitInfo
  # progenitor_state <- subset(pData(scds2), Pseudotime == 0)[, "State"]
  # branch_states <- setdiff(pData(scds2)$State, progenitor_state)
  newdata<-data.frame(Pseudotime = vectorPseudotime, Branch = as.factor(unique(as.character(pData(new_scds)$Branch))))
  # Update NewData with new values depending on new Pseudotime (specific to the lineage)
  Lineage_Gene_expr <- genSmoothCurves(new_scds[, ], cores = 1, trend_formula = "~sm.ns(Pseudotime, df=2)", relative_expr = T, new_data = newdata)
  if(DataMethodTransformation=="vst"){ # DataTransformation
    return(vstExprs(new_scds,expr_matrix=Lineage_Gene_expr)) 
  }
  if(DataMethodTransformation=="log2"){ 
    return(log2(Lineage_Gene_expr+1)) # log2(0) is impossible (-inf)
  }
  if(DataMethodTransformation=="log10"){
    return(log10(Lineage_Gene_expr+1)) # log2(0) is impossible (-inf)
  }
  return(Lineage_Gene_expr) # No transformation
}

# Return a VectorPseudotime lineage specific

calculPseudotime<-function(scds,first_lineage_cell,last_lineage_cell,begin=F,MinvectorPseudotime=14.36777){
  #begin : Var to control if we used first cell of the tree
  require(igraph) 
  n<-2.3
  if(begin) n<-1
  pr_graph_cell_proj_mst <- minSpanningTree(scds) # Tree caracterization
  path_to_ancestor <- shortest_paths(pr_graph_cell_proj_mst, closest_vertex[first_lineage_cell,], closest_vertex[last_lineage_cell,])
  path_to_ancestor <- names(unlist(path_to_ancestor$vpath)) # Path cells that describes the lineage
  vectorPseudotime <- seq(MinvectorPseudotime,(n*pData(scds)[last_lineage_cell,"Pseudotime"]),by=(n*(pData(scds)[last_lineage_cell,"Pseudotime"]-pData(scds)[first_lineage_cell,"Pseudotime"]))/(length(row.names(closest_vertex)[which(closest_vertex[,1]%in%path_to_ancestor)])/6))
  return(vectorPseudotime)
}


OldcalculPseudotime<-function(scds,first_lineage_cell,last_lineage_cell,MinvectorPseudotime=0){
  require(igraph) 
  pr_graph_cell_proj_mst <- minSpanningTree(scds) # Tree caracterization
  path_to_ancestor <- shortest_paths(pr_graph_cell_proj_mst, first_lineage_cell, last_lineage_cell)
  path_to_ancestor <- names(unlist(path_to_ancestor$vpath)) # Path cells that describes the lineage
  vectorPseudotime <- seq(MinvectorPseudotime,(MinvectorPseudotime+(length(row.names(closest_vertex)[which(closest_vertex[,1]%in%path_to_ancestor)]))/2),length.out=((MinvectorPseudotime+length(row.names(closest_vertex)[which(closest_vertex[,1]%in%path_to_ancestor)]))/2))
  return(vectorPseudotime)
}


NumberCellInLineage<-function(scds,first_lineage_cell,last_lineage_cell){
  require(igraph) 
  pr_graph_cell_proj_mst <- minSpanningTree(scds) # Tree caracterization
  path_to_ancestor <- shortest_paths(pr_graph_cell_proj_mst, first_lineage_cell, last_lineage_cell)
  path_to_ancestor <- names(unlist(path_to_ancestor$vpath)) # Path cells that describes the lineage
  return(length(row.names(closest_vertex)[which(closest_vertex[,1]%in%path_to_ancestor)])/6) # 1 Theorical cell = 4 cells
}
  
NewSCDSSmoothCurveExpressionLineage<-function(scds,first_lineage_cell,last_lineage_cell,vectorPseudotime,DataMethodTransformation="N"){
  # scds : a CellDataSet Object obtain with Monocle
  # first_lineage_cell and last_lineage cell : two cells/sample that delimitated the lineage of interest. They are rownames of scds expression matrix "exprs(scds)"
  # DataMethodTransformation : Transformation of data, log2, vst, log or nothing by default.
  # MinvectorPseudotime : Min value for the pseudotime. Scale and max depend of number of cells in the lineage
  # vectorPseudotime : Vector of values for the pseudotime. (seq(0,100,length.out=100) : 100 values betweens 0 and 100 with same delta between each values)
  require(igraph) 
  require(scater)
  require(monocle)
  require(Matrix)
  if(!(DataMethodTransformation %in%c("N","vst","log10","log2"))){
    stop('DataMethodTransformation must be "vst", "log10","log2 or "N"')
  }
  if((!(first_lineage_cell %in%paste("Y",V(scds@minSpanningTree),sep="_")))|(!(last_lineage_cell %in%paste("Y",V(scds@minSpanningTree),sep="_")))){
    stop("At least one Cell don't exist")
  }
  if((!(is.numeric(vectorPseudotime)))|(length(vectorPseudotime)<10)){ # No numeric values or null object
    stop("vectorPseudotime must contains at least 10 numeric continuous values")
  }
  pr_graph_cell_proj_mst <- minSpanningTree(scds) # Tree caracterization
  path_to_ancestor <- shortest_paths(pr_graph_cell_proj_mst, first_lineage_cell, last_lineage_cell)
  path_to_ancestor <- names(unlist(path_to_ancestor$vpath)) # Path cells that describes the lineage
  scds2 <- scds[, row.names(pData(scds[,row.names(closest_vertex)[which(closest_vertex[,1]%in%path_to_ancestor)]]))] #or just union(ancestor_cells, branch_cells)
  Pseudotime <- pData(scds2)$Pseudotime 
  pData <- pData(scds2)
  max_pseudotime_on_path <- max(pData[row.names(closest_vertex)[which(closest_vertex[,1]%in%path_to_ancestor)],]$Pseudotime) # Pseudotemps max branch
  pData$Pseudotime <- 100 * pData$Pseudotime / max_pseudotime_on_path # Adjusted Pseudotime
  pData$original_cell_id <- row.names(pData)
  exprs_data<-exprs(scds2) # Copy gene expression
  pData$Branch <- paste(first_lineage_cell,last_lineage_cell,sep="_to_") # Lineage 
  pData$State <- factor(pData$State)
  Size_Factor <- pData$Size_Factor
  fData <- fData(scds2)
  colnames(exprs_data) <- row.names(pData) # Create a pData annotation object from a cellDataSet Object
  options(warn=-1) ## Don't show warning (miss column gene_short_name in featureData that is necessary for others functions)
  new_scds <- newCellDataSet(as.matrix(exprs_data),
                             phenoData = new("AnnotatedDataFrame", data = pData),
                             featureData = new("AnnotatedDataFrame", data = fData),
                             expressionFamily=scds@expressionFamily,
                             lowerDetectionLimit=scds@lowerDetectionLimit)
  options(warn=0) 
  pData(new_scds)$State <- as.factor(pData(new_scds)$State)
  pData(new_scds)$Size_Factor <- Size_Factor
  new_scds@dispFitInfo <- scds@dispFitInfo
  # progenitor_state <- subset(pData(scds2), Pseudotime == 0)[, "State"]
  # branch_states <- setdiff(pData(scds2)$State, progenitor_state)
  # vectorPseudotime<-seq(MinvectorPseudotime,((MinvectorPseudotime+length(row.names(closest_vertex)[which(closest_vertex[,1]%in%path_to_ancestor)]))/2),length.out=((MinvectorPseudotime+length(row.names(closest_vertex)[which(closest_vertex[,1]%in%path_to_ancestor)]))/2))
  # 1 Pseudotime = 4 cells
  newdata<-data.frame(Pseudotime = vectorPseudotime, Branch = as.factor(unique(as.character(pData(new_scds)$Branch))))
  # Update NewData with new values depending on new Pseudotime (specific to the lineage)
  Lineage_Gene_expr <- genSmoothCurves(new_scds[, ], cores = 1, trend_formula = "~sm.ns(as.numeric(Pseudotime), df=2)", relative_expr = T, new_data = newdata)
  if(DataMethodTransformation=="vst"){ # DataTransformation
    return(vstExprs(new_scds,expr_matrix=Lineage_Gene_expr)) 
  }
  if(DataMethodTransformation=="log2"){ 
    return(log2(Lineage_Gene_expr+1)) # log2(0) is impossible (-inf)
  }
  if(DataMethodTransformation=="log10"){
    return(log10(Lineage_Gene_expr+1)) # log2(0) is impossible (-inf)
  }
  return(Lineage_Gene_expr) # No transformation
}

Cells_dataset<-function(scds,first_lineage_cell,last_lineage_cell){
  require(igraph) 
  require(scater)
  require(monocle)
  require(Matrix)
  pr_graph_cell_proj_mst <- minSpanningTree(scds) # Tree caracterization
  path_to_ancestor <- shortest_paths(pr_graph_cell_proj_mst, first_lineage_cell, last_lineage_cell)
  path_to_ancestor <- names(unlist(path_to_ancestor$vpath)) # Path cells that describes the lineage
}

NewSCDSS<-function(scds,all_cells_dataset){
  # scds : a CellDataSet Object obtain with Monocle
  # first_lineage_cell and last_lineage cell : two cells/sample that delimitated the lineage of interest. They are rownames of scds expression matrix "exprs(scds)"
  require(igraph) 
  require(scater)
  require(monocle)
  require(Matrix)
  scds2 <- scds[, row.names(pData(scds[,row.names(closest_vertex)[which(closest_vertex[,1]%in%all_cells_dataset)]]))] #or just union(ancestor_cells, branch_cells)
  Pseudotime <- pData(scds2)$Pseudotime 
  pData <- pData(scds2)
  max_pseudotime_on_path <- max(pData[row.names(closest_vertex)[which(closest_vertex[,1]%in%all_cells_dataset)],]$Pseudotime) # Pseudotemps max branch
  pData$Pseudotime <- 100 * pData$Pseudotime / max_pseudotime_on_path # Adjusted Pseudotime
  pData$original_cell_id <- row.names(pData)
  exprs_data<-exprs(scds2) # Copy gene expression
  pData$Branch <- pData$Branches # Lineage 
  pData$State <- factor(pData$State)
  Size_Factor <- pData$Size_Factor
  fData <- fData(scds2)
  colnames(exprs_data) <- row.names(pData) # Create a pData annotation object from a cellDataSet Object
  options(warn=-1) ## Don't show warning (miss column gene_short_name in featureData that is necessary for others functions)
  new_scds <- newCellDataSet(as.matrix(exprs_data),
                             phenoData = new("AnnotatedDataFrame", data = pData),
                             featureData = new("AnnotatedDataFrame", data = fData),
                             expressionFamily=scds@expressionFamily,
                             lowerDetectionLimit=scds@lowerDetectionLimit)
  options(warn=0) 
  pData(new_scds)$State <- as.factor(pData(new_scds)$State)
  pData(new_scds)$Size_Factor <- Size_Factor
  new_scds@dispFitInfo <- scds@dispFitInfo
  return(new_scds) # No transformation
}

SmoothRegLineage<-function(new_scds,vectorPseudotime,DataMethodTransformation="N"){
  Lineage_Gene_expr <- genSmoothCurves(new_scds[,], cores = 1, trend_formula = "~as.numeric(Pseudotime):as.factor(Branch)", relative_expr = T, new_data = vectorPseudotime)
  if(DataMethodTransformation=="vst"){ # DataTransformation
    return(vstExprs(new_scds,expr_matrix=Lineage_Gene_expr)) 
  }
  if(DataMethodTransformation=="log2"){ 
    return(log2(Lineage_Gene_expr+1)) # log2(0) is impossible (-inf)
  }
  if(DataMethodTransformation=="log10"){
    return(log10(Lineage_Gene_expr+1)) # log2(0) is impossible (-inf)
  }
  return(Lineage_Gene_expr) # No transformation
}
# 
# 
# f_expression<-t(round(exprs(scdsbis["POU5F1",])))
# modelFormulaStr <- paste("f_expression","~as.numeric(pData(scdsbis)$Pseudotime)*as.factor(as.character(pData(scdsbis)$Branches))", sep = "")
# FM_fit <- VGAM::vglm(as.formula(modelFormulaStr), family = scdsbis@expressionFamily, epsilon=1E-4)
# vectorPseudotime$Branch<-as.character(vectorPseudotime$Branch)
# vectorPseudotime$Branch[1:15]<-"2.Morula"
# vectorPseudotime$Branch[117:(117+floor(NumberCellInLineage(scds,"Y_68",closest_vertex["E6.6.714",])))]<-"5.Early trophectoderm"
# vectorPseudotime$Branch<-as.factor(vectorPseudotime$Branch)
# coeff<-FM_fit@coefficients
# names(coeff)<-c("Intercept",levels(vectorPseudotime$Branch))
# coeff<-c(coeff,0)
# coeff<-c(coeff,0)
# names(coeff)[2]<-"Pseudotemps"
# names(coeff)[15]<-"2.Morula"
# names(coeff)[16]<-paste("Inter","2.Morula","")
# names(coeff)[9:14]<-c(paste("Inter",levels(vectorPseudotime$Branch)[2:7],""))
# POU5F1<-(as.numeric(coeff["Intercept"])+as.numeric(coeff[as.character(vectorPseudotime$Branch)]))+((as.numeric(coeff["Pseudotemps"])+as.numeric(coeff[paste("Inter",as.character(vectorPseudotime$Branch),"")]))*as.numeric(as.character(vectorPseudotime$Pseudotime)))
# plot(vectorPseudotime$Pseudotime,POU5F1)
# # GATA6<-POU5F1
# 
# ExprSmooth<-rbind(POU5F1,GATA6)
# colnames(exprSmooth)<-rn(vectorPseudotime)
#### Test Function differents inputs

# vectorPseudotime<-seq(0,100,length.out=25)
# TE<-SmoothCurveExpressionLineage(scds,"M42C10","E7.13.881",vectorPseudotime,"log2")
# vectorPseudotime<-seq(0,100,length.out=100)
# EPI<-SmoothCurveExpressionLineage(scds,"M42C10","E6.10.1057",vectorPseudotime,"log10")
# vectorPseudotime<-seq(0,200,length.out=400)
# PE<-SmoothCurveExpressionLineage(scds,"M42C10","E7.14.905",vectorPseudotime,"vst")


# PlotSmoothCurveExpressionLineage<-function(exprLog,Lineage,col){
#   plot(x=colnames(exprLog),y=exprLog[,],col=col,pch=16,xlab="Pseudotime",ylab=paste("Expr Gene",row.names(exprLog),sep=" "),main=paste("Calculate Expression Pseudotime",Lineage,sep=" "))
# }


###---------------------###
### Mutual Information ####
###---------------------###


MutualInformation<-function(exprLog,methodMI="emp",disc="equalfreq",nbins=nrow(exprLog)^(1/3),Discretization=T,OneGene=TRUE,Gene=NULL){ 
  ## exprLog : genes expression dataframe 
  ## methodMI : Method estimator for the Mutual Information. 'emp' for empirical, 'mm' for Miller-Madow, 'shrink' for shrinkage and 'sg' for Schurmann-Grassberger 
  ## disc : Discretization method."equalwidth" and "equalfreq" discretizes each random variable (each column) of the data into nbins. "globalequalwidth" discretizes the range of the random vector data into nbins. 
  ## nbins : number of bins to be used for the discretization 
  ## Discretization : Discretize data is necessary to execute Mutual Information 
  ## OneByOne :  All genes are computed one by one with an output file or all in a dataframe in one time. 
  require(infotheo) 
  exprLog<-t(exprLog) # Genes are variables, we need genes in columns 
  if(!(methodMI%in%c("emp","mm","shrink","sg"))){ 
    stop("Method Mutual Information don't exist. They are only 'emp', 'mm', 'shrink' and 'sg' ") 
  } 
  if(Discretization){ 
    if(!(disc%in%c("equalfreq","equalwidth","globalequalwidth"))){ 
      stop("Method Discretization don't exist. They are only 'equalfreq', 'equalwidth' and 'globalequalwidth ") 
    } 
    dat<-discretize(exprLog,disc=disc,nbins=nbins) # discretes values needed 
    row.names(dat)<-row.names(exprLog) 
  } 
  print("Mutual information computing") 
  if(OneGene){ 
    if(is.null(Gene)) stop("You must choose a gene of interest")
    MI<-c()
    for(i in 1:ncol(dat)){
      res<-mutinformation(dat[,i],dat[,Gene], method=methodMI) ## One gene against others
      MI<-c(MI,res)
    }
    names(MI)<-colnames(dat)
    return(MI)
  } 
  else{ # One gene with all others genes 
    return(mutinformation(dat, method=methodMI)) ## Take a long time for a big dataset (~5mins 1000genes) 
  } 
} 


###---------------------###
### Network into json  ####
###---------------------###



exportGraph <- function(g,filename){
  #convert graph into a list
  graph <- list()
  graph$edges <- as_data_frame(g, what = "edges")
  graph$vertices <- as_data_frame(g, what="vertices")
  # polish vertices
  row.names(graph$vertices) <- NULL
  if(ncol(graph$vertices)==0) graph$vertices <- NULL #in case the 
  
  graph$directed <- is.directed(g)
  graph$name <- g$name
  #convert list into a json
  json.content <- toJSON(graph, pretty=TRUE)
  #write json into a file
  if(!missing(filename)) {
    sink(filename)
    cat(json.content)
    sink()
  }
  #return the json in case you need it
  return(json.content)
}

exportGML <- function(graph, filename) {
  file <- file(filename, "w")
  cat("Creator \"igraph exportGML\"\n", file=file)
  cat("Version 1.0\n", file=file)
  cat("graph\n[\n", file=file)
  cat("  directed", as.integer(is.directed(graph)), "\n", file=file)
  for (i in seq_len(vcount(graph))) {
    cat("  node\n  [\n", file=file)
    cat("    id", i-1, "\n", file=file)
    cat("    graphics\n    [\n", file=file)
    cat("      fill \"", V(graph)$color[i], "\"\n", sep="", file=file)
    cat("      type \"rectangle\"\n", file=file)
    cat("      outline \"#000000\"\n", file=file)
    cat("    ]\n", file=file)
    cat("  ]\n", file=file)
  }
  el <- get.edgelist(graph, names=FALSE)
  for (i in seq_len(nrow(el))) {
    cat("  edge\n  [\n", file=file)
    cat("    source", el[i,1], "\n", file=file)
    cat("    target", el[i,2], "\n", file=file)
    cat("  ]\n", file=file)
  }
  cat("]\n", file=file)
  close(file)
}



#' Allows to convert a json into an igraph object
#' It is essential that the json contains an igraph object
#' 
#' @param filename either the json or the filename of the file containing the json
importGraph <- function(filename){
  built.graph <- fromJSON(filename, flatten=TRUE)
  if("vertices" %in% names(built.graph)){
    built.g <- graph_from_data_frame(built.graph$edges, directed=built.graph$directed, vertices=built.graph$vertices)
  }else{
    built.g <- graph_from_data_frame(built.graph$edges, directed=built.graph$directed)
  }
  if("name" %in% names(built.graph)){
    built.g$name <- built.graph$name
  }
  return(built.g)
}
