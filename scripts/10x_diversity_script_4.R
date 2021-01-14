library(dplyr)
library(Seurat)
library(ggplot2)
library(treemap)
library(circlize)
library(gridExtra)

## there is more info in the comments (and things to add lower)

setwd("/SC_seq/Analysis")

expr.data_SVH1a <- Read10X(data.dir = "/filtered_feature_bc_matrix")
# Initialize the Seurat object with the raw (non-normalized data).

seurat_object_SVH1a = CreateSeuratObject(counts = expr.data_SVH1a$`Gene Expression`,project = "SeuratProject1a")
seurat_object_SVH1a[['ADT']] <- CreateAssayObject(counts = expr.data_SVH1a$`Antibody Capture`[,colnames(seurat_object_SVH1a)],min.cells=1)

# Run the standard workflow for visualization and clustering

seurat_object_SVH1a<-FindVariableFeatures(seurat_object_SVH1a)
seurat_object_SVH1a <- ScaleData(seurat_object_SVH1a, verbose = FALSE)
seurat_object_SVH1a <- RunPCA(seurat_object_SVH1a, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
seurat_object_SVH1a <- RunTSNE(seurat_object_SVH1a, reduction = "pca", dims = 1:20)
seurat_object_SVH1a <- FindNeighbors(seurat_object_SVH1a, reduction = "pca", dims = 1:20)
seurat_object_SVH1a <- FindClusters(seurat_object_SVH1a, resolution = 0.8)

p1 <- DimPlot(seurat_object_SVH1a, reduction = "tsne", group.by = "orig.ident")
p2 <- DimPlot(seurat_object_SVH1a, reduction = "tsne", label = TRUE)
CombinePlots(plots = list(p1, p2))

DefaultAssay(seurat_object_SVH1a) <- "RNA"
B2.markers <- FindConservedMarkers(seurat_object_SVH1a, ident.1 = 2, grouping.var = "orig.ident", verbose = FALSE)
head(B2.markers)

###files lung&BAL
seurat_object<-seurat_object_SVH1a
seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^mt-")
seurat_object <- subset(seurat_object, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5)
seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seurat_object), 10)

# plot variable features with and without labels
all.genes <- rownames(seurat_object)
seurat_object <- ScaleData(seurat_object, features = all.genes)
seurat_object <- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object))
seurat_object <- FindNeighbors(seurat_object, dims = 1:10)
seurat_object <- FindClusters(seurat_object, resolution = 0.5)
seurat_object <- RunTSNE(seurat_object, dims = 1:10)
DimPlot(seurat_object, reduction = "tsne")

markers.all=FindAllMarkers(seurat_object,thresh.test = 3,test.use = "roc", do.print = TRUE)


##loading BCR data
bcr.contig1a <- read.csv("SVH1a/BCR/filtered_contig_annotations.csv", sep=",")
head(bcr.contig1a)
bcr.contig1a$barcode <- gsub("-1", "", bcr.contig1a$barcode)

bcr.contig1a.productive<-subset(bcr.contig1a, bcr.contig1a$productive=="True")
bcr.contig1a.productive.heavy<-subset(bcr.contig1a.productive , grepl(glob2rx("IGH*"), bcr.contig1a.productive$c_gene))
bcr.contig1a.productive.light<-subset(bcr.contig1a.productive ,grepl(glob2rx("IGL*"), bcr.contig1a.productive$c_gene) | grepl(glob2rx("IGK*"), bcr.contig1a.productive$c_gene))
bcr.contig1a.productive.heavy <- bcr.contig1a.productive.heavy[!duplicated(bcr.contig1a.productive.heavy$barcode), ]
bcr.contig1a.productive.light <- bcr.contig1a.productive.light[!duplicated(bcr.contig1a.productive.light$barcode), ]
bcr1a.heavy<-bcr.contig1a.productive.heavy[,c("v_gene","d_gene","j_gene","c_gene", "cdr3","cdr3_nt")]
rownames(bcr1a.heavy)<-paste(bcr.contig1a.productive.heavy[,1],"-1",sep="")
colnames(bcr1a.heavy)<-c("h.v_gene","h.d_gene","h.j_gene","h.c_gene", "h.cdr3","h.cdr3_nt")
bcr1a.light<-bcr.contig1a.productive.light[,c("v_gene","d_gene","j_gene","c_gene", "cdr3","cdr3_nt")]
rownames(bcr1a.light)<-paste(bcr.contig1a.productive.light[,1],"-1", sep="")
colnames(bcr1a.light)<-c("l.v_gene","l.d_gene","l.j_gene","l.c_gene", "l.cdr3","l.cdr3_nt")

SVH1a<-seurat_object
SVH1a_BCR <- AddMetaData(object=SVH1a, metadata=bcr1a.heavy)
SVH1a_BCR <- AddMetaData(object=SVH1a_BCR, metadata=bcr1a.light)
clusters.n<-length(levels(SVH1a_BCR@meta.data$seurat_clusters))
## Generate a table with VH family data

Vhfamtable1a<-as.data.frame(matrix(0,length(SVH1a_BCR@meta.data$h.v_gene),1))
colnames(Vhfamtable1a)<-c("h.V.fam")
Vhfamtable1a[,1]<-c("no_IGHV_data")
rownames(Vhfamtable1a)<-rownames(SVH1a_BCR@meta.data)
for(n in 1:16){
  grx <- glob2rx(paste("IGHV",n,"-*", sep=""))
  for(r in 1:length(SVH1a_BCR@meta.data$h.v_gene)){
    if(grepl(grx, as.data.frame(SVH1a_BCR@meta.data$h.v_gene)[r,1])==TRUE){
      Vhfamtable1a[r,1]<-paste("IGHV", n, sep="-")
    }
  }
}


SVH1a_BCR <- AddMetaData(object=SVH1a_BCR, metadata=Vhfamtable1a)

## Generate a table with VL family data

Vlfamtable1a<-as.data.frame(matrix(0,length(SVH1a_BCR@meta.data$l.v_gene),1))
colnames(Vlfamtable1a)<-c("L.V.fam")
Vlfamtable1a[,1]<-c("no_IGHV_data")
rownames(Vlfamtable1a)<-rownames(SVH1a_BCR@meta.data)
for(r in 1:length(SVH1a_BCR@meta.data$l.v_gene)){
  if(grepl(glob2rx("IGKV*"), as.data.frame(SVH1a_BCR@meta.data$l.v_gene)[r,1])==TRUE){
    for(n in 1:19){
      grx <- glob2rx(paste("IGKV",n,"-*", sep=""))
        if(grepl(grx, as.data.frame(SVH1a_BCR@meta.data$l.v_gene)[r,1])==TRUE){
          Vlfamtable1a[r,1]<-paste("IGKV", n, sep="-")
      }
    }
  }
  else{
    for(n in 1:8){
      grx <- glob2rx(paste("IGLV",n,"-*", sep=""))
      if(grepl(grx, as.data.frame(SVH1a_BCR@meta.data$l.v_gene)[r,1])==TRUE){
        Vlfamtable1a[r,1]<-paste("IGLV", n, sep="-")
      }
    }
    
  }
}

SVH1a_BCR <- AddMetaData(object=SVH1a_BCR, metadata=Vlfamtable1a)


## barplots VH family usage

Vfamtabley<-as.data.frame(matrix(0,16,length(levels(SVH1a_BCR@meta.data$seurat_clusters))))
rownames(Vfamtabley)<-seq(1:16)
Vfamtabley[,1:length(levels(SVH1a_BCR@meta.data$seurat_clusters))]<-0
colnames(Vfamtabley)<-paste("Seurat_cluster", levels(SVH1a_BCR@meta.data$seurat_clusters), sep="_")
for(x in 1:length(levels(SVH1a_BCR@meta.data$seurat_clusters))){
    object.x<-subset(SVH1a_BCR, seurat_clusters==x)
    Vfamtablex<-as.data.frame(matrix(0,16,1))
    rownames(Vfamtablex)<-seq(1:16)
    Vfamtablex[,1]<-0
    colnames(Vfamtablex)<-c("Freq")
    for(n in 1:16){
      grx <- glob2rx(paste("IGHV",n,"-*", sep=""))
      for(r in 1:length(object.x@meta.data$h.v_gene)){
        if(grepl(grx, as.data.frame(object.x@meta.data$h.v_gene)[r,1])==TRUE){
          Vfamtablex[n,1]<-Vfamtablex[n,1]+1
        }
      }
      Vfamtabley[,as.numeric(x)+1]<-Vfamtablex
    }
    Vfamtabley
}
Vfamtabley

myplots <- list()
for(d in 1:length(levels(unique(SVH1a_BCR@meta.data$seurat_clusters)))){
  cluster.data <- as.data.frame(Vfamtabley[,d])
  colnames(cluster.data)<-c("Freq")
  cluster.data$IGHV<-paste("IGHV",seq(1:16),sep="-")
  nam<-paste("p",d-1,sep="")
  px<-ggplot(cluster.data, aes(x=IGHV, y=Freq, fill = as.factor(IGHV)), title="IGHV family usage")+
    geom_col()+
    labs(y="Cell number",title = paste("IGHV family usage Lung - cluster", d-1,sep=" "), fill='IGHV family')+
    scale_fill_manual(values = colorRampPalette(c("darkblue", "light blue"))(16), labels = c(paste("IGHV-",seq(1:16))))+
    geom_label(data=NULL,aes(label=Freq), fill="white",size=2)+
    theme(
      panel.background = element_rect(fill = "white"), # bg of the panel
      plot.background = element_rect(fill = "white"), # bg of the plot
      panel.grid.major = element_blank(), # get rid of major grid
      panel.grid.minor = element_blank(), # get rid of minor grid
      legend.background = element_rect(fill = "transparent"), # get rid of legend bg
      legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
      legend.position="none",
      legend.key.size = unit(0.1, "cm"),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  assign(nam, px)
  nam
  myplots[[d]] <- px
  myplots
}
pdf(file="test_IGHVusage_cluster.pdf", width=20, height=10)
grid.arrange(grobs = myplots, ncol = 3)
dev.off()


## barplots VL family usage

Vfamtabley<-as.data.frame(matrix(0,19+8,length(levels(SVH1a_BCR@meta.data$seurat_clusters))))
rownames(Vfamtabley)<-seq(1:27)
Vfamtabley[,1:length(levels(SVH1a_BCR@meta.data$seurat_clusters))]<-0
colnames(Vfamtabley)<-paste("Seurat_cluster", levels(SVH1a_BCR@meta.data$seurat_clusters), sep="_")
rownames(Vfamtabley)[1:19]<-paste("IGKV",seq(1:19),sep="-")
rownames(Vfamtabley)[20:27]<-paste("IGLV",seq(1:8),sep="-")
for(x in levels(SVH1a_BCR@meta.data$seurat_clusters)){
  object.x<-subset(SVH1a_BCR, seurat_clusters==x)
  Vfamtablex<-as.data.frame(matrix(0,27,1))
  rownames(Vfamtablex)[1:19]<-paste("IGKV",seq(1:19),sep="-")
  rownames(Vfamtablex)[20:27]<-paste("IGLV",seq(1:8),sep="-")
  Vfamtablex[,1]<-0
  colnames(Vfamtablex)<-c("Freq")
  for(r in 1:length(object.x@meta.data$l.v_gene)){
    if(grepl(glob2rx("IGKV*"), as.data.frame(object.x@meta.data$l.v_gene)[r,1])==TRUE){
      for(n in 1:19){
        grx <- glob2rx(paste("IGKV",n,"-*", sep=""))
        if(grepl(grx, as.data.frame(object.x@meta.data$l.v_gene)[r,1])==TRUE){
          Vfamtablex[n,1]<-Vfamtablex[n,1]+1
        }
      }
    }
    if(grepl(glob2rx("IGLV*"), as.data.frame(object.x@meta.data$l.v_gene)[r,1])==TRUE){
      for(n in 1:8){
        grx <- glob2rx(paste("IGLV",n,"*", sep=""))
        row<-n+19
        if(grepl(grx, as.data.frame(object.x@meta.data$l.v_gene)[r,1])==TRUE){
          Vfamtablex[row,1]<-Vfamtablex[row,1]+1
        }
      }
      
    }
  }
  Vfamtabley[,as.numeric(x)+1]<-Vfamtablex
  Vfamtabley
  }
Vfamtabley

  vl.fam.plots <- list()
for(d in 1:length(levels(unique(SVH1a_BCR@meta.data$seurat_clusters)))){
  cluster.data <- as.data.frame(Vfamtabley[,d])
  colnames(cluster.data)<-c("Freq")
  cluster.data$IGLV<-rownames(Vfamtabley)
  nam<-paste("p",d-1,sep="")
  px<-ggplot(cluster.data, aes(x=IGLV, y=Freq, fill = as.factor(IGLV)), title="IGLV family usage")+
    geom_col()+
    labs(y="Cell number",title = paste("IGLV family usage Lung - cluster", d-1,sep=" "), fill='IGLV family')+
    scale_fill_manual(values = c(colorRampPalette(c("darkblue", "light blue"))(19),colorRampPalette(c("darkred", "red"))(8)), labels = c(rownames(Vfamtabley)))+
    geom_label(data=NULL,aes(label=Freq), fill="white",size=2)+
    theme(
      panel.background = element_rect(fill = "white"), # bg of the panel
      plot.background = element_rect(fill = "white"), # bg of the plot
      panel.grid.major = element_blank(), # get rid of major grid
      panel.grid.minor = element_blank(), # get rid of minor grid
      legend.background = element_rect(fill = "transparent"), # get rid of legend bg
      legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
      legend.position="none",
      legend.key.size = unit(0.1, "cm"),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  assign(nam, px)
  nam
  vl.fam.plots[[d]] <- px
  vl.fam.plots
}
pdf(file="test_IGLVusage_cluster.pdf", width=20, height=10)
grid.arrange(grobs = vl.fam.plots, ncol = 3)
dev.off()



##isotype usage lc
iso.tabley<-as.data.frame(matrix(0,length(unique(SVH1a_BCR@meta.data$l.c_gene[!is.na(SVH1a_BCR@meta.data$l.c_gene)==TRUE])),length(levels(SVH1a_BCR@meta.data$seurat_clusters))))
rownames(iso.tabley)<-unique(SVH1a_BCR@meta.data$l.c_gene[!is.na(SVH1a_BCR@meta.data$l.c_gene)==TRUE])
iso.tabley[,1:length(levels(SVH1a_BCR@meta.data$seurat_clusters))]<-0
colnames(iso.tabley)<-paste("Seurat_cluster", levels(SVH1a_BCR@meta.data$seurat_clusters), sep="_")
clusters.n<-length(levels(SVH1a_BCR@meta.data$seurat_clusters))
isotype.names<-rownames(iso.tabley)
for(x in 0:clusters.n){
    object.x<-subset(SVH1a_BCR, seurat_clusters==x)
    iso.tablex<-as.data.frame(matrix(0,length(isotype.names),1))
    rownames(iso.tablex)<-isotype.names
    iso.tablex[,1]<-0
    colnames(iso.tablex)<-c("Freq")
    for(n in 1:nrow(iso.tablex)){
      grx <- rownames(iso.tablex)[n]
      for(r in 1:length(object.x@meta.data$l.c_gene)){
        if(grepl(grx, as.data.frame(object.x@meta.data$l.c_gene)[r,1])==TRUE){
          iso.tablex[n,1]<-iso.tablex[n,1]+1
        }
      }
    }
    iso.tabley[,x+1]<-iso.tablex
    iso.tabley
  }

iso.tabley



iso.plots <- list()
for(d in 1:length(isotype.names)){
  cluster.data <- as.data.frame(iso.tabley[,d])
  colnames(cluster.data)<-c("Freq")
  cluster.data$isotype<-isotype.names
  nam<-paste("p",d,sep="")
  px<-ggplot(cluster.data, aes(x=isotype, y=Freq, fill = as.factor(isotype)), title="IGH isotype usage")+
    geom_col()+
    labs(y="Cell number",title = paste("IGH isotype usage Lung - cluster", d,sep=" "), fill='Isotype')+
    scale_fill_manual(values = colorRampPalette(c("darkblue", "light blue"))(length(isotype.names)), labels = isotype.names)+
    geom_label(data=NULL,aes(label=Freq), fill="white",size=2)+
    theme(
      panel.background = element_rect(fill = "white"), # bg of the panel
      plot.background = element_rect(fill = "white"), # bg of the plot
      panel.grid.major = element_blank(), # get rid of major grid
      panel.grid.minor = element_blank(), # get rid of minor grid
      legend.background = element_rect(fill = "transparent"), # get rid of legend bg
      legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
      legend.position="none",
      legend.key.size = unit(0.1, "cm"),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  assign(nam, px)
  nam
  iso.plots[[d]] <- px
  iso.plots
}



pdf(file="test_light_chain_isotype_usage_cluster.pdf", width=20, height=10)
grid.arrange(grobs = iso.plots, ncol = 3)
dev.off()

##isotype usage hc
iso.tabley<-as.data.frame(matrix(0,length(unique(SVH1a_BCR@meta.data$h.c_gene[!is.na(SVH1a_BCR@meta.data$h.c_gene)==TRUE])),length(levels(SVH1a_BCR@meta.data$seurat_clusters))))
rownames(iso.tabley)<-unique(SVH1a_BCR@meta.data$h.c_gene[!is.na(SVH1a_BCR@meta.data$h.c_gene)==TRUE])
iso.tabley[,1:length(levels(SVH1a_BCR@meta.data$seurat_clusters))]<-0
colnames(iso.tabley)<-paste("Seurat_cluster", levels(SVH1a_BCR@meta.data$seurat_clusters), sep="_")
clusters.n<-length(levels(SVH1a_BCR@meta.data$seurat_clusters))
isotype.names<-rownames(iso.tabley)
for(x in 0:clusters.n){
  object.x<-subset(SVH1a_BCR, seurat_clusters==x)
  iso.tablex<-as.data.frame(matrix(0,length(isotype.names),1))
  rownames(iso.tablex)<-isotype.names
  iso.tablex[,1]<-0
  colnames(iso.tablex)<-c("Freq")
  for(n in 1:nrow(iso.tablex)){
    grx <- rownames(iso.tablex)[n]
    for(r in 1:length(object.x@meta.data$h.c_gene)){
      if(grepl(grx, as.data.frame(object.x@meta.data$h.c_gene)[r,1])==TRUE){
        iso.tablex[n,1]<-iso.tablex[n,1]+1
      }
    }
  }
  iso.tabley[,x+1]<-iso.tablex
  iso.tabley
}

iso.tabley



iso.plots <- list()
for(d in 1:length(isotype.names)){
  cluster.data <- as.data.frame(iso.tabley[,d])
  colnames(cluster.data)<-c("Freq")
  cluster.data$isotype<-isotype.names
  nam<-paste("p",d,sep="")
  px<-ggplot(cluster.data, aes(x=isotype, y=Freq, fill = as.factor(isotype)), title="IGH isotype usage")+
    geom_col()+
    labs(y="Cell number",title = paste("IGH isotype usage Lung - cluster", d,sep=" "), fill='Isotype')+
    scale_fill_manual(values = colorRampPalette(c("darkblue", "light blue"))(length(isotype.names)), labels = isotype.names)+
    geom_label(data=NULL,aes(label=Freq), fill="white",size=2)+
    theme(
      panel.background = element_rect(fill = "white"), # bg of the panel
      plot.background = element_rect(fill = "white"), # bg of the plot
      panel.grid.major = element_blank(), # get rid of major grid
      panel.grid.minor = element_blank(), # get rid of minor grid
      legend.background = element_rect(fill = "transparent"), # get rid of legend bg
      legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
      legend.position="none",
      legend.key.size = unit(0.1, "cm"),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  assign(nam, px)
  nam
  iso.plots[[d]] <- px
  iso.plots
}



pdf(file="test_isotype_usage_cluster.pdf", width=20, height=10)
grid.arrange(grobs = iso.plots, ncol = 3)
dev.off()


  
###circosplots_clone_VH
pdf(file="test_IGHV_IGLV_usage_cluster.pdf", width=20, height=5)

par(mfrow=c(1,clusters.n))
for(x in 0:clusters.n){
  object.x<-subset(SVH1a_BCR, seurat_clusters==x)
  V.heavy.light<-as.data.frame(object.x@meta.data$h.V.fam)
  colnames(V.heavy.light)<-c("heavy")
  V.heavy.light$light<-object.x@meta.data$l.v_gene
  V.heavy.light2<-as.data.frame(V.heavy.light[which(V.heavy.light!="no_IGHV_data",TRUE),])
  Vlink<-table(V.heavy.light2)
  cols<-colorRampPalette(c("darkblue", "grey"))(1)
  par(cex = 0.4, mar = c(0, 0, 0, 0))
  chordDiagram(as.matrix(Vlink), annotationTrack = "grid", preAllocateTracks = 1, grid.col = cols)
  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 25),cex=0.5)
    circos.axis(h = "top", labels.cex = 0.8, major.tick.length = 0.2, sector.index = sector.name, track.index = 2)
    title(paste("IGHV-IGLV usage Lung - cluster",x,sep=" "))
    for(n in 1:16){
      grx <- glob2rx(paste("Musmus IGHV",n,"-*", sep=""))
      if(grepl(grx, sector.name)==TRUE){
        sn<-sector.name
        circos.text(mean(get.cell.meta.data ("xlim", sector.index = sn)),get.cell.meta.data ("ylim", sector.index = sn)[1] +0.1, n, facing = "clockwise", niceFacing = TRUE, adj = c(0, 25))
        
      }
    }
  }, bg.border = NA)
  
}
dev.off()

##CDR3 length
cdr3.length.plots<-list()
for(x in 0:clusters.n){
  df<-as.data.frame(matrix(0,25,1))
  df$Var1<-as.numeric(seq(1:25))
  df$Freq<-1
  df<-df[,2:3]
  
  
  object.x<-subset(SVH1a_BCR, seurat_clusters==x)
  
  cdr3.heavy<-as.data.frame(object.x@meta.data$h.cdr3)
  cdr3.heavy<-cdr3.heavy[!is.na(cdr3.heavy)==TRUE,]
  cdr3.heavy<-as.data.frame(table(nchar(cdr3.heavy)))
  cdr3.heavy[,1]<-as.numeric(cdr3.heavy[,1])
  cdr3.heavy[,2]<-as.numeric(cdr3.heavy[,2])
  
  df2<-left_join(df,cdr3.heavy,by=c("Var1"))
  
  cdr3.light<-as.data.frame(object.x@meta.data$l.cdr3)
  cdr3.light<-cdr3.light[!is.na(cdr3.light)==TRUE,]
  cdr3.light<-as.data.frame(table(nchar(cdr3.light)))
  cdr3.light[,1]<-as.numeric(cdr3.light[,1])
  cdr3.light[,2]<-as.numeric(cdr3.light[,2])
  
  df3<-left_join(df2,cdr3.light,by=c("Var1"))
  
  df3$Freq.y<-(df3$Freq.y/sum(df3$Freq.y,na.rm=TRUE))*100
  df3$Freq<-(df3$Freq/sum(df3$Freq,na.rm=TRUE))*100
  
  px<-ggplot(df3, aes(Var1)) + 
    geom_line(aes(y = Freq.y, colour = "black")) + 
    geom_line(aes(y = Freq, colour = "red"))+
    labs(x="CDR3 length (AA)",y="Frequency of cells",title = paste("CDR3 length - cluster", x,sep=" "))+
    scale_x_continuous(limits = c(0, 26),breaks = c(seq(0:26)))+
    scale_color_manual(name=c("chain"),labels = c("heavy chain", "light chain"), values = c("black", "red"))+
    theme(
      panel.background = element_rect(fill = "white"), # bg of the panel
      plot.background = element_rect(fill = "white"), # bg of the plot
      axis.line = element_line(colour="black"),
      panel.grid.major = element_blank(), # get rid of major grid
      panel.grid.minor = element_blank(), # get rid of minor grid
      legend.background = element_rect(fill = "transparent"), # get rid of legend bg
      legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
      legend.key.size = unit(0.1, "cm"),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  nam<-paste("p",x,sep="")
  assign(nam, px)
  nam
  cdr3.length.plots[[x+1]] <- px
  cdr3.length.plots
}
pdf(file="test_cdr3lengths.pdf", width=20, height=10)
grid.arrange(grobs = cdr3.length.plots, ncol = 3)
dev.off()














