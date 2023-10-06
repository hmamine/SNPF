#require packages for the analysis the analysis
pkg=c("plyr", "ggplot2", "phyloseq", "data.table", "ape", "PMCMRplus","UpSetR", "venn","seqtime","metagenomeSeq","ggtern")
lapply(pkg, library, character.only = TRUE)

rm(list = ls()[grep("*.*", ls())])
set.seed(123)
Palette.phylum1 <-c("#EF5656","#47B3DA","#F7A415","#c1c5c1","#55555f","#2BB065")
Palette.phylum2 <-c("#EF5656","#47B3DA","#c1c5c1","#55555f","#2BB065")
ColPhylum<-c("Actinobacteriota","Bacteroidota","Firmicutes","Proteobacteria")
color_palette<-c("#00b2b2","#b02b76","#155832", "#000000","#ffa500")

lines <-data.frame(x=c(0.5,0,0.5),y=c(0.5,0.5,0),z=c(0,0.5,0.5), xend=c(1,1,1)/3, yend=c(1,1,1)/3, zend=c(1,1,1)/3)

minZ=3
maxZ=9
RAN=c(minZ,maxZ)

Shape.1=c(19,21,21)
DT1<-vector( "list" )
DT2<-vector( "list" )
DT3<-vector( "list" )
#log2 fold change
	folch=2

#Adjusted p-value
	alpha=0.01

theme_new <- theme (
	panel.grid.minor = element_blank(),
	panel.grid.major = element_blank(),
	plot.background = element_blank(),
	legend.position="none",
	)

# upload and prepare phyloseq objects***
	mat=read.table( "otu_table.txt", sep="\t", row.names=1, header=T)
	mat=as.matrix(mat)
	OTU=otu_table(mat, taxa_are_rows=T) 
	tax=read.table("taxonomy.txt", sep="\t", row.names=1, header=1)
	tax=as.matrix(tax)
	TAXA=tax_table(tax)
	sd=read.table("sample_data.txt", sep="\t", row.names=1, header=1)
	SD=sample_data(sd)

	physeq=phyloseq(OTU, TAXA, SD) 
	physeq1=subset_samples(physeq, Cond != "none")
	
	NTAXA=tax_table(physeq1)
	NSD=sample_data(physeq1)
	
	otu.I=otu_table(physeq1)
	taxa.I=tax_table(physeq1)
	sd.I=sample_data(physeq1)
	flt=filterTaxonMatrix(otu.I,minocc=2, keepSum = TRUE, return.filtered.indices = TRUE)
	otu.I.filtered=flt$mat
	taxa.I.filtered=taxa.I[setdiff(1:nrow(taxa.I),flt$filtered.indices),]
	dummyTaxonomy=c("k__dummy","p__","c__","o__","f__","g__","s__","sc__","ot_")
	TAXA.filt=rbind(taxa.I.filtered, dummyTaxonomy)

	rownames(TAXA.filt)[nrow(TAXA.filt)]="SUM"
	rownames(otu.I.filtered)[nrow(otu.I.filtered)]="SUM"
	NewOTU=otu_table(otu.I.filtered, taxa_are_rows = TRUE)
	NewTaxa=tax_table(TAXA.filt)
	
	DT.taxa=data.table( as(NewTaxa,"matrix"), keep.rownames=T, key="rn" )
#	DT.taxa=data.table( as(Taxa,"matrix"), keep.rownames=T, key="rn" )

	physeq.f=phyloseq(NewOTU, NewTaxa, sd.I)
#	physeq.Ctl=subset_samples(physeq.f, Cond == "H")
#	physeq.Fus=subset_samples(physeq.f, Cond == "H")
	physeq.f=tax_glom(physeq.f, "Genus") 
	LIST1<-list("H","F")
for ( j in LIST1 )
	{
	print( j )
	physeq.j=subset_samples( physeq.f, Cond == j)
	
#	computing differentially abundant ASVs	
	physeq.i=subset_samples( physeq.j, Treat %in% c("nS","cS") )	
	m=as( otu_table( physeq.i ), "matrix" ) + 1L
	t=data.frame( as( tax_table( physeq.i ), "matrix" ) )
	T=AnnotatedDataFrame( t )
	s=as( sample_data( physeq.i ), "data.frame" )
	S=AnnotatedDataFrame( s )
	obj=newMRexperiment( m, phenoData=S, featureData=T ) 
	p=cumNormStatFast( obj )
	objTrim=cumNorm( obj, p=p )
	Treatment = pData( obj )$Treat
	settings = zigControl( maxit=20, verbose=TRUE )
	dsg1=model.matrix( ~0+Treat, data=s )
	res1=fitZig( obj=objTrim, mod=dsg1, control=settings, useCSSoffset = TRUE)
	zigFit1=res1@fit
	finalMod1=res1@fit$design
	c.mat1 = makeContrasts ( TreatnS	-	TreatcS, levels = finalMod1)
	fit1 = contrasts.fit( zigFit1, c.mat1 )
	fit1 = eBayes( fit1 )
	DT_1=fit1$coefficients
	DT_1p=fit1$p.value
	DT.zig<-topTable(fit1, number="inf",  adjust.method="BH" )
	DT.zig$rn<-rownames(DT.zig)
	DT.zig<-data.table(DT.zig, key="rn")
	DT=DT.taxa[ DT.zig, ]	
	DT$Enr <- ifelse ( DT$logFC > folch & DT$adj.P.Val < alpha, "Enr.nS", ifelse (DT$logFC < -folch & DT$adj.P.Val < alpha, "Enr.cS", "NS"))
	DT1[[j]]<-DT
	print ( "########################################################################")
	
#	computing differentially abundant ASVs	
	physeq.i=subset_samples( physeq.j, Treat %in% c("nS","bS") )	
	m=as( otu_table( physeq.i ), "matrix" ) + 1L
	t=data.frame( as( tax_table( physeq.i ), "matrix" ) )
	T=AnnotatedDataFrame( t )
	s=as( sample_data( physeq.i ), "data.frame" )
	S=AnnotatedDataFrame( s )
	obj=newMRexperiment( m, phenoData=S, featureData=T ) 
	p=cumNormStatFast( obj )
	objTrim=cumNorm( obj, p=p )
	Treatment = pData( obj )$Treat
	settings = zigControl( maxit=20, verbose=TRUE )
	dsg1=model.matrix( ~0+Treat, data=s )
	res1=fitZig( obj=objTrim, mod=dsg1, control=settings, useCSSoffset = TRUE)
	zigFit1=res1@fit
	finalMod1=res1@fit$design
	c.mat1 = makeContrasts ( TreatnS	-	TreatbS, levels = finalMod1)
	fit1 = contrasts.fit( zigFit1, c.mat1 )
	fit1 = eBayes( fit1 )
	DT_1=fit1$coefficients
	DT_1p=fit1$p.value
	DT.zig<-topTable(fit1, number="inf",  adjust.method="BH" )
	DT.zig$rn<-rownames(DT.zig)
	DT.zig<-data.table(DT.zig, key="rn")
	DT=DT.taxa[ DT.zig, ]	
	DT$Enr <- ifelse ( DT$logFC > folch & DT$adj.P.Val < alpha, "Enr.nS", ifelse (DT$logFC < -folch & DT$adj.P.Val < alpha, "Enr.bS", "NS"))
	DT2[[j]]<-DT
	print ( "########################################################################")
	
#	computing differentially abundant ASVs	
	physeq.i=subset_samples( physeq.j, Treat %in% c("cS","bS") )	
	m=as( otu_table( physeq.i ), "matrix" ) + 1L
	t=data.frame( as( tax_table( physeq.i ), "matrix" ) )
	T=AnnotatedDataFrame( t )
	s=as( sample_data( physeq.i ), "data.frame" )
	S=AnnotatedDataFrame( s )
	obj=newMRexperiment( m, phenoData=S, featureData=T ) 
	p=cumNormStatFast( obj )
	objTrim=cumNorm( obj, p=p )
	Treatment = pData( obj )$Treat
	settings = zigControl( maxit=20, verbose=TRUE )
	dsg1=model.matrix( ~0+Treat, data=s )
	res1=fitZig( obj=objTrim, mod=dsg1, control=settings, useCSSoffset = TRUE)
	zigFit1=res1@fit
	finalMod1=res1@fit$design
	c.mat1 = makeContrasts ( TreatcS	-	TreatbS, levels = finalMod1)
	fit1 = contrasts.fit( zigFit1, c.mat1 )
	fit1 = eBayes( fit1 )
	DT_1=fit1$coefficients
	DT_1p=fit1$p.value
	DT.zig<-topTable(fit1, number="inf",  adjust.method="BH" )
	DT.zig$rn<-rownames(DT.zig)
	DT.zig<-data.table(DT.zig, key="rn")
	DT=DT.taxa[ DT.zig, ]	
	DT$Enr <- ifelse ( DT$logFC > folch & DT$adj.P.Val < alpha, "Enr.cS", ifelse (DT$logFC < -folch & DT$adj.P.Val < alpha, "Enr.bS", "NS"))
	DT3[[j]]<-DT

	k = ifelse ( j == "H", "Healthy", "Fusarium") 
	print ( "###################################################################")
	print ( paste0("###### analysis done for samples in condition: ", k," ######"))
	print ( "###################################################################")
	}
# transformt counts to rel. abundance
	physeq.ra = transform_sample_counts(physeq.f, function(x) x/sum(x))
	
# computing arithmetic mean of rel. abund. by genotype
	DT.taxa=data.table( as(tax_table(physeq.f),"matrix"), keep.rownames=T, key="rn" )
	DT.ra=data.table(otu_table(physeq.ra), keep.rownames=T)
	DT.melt=data.table(melt(otu_table(physeq.ra)), keep.rownames=F)
	DT_sd=data.table(data.frame(sample_data(physeq.ra)), keep.rownames=T, key="rn")
	DT.merge=merge(DT.melt, DT_sd, by.x="Var2", by.y="rn")
	

	Enr.HnS=c(unique(DT1$H[DT1$H$Enr == "Enr.nS"]$rn, DT2$H[DT2$H$Enr == "Enr.nS"]$rn))
	Enr.HcS=c(unique(DT1$H[DT1$H$Enr == "Enr.cS"]$rn, DT3$H[DT3$H$Enr == "Enr.cS"]$rn))
	Enr.HbS=c(unique(DT2$H[DT2$H$Enr == "Enr.bS"]$rn, DT3$H[DT3$H$Enr == "Enr.bS"]$rn))
	
	LIST.H<-c(unique(DT.merge[DT.merge$Cond=="H",]$Var2))
	DT.ra$MeanH<-DT.ra[,.(rowMeans(.SD,na.rm=TRUE)),.SDcols = LIST.H]
	LIST.HnS<-c(unique(DT.merge[DT.merge$Concat=="HnS",]$Var2))
	DT.ra$MeanHnS<-DT.ra[,.(rowMeans(.SD,na.rm=TRUE)),.SDcols = LIST.HnS]	
	LIST.HcS<-c(unique(DT.merge[DT.merge$Concat=="HcS",]$Var2))
	DT.ra$MeanHcS<-DT.ra[,.(rowMeans(.SD,na.rm=TRUE)),.SDcols = LIST.HcS]
	LIST.HbS<-c(unique(DT.merge[DT.merge$Concat=="HbS",]$Var2))
	DT.ra$MeanHbS<-DT.ra[,.(rowMeans(.SD,na.rm=TRUE)),.SDcols = LIST.HbS]
	
	Uniq.H=unique(c(Enr.HnS,Enr.HbS, Enr.HcS))
	VENN=venn(list(Enr.HnS,Enr.HcS,Enr.HbS), zcolor="style", snames=c("HnS","HcS","HbS"))
	
	Enr.FnS=c(unique(DT1$F[DT1$F$Enr == "Enr.nS"]$rn, DT2$F[DT2$F$Enr == "Enr.nS"]$rn))
	Enr.FcS=c(unique(DT1$F[DT1$F$Enr == "Enr.cS"]$rn, DT3$F[DT3$F$Enr == "Enr.cS"]$rn))
	Enr.FbS=c(unique(DT2$F[DT2$F$Enr == "Enr.bS"]$rn, DT3$F[DT3$F$Enr == "Enr.bS"]$rn))
	
	Uniq.F=unique(c(Enr.FnS,Enr.FbS, Enr.FcS))
	VENN=venn(list(Enr.FnS,Enr.FcS,Enr.FbS), zcolor="style", snames=c("FnS","FcS","FbS"))
	
	LIST.H<-c(unique(DT.merge[DT.merge$Cond=="H",]$Var2))
	DT.ra$MeanH<-DT.ra[,.(rowMeans(.SD,na.rm=TRUE)),.SDcols = LIST.H]
	LIST.HnS<-c(unique(DT.merge[DT.merge$Concat=="HnS",]$Var2))
	DT.ra$MeanHnS<-DT.ra[,.(rowMeans(.SD,na.rm=TRUE)),.SDcols = LIST.HnS]	
	LIST.HcS<-c(unique(DT.merge[DT.merge$Concat=="HcS",]$Var2))
	DT.ra$MeanHcS<-DT.ra[,.(rowMeans(.SD,na.rm=TRUE)),.SDcols = LIST.HcS]
	LIST.HbS<-c(unique(DT.merge[DT.merge$Concat=="HbS",]$Var2))
	DT.ra$MeanHbS<-DT.ra[,.(rowMeans(.SD,na.rm=TRUE)),.SDcols = LIST.HbS]
	
	LIST.F<-c(unique(DT.merge[DT.merge$Cond=="F",]$Var2))
	DT.ra$MeanF<-DT.ra[,.(rowMeans(.SD,na.rm=TRUE)),.SDcols = LIST.F]
	LIST.FnS<-c(unique(DT.merge[DT.merge$Concat=="FnS",]$Var2))
	DT.ra$MeanFnS<-DT.ra[,.(rowMeans(.SD,na.rm=TRUE)),.SDcols = LIST.FnS]
	LIST.FcS<-c(unique(DT.merge[DT.merge$Concat=="FcS",]$Var2))
	DT.ra$MeanFcS<-DT.ra[,.(rowMeans(.SD,na.rm=TRUE)),.SDcols = LIST.FcS]
	LIST.FbS<-c(unique(DT.merge[DT.merge$Concat=="FbS",]$Var2))
	DT.ra$MeanFbS<-DT.ra[,.(rowMeans(.SD,na.rm=TRUE)),.SDcols = LIST.FbS]
	
	DT.H=DT.ra[,c("rn","MeanH","MeanHnS","MeanHcS","MeanHbS")]
	DT.F=DT.ra[,c("rn","MeanF","MeanFnS","MeanFcS","MeanFbS")]
	
	DT.H=merge(DT.H, DT.taxa, by.x="rn", by.y="rn") 
	DT.F=merge(DT.F, DT.taxa, by.x="rn", by.y="rn") 
	
	DT.H$ColPhylum <- ifelse(!DT.H$Phylum %in% ColPhylum, "Other", ifelse(DT.H$Phylum == "Bacteroidota", "Bacteroidetes", ifelse(DT.H$Phylum == "Firmicutes", "Firmicutes", ifelse (DT.H$Phylum == "Actinobacteriota", "Actinobacteria", "Proteobacteria"))))
	DT.F$ColPhylum <- ifelse(!DT.F$Phylum %in% ColPhylum, "Other", ifelse(DT.F$Phylum == "Bacteroidota", "Bacteroidetes", ifelse(DT.F$Phylum == "Firmicutes", "Firmicutes", ifelse (DT.F$Phylum == "Actinobacteriota", "Actinobacteria", "Proteobacteria"))))

	
	DT.H$Col <- ifelse(!DT.H$rn %in% Uniq.H, "NS", DT.H$ColPhylum)
	DT.F$Col <- ifelse(!DT.F$rn %in% Uniq.F, "NS", DT.F$ColPhylum)
	
	DT.H$Sum=rowSums(DT.H[,c(2:5)])
	DT.F$Sum=rowSums(DT.F[,c(2:5)])
	
	DT.H=DT.H[DT.H$Sum !=0]
	DT.F=DT.F[DT.F$Sum !=0]
	
	DT.H$Name=ifelse(DT.H$Col != "NS", DT.H$Genus, "")
	
	minRA=ifelse(min(DT.H$MeanH)<min(DT.F$MeanF), min(DT.H$MeanH), min(DT.F$MeanF))
	maxRA=ifelse(max(DT.H$MeanH)>max(DT.F$MeanF), max(DT.H$MeanH), max(DT.F$MeanF))
	LIM=c(minRA,maxRA)
	SEQ=seq(0.01, 0.07, 0.02)
	
	pp1=ggtern(DT.H, mapping=aes(MeanHnS, MeanHcS, MeanHbS))

	p1=pp1+
	geom_segment(data=lines, aes(x,y,z, xend=xend, yend=yend, zend=zend),color="black", size=1,linetype="twodash")+
	geom_point(aes(size=MeanH, color=Col), alpha=.8)+  
	scale_color_manual(values=Palette.phylum1)+ 
	scale_size_continuous(range = RAN,limits=LIM, breaks=SEQ) +
	theme_bw()+ theme_new +  labs(x="nS", y="cS", z="bS")

#	pdf("ternary_noninfest.pdf", useDingbats=FALSE)
	print(p1+gtitle("fig3D"))
#	dev.off()
	
	DT.F$Name=ifelse(DT.F$Col != "NS", DT.F$Genus, "")
	
	pp2=ggtern(DT.F, mapping=aes(MeanFnS, MeanFcS, MeanFbS))	
	
	p2=pp2 + 
	geom_segment(data=lines, aes(x,y,z, xend=xend, yend=yend, zend=zend),color="black", size=1,linetype="twodash")+
	geom_point(aes(size=MeanF, color=Col), alpha=.8)+  
	scale_color_manual(values=Palette.phylum2)+ 
	scale_size_continuous(range = RAN,limits=LIM, breaks=SEQ) +
	theme_bw()+ theme_new +  labs(x="nS", y="cS", z="bS")

#	pdf("ternary_Fo-infested.pdf", useDingbats=FALSE)
	print(p2gtitle("fig3D"))
#	dev.off()
	
#	gridExtra::grid.arrange(print(p1), print(p2), nrow=2, ncol=2)
	



