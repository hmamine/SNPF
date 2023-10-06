#require packages for the analysis the analysis
pkg=c("plyr", "ggplot2", "phyloseq", "data.table", "ape", "PMCMRplus","UpSetR", "venn","seqtime","metagenomeSeq")
lapply(pkg, library, character.only = TRUE)

rm(list = ls()[grep("*.*", ls())])
set.seed(123)
Palette.phylum <-c("#F7A415","#88888a","#2BB065")
Palette.phylum <-c("#EF5656","#47B3DA","#F7A415","#c1c5c1","#55555f","#2BB065")
Palette.phylum1 <-c("#EF5656","#47B3DA","#c1c5c1","#55555f","#2BB065")
ColPhylum<-c("Actinobacteriota","Bacteroidota","Firmicutes","Proteobacteria")
color_palette<-c("#00b2b2","#b02b76","#155832", "#000000","#ffa500")

Shape.1=c(19,21,21)
DT<-vector( "list" )
#log2 fold change
	folch=2

#Adjusted p-value
	alpha=0.01

theme_new <- theme (
	panel.grid.minor = element_blank(),
	legend.position="none",
	axis.title.x=element_blank(),
	axis.title.y=element_blank(),
	axis.text.y=element_blank(),
#	axis.text.y=element_text(size=4, colour="black"),
#	axis.text.x=element_blank(),
	axis.ticks.x=element_blank(),
	panel.grid.major = element_blank()
	)
theme_new1 <- theme(
	panel.grid.major = element_blank(), 
	legend.position="none",
	axis.title.y=element_blank(),
	axis.title.x=element_blank(),
	panel.grid.minor = element_blank(),
	axis.text.x = element_blank(),
	axis.ticks.x=element_blank(),
#	axis.text.x = element_text(angle=90, vjust=1, size=9, color="black"),
	axis.text.y = element_text(size=8, color="black"),
	)
theme_new2 <- theme(
	legend.position="none",
	axis.title.y=element_blank(),
	axis.title.x=element_blank(),
	panel.grid.minor = element_blank(),
	panel.grid.major = element_blank(), 
	panel.border = element_blank(),
	)
Manhat<-function(mapping){
		p<- ggplot( DT.TAXA , mapping )
		p+geom_jitter( size=2 )+theme_bw( )+
		geom_hline(yintercept = -log2(alpha), linetype="longdash")+
		scale_color_manual(values=Palette.phylum)+
		scale_shape_manual(values=Shape.1)+
		scale_y_continuous(limits=c(0,70),breaks=c(0,10,20,30,40,50,60)) +
		theme_new1
	}	
Manhat1<-function(mapping){
		p<- ggplot( DT.TAXA , mapping )
		p+geom_jitter( size=2 )+theme_bw( )+
		geom_hline(yintercept = -log2(alpha), linetype="longdash")+
		scale_color_manual(values=Palette.phylum)+
		scale_shape_manual(values=Shape.1)+
#		scale_y_continuous(breaks=c(0,10,20,30,40,50,60)) +
		theme_new1
	}
Volcan<-function(mapping){
		p<- ggplot( DT.TAXA , mapping )
		p+geom_jitter( size=2 )+ theme_bw( )+
		scale_color_manual(values=Palette.phylum)+
		geom_hline(yintercept = -log2(alpha), linetype="longdash")+
		geom_vline(xintercept = c(-folch, folch), linetype="longdash")+
		scale_x_continuous( limits=c(-13, 13 ), breaks=c(-11,-9,-7,-5,-3,0,3,5,7,9,11)) +
		scale_y_continuous(limits=c(0,30),breaks=c(0,5,10,15,20,25)) +
		theme_new2
	}
Volcan1<-function(mapping){
		p<- ggplot( DT.TAXA , mapping )
		p+geom_jitter( size=2 )+ theme_bw( )+
		scale_color_manual(values=Palette.phylum1)+
		geom_hline(yintercept = -log2(alpha), linetype="longdash")+
		geom_vline(xintercept = c(-folch, folch), linetype="longdash")+
		scale_x_continuous( limits=c(-13, 13 ), breaks=c(-11,-9,-7,-5,-3,0,3,5,7,9,11)) +
		scale_y_continuous(limits=c(0,30),breaks=c(0,5,10,15,20,25)) +
		theme_new2
	}
	
BiPlot<-function(mapping){
		p<- ggplot( DT.TAXA , mapping )
		p+
		geom_hline(yintercept = 0, color="black")+
		geom_vline(xintercept = 0, color="black")+
		geom_point( size=3 )+ theme_bw( )+
		scale_color_manual(values=Palette.phylum)+
		scale_x_continuous(limits=c(-13, 13 ), breaks=c(-11,-9,-7,-5,-3,-1,0,1,3,5,7,9,11)) +
		scale_y_continuous(limits=c(-13, 13 ), breaks=c(-11,-9,-7,-5,-3,-1,0,1,3,5,7,9,11)) +
		theme_new2
	}

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
	LIST2<-list ("nS","cS","bS")
for ( j in LIST1 )
	{
	print( j )
	physeq.j=subset_samples( physeq.f, Cond == j)
	
		for( i in LIST2 ) 
		{
			print( i )
			physeq.i=subset_samples( physeq.j, Treat %in% c("CTL",i) )
		
#Computing differentially abundant ASVs
			m=as( otu_table( physeq.i ), "matrix" ) + 1L
			t=data.frame( as( tax_table( physeq.i ), "matrix" ) )
			T=AnnotatedDataFrame( t )
			s=as( sample_data( physeq.i ), "data.frame" )
			S=AnnotatedDataFrame( s )
			obj=newMRexperiment( m, phenoData=S, featureData=T ) 
			p=cumNormStatFast( obj )
			objTrim=cumNorm( obj, p=p )
			Treatment = pData( obj )$Cond1
			settings = zigControl( maxit=20, verbose=TRUE )
			dsg1=model.matrix( ~0+Cond1, data=s )
			res1=fitZig( obj=objTrim, mod=dsg1, control=settings, useCSSoffset = TRUE)
			zigFit1=res1@fit
			finalMod1=res1@fit$design
			c.mat1 = makeContrasts ( Cond1Ctl	-	Cond1Treat, levels = finalMod1)
			fit1 = contrasts.fit( zigFit1, c.mat1 )
			fit1 = eBayes( fit1 )
			DT_1=fit1$coefficients
			DT_1p=fit1$p.value

			DT.zig<-topTable(fit1, number="inf",  adjust.method="BH" )
			DT.zig$rn<-rownames(DT.zig)
			DT.zig<-data.table(DT.zig, key="rn")
			DT1=DT.taxa[ DT.zig, ]
		
			DT1$ColPhylum <- ifelse(!DT1$Phylum %in% ColPhylum, "Other", ifelse(DT1$Phylum == "Bacteroidota", "Bacteroidetes", ifelse(DT1$Phylum == "Firmicutes", "Firmicutes", ifelse (DT1$Phylum == "Actinobacteriota", "Actinobacteria", "Proteobacteria"))))

			DT1$FC <- ifelse ( DT1$logFC > folch & DT1$adj.P.Val < alpha, "Dep", ifelse (DT1$logFC < -folch & DT1$adj.P.Val < alpha, "Enr", "NS"))
			DT1$Col<- ifelse(DT1$FC =="NS", "NS", DT1$ColPhylum)

	DT[[j]][[i]]<-DT1
	}
	k = ifelse ( j == "H", "Healthy", "Fusarium") 
	print ( "############################################################")
	print ( paste0("###### analysis done for samples in condition: ", k," ######"))
	print ( "############################################################")
	}

	DT.HnS=DT$H$nS
	nS.Hpos=DT.HnS[DT.HnS$FC=="Enr"]$rn
	nS.Hneg=DT.HnS[DT.HnS$FC=="Dep"]$rn
	
	DT.FnS=DT$F$nS
	nS.Fpos=DT.FnS[DT.FnS$FC=="Enr"]$rn
	nS.Fneg=DT.FnS[DT.FnS$FC=="Dep"]$rn
	
	DT.HcS=DT$H$cS
	cS.Hpos=DT.HcS[DT.HcS$FC=="Enr"]$rn
	cS.Hneg=DT.HcS[DT.HcS$FC=="Dep"]$rn
	
	DT.FcS=DT$F$cS
	cS.Fpos=DT.FcS[DT.FcS$FC=="Enr"]$rn
	cS.Fneg=DT.FcS[DT.FcS$FC=="Dep"]$rn
	
	DT.HbS=DT$H$bS
	bS.Hpos=DT.HbS[DT.HbS$FC=="Enr"]$rn
	bS.Hneg=DT.HbS[DT.HbS$FC=="Dep"]$rn
	
	DT.FbS=DT$F$bS
	bS.Fpos=DT.FbS[DT.FbS$FC=="Enr"]$rn
	bS.Fneg=DT.FbS[DT.FbS$FC=="Dep"]$rn
	
	DT.taxa=data.table( as(tax_table(physeq.f),"matrix"), keep.rownames=T, key="rn" )
	DT.taxa$ColPhylum <- ifelse(!DT.taxa$Phylum %in% ColPhylum, "Other", ifelse(DT.taxa$Phylum == "Bacteroidota", "Bacteroidetes", ifelse(DT.taxa$Phylum == "Firmicutes", "Firmicutes", ifelse (DT.taxa$Phylum == "Actinobacteriota", "Actinobacteria", "Proteobacteria"))))
		
	nS.neg=venn(list(nS.Hneg,nS.Fneg), zcolor="style", snames=c("HN","FN"))
	attr.nSneg<-attr(nS.neg,"intersection")
	
	nS.pos=venn(list(nS.Hpos,nS.Fpos), zcolor="style", snames=c("HP","FP"))
	attr.nSpos<-attr(nS.pos,"intersection")
	
	cS.neg=venn(list(cS.Hneg,cS.Fneg), zcolor="style", snames=c("HN","FN"))
	attr.cSneg<-attr(cS.neg,"intersection")
	
	cS.pos=venn(list(cS.Hpos,cS.Fpos), zcolor="style", snames=c("HP","FP"))
	attr.cSpos<-attr(cS.pos,"intersection")
	
	bS.neg=venn(list(bS.Hneg,bS.Fneg), zcolor="style", snames=c("HN","FN"))
	attr.bSneg<-attr(bS.neg,"intersection")
	
	bS.pos=venn(list(bS.Hpos,bS.Fpos), zcolor="style", snames=c("HP","FP"))
	attr.bSpos<-attr(bS.pos,"intersection")
	
	H.pos=venn(list(nS.Hpos,cS.Hpos,bS.Hpos), zcolor="style", snames=c("nS","cS","bS"))
	attr.Hpos<-attr(H.pos,"intersection")

	H.neg=venn(list(nS.Hneg,cS.Hneg,bS.Hneg), zcolor="style", snames=c("nS","cS","bS"))
	attr.Hneg<-attr(H.neg,"intersection")
	
	F.pos=venn(list(nS.Fpos,cS.Fpos,bS.Fpos), zcolor="style", snames=c("nS","cS","bS"))
	attr.Fpos<-attr(F.pos,"intersection")

	F.neg=venn(list(nS.Fneg,cS.Fneg,bS.Fneg), zcolor="style", snames=c("nS","cS","bS"))
	attr.Fneg<-attr(F.neg,"intersection")

	V1=venn(list(attr.Hneg$`nS:cS:bS`,attr.Fneg$`nS:cS:bS`), zcolor="style", snames=c("Dep_ctl","Dep_Fus"))
	attr.V1<-attr(V1,"intersection")
	DT.tmp=DT.taxa[DT.taxa$rn %in% unique(c(attr.Hneg$`nS:cS:bS`,attr.Fneg$`nS:cS:bS`)),]
	write.table(DT.tmp, "Dif_Dep.txt", sep="\t")
	
	V2=venn(list(attr.Hpos$`nS:cS:bS`,attr.Fpos$`nS:cS:bS`), zcolor="style", snames=c("Enr_ctl","Enr_Fus"))
	attr.V2<-attr(V2,"intersection")
	DT.tmp=DT.taxa[DT.taxa$rn %in% unique(c(attr.Hpos$`nS:cS:bS`,attr.Fpos$`nS:cS:bS`)),]

	LIST.nS<-unique(c(attr.V1$'Dep_ctl:Dep_Fus', attr.V2$'Enr_ctl:Enr_Fus',attr.nSneg$'HN:FN', attr.nSpos$'HP:FP'))
	
	LIST.cS<-unique(c(attr.V1$'Dep_ctl:Dep_Fus', attr.V2$'Enr_ctl:Enr_Fus',attr.cSneg$'HN:FN', attr.cSpos$'HP:FP'))
	
	LIST.bS<-unique(c(attr.V1$'Dep_ctl:Dep_Fus', attr.V2$'Enr_ctl:Enr_Fus',attr.bSneg$'HN:FN', attr.bSpos$'HP:FP'))
	

	DT.HnS$Gname<-ifelse(DT.HnS$rn %in% LIST.nS, DT.HnS$Genus, "")	
	DT.HcS$Gname<-ifelse(DT.HcS$rn %in% LIST.cS, DT.HcS$Genus, "")	
	DT.HbS$Gname<-ifelse(DT.HbS$rn %in% LIST.bS, DT.HcS$Genus, "")	
	DT.FnS$Gname<-ifelse(DT.FnS$rn %in% LIST.nS, DT.FnS$Genus, "")	
	DT.FcS$Gname<-ifelse(DT.FcS$rn %in% LIST.cS, DT.FcS$Genus, "")	
	DT.FbS$Gname<-ifelse(DT.FbS$rn %in% LIST.bS, DT.FbS$Genus, "")
	
	DT.nS<-DT.taxa[DT.HnS[,c("rn","logFC","adj.P.Val")],] [DT.FnS[,c("rn","logFC","adj.P.Val")],]
	setnames(DT.nS, c("logFC","adj.P.Val","i.logFC","i.adj.P.Val"), c("FCH","pvalH","FCF","pvalF"))
	DT.nS$Gname<-ifelse(DT.nS$rn %in% LIST.nS, DT.nS$Genus, "")
	DT.nS$Col<-ifelse(DT.nS$pvalH > alpha & DT.nS$pvalF > alpha, "NS",DT.nS$ColPhylum)
	DT.TAXA<-DT.nS
	bp.nS=BiPlot(aes(x=FCH, y=FCF, color=Col)) + geom_text(aes(label=Gname),size=3, color="black")
	
	
	DT.cS<-DT.taxa[DT.HcS[,c("rn","logFC","adj.P.Val")],] [DT.FcS[,c("rn","logFC","adj.P.Val")],]
	setnames(DT.cS, c("logFC","adj.P.Val","i.logFC","i.adj.P.Val"), c("FCH","pvalH","FCF","pvalF"))
	DT.cS$Gname<-ifelse(DT.cS$rn %in% LIST.cS, DT.cS$Genus, "")
	DT.cS$Col<-ifelse(DT.cS$pvalH > alpha & DT.cS$pvalF > alpha, "NS",DT.cS$ColPhylum)
	DT.TAXA<-DT.cS
	bp.cS=BiPlot(aes(x=FCH, y=FCF, color=Col)) + geom_text(aes(label=Gname),size=3, color="black")
	
	
	DT.bS<-DT.taxa[DT.HbS[,c("rn","logFC","adj.P.Val")],] [DT.FbS[,c("rn","logFC","adj.P.Val")],]
	setnames(DT.bS, c("logFC","adj.P.Val","i.logFC","i.adj.P.Val"), c("FCH","pvalH","FCF","pvalF"))
	DT.bS$Gname<-ifelse(DT.bS$rn %in% LIST.bS, DT.bS$Genus, "")
	DT.bS$Col<-ifelse(DT.bS$pvalH > alpha & DT.bS$pvalF > alpha, "NS",DT.bS$ColPhylum)
	DT.TAXA<-DT.bS
	bp.bS=BiPlot(aes(x=FCH, y=FCF, color=Col)) + geom_text(aes(label=Gname),size=3, color="black")
	
	gridExtra::grid.arrange(bp.nS, bp.cS, bp.bS, ncol=2, nrow=2)
	
#	pdf("biplot_nS.pdf", useDingbats=FALSE)
#	print(bp.nS)
#	dev.off()

break()

	DT.TAXA<-DT.HnS
#	p1<-Manhat( aes( y=-log2(adj.P.Val), x=ColPhylum, color = Col, shape=FC )) + ggtitle("nS in CTL") + geom_text(aes(label=Gname))
#	p2<-Volcan( aes( y=-log2(adj.P.Val), x=logFC, color = Col )) + ggtitle("nS in CTL") + geom_text(aes(label=Gname),size=3, color="black")
	
	DT.TAXA<-DT.HcS
#	p3<-Manhat( aes( y=-log2(adj.P.Val), x=ColPhylum, color = Col, shape=FC )) + ggtitle("cS in CTL") + geom_text(aes(label=Gname))
#	p4<-Volcan( aes( y=-log2(adj.P.Val), x=logFC, color = Col )) + ggtitle("cS in CTL") + geom_text(aes(label=Gname), size=3, color="black")
	
	DT.TAXA<-DT.HbS
#	p5<-Manhat( aes( y=-log2(adj.P.Val), x=ColPhylum, color = Col, shape=FC )) + ggtitle("bS in CTL") + geom_text(aes(label=Gname))
#	p6<-Volcan1( aes( y=-log2(adj.P.Val), x=logFC, color = Col )) + ggtitle("bS in CTL") + geom_text(aes(label=Gname), size=3, color="black")
	
	DT.TAXA<-DT.FnS
#	p7<-Manhat( aes( y=-log2(adj.P.Val), x=ColPhylum, color = Col, shape=FC ))  + ggtitle("nS in Fus") + geom_text(aes(label=Gname))
#	p8<-Volcan( aes( y=-log2(adj.P.Val), x=logFC, color = Col )) + ggtitle("nS in Fus") + geom_text(aes(label=Gname), size=3, color="black")
	 
	DT.TAXA<-DT.FcS
#	p9<-Manhat( aes( y=-log2(adj.P.Val), x=ColPhylum, color = Col, shape=FC )) + ggtitle("cS in Fus") + geom_text(aes(label=Gname))
#	p10<-Volcan1( aes( y=-log2(adj.P.Val), x=logFC, color = Col )) + ggtitle("cS in Fus") + geom_text(aes(label=Gname), size=3, color="black")
	
	DT.TAXA<-DT.FbS
#	p11<-Manhat( aes( y=-log2(adj.P.Val), x=ColPhylum, color = Col, shape=FC )) + ggtitle("bS in Fus") + geom_text(aes(label=Gname))
#	p12<-Volcan( aes( y=-log2(adj.P.Val), x=logFC, color = Col )) + ggtitle("bS in Fus") + geom_text(aes(label=Gname), size=3, color="black")
	

#	pdf("manhattan_plots.pdf",useDingbats=FALSE)
#	gridExtra::grid.arrange(p1,p7,p3,p9,p5,p11,  ncol=2, nrow=3)
#	dev.off()

#	pdf("volcano_plots.pdf",useDingbats=FALSE)
#	gridExtra::grid.arrange(p2,p8,p4,p10,p6,p12,  ncol=2, nrow=3)
#	dev.off()



