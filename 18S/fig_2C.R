#***require packages for the analysis the analysis
pkg=c("plyr", "ggplot2", "phyloseq", "data.table", "reshape2", "ape", "phytools", "metagenomeSeq","PMCMRplus")
lapply(pkg, library, character.only = TRUE)

rm(list = ls()[grep("*.*", ls())])

color_palette<-c("#b02b76","#ffa500")
shape_palette=c(17,15,19,18)
shape_palette=c(0,1,2,4)

theme_new <- theme (
	panel.grid.major = element_blank(),
	legend.position="none",
	panel.grid.minor = element_blank(),
	)
meth1=("PCoA")
meth2=("CAP")
dist1=("bray")
																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																											
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
	
#normalization of count reads using CSS 
	otumat=as( otu_table( physeq1 ), "matrix" )
	mp=newMRexperiment( ( otumat ) )
	physeq.norm=phyloseq( otu_table( MRcounts( cumNorm( mp, p=cumNormStat( mp ) ),
	norm=TRUE, log=TRUE ),taxa_are_rows=TRUE ), NTAXA, NSD )

#computing Bray-Curtis distances
	dist.BC=distance( physeq.norm, dist1 )

#computing unconstrained PCoA
	pcoa.BC=ordinate( physeq.norm, meth1 , dist.BC)
	
	pp<-plot_ordination(physeq.norm, pcoa.BC, shape="Treat", color="Cond")
	pp$layers<-pp$layers[-1]
	
	p1=pp+geom_point(size=4)+theme_bw()+theme_new+
	scale_colour_manual(values=color_palette)+
	scale_shape_manual(values=shape_palette)
	
#	print(p1)
	
#	pdf("figure_pcoa_18S.pdf",useDingbats=FALSE)
	gridExtra::grid.arrange(p1, nrow=2, ncol=2)
#	dev.off()

	


