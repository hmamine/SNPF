#***require packages for the analysis the analysis
pkg=c("plyr", "ggplot2", "phyloseq", "data.table","agricolae", "metagenomeSeq", "ape", "vegan", "dplyr", "seqtime","picante","PMCMRplus")
lapply(pkg, library, character.only = TRUE)

rm(list = ls()[grep("*.*", ls())])
set.seed(131)
color_palette<-c("#b02b76","#ffa500")
color_palette1<-c("#ffa500")
color_palette2<-c("#b02b76")


theme_new <- theme (
	axis.title.y=element_blank(),
	axis.title.x=element_blank(),
	axis.ticks.x=element_blank(),
	axis.text.x = element_blank(),
#	axis.text.y = element_blank(),
#	legend.position="none",
#	axis.text.x = element_text(angle=90, vjust=1, size=9, color="black"),
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(),
	)
dist=( "bray" )

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

#computing BC distances
	DT.sd=data.table(as(sample_data(physeq.norm), "data.frame"),keep.rownames=T, key="rn")
	dist.BC=distance(physeq.norm, dist )
	melt.dist.BC=reshape2::melt(as.matrix( dist.BC ))

	DT.dist = melt.dist.BC %>%
   	filter(as.character(Var1) != as.character(Var2)) %>%
   	mutate_if(is.factor,as.character)
 	
   	DT.dist=data.table(DT.dist ,  key=c("Var1","Var2"))
   	DT.merge=merge( DT.dist, DT.sd, by.x="Var2", by.y="rn" )
	DT.merge=merge( DT.merge, DT.sd, by.x="Var1", by.y="rn" )
	
	DT.H=DT.merge[DT.merge$Cond.x=="H" & DT.merge$Cond.y=="H"  ]
	DT.H=DT.H[DT.H$Treat.x=="CTL" & DT.H$Treat.y != "CTL" ]
	pp1=ggplot(DT.H, aes(x=Concat.y, y=value))
	Ord1=c("HnS","HcS","HbS")
	pp1$data$Concat.y <- ordered(pp1$data$Concat.y, levels=Ord1 )
	p1=pp1+geom_violin(aes(fill=Cond.x))+geom_boxplot(width=0.1,outlier.size =0.5)+
	scale_fill_manual(values=color_palette1)+
	theme_bw()+ theme_new + ylab("BC distances")
	
	DT.F=DT.merge[DT.merge$Cond.x=="F" & DT.merge$Cond.y=="F"  ]
	DT.F=DT.F[DT.F$Treat.x=="CTL" & DT.F$Treat.y != "CTL" ]
	pp2=ggplot(DT.F, aes(x=Concat.y, y=value))
	Ord2=c("FnS","FcS","FbS")
	pp2$data$Concat.y <- ordered(pp2$data$Concat.y, levels=Ord2 )
	p2=pp2+geom_violin(aes(fill=Cond.x))+geom_boxplot(width=0.1,outlier.size =0.5)+
	scale_fill_manual(values=color_palette2)+
	theme_bw()+ theme_new + ylab("BC distances")
	
	DT.comb=rbind(DT.H,DT.F)
	pp3=ggplot(DT.comb, aes(x=Concat.y, y=value))
	Ord3=c("HnS","FnS","HcS","FcS","HbS","FbS")
	pp3$data$Concat.y <- ordered(pp3$data$Concat.y, levels=Ord3 )
	p3=pp3+geom_violin(aes(fill=Cond.x))+geom_boxplot(width=0.1,outlier.size =0.5)+
	scale_fill_manual(values=color_palette)+
	theme_bw()+ theme_new + ylab("BC distances")

#	pdf("BC_distances.pdf",paper="A4" ,useDingbats=FALSE)
	gridExtra::grid.arrange(p3, nrow=2, ncol=2)
#	dev.off()
	
#	aov.out<-aov(data=p3$data, value ~ Concat.y)
#	print(summary(aov.out))
#	print(TukeyHSD(aov.out))
#	print(HSD.test(aov.out,trt="Concat"))
	
#	print(kruskal.test(data=p3$data, value ~ Concat.y))
#	print(kwAllPairsConoverTest(value~as.factor(Concat.y),data=p3$data,  p.adjust.method="BH" ))











