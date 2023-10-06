#***require packages for the analysis the analysis
pkg=c("plyr", "ggplot2", "phyloseq", "data.table", "reshape2", "ape", "metagenomeSeq", "magrittr","RColorBrewer", "PMCMRplus", "agricolae")
lapply(pkg, library, character.only = TRUE)

rm(list = ls()[grep("*.*", ls())])

color_palette<-c("#b02b76","#ffa500","#00b2b2","#155832","#4f7b95")
list_taxa<-c("Fungi","Ciliophora","Chlorophyta","Streptophyta")
p<-vector("list")
DT.K<-vector("list")
Ord1<-c("HCTL","FCTL","HnS","FnS","HcS","FcS","HbS","FbS")
g_legend <- function(a.gplot){ 
    tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
    legend <- tmp$grobs[[leg]] 
    legend
} 
theme_new <- theme (
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(),
	legend.position="none",
	axis.text.y=element_text(size=9, colour="black"),
	axis.title.y=element_text(size=10, colour="black"),
#	axis.title.y=element_blank(),
#	axis.text.x = element_blank(),
#	axis.title.x = element_blank(),
	axis.ticks.x=element_blank(),
	)
	
# upload and prepare phyloseq objects***
	mat=read.table( "otu_table.txt", sep="\t", row.names=1, header=T)
	mat=as.matrix(mat)
	OTU=otu_table(mat, taxa_are_rows=T) 
	tax=read.table("taxonomy.txt", sep="\t", row.names=1, header=1)
	tax=as.matrix(tax)
	Taxa=tax_table(tax)
	sd=read.table("sample_data.txt", sep="\t", row.names=1, header=1)
	SD=sample_data(sd)
	
	physeq=phyloseq(OTU, Taxa, SD) 
	physeq1=subset_samples(physeq, Cond != "none")
###transform counts to rel. abund.
	physeq.RA = transform_sample_counts(physeq1, function(x) x/sum(x)) 
###agglomerate taxa to Kingdom rank	
	physeq.agg=tax_glom(physeq.RA, "Kingdom") 	
	m=as( otu_table( physeq.agg ), "matrix" )
	DT_m <- m %>% 
	melt(id.vars = "variable", value.name="relabund")
	DT.m=data.table(DT_m, keep.rownames=F, key="Var2") 
	DT.taxa=data.table( as(Taxa,"matrix"), keep.rownames=T, key="rn")
	DT.SD=data.table( as(SD,"data.frame"), keep.rownames=T, key="rn")
	
	DT=merge(DT.m, DT.taxa, by.x="Var1",by.y="rn")
	DT=merge(DT, DT.SD, by.x="Var2",by.y="rn")
	
	getPalette = colorRampPalette(brewer.pal(9, "Set1"))
	Palette.Genus<-getPalette(40)

	for (i in list_taxa)
	{  
		print (paste0("subset Kingdom: ",i))
 		DT1=DT[DT$Kingdom == i]
		pp1=ggplot(DT1, aes(x=Concat, y=relabund))
		pp1$data$Concat <- ordered(pp1$data$Concat, levels=Ord1 )
		p1=pp1+ 
		geom_boxplot(color="black")+ theme_bw()+ theme_new +
		ylab(paste0("Relative Abundance of ",i))
		print(p1)
		p[[i]]<-p1
		DT.K[[i]]<-DT1		
#	print(p1)
	}

#	pdf("RA_plot_1.pdf", useDingbats=FALSE)
	gridExtra::grid.arrange(p$Fungi,p$Ciliophora,p$Chlorophyta,nrow=2, ncol=2)
#	dev.off()
	
	
#	aov.out<-aov(data=p$Fungi$data, relabund ~ Concat)
#	print(summary(aov.out))
#	print(TukeyHSD(aov.out))
#	print(HSD.test(aov.out,trt="Concat"))
	print(kruskal.test(data=p$Fungi$data, relabund ~ Concat))
#	print(kwAllPairsConoverTest(value~as.factor(Concat),data=subp2$data, p.adjust.method="BH")) 

#	aov.out<-aov(data=p$Ciliophora$data, relabund ~ Concat)
#	print(summary(aov.out))
#	print(TukeyHSD(aov.out))
#	print(HSD.test(aov.out,trt="Concat"))
	print(kruskal.test(data=p$Ciliophora$data, relabund ~ Concat))
#	print(kwAllPairsConoverTest(relabund ~as.factor(Concat),data=p$Ciliophora$data, p.adjust.method="BH"))

#	aov.out<-aov(data=p$Chlorophyta$data, relabund ~ Concat)
#	print(summary(aov.out))
#	print(TukeyHSD(aov.out))
#	print(HSD.test(aov.out,trt="Concat"))
	print(kruskal.test(data=p$Chlorophyta$data, relabund ~ Concat))
#	print(kwAllPairsConoverTest(relabund ~as.factor(Concat),data=p$Chlorophyta$data, p.adjust.method="BH"))

#	aov.out<-aov(data=p$Streptophyta$data, relabund ~ Concat)
#	print(summary(aov.out))
#	print(TukeyHSD(aov.out))
#	print(HSD.test(aov.out,trt="Concat"))
	print(kruskal.test(data=p$Streptophyta$data, relabund ~ Concat))
#	print(kwAllPairsConoverTest(relabund ~as.factor(Concat),data=p$Chlorophyta$data, p.adjust.method="BH"))




























	
	
	

