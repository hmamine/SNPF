#***require packages for the analysis the analysis
pkg=c("ggplot2", "phyloseq", "ape", "picante", "data.table","PMCMRplus", "agricolae")
lapply(pkg, library, character.only = TRUE)

rm(list = ls()[grep("*.*", ls())])

#tabOTU <- import_mothur(mothur_shared_file="asv.txt")
#write.table(tabOTU, "otu_table.txt", sep="\t")

alphaInd = c("Shannon", "Observed")
color_palette<-c("#b02b76","#ffa500")

set.seed(123)
theme_new <- theme (
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(),
	axis.title.y=element_blank(),
	axis.title.x=element_blank(),
	legend.position="none",
	axis.text.x = element_text(angle=90, vjust=1, size=10, color="black")
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
	
#Rarefication to even depth based on smallest sample size
rf=min(sample_sums(physeq1))
physeq.rrf=rarefy_even_depth(physeq1, rf, replace=TRUE, rngseed = 131)

#Ploting species richness	
	p_nat=plot_richness(physeq.rrf,"Concat","Cond" ,measures=alphaInd)
	p_nat$layers <- p_nat$layers[-1]
	
	Ord1=c("HCTL","FCTL","HnS","FnS","HcS","FcS","HbS","FbS")
	p_nat$data$Concat <- ordered(p_nat$data$Concat, levels=Ord1 )

	p1=p_nat+geom_boxplot(data=p_nat$data, aes(x=Concat, y=value, color=Cond))+
	geom_point(size=.5)+theme_bw()+
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	scale_colour_manual(values=color_palette)
	print(p1)

	subp2 = ggplot(data=p_nat$data[p_nat$data$variable=="Observed",], aes(x=Concat, y=value))+
	geom_boxplot(aes(fill=Cond), outlier.shape = 1, outlier.size= 1)+
	geom_point(size=.5)+theme_bw()+ theme_new +
	scale_fill_manual(values=color_palette) + 
	facet_wrap(~variable) 
	ylab("Rhichness")

#	aov.out<-aov(data=subp2$data, value ~ Concat)
#	print(summary(aov.out))
#	print(TukeyHSD(aov.out))
#	print(HSD.test(aov.out,trt="Concat"))
#	print(kruskal.test(data=subp2$data, value ~ Concat))
#	print(kwAllPairsConoverTest(value~as.factor(Concat),data=subp2$data, p.adjust.method="BH")) 

#Pairwise significance testing - observed CTUs
#	wilcox.test(data=subp2$data[subp2$data$Treat=="CTL",], value ~ Cond)
#	wilcox.test(data=subp2$data[subp2$data$Treat=="nS",], value ~ Cond)
#	wilcox.test(data=subp2$data[subp2$data$Treat=="cS",], value ~ Cond)
#	wilcox.test(data=subp2$data[subp2$data$Treat=="bS",], value ~ Cond)

	subp1 = ggplot(data=p_nat$data[p_nat$data$variable=="Shannon",], aes(x=Concat, y=value))+
	geom_boxplot(aes(fill=Cond), outlier.shape = 1, outlier.size= 1)+
	geom_point(size=.5)+theme_bw()+ theme_new +
	scale_fill_manual(values=color_palette) + 
	facet_wrap(~variable)+
	ylab("Shannon index")
	
#	aov.out<-aov(data=subp1$data, value ~ Concat)
#	print(summary(aov.out))
#	print(TukeyHSD(aov.out))
#	print(HSD.test(aov.out,trt="Concat"))
#	print(kruskal.test(data=subp2$data, value ~ Concat))
#	print(kwAllPairsConoverTest(value~as.factor(Concat),data=subp2$data, p.adjust.method="BH")) 
	
#Pairwise significance testing - Shannon index
#	wilcox.test(data=subp1$data[subp1$data$Treat=="CTL",], value ~ Cond)
#	wilcox.test(data=subp1$data[subp1$data$Treat=="nS",], value ~ Cond)
#	wilcox.test(data=subp1$data[subp1$data$Treat=="cS",], value ~ Cond)
#	wilcox.test(data=subp1$data[subp1$data$Treat=="bS",], value ~ Cond)


#	pdf("alpha_diversity.pdf",paper="A4" ,useDingbats=FALSE)
	gridExtra::grid.arrange(subp1, subp2, nrow=2, ncol=2)
#	dev.off()

	
 
            
         

