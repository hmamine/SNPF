#require packages for the analysis the analysis
pkg=c("ggplot2", "data.table", "reshape2","PMCMRplus","venn","dplyr","packcircles","scales","tidyr","igraph","RColorBrewer")
lapply(pkg, library, character.only = TRUE)

rm(list = ls()[grep("*.*", ls())])


list.color<-c("#b02b76","#ffa500")
list.color1<-c("#ffa500","#6c6b20","#b36200")
list.color2<-c("#0072ff","#398000","#800039")
theme_new <- theme (
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(),
	legend.position="none",
	axis.title.x=element_blank(),
	axis.title.y=element_blank(),
#	axis.text.y=element_blank(),
	axis.text.y=element_text(size=8, color="black"),
	axis.ticks.x=element_blank(),
	)

theme_new1 <- theme(
	panel.grid.major = element_blank(), 
	panel.grid.minor = element_blank(),
	axis.title.y=element_blank(),
	axis.title.x=element_blank(),
#	axis.text.y=element_blank(),
	axis.text.y=element_text(color="black", size=8)
	)

###################################################################################

	tab=read.table("table_16S18S.txt", sep="\t", row.names=1, header=1)
	

	
	Order1<-c("HCTL","FCTL","HnS","FnS","HcS","FcS","HbS","FbS")

	pp1<-ggplot(data=tab, aes(x=Concat, y=Value))

	pp1$data$Concat <- ordered(pp1$data$Concat, levels= Order1)

	p1=pp1+geom_boxplot(aes(color=Cond))+ geom_point(color="black", size=2)+
	geom_hline(yintercept = 0, linetype="longdash")+
	scale_color_manual(values=list.color)+ 
	theme_bw()+ theme_new 

	print(p1)
#	pdf("qPCR_ratio_16S_18S.pdf",useDingbats=FALSE)
#	gridExtra::grid.arrange(p1,nrow=2, ncol=2)
#	dev.off()

#	aov.out<-aov(data=p1$data, Value ~ Concat)
#	print(summary(aov.out))

#	print(kruskal.test(data=p1$data, Value ~ Concat))
	
	
	
	
	
	
	
