
#Comment next two lines out to install the EdgeR
#source("http://bioconductor.org/biocLite.R")
#biocLite("edgeR")
library("edgeR")
#citation("edgeR")

setwd("/home/ruolin/Research/WaspsProject/HTSeq-Count-BFAST-raw-reads-count")
G1=read.table("BFAST-G1-r1.0.txt")
G2=read.table("BFAST-G2-r1.0.txt")
G3=read.table("BFAST-G3-r1.0.txt")
G4=read.table("BFAST-G4-r1.0.txt")
G5=read.table("BFAST-G5-r1.0.txt")
S1=read.table("BFAST-S1-r1.0.txt")
S2=read.table("BFAST-S2-r1.0.txt")
S3=read.table("BFAST-S3-r1.0.txt")
S4=read.table("BFAST-S4-r1.0.txt")
NW1=read.table("BFAST-NW1-r1.0.txt")
NW2=read.table("BFAST-NW2-r1.0.txt")
NW3=read.table("BFAST-NW3-r1.0.txt")
NW4=read.table("BFAST-NW4-r1.0.txt")
NW5=read.table("BFAST-NW5-r1.0.txt")

print( "Reading input files finished!")

### creating input for edgeR
wasps_table=cbind(G1,G2,G3,G4,G5,NW1,NW2,NW3,NW4,NW5,S1,S2,S3,S4)
wasps_table=wasps_table[,c(1,2,4,6,8,10,12,14,16,18,20,22,24,26,28)]
length(unique(rownames(wasps_table)))
colnames(wasps_table)=c("gene_id","G1","G2","G3","G4","G5","NW1","NW2","NW3","NW4","NW5","S1","S2","S3","S4") ## set colnames
replicates=c(rep("gyne",5),rep("normal",5),rep("stylopized",4))
rownames(wasps_table)=wasps_table[,1]
wasps_table=wasps_table[,-1]
wasps_table=wasps_table[1:11506,]
sum(wasps_table[,1])
#filter<-function(df, min_count){
  #df = df[apply(df[, -1], MARGIN = 1, function(x) all(x >= min_count)), ]
  #return (df)
#}

#wasps_table=filter(wasps_table, 1)
### Runint edge for DE genes
wasps_d=DGEList(count=wasps_table,group=replicates)
cpm.wasps_d <- cpm(wasps_d)
#cpm.wasps_d
print(paste("number of transcripts before filtering", dim(wasps_table)[1]))
###Filtering gene with low read counts. 
wasps_d <- wasps_d[ rowSums(cpm.wasps_d>3) >=3,]
wasps_d <- calcNormFactors(wasps_d,)
wasps_table=wasps_d$count
print(paste("number of transcripts after filtering", dim(wasps_table)[1]))

wasps_d <- estimateCommonDisp(wasps_d, verbose=TRUE)
wasps_d <- estimateTagwiseDisp(wasps_d,verbose=TRUE)

#####DE number between Normal workers and queens
wasps_et <- exactTest(wasps_d,pair=c(1,2))
summary(decideTestsDGE(wasps_et, p.value=0.05))
NWandG=topTags(wasps_et, n=305)

####DE number between number workers and stylopized workers
wasps_et <- exactTest(wasps_d,pair=c(2,3))
summary(decideTestsDGE(wasps_et, p.value=0.05))
SandNW=topTags(wasps_et, n=51)

####DE number between queens and stylopized workers
wasps_et <- exactTest(wasps_d,pair=c(1,3))
summary(decideTestsDGE(wasps_et, p.value=0.05))
SandG=topTags(wasps_et, n=90)



#########
##### PCA 
#########
library(RColorBrewer)
two_inter = intersect( rownames(NWandG), rownames(SandNW))
three_inter = intersect(rownames(SandG), two_inter)

wasps_table[three_inter,]
pca= prcomp(wasps_table[three_inter,], scale.=T)$rotation
col4 = as.factor(c("G","G","G","G","G","NW","NW","NW","NW","NW","S","S","S","S"))
pca_wg = cbind(pca, col4)
pca_wg
plot(pca_wg[,1], pca_wg[,2], col=pca_wg[,4], pch=1)
legend('topright', legend=c("Gynes, Normal workers, Stylopized workers"), col=1:3)
warnings
#tab-delimited gene numbers
#write.table(NWandG,file="SandNW.txt", quote=F, sep="\t")
#write.table(SandG,file="SandG.txt", quote=F, sep="\t")
#write.table(SandNW,file="SandNW.txt", quote=F, sep="\t")



##################
### Venn Diagram############
################
require(VennDiagram)
require(RColorBrewer)

colors <- brewer.pal(3, "BuPu")
cat.colors = brewer.pal(3,"Set1")
venn.plot = draw.triple.venn(area1=305, 
                             area2=51, 
                             area3=90, 
                             n12=length(intersect(rownames(SandNW), rownames(NWandG))), 
                             n23=length(intersect(rownames(SandNW), rownames(SandG))), 
                             n13=length(intersect(rownames(SandG), rownames(NWandG))), 
                             n123=length(intersect(intersect(rownames(SandG), rownames(NWandG)),rownames(SandNW))),
                             category=c("NWvsG", "SWvsNW", "SWvsG"), fill=colors, cex = 2,cat.default.pos="text", cat.cex=2, 
                             cat.dist = -0.05, cat.col=cat.colors                
)
grid.draw(venn.plot)


### The rest of the code generated two excel files. 

### expressed_gene_lists.xlsx shows the everage reads counts for three groups: 
### Normal Workers (NW),  Stylopized Workers (SW) and Gynes(G)

### overall_DEG_list_addFDR.xlsx lists all differential expressed genes between the 3 comparisons
### NW <->SW,   SW <->G,  G<-> NW.

colnames(wasps_table)=c("G","G","G","G","G","NW","NW","NW","NW","NW","S","S","S","S")
t_table=t(wasps_table)### talbe transpose for using aggregate
agg <- aggregate(x = t_table, by = list(rownames(t_table)), FUN = "mean", na.rm = T)#### average
t_agg <- t(agg)##transpose back
t_agg <-(t_agg[-1,])

class(t_agg) <- "numeric"
colnames(t_agg)=c("read_count(G)","read_count(NW)","read_count(S)") ### rename each column  
### need to change 
SandG.list=t_agg[rownames(SandG),] #### in t_agg table pull out the DE genes between SandG
SandG.table=SandG$table
SandG.bind= cbind(SandG.list,SandG.table)
SandG.new=SandG.bind

colnames(SandG.new)[4:6]=c("logFD:NW vs G","logFD:S vs G","logFD:S vs NW")
SandG.new["FDR:S vs G"]=NA
SandG.new["FDR:S vs NW"]=NA
colnames(SandG.new)[7]="FDR:NW vs G"
SandG.new[,5]=SandG.new[,4]
SandG.new[,c(4,6)]=NA
SandG.new[,8]=SandG.new[,7]
SandG.new[,c(7,9)]=NA


NWandG.list=t_agg[rownames(NWandG),] #### in t_agg table pull out the DE genes between SandG
NWandG.table=NWandG$table
NWandG.bind= cbind(NWandG.list,NWandG.table)
NWandG.new=NWandG.bind
colnames(NWandG.new)[4:6]=c("logFD:NW vs G","logFD:S vs G","logFD:S vs NW")
NWandG.new["FDR:S vs G"]=NA
NWandG.new["FDR:S vs NW"]=NA
colnames(NWandG.new)[7]="FDR:NW vs G"
NWandG.new[,c(5,6)]=NA
NWandG.new[,c(8,9)]=NA

overlap1 = intersect(rownames(SandG), rownames(NWandG))
overlap2 = intersect(rownames(SandNW), rownames(NWandG))

NWandG.nonoverlap=NWandG.new[!rownames(NWandG.bind)%in%overlap1,]
NWandG.overlap=NWandG.new[rownames(NWandG.bind)%in%overlap1,]
names.overlap=rownames(NWandG.overlap)

overall.table=SandG.new
overall.table[names.overlap,4]=NWandG.overlap[,4]
overall.table[names.overlap,7]=NWandG.overlap[,7]
overall.table=rbind(overall.table,NWandG.nonoverlap)

SandNW.list=t_agg[rownames(SandNW),] #### in t_agg table pull out the DE genes between SandG
SandNW.table=SandNW$table
SandNW.bind= cbind(SandNW.list,SandNW.table)
SandNW.new=SandNW.bind
colnames(SandNW.new)[4:6]=c("logFD:NW vs G","logFD:S vs G","logFD:S vs NW")
SandNW.new["FDR:S vs G"]="NA"
SandNW.new["FDR:S vs NW"]="NA"
colnames(SandNW.new)[7]="FDR:NW vs G"

SandNW.new[,6]=SandNW.new[,4]
SandNW.new[,4:5]=NA
SandNW.new[,9]=SandNW.new[,7]
SandNW.new[,7:8]=NA

SandNW.nonoverlap=SandNW.new[!rownames(SandNW.new)%in%overlap2,]
SandNW.overlap=SandNW.new[rownames(SandNW.new)%in%overlap2,]
names.overlap2=rownames(SandNW.overlap)
overall.table[names.overlap2,6]=SandNW.overlap[,6]
overall.table[names.overlap2,9]=SandNW.overlap[,9]
overall.table=rbind(overall.table,SandNW.nonoverlap)

#install.packages("xlsx")
library(xlsx)
write.xlsx(overall.table,file="overall_DEG_list_addFDR.xlsx")
write.xlsx(t_agg,file="expressed_gene_lists.xlsx")


