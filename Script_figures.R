library(VariantAnnotation)
library(StructuralVariantAnnotation)
library(stringr)
library(ggplot2)
library(reshape2)

#### function
cbind.fill <- function(...){
  nm <- list(...) 
  nm <- lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow)) 
  do.call(cbind, lapply(nm, function (x) 
    rbind(x, matrix(, n-nrow(x), ncol(x))))) 
}
# simpleEventType <- function(gr) {
#   return(ifelse(seqnames(gr) != seqnames(partner(gr)), "ITX", # inter-chromosomosal
#                 ifelse(gr$insLen >= abs(gr$svLen) * 0.7, "INS", # TODO: improve classification of complex events
#                        ifelse(strand(gr) == strand(partner(gr)), "INV",
#                               ifelse(xor(start(gr) < start(partner(gr)), strand(gr) == "-"), "DEL",
#                                      "DUP")))))
# }

### load data 
grlist <- paste("/DATA/",list.files("/DATA/"),sep = "")

grlist <- c("WT 1","WT 2","WT 3","WT 4","DM 1","DM 2","DM 3","DM 4",
            "WT UV 1","WT UV 2","DM UV 1","DM UV 2",
            "KR 1","KR 2","KR UV 1","KR UV 2",
            "Rev1 1","Rev1 2","Rev1 UV 1", "Rev1 UV 2")

all.chr <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chrY","chrX")

### load SV vcfs 
for(i in 1:length(grlist)){
  vcf <- VariantAnnotation::readVcf(gridsslist[i], "mm10")
  assign(grlist[i],vcf)
  print(i)
}

##################### Figure 4A ##################### 

### make DF with SVs
alltypes <- data.frame(Del=numeric(),Dup=numeric(),Ins=numeric(),Inv=numeric(),Itx=numeric())
for(i in 1:length(grlist)){
  gr <- get(grlist[i])
  gr <- breakpointRanges(gr)
  gr <- gr[gr$FILTER == "PASS" & partner(gr)$FILTER == "PASS"]
  gr <- gr[gr@strand == "+",]
  gr <- gr[gr@seqnames %in% all.chr]
  gr <- gr[gr$QUAL >= 1000]
  gr <- gr[gr$partner %in% names(gr)]
  assign(paste(grlist[i],"gr",sep=""),gr)
  svtype <- simpleEventType(gr)
  alltypes[i,] <- table(svtype)
  print(i)
}
rownames(alltypes) <- grlist

ggplot(data=all, aes(x=Var1, y=value, fill=Var2)) + geom_bar(stat="identity", position=position_dodge())



##################### Figure 4B-E ##################### 

### get deletion sizes, deletion insertions size and homology length
all_del_info <- size <- delins <- homlen <- delinseq <- NULL
for(i in 1:length(grlist)){
  gr <- get(paste(grlist[i],"gr",sep=""))
  svtype <- simpleEventType(gr)
  del <- gr[svtype == "DEL",]
  delp <- del[del@strand == "+",]
  size[[i]] <- delp$svLen
  delins[[i]] <- delp$insLen
  homlen[[i]] <- delp$HOMLEN
  delinseq[[i]] <- delp$insSeq
  geno <- rep(grlist[i],length(delp))
  all_del_info[[i]] <- deletioninfo
}

### For the density plots (use different replicates for the different plots)
WT_comb <- as.data.frame(unlist(Del_sizes[1:4]))
DM_comb <- as.data.frame(unlist(Del_sizes[5:8]))
WT_UV <- as.data.frame(unlist(Del_sizes[9:10]))
DM_UV <- as.data.frame(unlist(Del_sizes[11:12]))
both_comb <- cbind.fill(WT_comb,DM_comb,DM_CST,DM_BCL)
colnames(both_comb) <- c("WT","DM","WT UV","DM UV")
meltboth <- melt(both_comb)

ggplot(data = meltboth, aes(x=value, fill=Var2)) + 
  geom_density(alpha = 0.5) +
  labs(x="Log2 deletion size") +
  theme_minimal()


### get the numbers of deletions in specific segments of sizes. (Bar graphs)
Del_sizes <- do.call(cbind.fill,size)
colnames(Del_sizes) <- grlist
# both_real <- both_real[,c(1,2,5,6,9,3,4,7,8,10)]
Del_sizes <- Del_sizes*-1

numdel <- data.frame(matrix(0, nrow = 3, ncol = 20, dimnames = list(NULL, paste0(grlist))) )
for(t in 1:20){
  getit <- Del_sizes[,t][Del_sizes[,t] < 362]
  numdel[1,t] <- length(na.omit(getit))
  getit <- Del_sizes[,t][Del_sizes[,t] > 362 & Del_sizes[,t] < 4096] # type 3 deletions
  numdel[2,t] <- length(na.omit(getit))
  getit <- Del_sizes[,t][Del_sizes[,t] > 4096]
  numdel[3,t] <- length(na.omit(getit))
}

##################### Figure 4F ##################### 

### Creates dataframe with the number of breakpoints per number of homologous basepairs per genotype
mh_df <- data.frame(Var1 = 0:100)
for(i in 1:length(all_del_info)){
  temp <- all_del_info[[i]]
  temp_mh <- as.data.frame(table(temp$Homology_length)/nrow(temp))
  mh_df <- merge(mh_df,temp_mh, by="Var1", all =T)
}
colnames(mh_df) <- c("bp MH",grlist)
write.table(mh_df, "Percentage MH per base SM samples.txt", quote = F, sep="\t")


##################### Figure 5A ##################### 

### First move all SV data files to a seperate dir

allgrids <- list.files("/DATA/HW/all_SV/")
tumortype <- read.table(file = "/DATA/HW/metadata/metadata.tsv", sep = '\t', header = T)
tissueID <- tumortype$sampleId
tissuegrid <- unique(grep(paste(tissueID,collapse="|"), allgrids, value=TRUE))
tissuegridpath <- paste("/DATA/HW/all_SV/", tissuegrid, sep="/")

### Get, and make an object of, the SV deletion information of each tumor type
alltype <- unique(tumortype$primaryTumorLocation)
for(ai in 1:length(alltype)){
  type <- alltype[ai]
  tissueID <- tumortype$sampleId[tumortype$primaryTumorLocation %in% type]
  tissuegrid <- unique(grep(paste(tissueID,collapse="|"), allgrids, value=TRUE))
  tissuegridpath <- paste("/DATA/HW/all_SV/", tissuegrid, sep="/")
  all_del_info <- size <- delins <- homlen <- delinseq <- NULL
  for(i in 1:length(tissuegrid)){
    # Load vcf
    gr <- VariantAnnotation::readVcf(tissuegridpath[i], "mm10")
    # Make Granges
    gr <- breakpointRanges(gr)
    gr <- gr[gr$FILTER == "PASS" & partner(gr)$FILTER == "PASS"]
    # sometimes NA have to remove
    temp <- gr$QUAL >= 1000
    temp <- replace(temp, is.na(temp), FALSE)
    gr <- gr[temp]
    gr <- gr[gr$partner %in% names(gr)] #filter orphaned break partners
    svtype <- simpleEventType(gr)
    del <- gr[svtype == "DEL",]
    delp <- del[del@strand == "+",] # Only take one side of the deletion (for the distribution)
    size[[i]] <- delp$svLen
    delins[[i]] <- delp$insLen
    homlen[[i]] <- delp$HOMLEN
    delinseq[[i]] <- delp$insSeq
    geno <- rep(tissuegrid[i],length(delp))
    deletioninfo <- data.frame(Deletion_size=size[[i]],Deletion_insertion=delins[[i]], Homology_length=homlen[[i]], Deletion_insertion_seq=delinseq[[i]], Genotype=geno)
    all_del_info[[i]] <- deletioninfo
  }
  assign(paste("all_del_info_",type,sep=""),all_del_info)
}

### Get the deletion sizes, handy for later
# assign("tumor_del", get(paste("all_del_info_",type,sep="")))
# size <- NULL
# for(i in 1:length(tumor_del)){
#   size[[i]] <- tumor_del[[i]]$Deletion_size
# }
# allsize <- do.call(cbind.fill,size)
# colnames(allsize) <- tissuegrid
# allsize <- as.data.frame(allsize)
# allsize <- allsize * -1

### replace type with the types in figure
type <- "Lung"
tissueID <- tumortype$sampleId[tumortype$primaryTumorLocation %in% type]
tissuegrid <- unique(grep(paste(tissueID,collapse="|"), allgrids, value=TRUE))

assign("tumor_del", get(paste("all_del_info_",type,sep="")))
size <- NULL
for(i in 1:length(tumor_del)){
  size[[i]] <- tumor_del[[i]]$Deletion_size
}

allsize <- do.call(cbind.fill,size)
colnames(allsize) <- tissuegrid
allsize <- as.data.frame(allsize)
allsize <- allsize * -1

### get randoms
random_lung <- sample(allsize,50)
HW <- log2(random_lung)
HW_comb <- as.data.frame(unlist(HW[1:50]))

both_comb <- cbind.fill(HW_comb,WT_comb,DM_comb)
colnames(both_comb) <- c("HW lung random","WT","DM")
meltboth <- melt(both_comb)

ggplot(data = meltboth, aes(x=value,fill=Var2)) + 
  geom_density(alpha = 0.5) +
  labs(x="Log2 deletion size") +
  theme_minimal() +
  geom_vline(xintercept = c(11),colour=c("red"), linetype = c("dashed"))


### for plotting a selection type 3 high tumors

relative <- NULL
for(ii in 1:ncol(allsize)){
  temp <- allsize[,ii]
  relative[[ii]] <- length(na.omit(temp[temp > 362 & temp < 4096]))/length(na.omit(temp))
}

relative <- as.character(relative)
names(relative) <- colnames(allsize)
ordered <- sort(relative,decreasing = TRUE)
allsize_peak3 <- allsize[names(ordered)]

### only take samples that have more then 10 deletions
b <- NULL
for(a in 1:ncol(allsize_peak3)){
  if(length(na.omit(allsize_peak3[,a])) < 10) b <- append(b,a)
}
allsize_peak3_fil <- allsize_peak3[-b]

### select tumors with high type 3 deletions for plotting
lung_peak3 <- allsize_peak3_fil[c(1:5)]
HW <- log2(lung_peak3)
HW_select_lung <- as.data.frame(unlist(HW[1:5]))

both_comb <- cbind.fill(HW_select_lung,WT_comb,DM_comb)
colnames(both_comb) <- c("HW lung selected","WT","DM")
meltboth <- melt(both_comb)

ggplot(data = meltboth, aes(x=value,fill=Var2)) + 
  geom_density(alpha = 0.5) +
  labs(x="Log2 deletion size") +
  theme_minimal() +
  geom_vline(xintercept = c(11),colour=c("red"), linetype = c("dashed"))


##################### Figure 5B ##################### 

total_dels <- NULL
numbers <- data.frame(matrix(nrow = 4, ncol = 5))
# Lung
# Esophagus
# Breast
# Colon/Rectum
# Biliary
all_types <- c("Breast","Lung","Esophagus","Colon/Rectum","Biliary")
for(i in 1:length(all_types)){
  type <- all_types[i]
  tissueID <- tumortype$sampleId[tumortype$primaryTumorLocation %in% type]
  tissuegrid <- unique(grep(paste(tissueID,collapse="|"), allgrids, value=TRUE))
  
  if(type == "Colon/Rectum"){
    tumor_del <- `all_del_info_colon rectum`
  } else {
    assign("tumor_del", get(paste("all_del_info_",type,sep="")))
  }
  size <- NULL
  for(iii in 1:length(tumor_del)){
    size[[iii]] <- tumor_del[[iii]]$Deletion_size
  }
  
  allsize <- do.call(cbind.fill,size)
  colnames(allsize) <- tissuegrid
  allsize <- as.data.frame(allsize)
  allsize <- allsize * -1
  
  relative <- NULL
  for(j in 1:ncol(allsize)){
    temp <- allsize[,j]
    relative[j] <- length(na.omit(temp[temp > 362 & temp < 4096]))/length(na.omit(temp))
  }
  names(relative) <- colnames(allsize)
  ordered <- sort(relative,decreasing = TRUE)
  allsize_peak3 <- allsize[names(ordered)]
  # }
  
  b <- NULL
  for(a in 1:ncol(allsize_peak3)){
    if(length(na.omit(allsize_peak3[,a])) < 10) b <- append(b,a)
  }
  if(length(b) > 0){
    allsize_peak3_fil <- allsize_peak3[-b]
  } else {
    allsize_peak3_fil <- allsize_peak3
  }
  ### makes DF with % and number of peak 3 deletions per tumor
  number3 <- total <- NULL
  for(jj in 1:ncol(allsize)){
    temp <- allsize[,jj]
    number3[jj] <- length(na.omit(temp[temp > 362 & temp < 4096]))
    total[jj] <- length(na.omit(temp))
  }
  names(number3) <- colnames(allsize)
  names(total) <- colnames(allsize)
  sum_peak3 <- merge(as.data.frame(number3),as.data.frame(relative),by="row.names")
  rownames(sum_peak3) <- sum_peak3$Row.names
  sum_peak3 <- merge(as.data.frame(total),sum_peak3,by="row.names")
  sum_peak3 <- sum_peak3[,-3]
  sum_peak3_order <-  sum_peak3[order(sum_peak3$relative, decreasing = T),]
  
  if(type == "Colon/Rectum"){
    write.table(sum_peak3_order, paste("number and percentage peak3 colon rectum.txt", sep=""),row.names = F, quote = F, sep="\t")
  } else {
    write.table(sum_peak3_order, paste("number and percentage peak3 ",type ,".txt", sep=""),row.names = F, quote = F, sep="\t")
  }
  

##################### Figure 5C ##################### 
  ### get tumors with % cutoff of peak 3 deletions
  percent <- c(.1,.2,.25,.3)
  for(ii in 1:length(percent)){
    select_x <- sum_peak3[sum_peak3$relative > percent[ii],1]
    
    HW <- log2(allsize_peak3_fil[,colnames(allsize_peak3_fil) %in% select_x])
    
    HW_comb <- as.data.frame(unlist(HW))
    
    both_comb <- cbind.fill(HW_comb,WT_comb,DM_comb)
    colnames(both_comb) <- c(paste("HW",type,percent[ii],"percent"),"WT","DM")
    # both_comb <- cbind.fill(HW_comb)
    # colnames(both_comb) <- c(paste("HW",type,percent[ii],"percent"))
    meltboth <- melt(both_comb)
    
    p <- ggplot(data = meltboth, aes(x=value,fill=Var2)) + 
      geom_density(alpha = 0.5) +
      labs(x="Log2 deletion size") +
      scale_fill_manual(values=c("#C77CFF", "#F8766D", "#00BFC4")) +
      theme_minimal()

    numbers[ii,i] <- length(select_x)
    
    
    if(type == "Colon/Rectum"){
      png(paste("Density HW only Colon rectum cutoff at ", percent[ii] ,".png", sep=""), width = 1500 ,height = 830)
    } else {
      png(paste("Density HW only ",type, " cutoff at ", percent[ii] ,".png", sep=""), width = 1500 ,height = 830)
    }
    
    print(p)
    dev.off()
  }
  total_dels[i] <- length(allsize_peak3_fil)
}


##################### Figure 5D ##################### 
### Get list of dna repair genes from MD Anderson (https://www.mdanderson.org/documents/Labs/Wood-Laboratory/human-dna-repair-genes.html)
allrepair <- read_xlsx(path = "DNA repair genes.xlsx", col_names =F)
allrepair <- allrepair$`allrepair$...1`
genes_selected <- gsub("\\s*\\([^\\)]+\\)","",allrepair)

allfiles <- list.files("/DATA/HW/somatics/")
CNV <- paste("/DATA/HW/somatics/",allfiles,"/purple",sep="")

total_del_amp <- NULL
for(i in 1:length(CNV)){
  amp <- del <- NULL
  ### grab the CNV gene file
  files <- list.files(CNV[i])
  temp <- files[!is.na(str_extract(files,"gene"))]
  temp <- read.table(paste(CNV[i],temp,sep="/"),header = T)
  temp_gene <- temp[temp$gene %in% genes_selected,]
  if(!any(temp_gene$minCopyNumber < 0.5 | temp_gene$maxCopyNumber > 5)){
    next
  }
  if(any(temp_gene$minCopyNumber < 0.5)){
    temp_del <- temp_gene[temp_gene$minCopyNumber < 0.5,c(4,5,6)]
    del <- cbind(temp_del,rep("DEL",nrow(temp_del)))
    colnames(del) <- c("Gene","minCopyNumber","maxCopyNumber","CNV")
  }
  if(any(temp_gene$maxCopyNumber > 5)){
    temp_amp <- temp_gene[temp_gene$maxCopyNumber > 20,c(4,5,6)]
    amp <- cbind(temp_amp,rep("AMP",nrow(temp_amp)))
    colnames(amp) <- c("Gene","minCopyNumber","maxCopyNumber","CNV")
  }
  del_amp <- as.data.frame(rbind(amp,del))
  del_amp$patient <- str_extract(files[!is.na(str_extract(files,"gene"))], "[^.]+")
  colnames(del_amp) <- c("Gene","minCopyNumber","maxCopyNumber","CNV","patient")
  total_del_amp <- rbind(total_del_amp,del_amp)
}
write.table(total_del_amp, "MD anderson genes AMP or DEL HW2.txt",row.names = F, quote = F, sep="\t")



##################### Figure 5E ##################### 

# top 6 AMP : BRIP1, RRM2B, SP011, NEIL2, POLB, EME1
# DEL: GTF2H5, BRCA2, PDS5B, ATRX, RAD51B
DEL <- table(total_del_amp$Gene[total_del_amp$CNV == "DEL"])
DEL <- DEL[order(-as.numeric(DEL))]

AMP <- table(moreAMP$Gene)
AMP <- AMP[order(-as.numeric(AMP))]

### Top 6 most frequent AMP/DEL genes, combine these 
### have to change the ATM variable from AMP/DEL inside the loop

AMP_deletions <- NULL
genes_final_selected <- names(AMP)[1:6]
genes_final_selected <- names(DEL)[2:7] # Dont take P53 bit different so take top 2-7

### run twice, once for DEL and for AMP
del_info_ATM_all <- NULL
for(i in 1:length(genes_final_selected)){
  gene <- genes_final_selected[i]
  
  # for AMP
  # ATM <- unique(moreAMP$patient[moreAMP$Gene == gene])
  # For DEL
  ATM <- unique(total_del_amp$patient[total_del_amp$Gene == gene & total_del_amp$CNV == "DEL"])
  
  del_info_ATM <- NULL
  for(ai in 1:length(alltype)){
    type <- alltype[ai]
    temp <- get(paste("all_del_info_",type,sep=""))
    tissueID <- tumortype$sampleId[tumortype$primaryTumorLocation %in% type]
    tissuegrid <- unique(grep(paste(tissueID,collapse="|"), allgrids, value=TRUE))
    names(temp) <- tissuegrid
    pattern <- paste(ATM, collapse = "|") 
    temp <- temp[grepl(pattern, names(temp))]
    del_info_ATM <- append(del_info_ATM,temp)
  }
  del_info_ATM_all[[i]] <- del_info_ATM
  # if(length(del_info_ATM) == 1){next}
  size <- hom_ATM <- list()
  for(ii in 1:length(del_info_ATM)){
    size[[ii]] <- del_info_ATM[[ii]]$Deletion_size
    temp1 <- del_info_ATM[[ii]]
    temp1$Deletion_size <- temp1$Deletion_size * -1
    temp_hom <- temp1$Homology_length[temp1$Deletion_size >= 724 & temp1$Deletion_size <= 3444]
    hom_ATM[[ii]] <- temp_hom
  }
  ATM_size <- do.call(cbind.fill,size)
  ATM_size_comb <- do.call(c,size)
  ATM_size_comb <- as.data.frame(ATM_size_comb)
  ATM_size_comb <- ATM_size_comb * -1
  
  ### For making taking the random samples in the same ratio (10x more)
  random_types <- table(tumortype$primaryTumorLocation[tumortype$sampleId %in% ATM])
  random_types <- as.data.frame(random_types*10)
  all_random_del <- NULL
  for(ii in 1:nrow(random_types)){
    type_temp <- random_types$Var1[ii]
    freq_temp <- random_types$Freq[ii]
    temp_del <- get(paste("all_del_info_",type_temp,sep=""))
    tissueID <- tumortype$sampleId[tumortype$primaryTumorLocation %in% type_temp]
    tissuegrid <- unique(grep(paste(tissueID,collapse="|"), allgrids, value=TRUE))
    names(temp_del) <- tissuegrid
    random_del <- temp_del[!names(temp_del) %in% names(del_info_ATM)]
    if(freq_temp > length(random_del)){
      temp_del_random <- random_del
      # print(paste("For",gene,freq_temp,type_temp,"should be picked from random but only",length(random_del),"avaiable",sep=" "))
    }else{
      temp_del_random <- sample(random_del,freq_temp)
    }
    all_random_del <- append(all_random_del,temp_del_random)
  }
  random_size <- list()
  for(iii in 1:length(all_random_del)){
    random_size[[iii]] <- all_random_del[[iii]]$Deletion_size
  }
  allrandom_size_comb <- do.call(c,random_size)
  allrandom_size_comb <- as.data.frame(allrandom_size_comb)
  allrandom_size_comb <- allrandom_size_comb * -1
  
  
  both <- cbind.fill(ATM_size_comb,allrandom_size_comb)
  colnames(both) <- c(paste(gene,status),paste(gene,"no",status))
  AMP_deletions[[i]] <- both
}
# Combine all data from AMP genes 
AMP_deletions_comb <- do.call("rbind", AMP_deletions)
colnames(AMP_deletions_comb) <- c("Repair gene AMP", "Random selected")

DEL_deletions_comb <- do.call("rbind", AMP_deletions)
colnames(DEL_deletions_comb) <- c("Repair gene DEL", "Random selected")

### combine DEL/AMP density
DEL_AMP_deletions_comb <- cbind.fill(DEL_deletions_comb, AMP_deletions_comb)


### make combined density plot
ggplot(data = meltsize, aes(x=value,fill=`CN status`)) +
  geom_density(alpha = 0.5) +
  labs(x="Log2 deletion size") +
  theme_classic() +
  theme_minimal() +
  geom_vline(xintercept = c(11),colour=c("red"), linetype = c("dashed"))


##################### Figure 5F (and supplementary figure 9) ##################### 
treatment <- read.delim("/DATA/HW/metadata/metadata.tsv",sep="\t")

types_t <- unique(treatment$treatmentType)

SV_names <- str_extract(list.files("/DATA/HW/all_SV/"), "[^.]+")
names(all_del_info) <- SV_names


all_types <- c("Breast","Lung","Esophagus","Colon/Rectum","Biliary")

treatment_type <- treatment[treatment$treatmentType %in% types_t[1],]
HW <- all_del_info[names(all_del_info) %in% treatment_type$sampleId]
HW <- do.call(rbind,HW)
HW$Deletion_size <- HW$Deletion_size * -1
control <- HW

p <- NULL
for(i in 1:length(types_t)){
  treatment_type <- treatment[treatment$treatmentType %in% types_t[i],]
  
  HW <- all_del_info[names(all_del_info) %in% treatment_type$sampleId]
  HW <- do.call(rbind,HW)
  HW$Deletion_size <- HW$Deletion_size * -1
  
  both <- cbind.fill(HW$Deletion_size,control$Deletion_size)
  
  colnames(both) <- c(paste("HW",types_t[i]),"NULL control")
  meltboth <- melt(log2(both))
  
### individual treatment plots
  p[[i]] <- ggplot(data = meltboth, aes(x=value,fill=Var2)) +
    geom_density(alpha = 0.5) +
    labs(x="Log2 deletion size") +
    scale_fill_manual(values=c("#00BA38", "#F8766D", "#00BFC4")) +
    theme_minimal() +
    geom_vline(xintercept = c(11),colour=c("red"), linetype = c("dashed"))
  }
types_t[9] <-  "Androgen_estrogen deprivation therapy"
for (i in 1:length(p)) {
  file_name = paste("Density plot treatment type ",types_t[i], ".jpeg", sep="")
  jpeg(file_name,width = 1200, height = 800, units = "px")
  print(p[[i]])
  dev.off()
}

### Combine chemotherapy plot with experimental and target + chemo with immuno and hormonal
treatment_type <- treatment[treatment$treatmentType %in% types_t[1],] ### change the types_t for different treatment

HW <- all_del_info[names(all_del_info) %in% treatment_type$sampleId]
HW <- do.call(rbind,HW)
HW$Deletion_size <- HW$Deletion_size * -1

# HW_chemo <- HW
# HW_immun <- HW
# HW_target <- HW
# HW_exp <- HW
# HW_horm <- HW


both <- cbind.fill(HW_chemo$Deletion_size,HW_immun$Deletion_size,HW_horm$Deletion_size )

colnames(both) <- c("Chemotherapy","Immunotherapy","Hormonal therapy")

both <- cbind.fill(HW_chemo$Deletion_size,HW_target$Deletion_size,HW_exp$Deletion_size )

colnames(both) <- c("Chemotherapy","Targeted therapy","Experimental therapy")
meltboth <- melt(log2(both))

ggplot(data = meltboth, aes(x=value,fill=Var2)) +
  geom_density(alpha = 0.5) +
  labs(x="Log2 deletion size") +
  scale_fill_manual(values=c("#F8766D", "#00BA38", "#619CFF")) +
  theme_minimal() +
  geom_vline(xintercept = c(10.5),colour=c("red"), linetype = c("dashed"))



##################### Supplementary figure 8 ##################### 
### continue from Figure 5C/D
total_type <- tumortype$primaryTumorLocation[tumortype$sampleId %in% total_del_amp$patient]
total_del_amp  <- total_del_amp[order(total_del_amp$patient),]
total_del_amp$type <- rep(total_type,table(total_del_amp$patient))
### Generate the plots
genes_final_selected <- unique(total_del_amp$Gene[total_del_amp$CNV == status])
for(i in 1:length(genes_final_selected)){
  gene <- genes_final_selected[i]
  ATM <- unique(total_del_amp$patient[total_del_amp$Gene == gene & total_del_amp$CNV == status])
  del_info_ATM <- NULL
  for(ai in 1:length(alltype)){
    type <- alltype[ai]
    temp <- get(paste("all_del_info_",type,sep=""))
    tissueID <- tumortype$sampleId[tumortype$primaryTumorLocation %in% type]
    tissuegrid <- unique(grep(paste(tissueID,collapse="|"), allgrids, value=TRUE))
    names(temp) <- tissuegrid
    pattern <- paste(ATM, collapse = "|") 
    temp <- temp[grepl(pattern, names(temp))]
    del_info_ATM <- append(del_info_ATM,temp)
  }
  # if(length(del_info_ATM) == 1){next}
  size <- hom_ATM <- list()
  for(ii in 1:length(del_info_ATM)){
    size[[ii]] <- del_info_ATM[[ii]]$Deletion_size
    temp1 <- del_info_ATM[[ii]]
    temp1$Deletion_size <- temp1$Deletion_size * -1
    temp_hom <- temp1$Homology_length[temp1$Deletion_size >= 724 & temp1$Deletion_size <= 3444]
    hom_ATM[[ii]] <- temp_hom
  }
  ATM_size <- do.call(cbind.fill,size)
  ATM_size_comb <- do.call(c,size)
  ATM_size_comb <- as.data.frame(ATM_size_comb)
  ATM_size_comb <- ATM_size_comb * -1
  
  ### For making taking the random samples in the same ratio (10x more)
  random_types <- table(tumortype$primaryTumorLocation[tumortype$sampleId %in% ATM])
  random_types <- as.data.frame(random_types*10)
  all_random_del <- NULL
  for(ii in 1:nrow(random_types)){
    type_temp <- random_types$Var1[ii]
    freq_temp <- random_types$Freq[ii]
    temp_del <- get(paste("all_del_info_",type_temp,sep=""))
    tissueID <- tumortype$sampleId[tumortype$primaryTumorLocation %in% type_temp]
    tissuegrid <- unique(grep(paste(tissueID,collapse="|"), allgrids, value=TRUE))
    names(temp_del) <- tissuegrid
    random_del <- temp_del[!names(temp_del) %in% names(del_info_ATM)]
    if(freq_temp > length(random_del)){
      temp_del_random <- random_del
      # print(paste("For",gene,freq_temp,type_temp,"should be picked from random but only",length(random_del),"avaiable",sep=" "))
    }else{
      temp_del_random <- sample(random_del,freq_temp)
    }
    all_random_del <- append(all_random_del,temp_del_random)
  }
  random_size <- list()
  for(iii in 1:length(all_random_del)){
    random_size[[iii]] <- all_random_del[[iii]]$Deletion_size
  }
  allrandom_size_comb <- do.call(c,random_size)
  allrandom_size_comb <- as.data.frame(allrandom_size_comb)
  allrandom_size_comb <- allrandom_size_comb * -1
  
  
  both <- cbind.fill(ATM_size_comb,allrandom_size_comb)
  colnames(both) <- c(paste(gene,status),paste(gene,"no",status))
  meltsize <- melt(log2(both))
  colnames(meltsize)[2] <- c("CN status")
  plot <- ggplot(data = meltsize, aes(x=value,fill=`CN status`)) +
    geom_density(alpha = 0.5) +
    labs(x="Log2 deletion size") +
    theme_classic() +
    theme_minimal() +
    geom_vline(xintercept = c(5,6,7,8.5,9.5,11.75,10.5),colour=c("blue","blue","black","black","green","green","red"), linetype = c("solid","solid","solid","solid","solid","solid","dashed"))
  ggsave(plot, file=paste("AMP over 20/Density deletion size CNV ",gene," ",status," all tumors vs random all tumors HW2",".png",sep=""), width = 10, height = 6.1, units = "in")
}

##################### Supplementary figure 7 ##################### 


ref_genome <- "BSgenome.Mmusculus.UCSC.mm10"
library(ref_genome, character.only = TRUE)
library(MutationalPatterns)
library(gridExtra)

### VCF needs chr in contig
listvcf <- list.files("/path_to_SNV_VCF/")
sample_names <-  sub("\\_.*", "", listvcf)
listvcf <- paste("/path_to_SNV_VCF/",list.files("/path_to_SNV_VCF/"),sep = "")
vcfs <- read_vcfs_as_granges(listvcf, sample_names, ref_genome)

seqlevelsStyle(vcfs_all) = 'UCSC'

muts = mutations_from_vcf(vcfs[[1]])
types = mut_type(vcfs[[1]])
context = mut_context(vcfs[[1]], ref_genome)
type_context = type_context(vcfs[[1]], ref_genome)

### Base substitutions (point mutations)
type_occurrences <- mut_type_occurrences(vcfs_all, ref_genome)
plot_spectrum(type_occurrences, CT = TRUE)
p1 <- plot_spectrum(type_occurrences, by = sample_names, CT = TRUE, legend = TRUE)

grid.arrange(p1, ncol = 3, widths = c(3, 3, 1.75))

### 96 mutational signature 
mut_mat <- mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)
plot_96_profile(mut_mat,condensed = TRUE)
# comparison plot
plot_compare_profiles(mut_mat[,1],mut_mat[,2],profile_names = c("WT", "DM"),condensed = TRUE)
