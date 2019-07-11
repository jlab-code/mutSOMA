makeVCFpedigreeTEMP<-function(genome.size, input.dir)
{
  vcf<-read.vcfR(paste(input.dir, "tree14_step3.vcf.gz", sep=""), verbose=TRUE)
  fix<- vcf@fix
  gt <- vcf@gt[,-1]
  #gt[,1]<-sub("\t", "", gt[,1])
  
  for (a in 1:ncol(gt))
  {
    gt[,a]<-vapply(strsplit(gt[,a], ":"), '[', 1, FUN.VALUE = character(1))
  }
  
  vcf<-cbind(fix[,c(1:5)], gt)
  
  
  ### Keeping only clear SNP substitutions, with a single alternative allele
  vcf<-vcf[which(vcf[,"REF"] %in% c("A", "C", "G", "T")),]
  vcf<-vcf[which(vcf[,"ALT"] %in% c("A", "C", "G", "T")),]
  
  ### Removing possible non-unique IDs
  vcf14<-vcf[which(duplicated(vcf[,"ID"]) == FALSE),]
  
  vcf<-read.vcfR(paste(input.dir, "tree13_step3.vcf.gz", sep=""), verbose=TRUE)
  fix<- vcf@fix
  gt <- vcf@gt[,-1]
  #gt[,1]<-sub("\t", "", gt[,1])
  
  for (a in 1:ncol(gt))
  {
    gt[,a]<-vapply(strsplit(gt[,a], ":"), '[', 1, FUN.VALUE = character(1))
  }
  
  vcf<-cbind(fix[,c(1:5)], gt)
  
  
  ### Keeping only clear SNP substitutions, with a single alternative allele
  vcf<-vcf[which(vcf[,"REF"] %in% c("A", "C", "G", "T")),]
  vcf<-vcf[which(vcf[,"ALT"] %in% c("A", "C", "G", "T")),]
  
  ### Removing possible non-unique IDs
  vcf13<-vcf[which(duplicated(vcf[,"ID"]) == FALSE),]
  
  ### Combining trees
  inter_id<-intersect(vcf14[,"ID"], vcf13[,"ID"])
  inter_vcf13<-vcf13[which(is.element(vcf13[,"ID"], inter_id) == TRUE),]
  no_inter_vcf13<-vcf13[which(is.element(vcf13[,"ID"], inter_id) == FALSE),]
  inter_vcf14<-vcf14[which(is.element(vcf14[,"ID"], inter_id) == TRUE),]
  no_inter_vcf14<-vcf14[which(is.element(vcf14[,"ID"], inter_id) == FALSE),]
  vcf1<-cbind(inter_vcf13, inter_vcf14[,6:ncol(inter_vcf14)])
  
  no_vcf14_fill<-matrix("0/0", ncol=length(6:ncol(inter_vcf13)), nrow=nrow(no_inter_vcf14))
  colnames(no_vcf14_fill)<-colnames(inter_vcf13)[6:ncol(inter_vcf13)]
  vcf2<-cbind(no_inter_vcf14[,1:5], no_vcf14_fill, no_inter_vcf14[, 6:ncol(no_inter_vcf14)])
  
  no_vcf13_fill<-matrix("0/0", ncol=length(6:ncol(inter_vcf14)), nrow=nrow(no_inter_vcf13))
  colnames(no_vcf13_fill)<-colnames(inter_vcf14)[6:ncol(inter_vcf13)]
  vcf3<-cbind(no_inter_vcf13, no_vcf13_fill)
  
  ### Final dataset including SNPs locations unique to tree 13 and tree 14 (aquired after split or unique false positives/false negatives ?) 
  ### AND SNP locations shared by both (aquired before split or recurrent or shared false positives?)
  vcf<-rbind(vcf1, vcf2, vcf3)
  
  # Let's only use high confidence ones
  #vcf<-vcf1
  
  gt<-vcf[,6:ncol(vcf)]
  out3<-apply(gt, 1, unique)
  out3<-sapply(out3, length)
  index<-which(out3 > 1)
  vcf<-vcf[index,]
  vcf<-vcf[,which(is.element(colnames(vcf), c("tree13.4", "tree14.1")) == FALSE)]
  
  for (a in 1:length(5:ncol(vcf)))
  {
    vcf<-vcf[!grepl("[.]/[.]", vcf[,a]),]
    
  }
  
  #### test
  calls.new<-matrix(NA, ncol=c(ncol(vcf)-5), nrow=nrow(vcf))
  
      for (i in 1:c(ncol(vcf)-5))
      {
        
          for (j in 1:nrow(vcf))
          {
            
            if (vcf[j, i+5] == "0/0"){calls.new[j, i]<-paste(vcf[j,"REF"], vcf[j, "REF"], sep="")}
            if (vcf[j, i+5] == "0/1"){calls.new[j, i]<-paste(vcf[j,"REF"], vcf[j, "ALT"], sep="")}
            if (vcf[j, i+5] == "1/1"){calls.new[j, i]<-paste(vcf[j,"ALT"], vcf[j, "ALT"], sep="")}
            if (is.element(vcf[j, i+5], c("0/0", "0/1", "1/1")) == FALSE){calls.new[j, i]<-paste(vcf[j,"REF"], vcf[j, "REF"], sep="")}
          }
      }
  
  
  colnames(calls.new)<-colnames(vcf[,6:ncol(vcf)])
  
  
  pairs<-colnames(vcf[,6:ncol(vcf)])
  #pairs<-combn(pairs, m=2)
  pairs<-expand.grid(pairs, pairs)
  pairs<-paste(pairs[,1], pairs[,2], sep="_")
  pairs<-matrix(pairs, ncol=8, nrow=8, byrow=F)
  pairs<-pairs[which(upper.tri(pairs) == TRUE)]
  pairs<-strsplit(pairs, "_")
  
  # c("AA", "CC", "TT", "GG", "AC", "AT", "AG", "CA", "CT", "CG", "TA", "TC", "TG", "GA", "GC", "GT")
  ## Coding the divergence between samples i and j
  Deffects<-matrix(as.numeric(c("0", "1", "1", "1", "0.5", "0.5", "0.5", "0.5", "1", "1", "0.5", "1", "1", "0.5", "1", "1", "1", "0", "1", "1", "0.5", "1", "1", 
                                "0.5", "0.5", "0.5", "1", "0.5", "1", "1", "0.5", "1", "1", "1", "0", "1", "1", "0.5", "1", "1", "0.5", "1", "0.5", "0.5", "0.5", 
                                "1", "1", "0.5", "1", "1", "1", "0", "1", "1", "0.5", "1", "1", "0.5", "1", "1", "0.5", "0.5", "0.5", "0.5", "0.5", "0.5", "1", 
                                "1", "0", "0.5", "0.5", "1", "1", "1", "1", "0.5", "1", "1", "0.5", "1", "0.5", "1", "0.5", "1", "0.5", "0", "0.5", "1", "0.5", 
                                "1", "1", "1", "1", "1", "1", "0.5", "0.5", "1", "1", "0.5", "0.5", "0.5", "0", "1", "1", "0.5", "1", "1", "0.5", "1", "1", 
                                "1", "0.5", "0.5", "1", "1", "1", "1", "1", "0", "0.5", "0.5", "0.5", "1", "1", "0.5", "1", "1", "1", "0.5", "0.5", "1", "1", 
                                "0.5", "1", "0.5", "0", "0.5", "1", "1", "1", "1", "1", "0.5", "1", "0.5", "1", "0.5", "1", "1", "0.5", "0.5", "0.5", "0", "1", 
                                "1", "0.5", "1", "1", "1", "0.5", "1", "0.5", "1", "1", "1", "1", "0.5", "1", "1", "0", "0.5", "0.5", "0.5", "1", "1", "1", 
                                "0.5", "0.5", "1", "0.5", "1", "1", "1", "1", "1", "0.5", "0", "0.5", "1", "0.5", "1", "1", "1", "0.5", "0.5", "1", "1", "0.5", 
                                "1", "1", "0.5", "0.5", "0.5", "0", "1", "1", "1", "0.5", "1", "1", "0.5", "1", "1", "1", "0.5", "1", "1", "0.5", "1", "1", 
                                "0", "0.5", "0.5", "1", "0.5", "1", "0.5", "0.5", "1", "1", "1", "1", "1", "1", "0.5", "1", "0.5", "0", "0.5", "1", "1", "0.5", 
                                "0.5", "1", "0.5", "1", "1", "0.5", "1", "1", "1", "1", "0.5", "0.5", "0")), nrow=16, ncol=16, byrow=TRUE)
  
  colnames(Deffects)<-c("AA", "CC", "TT", "GG", "AC", "AT", "AG", "CA", "CT", "CG", "TA", "TC", "TG", "GA", "GC", "GT")
  rownames(Deffects)<-c("AA", "CC", "TT", "GG", "AC", "AT", "AG", "CA", "CT", "CG", "TA", "TC", "TG", "GA", "GC", "GT")
  
  div<-NULL
  pair1<-NULL
  pair2<-NULL
  
  for (a in 1:length(pairs))
  {
    
    vec1<-calls.new[,which(colnames(calls.new) == pairs[[a]][1])]
    vec2<-calls.new[,which(colnames(calls.new) == pairs[[a]][2])]
    pair1[a]<-pairs[[a]][1]
    pair2[a]<-pairs[[a]][2]
    
    div.temp<-NULL
    
      for (b in 1:length(vec1))
      {
        div.temp[b]<-Deffects[which(rownames(Deffects) == vec1[b]), which(colnames(Deffects) == vec2[b])]
      }
    
    #zero<-length(which(vec1 == "0/0" & vec2 == "0/0" | vec1 == "0/1" & vec2 == "0/1" | vec1 == "1/1" & vec2 == "1/1"))
    #half<-length(which(vec1 == "0/1" & vec2 == "0/0" | vec1 == "0/1" & vec2 == "1/1" | vec1 == "0/0" & vec2 == "0/1" | vec1 == "1/1" & vec2 == "0/1"))
    #one<-length(which(vec1 == "0/0" & vec2 == "1/1" | vec1 == "1/1" & vec2 == "0/0"))
    
    #div[a]<-(zero* 0 + half * 1/2 + one)/genome.size
    div[a]<-sum(div.temp)/genome.size
    
  }
  
  div<-data.frame(pair1, pair2, div)
  
  ## Get Dmatrix
  pedigree<-read.table(paste(input.dir, "D-matrix_SingleCGfiltered.txt", sep=""), header=T)
  pedigree<-pedigree[,c(4,5,6)]
  sample.info<-read.csv(paste(input.dir, "sample_info.csv", sep=""), header=T)
  
  pedigree.out<-makePHYLO(tall=330, pedigree = pedigree, sample.info = sample.info)
  pedigree<-pedigree.out[[2]]
  
  
  ## Getting the data
  #pedigree<-read.table(paste(input.data.dir, "D-matrix_SingleCGfiltered.txt", sep=""), header=T)
  pedigree[,1]<-sub("_", ".", pedigree[,1])
  pedigree[,1]<-paste("tree", pedigree[,1], sep="")
  pedigree[,2]<-sub("_", ".", pedigree[,2])
  pedigree[,2]<-paste("tree", pedigree[,2], sep="")
  
  pedigree.new<-pedigree[,1:2]
  
  div.collect<-NULL
  for (a in 1:nrow(pedigree.new))
  {
    div.collect[a]<-div[which(paste(div[,1], div[,2]) == paste(pedigree.new[a,1], pedigree.new[a,2])),][,3]
  }
  
  pedigree.new<-cbind(pedigree.new, div.collect)
  
  
  pedigree[,3]<-div.collect
  pedigree2<-pedigree[,c(6,4,5,3)]
  
  #dthis<-as.numeric(pedigree2[,2]) + as.numeric(pedigree2[,3]) - as.numeric(2*pedigree2[,1])
  #plot(dthis, pedigree2[,4])
  
  output<-list(pedigree2, pedigree)
  names(output)<-c("pedigree", "pedigree.augmented")
  
  output
  
}

