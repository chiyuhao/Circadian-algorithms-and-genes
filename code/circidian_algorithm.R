#############################################meta2d
library(MetaCycle)
##ADM
meta2d(infile="D:/cyclingdata/Papio anubis (baboon)/64/adm/ADM2.csv", filestyle="csv", outdir="D:/cyclingdata/Papio anubis (baboon)/64/adm",timepoints="Line1", outIntegration="onlyIntegration")
##LUN
meta2d(infile="D:/cyclingdata/Papio anubis (baboon)/64/lun/LUN2.csv", filestyle="csv", outdir="D:/cyclingdata/Papio anubis (baboon)/64/lun",timepoints="Line1", outIntegration="onlyIntegration")
##LIV
meta2d(infile="D:/cyclingdata/Papio anubis (baboon)/64/LIV/LIV2.csv", filestyle="csv", outdir="D:/cyclingdata/Papio anubis (baboon)/64/LIV",timepoints="Line1", outIntegration="onlyIntegration")
##WAT
meta2d(infile="D:/cyclingdata/Papio anubis (baboon)/64/wat/WAT2.csv", filestyle="csv", outdir="D:/cyclingdata/Papio anubis (baboon)/64/wat",timepoints="Line1", outIntegration="onlyIntegration")
##KIC
meta2d(infile="D:/cyclingdata/Papio anubis (baboon)/64/KIC/KIC2.csv", filestyle="csv", outdir="D:/cyclingdata/Papio anubis (baboon)/64/KIC",timepoints="Line1", outIntegration="onlyIntegration")
##KIM
meta2d(infile="D:/cyclingdata/Papio anubis (baboon)/64/KIM/KIM2.csv", filestyle="csv", outdir="D:/cyclingdata/Papio anubis (baboon)/64/KIM",timepoints="Line1", outIntegration="onlyIntegration")
##CER
meta2d(infile="D:/cyclingdata/Papio anubis (baboon)/64/CER/CER2.csv", filestyle="csv", outdir="D:/cyclingdata/Papio anubis (baboon)/64/CER",timepoints="Line1", outIntegration="onlyIntegration")
##PON
meta2d(infile="D:/cyclingdata/Papio anubis (baboon)/64/pon/PON2.csv", filestyle="csv", outdir="D:/cyclingdata/Papio anubis (baboon)/64/pon",timepoints="Line1", outIntegration="onlyIntegration")
##AOR
meta2d(infile="D:/cyclingdata/Papio anubis (baboon)/64/AOR/AOR2.csv", filestyle="csv", outdir="D:/cyclingdata/Papio anubis (baboon)/64/AOR",timepoints="Line1", outIntegration="onlyIntegration")
##HEA
meta2d(infile="D:/cyclingdata/Papio anubis (baboon)/64/HEA/HEA2.csv", filestyle="csv", outdir="D:/cyclingdata/Papio anubis (baboon)/64/HEA",timepoints="Line1", outIntegration="onlyIntegration")
##MUA
meta2d(infile="D:/cyclingdata/Papio anubis (baboon)/64/MUA/MUA2.csv", filestyle="csv", outdir="D:/cyclingdata/Papio anubis (baboon)/64/MUA",timepoints="Line1", outIntegration="onlyIntegration")
##MUG
meta2d(infile="D:/cyclingdata/Papio anubis (baboon)/64/MUG/MUG2.csv", filestyle="csv", outdir="D:/cyclingdata/Papio anubis (baboon)/64/MUG",timepoints="Line1", outIntegration="onlyIntegration")






##############################################meta2darser
library(MetaCycle)
##ADR
meta2d(infile="D:/cyclingdata/mouse/gse54650/adr/adr.csv", filestyle="csv", outdir="D:/cyclingdata/mouse/gse54650/adr",timepoints="Line1",outRawData=TRUE)
##aor
meta2d(infile="D:/cyclingdata/mouse/gse54650/AOR/AOR.csv", filestyle="csv", outdir="D:/cyclingdata/mouse/gse54650/AOR",timepoints="Line1",outRawData=TRUE)
##BFAT
meta2d(infile="D:/cyclingdata/mouse/gse54650/BFAT/BFAT.csv", filestyle="csv", outdir="D:/cyclingdata/mouse/gse54650/BFAT",timepoints="Line1",outRawData=TRUE)
##BST
meta2d(infile="D:/cyclingdata/mouse/gse54650/BST/bst.csv", filestyle="csv", outdir="D:/cyclingdata/mouse/gse54650/BST",timepoints="Line1",outRawData=TRUE)
##CER
meta2d(infile="D:/cyclingdata/mouse/gse54650/CER/CER.csv", filestyle="csv", outdir="D:/cyclingdata/mouse/gse54650/CER",timepoints="Line1",outRawData=TRUE)
##HEA
meta2d(infile="D:/cyclingdata/mouse/gse54650/HEA/HEA.csv", filestyle="csv", outdir="D:/cyclingdata/mouse/gse54650/HEA",timepoints="Line1",outRawData=TRUE)
##HYP
meta2d(infile="D:/cyclingdata/mouse/gse54650/HYP/HYP.csv", filestyle="csv", outdir="D:/cyclingdata/mouse/gse54650/HYP",timepoints="Line1",outRawData=TRUE)
##KID
meta2d(infile="D:/cyclingdata/mouse/gse54650/KID/KID.csv", filestyle="csv", outdir="D:/cyclingdata/mouse/gse54650/KID",timepoints="Line1",outRawData=TRUE)
##LIV
meta2d(infile="D:/cyclingdata/mouse/gse54650/LIV/LIV.csv", filestyle="csv", outdir="D:/cyclingdata/mouse/gse54650/LIV",timepoints="Line1",outRawData=TRUE)
##LUN
meta2d(infile="D:/cyclingdata/mouse/gse54650/LUN/LUN.csv", filestyle="csv", outdir="D:/cyclingdata/mouse/gse54650/LUN",timepoints="Line1",outRawData=TRUE)
##MUS
meta2d(infile="D:/cyclingdata/mouse/gse54650/MUS/MUS.csv", filestyle="csv", outdir="D:/cyclingdata/mouse/gse54650/MUS",timepoints="Line1",outRawData=TRUE)
##WAT
meta2d(infile="D:/cyclingdata/mouse/gse54650/WAT/WAT.csv", filestyle="csv", outdir="D:/cyclingdata/mouse/gse54650/WAT",timepoints="Line1",outRawData=TRUE)
##
meta2d(infile="D:/cyclingdata/mouse/gse54650/adr/adr.csv", filestyle="csv", outdir="D:/cyclingdata/mouse/gse54650/adr",timepoints="Line1",outRawData=TRUE)


########################################################COSOPT
install.packages("E:/circidian_algorithm/cosopt-master/cosopt_0.3.1.tar.gz",repos = NULL, type = "source")
library(cosopt)
tissuename <- c("adr","AOR","BFAT","bst","CER","HEA","HYP","KID","LIV","LUN","MUS","WAT")
for(i in 1:length(tissuename))
{
  file_name <- paste("E:/circidian_algorithm/data/mouse/data/",tissuename[i],".csv",sep="")
  data <- read.csv(file_name)
  timepoints <- c(18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,62,64)
  sigma <- apply(data[,2:25],2,sd)
  data1 <- as.matrix(data[,2:25])
  cic_list <- c()
  for(j in 1:length(data[,1]))
  {
    temp_data <- as.vector(data1[j,])
    freq <- cosopt(temp_data,as.vector(sigma),as.vector(timepoints),plotting=FALSE)$freq
    period <- 1/freq
    if(period < 25 && period >23)
    {
      cic_list[(length(cic_list)+1)] <- as.character(data[j,1])
    }
  }
  file_out_name <- paste("E:/circidian_algorithm/result/mouse/cosopt/",tissuename[i],".csv",sep="")
  write.csv(cir_lsit, file_out_name)
}

#######################################################lomb
library(lomb)
tissuename <- c("adr","AOR","BFAT","bst","CER","HEA","HYP","KID","LIV","LUN","MUS","WAT")
for(j in 1:length(tissuename))
{
  file_name <- paste("E:/circidian_algorithm/data/mouse/data/",tissuename[j],".csv",sep="")
  data <- read.csv(file_name)
  data1 <- as.matrix(data[,2:25])
  pvalue_list <- c()
  period_list <- c()
  for(i in 1:length(data[,1]))
  {
    temp <- as.vector(data1[i,])
    lynx.spec <- lsp(temp,type='period',from=18,to=64,ofac=5,plot=FALSE)
    pvalue_list[i] <- lynx.spec$p.value
    period_list[i] <- lynx.spec$peak.at[1]
  }
  pval <- pvalue_list
  bhq <- p.adjust(pval,method="fdr")
  temp <- data.frame(data[,1], pvalue_list, bhq,period_list)
  file_out_name <- paste("E:/circidian_algorithm/result/mouse/lomb/",tissuename[j],".csv",sep="")
  print(tissuename[j])
  print(sum(pvalue_list < 0.05))
  print(sum(bhq < 0.05))
  write.csv(temp, file_out_name)

}

#fly
data <- read.csv("E:/circidian_algorithm/data/fly/JTK.flyold.csv")
data1 <- as.matrix(data[,2:13])
pvalue_list <- c()
period_list <- c()
for(i in 1:length(data[,1]))
{
  temp <- as.vector(data1[i,])
  if(sum(temp) == 0)
  {
    pvalue_list[i] <- NA
    period_list[i] <- NA
  }else
  {
    lynx.spec <- lsp(temp,type='period',from=4,to=48,ofac=5,plot=FALSE)
    pvalue_list[i] <- lynx.spec$p.value
    period_list[i] <- lynx.spec$peak.at[1]
  }
}
temp <- data.frame(pvalue_list, period_list)
write.csv(temp, "E:/circidian_algorithm/result/fly/lomb.csv")
result <- read.csv("E:/circidian_algorithm/result/fly/lomb.csv")
pval <- na.omit(result$pvalue_list)
bhq <- p.adjust(pval,method="fdr")
sum(bhq<0.05)
sum(pval<0.05)

#archaea
data <- read.csv("E:/circidian_algorithm/data/archaea/meta2d_AR2_achaea.csv")
data1 <- as.matrix(data[,2:37])
pvalue_list <- c()
period_list <- c()
for(i in 1:length(data[,1]))
{
  temp <- as.vector(data1[i,])
  if(sum(temp) == 0)
  {
    pvalue_list[i] <- NA
    period_list[i] <- NA
  }else
  {
    lynx.spec <- lsp(temp,type='period',from=4,to=72,ofac=5,plot=FALSE)
    pvalue_list[i] <- lynx.spec$p.value
    period_list[i] <- lynx.spec$peak.at[1]
  }
}
temp <- data.frame(pvalue_list, period_list)
write.csv(temp, "E:/circidian_algorithm/result/archaea/lomb.csv")
result <- read.csv("E:/circidian_algorithm/result/archaea/lomb.csv")
pval <- na.omit(result$pvalue_list)
bhq <- p.adjust(pval,method="fdr")
sum(bhq<0.05)
sum(pval<0.05)




######################################################RAIN
library(rain)

set.seed(123)

tissuename <- c("adr","AOR","BFAT","bst","CER","HEA","HYP","KID","LIV","LUN","MUS","WAT")
for(j in 1:length(tissuename))
{
  file_name <- paste("E:/circidian_algorithm/data/mouse/data/",tissuename[j],".csv",sep="")
  data <- read.csv(file_name)
  input_data <- t(data[,2:25])
  rainresult <- rain(input_data, period=24, deltat=2, peak.border=c(0.3,0.7),verbose=FALSE) 
  #如果有重复nr.series=n
  #results <- rain(t(menetRNASeqMouseLiver), deltat=4, period=24, nr.series=2,peak.border=c(0.3, 0.7), verbose=FALSE)
  pval <- rainresult$pVal
  bhq <- p.adjust(pval, method="fdr")
  temp <- data.frame(data[,1], pval, bhq, rainresult$phase)
  file_out_name <- paste("E:/circidian_algorithm/result/mouse/rain/",tissuename[j],".csv",sep="")
  print(tissuename[j])
  print(sum(pval < 0.05))
  print(sum(bhq < 0.05))
  write.csv(temp, file_out_name)
}

#fly
data <- read.csv("E:/circidian_algorithm/data/fly/JTK.flyold.csv")
input_data <- t(data[,2:13])
rainresult <- rain(input_data, period=24, deltat=4, peak.border=c(0.3,0.7),verbose=FALSE) 
#如果有重复nr.series=n
#results <- rain(t(menetRNASeqMouseLiver), deltat=4, period=24, nr.series=2,peak.border=c(0.3, 0.7), verbose=FALSE) 
write.csv(rainresult, "E:/circidian_algorithm/result/fly/rain.csv")
result <- read.csv("E:/circidian_algorithm/result/fly/rain.csv")
pval <- result$pVal
sum(pval<0.05)
bhq <- p.adjust(pval, method="fdr")
sum(bhq<0.05)

#archaea
data <- read.csv("E:/circidian_algorithm/data/archaea/meta2d_AR2_achaea.csv")
input_data <- t(data[,2:37])
rainresult <- rain(input_data, period=24, deltat=4, nr.series=2, peak.border=c(0.3,0.7),verbose=FALSE) 
#如果有重复nr.series=n
#results <- rain(t(menetRNASeqMouseLiver), deltat=4, period=24, nr.series=2,peak.border=c(0.3, 0.7), verbose=FALSE) 
write.csv(rainresult, "E:/circidian_algorithm/result/archaea/rain.csv")
result <- read.csv("E:/circidian_algorithm/result/archaea/rain.csv")
pval <- result$pVal
sum(pval<0.05)
bhq <- p.adjust(pval, method="fdr")
sum(bhq<0.05)


######################################################ABSR
library()

#######################################################FIsher G test
library(GeneCycle)


tissuename <- c("adr","AOR","BFAT","bst","CER","HEA","HYP","KID","LIV","LUN","MUS","WAT")
for(j in 1:length(tissuename))
{
  file_name <- paste("E:/circidian_algorithm/data/mouse/data/",tissuename[j],".csv",sep="")
  data <- read.csv(file_name)
  input_data <- t(data[,2:25])
  # p-values from Fisher's g test
  pval.caulobacter <- fisher.g.test(input_data)
  
  fdr.out <- fdrtool(pval.caulobacter, statistic="pvalue")
  # how many significant?
  temp <- data.frame(data[,1], fdr.out$pval, fdr.out$qval)
  print(tissuename[j])
  print(sum(fdr.out$pval < 0.05))
  print(sum(fdr.out$qval < 0.05))
  file_out_name <- paste("E:/circidian_algorithm/result/mouse/gtest/",tissuename[j],".csv",sep="")
  write.csv(temp, file_out_name)
  # sum(fdr.out$qval < 0.05) # tail area-based Fdr
  # sum(fdr.out$lfdr < 0.2) # density-based local fdr
  # sum(fdr.out$lfdr < 0.05) # density-based local fdr
}
#fly
data <- read.csv("E:/circidian_algorithm/data/fly/JTK.flyold.csv")
input_data <- t(data[,2:13])
# p-values from Fisher's g test
pval.caulobacter <- fisher.g.test(input_data)

fdr.out <- fdrtool(pval.caulobacter, statistic="pvalue")
# how many significant?
sum(fdr.out$qval < 0.05) # tail area-based Fdr
sum(fdr.out$lfdr < 0.2) # density-based local fdr
sum(fdr.out$lfdr < 0.05) # density-based local fdr

#archaea
data <- read.csv("E:/circidian_algorithm/data/archaea/meta2d_AR2_achaea.csv")
input_data <- t(data[,2:37])
# p-values from Fisher's g test
pval.caulobacter <- fisher.g.test(input_data)

fdr.out <- fdrtool(pval.caulobacter, statistic="pvalue")
# how many significant?
sum(fdr.out$qval < 0.05) # tail area-based Fdr
sum(fdr.out$lfdr < 0.2) # density-based local fdr
sum(fdr.out$lfdr < 0.05) # density-based local fdr






tissuename <- c("adr","AOR","BFAT","bst","CER","HEA","HYP","KID","LIV","LUN","MUS","WAT")
for(j in 1:length(tissuename))
{
  file_name <- paste("E:/circidian_algorithm/data/mouse/data/",tissuename[j],".csv",sep="")
  data <- read.csv(file_name,header=FALSE)
  new_data <- t(data)
  col_name <- data[,1]
  new_data <- new_data[-1,]
  colnames(new_data) <- col_name
  row_name <- data[1,2:25]
  rownames(new_data) <- row_name
  file_out_name <- paste("E:/circidian_algorithm/data/mouse/data/",tissuename[j],".txt",sep="")
  write.table(new_data[-1,], file_out_name, sep='/t',row.names = FALSE)
}




#################eJTK
tissuename <- c("adr","AOR","BFAT","bst","CER","HEA","HYP","KID","LIV","LUN","MUS","WAT")
for(j in 1:length(tissuename))
{
  file_name <- paste("E:/circidian_algorithm/result/mouse/eJTK/",tissuename[j],"_cos24_ph00-22_by2_a02-22_by2_OTHERTEXT_jtkout_GammaP.txt",sep="")
  data <- read.table(file_name, header = T)
  pval <- data$P
  bhq <- p.adjust(pval, method="fdr")
  print(tissuename[j])
  print(sum(pval < 0.05))
  print(sum(bhq < 0.05))
}























###########################################################4h
#######################################################lomb
library(lomb)
tissuename <- c("adr","AOR","BFAT","bst","CER","HEA","HYP","KID","LIV","LUN","MUS","WAT")
for(j in 1:length(tissuename))
{
  file_name <- paste("E:/circidian_algorithm/data/mouse/mouse_4h/",tissuename[j],".csv",sep="")
  data <- read.csv(file_name)
  data1 <- as.matrix(data[,2:13])
  pvalue_list <- c()
  period_list <- c()
  for(i in 1:length(data[,1]))
  {
    temp <- as.vector(data1[i,])
    lynx.spec <- lsp(temp,type='period',from=4,to=48,ofac=5,plot=FALSE)
    pvalue_list[i] <- lynx.spec$p.value
    period_list[i] <- lynx.spec$peak.at[1]
  }
  pval <- pvalue_list
  bhq <- p.adjust(pval,method="fdr")
  temp <- data.frame(data[,1], pvalue_list, bhq,period_list)
  file_out_name <- paste("E:/circidian_algorithm/result/mouse/4h/lomb/",tissuename[j],".csv",sep="")
  print(tissuename[j])
  print(sum(pvalue_list < 0.05))
  print(sum(bhq < 0.05))
  write.csv(temp, file_out_name)
  
}




######################################################RAIN
library(rain)

set.seed(123)

tissuename <- c("adr","AOR","BFAT","bst","CER","HEA","HYP","KID","LIV","LUN","MUS","WAT")
for(j in 1:length(tissuename))
{
  file_name <- paste("E:/circidian_algorithm/data/mouse/mouse_4h/",tissuename[j],".csv",sep="")
  data <- read.csv(file_name)
  input_data <- t(data[,2:13])
  rainresult <- rain(input_data, period=24, deltat=4, peak.border=c(0.3,0.7),verbose=FALSE) 
  #如果有重复nr.series=n
  #results <- rain(t(menetRNASeqMouseLiver), deltat=4, period=24, nr.series=2,peak.border=c(0.3, 0.7), verbose=FALSE)
  pval <- rainresult$pVal
  bhq <- p.adjust(pval, method="fdr")
  temp <- data.frame(data[,1], pval, bhq, rainresult$phase)
  file_out_name <- paste("E:/circidian_algorithm/result/mouse/4h/rain/",tissuename[j],".csv",sep="")
  print(tissuename[j])
  print(sum(pval < 0.05))
  print(sum(bhq < 0.05))
  write.csv(temp, file_out_name)
}




######################################################ABSR
library()

#######################################################FIsher G test
library(GeneCycle)


tissuename <- c("adr","AOR","BFAT","bst","CER","HEA","HYP","KID","LIV","LUN","MUS","WAT")
for(j in 1:length(tissuename))
{
  file_name <- paste("E:/circidian_algorithm/data/mouse/mouse_4h/",tissuename[j],".csv",sep="")
  data <- read.csv(file_name)
  input_data <- t(data[,2:13])
  # p-values from Fisher's g test
  pval.caulobacter <- fisher.g.test(input_data)
  
  fdr.out <- fdrtool(pval.caulobacter, statistic="pvalue")
  # how many significant?
  temp <- data.frame(data[,1], fdr.out$pval, fdr.out$qval)
  print(tissuename[j])
  print(sum(fdr.out$pval < 0.05))
  print(sum(fdr.out$qval < 0.05))
  file_out_name <- paste("E:/circidian_algorithm/result/mouse/4h/gtest/",tissuename[j],".csv",sep="")
  write.csv(temp, file_out_name)
  # sum(fdr.out$qval < 0.05) # tail area-based Fdr
  # sum(fdr.out$lfdr < 0.2) # density-based local fdr
  # sum(fdr.out$lfdr < 0.05) # density-based local fdr
}





#################eJTK
tissuename <- c("adr","AOR","BFAT","bst","CER","HEA","HYP","KID","LIV","LUN","MUS","WAT")
for(j in 1:length(tissuename))
{
  file_name <- paste("E:/circidian_algorithm/result/mouse/eJTK/",tissuename[j],"_cos24_ph00-22_by2_a02-22_by2_OTHERTEXT_jtkout_GammaP.txt",sep="")
  data <- read.table(file_name, header = T)
  pval <- data$P
  bhq <- p.adjust(pval, method="fdr")
  print(tissuename[j])
  print(sum(pval < 0.05))
  print(sum(bhq < 0.05))
}

############JTK
tissuename <- c("adr","AOR","BFAT","bst","CER","HEA","HYP","KID","LIV","LUN","MUS","WAT")
for(j in 1:length(tissuename))
{
  file_name <- paste("E:/circidian_algorithm/result/mouse/4h/JTK/JTKresult_",tissuename[j],".csv",sep="")
  data <- read.csv(file_name)
  print(tissuename[j])
  print(sum(data$ADJ.P < 0.05))
  print(sum(data$BH.Q < 0.05))
}

tissuename <- c("adr","AOR","BFAT","bst","CER","HEA","HYP","KID","LIV","LUN","MUS","WAT")
for(j in 1:length(tissuename))
{
  file_name <- paste("E:/circidian_algorithm/result/mouse/4h/JTK/ARSresult_",tissuename[j],".csv",sep="")
  data <- read.csv(file_name)
  print(tissuename[j])
  print(sum(data$pvalue < 0.05))
  print(sum(data$fdr_BH < 0.05))
}


tissuename <- c("adr","AOR","BFAT","bst","CER","HEA","HYP","KID","LIV","LUN","MUS","WAT")
for(j in 1:length(tissuename))
{
  file_name <- paste("E:/circidian_algorithm/result/mouse/4h/JTK/LSresult_",tissuename[j],".csv",sep="")
  data <- read.csv(file_name)
  print(tissuename[j])
  print(sum(data$p< 0.05))
  print(sum(data$BH.Q < 0.05))
}


tissuename <- c("adr","AOR","BFAT","bst","CER","HEA","HYP","KID","LIV","LUN","MUS","WAT")
for(j in 1:length(tissuename))
{
  file_name <- paste("E:/circidian_algorithm/result/mouse/4h/JTK/meta2d_",tissuename[j],".csv",sep="")
  data <- read.csv(file_name)
  print(tissuename[j])
  print(sum(data$meta2d_pvalue< 0.05))
  print(sum(data$meta2d_BH.Q < 0.05))
}



tissuename <- c("adr","AOR","BFAT","bst","CER","HEA","HYP","KID","LIV","LUN","MUS","WAT")
for(j in 1:length(tissuename))
{
  file_name <- paste("E:/circidian_algorithm/result/mouse/4h/eJTK/",tissuename[j],"_cos24_ph00-20_by4_a04-20_by4_OTHERTEXT_jtkout_GammaP.txt",sep="")
  data <- read.table(file_name, header = T)
  pval <- data$GammaP
  bhq <- p.adjust(pval, method="fdr")
  print(tissuename[j])
  print(sum(pval < 0.05))
  print(sum(bhq < 0.05))
  print(sum(data$GammaBH<0.05))
}

data <- read.csv("E:/circidian_algorithm/data/input.csv",row.names=1)
data <- t(data)
boxplot(data)




tissuename <- c("adr","AOR","BFAT","bst","CER","HEA","HYP","KID","LIV","LUN","MUS","WAT")
for(j in 1:length(tissuename))
{
  file_name <- paste("E:/circidian_algorithm/data/mouse/mouse_4h/",tissuename[j],".out.csv",sep="")
  data <- read.csv(file_name)
  pval <- data$p1
  bhq <- p.adjust(pval, method="fdr")
  temp <- data.frame(data$Code, pval, bhq)
  print(tissuename[j])
  print(sum(pval < 0.05))
  print(sum(bhq < 0.05))
}



file_name <- paste("E:/circidian_algorithm/data/mouse/mouse_4h/","LIV",".out.csv",sep="")
data <- read.csv(file_name)
pval <- data$p1
bhq <- p.adjust(pval, method="fdr")
temp <- data.frame(data$Code, pval, bhq)
output_circ <- temp[temp$bhq<0.05, ]$data.Code


#zeitzeiger
library(zeitzeiger)




###########################################################venn plot
#######################################################lomb
library(lomb)

file_name <- paste("E:/circidian_algorithm/data/mouse/mouse_4h/","LIV",".csv",sep="")
data <- read.csv(file_name)
data1 <- as.matrix(data[,2:13])
pvalue_list <- c()
period_list <- c()
for(i in 1:length(data[,1]))
{
  temp <- as.vector(data1[i,])
  lynx.spec <- lsp(temp,type='period',from=4,to=48,ofac=5,plot=FALSE)
  pvalue_list[i] <- lynx.spec$p.value
  period_list[i] <- lynx.spec$peak.at[1]
}
pval <- pvalue_list
bhq <- p.adjust(pval,method="fdr")
temp <- data.frame(data[,1], pvalue_list, bhq,period_list)
output <- temp[temp$bhq<0.05]
file_out_name <- paste("E:/circidian_algorithm/result/mouse/4h/lomb/",tissuename[j],".csv",sep="")
print(tissuename[j])
print(sum(pvalue_list < 0.05))
print(sum(bhq < 0.05))
write.csv(temp, file_out_name)





######################################################RAIN
library(rain)

set.seed(123)


file_name <- paste("E:/circidian_algorithm/data/mouse/mouse_4h/","LIV",".csv",sep="")
data <- read.csv(file_name)
input_data <- t(data[,2:13])
rainresult <- rain(input_data, period=24, deltat=4, peak.border=c(0.3,0.7),verbose=FALSE) 
#如果有重复nr.series=n
#results <- rain(t(menetRNASeqMouseLiver), deltat=4, period=24, nr.series=2,peak.border=c(0.3, 0.7), verbose=FALSE)
pval <- rainresult$pVal
bhq <- p.adjust(pval, method="fdr")
temp <- data.frame(data[,1], pval, bhq, rainresult$phase)
output_rain <- temp[temp$bhq<0.05,]$data...1.

#######################################################FIsher G test
library(GeneCycle)


file_name <- paste("E:/circidian_algorithm/data/mouse/mouse_4h/","LIV",".csv",sep="")
data <- read.csv(file_name)
input_data <- t(data[,2:13])
# p-values from Fisher's g test
pval.caulobacter <- fisher.g.test(input_data)

fdr.out <- fdrtool(pval.caulobacter, statistic="pvalue")
# how many significant?
temp <- data.frame(data[,1], fdr.out$pval, fdr.out$qval)
output_gtest <- temp[temp$fdr.out.qval<0.05,]$data...1.






#################eJTK

file_name <- paste("E:/circidian_algorithm/result/mouse/eJTK/","LIV","_cos24_ph00-22_by2_a02-22_by2_OTHERTEXT_jtkout_GammaP.txt",sep="")
data <- read.table(file_name, header = T)
pval <- data$P
bhq <- p.adjust(pval, method="fdr")
temp <- data.frame(data$ID,pval,bhq)
otput_ejtk <- temp[temp$bhq<0.05,]$data.ID

############JTK

file_name <- paste("E:/circidian_algorithm/result/mouse/4h/JTK/JTKresult_","LIV",".csv",sep="")
data <- read.csv(file_name)
temp <- data
output_JTK <- temp[temp$BH.Q<0.05,]$CycID




file_name <- paste("E:/circidian_algorithm/result/mouse/4h/JTK/ARSresult_","LIV",".csv",sep="")
data <- read.csv(file_name)
temp <- data
output_arser <- temp[temp$fdr_BH <0.05,]$CycID




file_name <- paste("E:/circidian_algorithm/result/mouse/4h/JTK/LSresult_","LIV",".csv",sep="")
data <- read.csv(file_name)
temp <- data
output_meta2dLS <- temp[temp$BH.Q <0.05,]$CycID




file_name <- paste("E:/circidian_algorithm/result/mouse/4h/JTK/meta2d_","LIV",".csv",sep="")
data <- read.csv(file_name)
temp <- data
fdr <- p.adjust(temp$meta2d_pvalue, method="fdr")
sum(fdr<0.05)
sum(temp$meta2d_BH.Q<0.05)
output_meta2d <- temp[temp$meta2d_BH.Q <0.05,]$CycID

write.csv(fdr, "E:/circidian_algorithm/result/mouse/4h/JTK/meta2dbhq.csv")




data <- read.csv("E:/circidian_algorithm/data/input.csv",row.names=1)
data <- t(data)
boxplot(data)



all_data <- list(output_rain, output_gtest, otput_ejtk, output_JTK, output_arser, output_meta2dLS, output_meta2d, output_circ)
library(VennDiagram)
venn.plot <- venn.diagram(list(rain=output_rain, jtk=output_JTK,arser=output_arser,circwave=output_circ,meta2d=output_meta2d), filename=NULL,lty = "dotted",
                          lwd = 2,
                          col = "black",  #"transparent",
                          fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
                          alpha = 0.60,
                          cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
                          cat.cex = 0.8,
                          cat.fontface = "bold",
                          margin = 0.07,
                          cex = 0.8)
grid.draw(venn.plot)



##########haystack

tissuename <- c("adr","AOR","BFAT","bst","CER","HEA","HYP","KID","LIV","LUN","MUS","WAT")
for(j in 1:length(tissuename))
{
  file_name <- paste("E:/circidian_algorithm/result/mouse/4h/haystack/",tissuename[j],".tsv",sep="")
  data <- read.table(file_name)
  pval <- data[,6]
  new_pval <- c(data[,6],rep(1,1000-length(pval)))
  fdr <- p.adjust(new_pval, method="fdr")
  print(sum(fdr<0.05))
}




























##########################################################################8 species
####################################meta2d
#mouse
##ADR
meta2d(infile="E:/circidian_algorithm/data/mouse/gse54650 4h/adr/adr.csv", filestyle="csv", outdir="E:/circidian_algorithm/data/mouse/4h_24_result/adr",timepoints="Line1",outRawData=TRUE,minper=20,maxper=28, cycMethod = c("ARS", "JTK", "LS"))
##aor
meta2d(infile="E:/circidian_algorithm/data/mouse/gse54650 4h/AOR/AOR.csv", filestyle="csv", outdir="E:/circidian_algorithm/data/mouse/4h_24_result/AOR",timepoints="Line1",outRawData=TRUE,minper=20,maxper=28, cycMethod = c("ARS", "JTK", "LS"))
##BFAT
meta2d(infile="E:/circidian_algorithm/data/mouse/gse54650 4h/BFAT/BFAT.csv", filestyle="csv", outdir="E:/circidian_algorithm/data/mouse/4h_24_result/BFAT",timepoints="Line1",outRawData=TRUE,minper=20,maxper=28, cycMethod = c("ARS", "JTK", "LS"))
##BST
meta2d(infile="E:/circidian_algorithm/data/mouse/gse54650 4h/BST/bst.csv", filestyle="csv", outdir="E:/circidian_algorithm/data/mouse/4h_24_result/BST",timepoints="Line1",outRawData=TRUE,minper=20,maxper=28, cycMethod = c("ARS", "JTK", "LS"))
##CER
meta2d(infile="E:/circidian_algorithm/data/mouse/gse54650 4h/CER/CER.csv", filestyle="csv", outdir="E:/circidian_algorithm/data/mouse/4h_24_result/CER",timepoints="Line1",outRawData=TRUE,minper=20,maxper=28, cycMethod = c("ARS", "JTK", "LS"))
##HEA
meta2d(infile="E:/circidian_algorithm/data/mouse/gse54650 4h/HEA/HEA.csv", filestyle="csv", outdir="E:/circidian_algorithm/data/mouse/4h_24_result/HEA",timepoints="Line1",outRawData=TRUE,minper=20,maxper=28, cycMethod = c("ARS", "JTK", "LS"))
##HYP
meta2d(infile="E:/circidian_algorithm/data/mouse/gse54650 4h/HYP/HYP.csv", filestyle="csv", outdir="E:/circidian_algorithm/data/mouse/4h_24_result/HYP",timepoints="Line1",outRawData=TRUE,minper=20,maxper=28, cycMethod = c("ARS", "JTK", "LS"))
##KID
meta2d(infile="E:/circidian_algorithm/data/mouse/gse54650 4h/KID/KID.csv", filestyle="csv", outdir="E:/circidian_algorithm/data/mouse/4h_24_result/KID",timepoints="Line1",outRawData=TRUE,minper=20,maxper=28, cycMethod = c("ARS", "JTK", "LS"))
##LIV
meta2d(infile="E:/circidian_algorithm/data/mouse/gse54650 4h/LIV/LIV.csv", filestyle="csv", outdir="E:/circidian_algorithm/data/mouse/4h_24_result/LIV",timepoints="Line1",outRawData=TRUE,minper=20,maxper=28, cycMethod = c("ARS", "JTK", "LS"))
##LUN
meta2d(infile="E:/circidian_algorithm/data/mouse/gse54650 4h/LUN/LUN.csv", filestyle="csv", outdir="E:/circidian_algorithm/data/mouse/4h_24_result/LUN",timepoints="Line1",outRawData=TRUE,minper=20,maxper=28, cycMethod = c("ARS", "JTK", "LS"))
##MUS
meta2d(infile="E:/circidian_algorithm/data/mouse/gse54650 4h/MUS/MUS.csv", filestyle="csv", outdir="E:/circidian_algorithm/data/mouse/4h_24_result/MUS",timepoints="Line1",outRawData=TRUE,minper=20,maxper=28, cycMethod = c("ARS", "JTK", "LS"))
##WAT
meta2d(infile="E:/circidian_algorithm/data/mouse/gse54650 4h/WAT/WAT.csv", filestyle="csv", outdir="E:/circidian_algorithm/data/mouse/4h_24_result/WAT",timepoints="Line1",outRawData=TRUE,minper=20,maxper=28, cycMethod = c("ARS", "JTK", "LS"))

#fly
meta2d(infile="E:/circidian_algorithm/data/fly/JTK.flyold.csv", filestyle="csv", outdir="E:/circidian_algorithm/result/fly/meta2d",timepoints="Line1",outRawData=TRUE,minper=20,maxper=28, cycMethod = c("ARS", "JTK", "LS"))
meta2d(infile="E:/circidian_algorithm/data/fly/JTK.flyyoung.csv", filestyle="csv", outdir="E:/circidian_algorithm/result/fly/meta2d",timepoints="Line1",outRawData=TRUE,minper=20,maxper=28, cycMethod = c("ARS", "JTK", "LS"))

#yeast
meta2d(infile="E:/circidian_algorithm/data/yeast/JTK.Sample2.csv", filestyle="csv", outdir="E:/circidian_algorithm/result/yeast/meta2d",timepoints="Line1",outRawData=TRUE,minper=60,maxper=180, cycMethod = c("ARS", "JTK", "LS"),ARSdefaultPer=120,ARSmle="nomle")
meta2d(infile="E:/circidian_algorithm/data/yeast/JTK.Sample6.csv", filestyle="csv", outdir="E:/circidian_algorithm/result/yeast/meta2d",timepoints="Line1",outRawData=TRUE,minper=300,maxper=420, cycMethod = c("ARS", "JTK", "LS"),ARSdefaultPer=360,ARSmle="nomle")


#arabidopsis
meta2d(infile="E:/circidian_algorithm/data/arabidopsis/flower.csv", filestyle="csv", outdir="E:/circidian_algorithm/result/arabidopsis/meta2d",timepoints="Line1",outRawData=TRUE,minper=20,maxper=28, cycMethod = c("ARS", "JTK", "LS"))
meta2d(infile="E:/circidian_algorithm/data/arabidopsis/JTK.Agse5612.csv", filestyle="csv", outdir="E:/circidian_algorithm/result/arabidopsis/meta2d",timepoints="Line1",outRawData=TRUE,minper=20,maxper=28, cycMethod = c("ARS", "JTK", "LS"))
meta2d(infile="E:/circidian_algorithm/data/arabidopsis/flower_mean.csv", filestyle="csv", outdir="E:/circidian_algorithm/result/arabidopsis/meta2d",timepoints="Line1",outRawData=TRUE,minper=20,maxper=28, cycMethod = c("ARS", "JTK", "LS"))


#neurospora
meta2d(infile="E:/circidian_algorithm/data/neurospora/JTKresult_RNA.csv", filestyle="csv", outdir="E:/circidian_algorithm/result/neurospora/meta2d",timepoints="Line1",outRawData=TRUE,minper=20,maxper=28, cycMethod = c("ARS", "JTK", "LS"))

#cyanobacteria
meta2d(infile="E:/circidian_algorithm/data/cyanobacteria/JTK.LLREPLICATE2.csv", filestyle="csv", outdir="E:/circidian_algorithm/result/cyanobacteria/meta2d",timepoints="Line1",outRawData=TRUE,minper=20,maxper=28, cycMethod = c("ARS", "JTK", "LS"))
meta2d(infile="E:/circidian_algorithm/data/cyanobacteria/cyanobacteria_mean.csv", filestyle="csv", outdir="E:/circidian_algorithm/result/cyanobacteria/meta2d",timepoints="Line1",outRawData=TRUE,minper=20,maxper=28, cycMethod = c("ARS", "JTK", "LS"))


#archaea
meta2d(infile="E:/circidian_algorithm/data/archaea/hal_1.csv", filestyle="csv", outdir="E:/circidian_algorithm/result/archaea/meta2d",timepoints="Line1",outRawData=TRUE,minper=20,maxper=28, cycMethod = c("ARS", "JTK", "LS"))
meta2d(infile="E:/circidian_algorithm/data/archaea/hal_2.csv", filestyle="csv", outdir="E:/circidian_algorithm/result/archaea/meta2d",timepoints="Line1",outRawData=TRUE,minper=20,maxper=28, cycMethod = c("ARS", "JTK", "LS"))
meta2d(infile="E:/circidian_algorithm/data/archaea/hal_3.csv", filestyle="csv", outdir="E:/circidian_algorithm/result/archaea/meta2d",timepoints="Line1",outRawData=TRUE,minper=20,maxper=28, cycMethod = c("ARS", "JTK", "LS"))
meta2d(infile="E:/circidian_algorithm/data/archaea/hal_4.csv", filestyle="csv", outdir="E:/circidian_algorithm/result/archaea/meta2d",timepoints="Line1",outRawData=TRUE,minper=20,maxper=28, cycMethod = c("ARS", "JTK", "LS"))








#############################################################rain
#fly
data <- read.csv("E:/circidian_algorithm/data/fly/JTK.flyold.csv")
input_data <- t(data[,2:13])
rainresult <- rain(input_data, period=24, deltat=4, peak.border=c(0.3,0.7),verbose=FALSE) 
pval <- rainresult$pVal
bhq <- p.adjust(pval, method="fdr")
temp <- data.frame(data[,1], pval, bhq, rainresult$phase)
print(sum(pval < 0.05))
print(sum(bhq < 0.05))
write.csv(temp, "E:/circidian_algorithm/result/fly/rain/old.csv")

data <- read.csv("E:/circidian_algorithm/data/fly/JTK.flyyoung.csv")
input_data <- t(data[,2:13])
rainresult <- rain(input_data, period=24, deltat=4, peak.border=c(0.3,0.7),verbose=FALSE) 
pval <- rainresult$pVal
bhq <- p.adjust(pval, method="fdr")
temp <- data.frame(data[,1], pval, bhq, rainresult$phase)
print(sum(pval < 0.05))
print(sum(bhq < 0.05))
write.csv(temp, "E:/circidian_algorithm/result/fly/rain/young.csv")

#arabidopsis
data <- read.csv("E:/circidian_algorithm/data/arabidopsis/flower.csv")
input_data <- t(data[,2:19])
rainresult <- rain(input_data, deltat=4, period=24, nr.series=3,peak.border=c(0.3, 0.7), verbose=FALSE)
pval <- rainresult$pVal
bhq <- p.adjust(pval, method="fdr")
temp <- data.frame(data[,1], pval, bhq, rainresult$phase)
print(sum(pval < 0.05))
print(sum(bhq < 0.05))
write.csv(temp, "E:/circidian_algorithm/result/arabidopsis/rain/flower.csv")

data <- read.csv("E:/circidian_algorithm/data/arabidopsis/JTK.Agse5612.csv")
input_data <- t(data[,2:14])
rainresult <- rain(input_data, deltat=4, period=24,peak.border=c(0.3, 0.7), verbose=FALSE)
pval <- rainresult$pVal
bhq <- p.adjust(pval, method="fdr")
temp <- data.frame(data[,1], pval, bhq, rainresult$phase)
print(sum(pval < 0.05))
print(sum(bhq < 0.05))
write.csv(temp, "E:/circidian_algorithm/result/arabidopsis/rain/seed.csv")

#yeast
data <- read.csv("E:/circidian_algorithm/data/yeast/JTK.Sample2.csv")
input_data <- t(data[,2:21])
rainresult <- rain(input_data, deltat=13, period=120,peak.border=c(0.3, 0.7), verbose=FALSE)
pval <- rainresult$pVal
bhq <- p.adjust(pval, method="fdr")
temp <- data.frame(data[,1], pval, bhq, rainresult$phase)
print(sum(pval < 0.05))
print(sum(bhq < 0.05))
write.csv(temp, "E:/circidian_algorithm/result/yeast/rain/sample2.csv")

data <- read.csv("E:/circidian_algorithm/data/yeast/JTK.Sample6.csv")
input_data <- t(data[,2:24])
rainresult <- rain(input_data, deltat=36, period=360,peak.border=c(0.3, 0.7), verbose=FALSE)
pval <- rainresult$pVal
bhq <- p.adjust(pval, method="fdr")
temp <- data.frame(data[,1], pval, bhq, rainresult$phase)
print(sum(pval < 0.05))
print(sum(bhq < 0.05))
write.csv(temp, "E:/circidian_algorithm/result/yeast/rain/sample6.csv")


#neurospora
data <- read.csv("E:/circidian_algorithm/data/neurospora/JTKresult_RNA.csv")
input_data <- t(data[,2:13])
rainresult <- rain(input_data, deltat=2, period=24,peak.border=c(0.3, 0.7), verbose=FALSE)
pval <- rainresult$pVal
bhq <- p.adjust(pval, method="fdr")
temp <- data.frame(data[,1], pval, bhq, rainresult$phase)
print(sum(pval < 0.05))
print(sum(bhq < 0.05))
write.csv(temp, "E:/circidian_algorithm/result/neurospora/rain/neurospora.csv")

#cyanobacteria
data <- read.csv("E:/circidian_algorithm/data/cyanobacteria/JTK.LLREPLICATE2.csv")
input_data <- t(data[,2:21])
rainresult <- rain(input_data, deltat=4, period=24, nr.series=2,peak.border=c(0.3, 0.7), verbose=FALSE)
pval <- rainresult$pVal
bhq <- p.adjust(pval, method="fdr")
temp <- data.frame(data[,1], pval, bhq, rainresult$phase)
print(sum(pval < 0.05))
print(sum(bhq < 0.05))
write.csv(temp, "E:/circidian_algorithm/result/neurospora/rain/flower.csv")

#archaea
data <- read.csv("E:/circidian_algorithm/data/archaea/hal_1.csv")
input_data <- t(data[,2:20])
rainresult <- rain(input_data, deltat=3, period=24,peak.border=c(0.3, 0.7), verbose=FALSE)
pval <- rainresult$pVal
bhq <- p.adjust(pval, method="fdr")
temp <- data.frame(data[,1], pval, bhq, rainresult$phase)
print(sum(pval < 0.05))
print(sum(bhq < 0.05))
write.csv(temp, "E:/circidian_algorithm/result/archaea/rain/hal1.csv")

data <- read.csv("E:/circidian_algorithm/data/archaea/hal_2.csv")
input_data <- t(data[,2:20])
rainresult <- rain(input_data, deltat=3, period=24,peak.border=c(0.3, 0.7), verbose=FALSE)
pval <- rainresult$pVal
bhq <- p.adjust(pval, method="fdr")
temp <- data.frame(data[,1], pval, bhq, rainresult$phase)
print(sum(pval < 0.05))
print(sum(bhq < 0.05))
write.csv(temp, "E:/circidian_algorithm/result/archaea/rain/hal2.csv")

data <- read.csv("E:/circidian_algorithm/data/archaea/hal_3.csv")
input_data <- t(data[,2:9])
rainresult <- rain(input_data, deltat=3, period=24,peak.border=c(0.3, 0.7), verbose=FALSE)
pval <- rainresult$pVal
bhq <- p.adjust(pval, method="fdr")
temp <- data.frame(data[,1], pval, bhq, rainresult$phase)
print(sum(pval < 0.05))
print(sum(bhq < 0.05))
write.csv(temp, "E:/circidian_algorithm/result/archaea/rain/hal3.csv")

data <- read.csv("E:/circidian_algorithm/data/archaea/hal_4.csv")
input_data <- t(data[,2:20])
rainresult <- rain(input_data, deltat=3, period=24,peak.border=c(0.3, 0.7), verbose=FALSE)
pval <- rainresult$pVal
bhq <- p.adjust(pval, method="fdr")
temp <- data.frame(data[,1], pval, bhq, rainresult$phase)
print(sum(pval < 0.05))
print(sum(bhq < 0.05))
write.csv(temp, "E:/circidian_algorithm/result/archaea/rain/hal4.csv")



################################################################Fisher-G test
#fly
data <- read.csv("E:/circidian_algorithm/data/fly/JTK.flyold.csv")
input_data <- t(data[,2:13])
pval.caulobacter <- fisher.g.test(input_data)
fdr.out <- fdrtool(pval.caulobacter, statistic="pvalue")
fdr <- p.adjust(fdr.out$pval,method="fdr")
temp <- data.frame(data[,1], fdr.out$pval, fdr)
print(sum(fdr.out$pval < 0.05))
print(sum(fdr.out$qval < 0.05))
write.csv(temp, "E:/circidian_algorithm/result/fly/fihsergtest/old.csv")

data <- read.csv("E:/circidian_algorithm/data/fly/JTK.flyyoung.csv")
input_data <- t(data[,2:13])
pval.caulobacter <- fisher.g.test(input_data)
fdr.out <- fdrtool(pval.caulobacter, statistic="pvalue")
fdr <- p.adjust(fdr.out$pval,method="fdr")
temp <- data.frame(data[,1], fdr.out$pval, fdr)
print(sum(fdr.out$pval < 0.05))
print(sum(fdr.out$qval < 0.05))
write.csv(temp, "E:/circidian_algorithm/result/fly/fihsergtest/young.csv")

#arabidopsis
data <- read.csv("E:/circidian_algorithm/data/arabidopsis/flower_mean.csv")
input_data <- t(data[,2:7])
pval.caulobacter <- fisher.g.test(input_data)
fdr.out <- fdrtool(pval.caulobacter, statistic="pvalue")
fdr <- p.adjust(fdr.out$pval,method="fdr")
temp <- data.frame(data[,1], fdr.out$pval, fdr)
print(sum(fdr.out$pval < 0.05))
print(sum(fdr.out$qval < 0.05))
write.csv(temp, "E:/circidian_algorithm/result/arabidopsis/fihsergtest/flower.csv")

data <- read.csv("E:/circidian_algorithm/data/arabidopsis/JTK.Agse5612.csv")
input_data <- t(data[,2:14])
pval.caulobacter <- fisher.g.test(input_data)
fdr.out <- fdrtool(pval.caulobacter, statistic="pvalue")
fdr <- p.adjust(fdr.out$pval,method="fdr")
temp <- data.frame(data[,1], fdr.out$pval, fdr)
print(sum(fdr.out$pval < 0.05))
print(sum(fdr < 0.05))
write.csv(temp, "E:/circidian_algorithm/result/arabidopsis/fihsergtest/seed.csv")

#yeast
data <- read.csv("E:/circidian_algorithm/data/yeast/JTK.Sample2.csv")
input_data <- t(data[,2:21])
pval.caulobacter <- fisher.g.test(input_data)
fdr.out <- fdrtool(pval.caulobacter, statistic="pvalue")
fdr <- p.adjust(fdr.out$pval,method="fdr")
temp <- data.frame(data[,1], fdr.out$pval, fdr)
print(sum(fdr.out$pval < 0.05))
print(sum(fdr < 0.05))
write.csv(temp, "E:/circidian_algorithm/result/yeast/fihsergtest/sample2.csv")

data <- read.csv("E:/circidian_algorithm/data/yeast/JTK.Sample6.csv")
input_data <- t(data[,2:25])
pval.caulobacter <- fisher.g.test(input_data)
fdr.out <- fdrtool(pval.caulobacter, statistic="pvalue")
fdr <- p.adjust(fdr.out$pval,method="fdr")
temp <- data.frame(data[,1], fdr.out$pval, fdr)
print(sum(fdr.out$pval < 0.05))
print(sum(fdr < 0.05))
write.csv(temp, "E:/circidian_algorithm/result/yeast/fihsergtest/sample6.csv")

#neurospora
data <- read.csv("E:/circidian_algorithm/data/neurospora/JTKresult_RNA.csv")
input_data <- t(data[,2:13])
pval.caulobacter <- fisher.g.test(input_data)
fdr.out <- fdrtool(pval.caulobacter, statistic="pvalue")
fdr <- p.adjust(fdr.out$pval,method="fdr")
temp <- data.frame(data[,1], fdr.out$pval, fdr)
print(sum(fdr.out$pval < 0.05))
print(sum(fdr < 0.05))
write.csv(temp, "E:/circidian_algorithm/result/neurospora/fihsergtest/neurospora.csv")

#cyanobacteria
data <- read.csv("E:/circidian_algorithm/data/cyanobacteria/cyanobacteria_mean.csv")
input_data <- t(data[,2:11])
pval.caulobacter <- fisher.g.test(input_data)
fdr.out <- fdrtool(pval.caulobacter, statistic="pvalue")
fdr <- p.adjust(fdr.out$pval,method="fdr")
temp <- data.frame(data[,1], fdr.out$pval, fdr)
print(sum(fdr.out$pval < 0.05))
print(sum(fdr < 0.05))
write.csv(temp, "E:/circidian_algorithm/result/cyanobacteria/fihsergtest/cyanobacteria.csv")

#archaea
data <- read.csv("E:/circidian_algorithm/data/archaea/hal_1.csv")
input_data <- t(data[,2:20])
pval.caulobacter <- fisher.g.test(input_data)
fdr.out <- fdrtool(pval.caulobacter, statistic="pvalue")
fdr <- p.adjust(fdr.out$pval,method="fdr")
temp <- data.frame(data[,1], fdr.out$pval, fdr)
print(sum(fdr.out$pval < 0.05))
print(sum(fdr < 0.05))
write.csv(temp, "E:/circidian_algorithm/result/archaea/fihsergtest/hal1.csv")

data <- read.csv("E:/circidian_algorithm/data/archaea/hal_2.csv")
input_data <- t(data[,2:20])
pval.caulobacter <- fisher.g.test(input_data)
fdr.out <- fdrtool(pval.caulobacter, statistic="pvalue")
fdr <- p.adjust(fdr.out$pval,method="fdr")
temp <- data.frame(data[,1], fdr.out$pval, fdr)
print(sum(fdr.out$pval < 0.05))
print(sum(fdr < 0.05))
write.csv(temp, "E:/circidian_algorithm/result/archaea/fihsergtest/hal2.csv")

data <- read.csv("E:/circidian_algorithm/data/archaea/hal_3.csv")
input_data <- t(data[,2:10])
pval.caulobacter <- fisher.g.test(input_data)
fdr.out <- fdrtool(pval.caulobacter, statistic="pvalue")
fdr <- p.adjust(fdr.out$pval,method="fdr")
temp <- data.frame(data[,1], fdr.out$pval, fdr)
print(sum(fdr.out$pval < 0.05))
print(sum(fdr < 0.05))
write.csv(temp, "E:/circidian_algorithm/result/archaea/fihsergtest/hal3.csv")

data <- read.csv("E:/circidian_algorithm/data/archaea/hal_4.csv")
input_data <- t(data[,2:20])
pval.caulobacter <- fisher.g.test(input_data)
fdr.out <- fdrtool(pval.caulobacter, statistic="pvalue")
fdr <- p.adjust(fdr.out$pval,method="fdr")
temp <- data.frame(data[,1], fdr.out$pval, fdr)
print(sum(fdr.out$pval < 0.05))
print(sum(fdr < 0.05))
write.csv(temp, "E:/circidian_algorithm/result/archaea/fihsergtest/hal4.csv")


















###############################################compare 2h and 4h
#meta2d
tissuename <- c("adr","AOR","BFAT","bst","CER","HEA","HYP","KID","LIV","LUN","MUS","WAT")
jtk_table <- data.frame(0,0,0,0,0)
ls_table <- data.frame(0,0,0,0,0)
arser_table <- data.frame(0,0,0,0,0)
meta2d_table <- data.frame(0,0,0,0,0)
for(j in 1:length(tissuename))
{
  file_name1 <- paste("E:/circidian_algorithm/result/mouse/2h/JTK/meta2d_",tissuename[j],".csv",sep="")
  file_name2 <- paste("E:/circidian_algorithm/result/mouse/4h/JTK/meta2d_",tissuename[j],".csv",sep="")
  data1 <- read.csv(file_name1)
  data2 <- read.csv(file_name2)
  
  gene1 <- data1[data1$JTK_BH.Q < 0.05, 1]
  gene2 <- data2[data2$JTK_BH.Q < 0.05, 1]
  temp <- c(length(gene1), length(gene2), length(intersect(gene1, gene2)),(length(intersect(gene1, gene2)) / (length(gene1) + length(gene2) - length(intersect(gene1, gene2)))),(length(intersect(gene1, gene2))/length(gene1)))
  jtk_table <- rbind(jtk_table, temp)
  
  gene1 <- data1[data1$LS_BH.Q < 0.05, 1]
  gene2 <- data2[data2$LS_BH.Q < 0.05, 1]
  temp <- c(length(gene1), length(gene2), length(intersect(gene1, gene2)),(length(intersect(gene1, gene2)) / (length(gene1) + length(gene2) - length(intersect(gene1, gene2)))),(length(intersect(gene1, gene2))/length(gene1) ))
  ls_table <- rbind(ls_table, temp)
  
  gene1 <- data1[data1$ARS_BH.Q < 0.05, 1]
  gene2 <- data2[data2$ARS_BH.Q < 0.05, 1]
  temp <- c(length(gene1), length(gene2), length(intersect(gene1, gene2)),(length(intersect(gene1, gene2)) / (length(gene1) + length(gene2) - length(intersect(gene1, gene2)))),(length(intersect(gene1, gene2))/length(gene1) ))
  arser_table <- rbind(arser_table, temp)
  
  gene1 <- data1[data1$meta2d_BH.Q < 0.05, 1]
  gene2 <- data2[data2$meta2d_BH.Q < 0.05, 1]
  temp <- c(length(gene1), length(gene2), length(intersect(gene1, gene2)),(length(intersect(gene1, gene2)) / (length(gene1) + length(gene2) - length(intersect(gene1, gene2)))),(length(intersect(gene1, gene2))/length(gene1) ))
  meta2d_table <- rbind(meta2d_table, temp)
}
write.csv(jtk_table,"E:/circidian_algorithm/result/mouse/jtk_table.csv")
write.csv(ls_table,"E:/circidian_algorithm/result/mouse/ls_table.csv")
write.csv(arser_table,"E:/circidian_algorithm/result/mouse/arser_table.csv")
write.csv(meta2d_table,"E:/circidian_algorithm/result/mouse/meta2d_table.csv")



#biocycle
tissuename <- c("adr","AOR","BFAT","bst","CER","HEA","HYP","KID","LIV","LUN","MUS","WAT")
biocycle_table <- data.frame(0,0,0,0,0)
for(j in 1:length(tissuename))
{
  file_name1 <- paste("E:/circidian_algorithm/result/mouse/2h/biocycle/",tissuename[j],".tsv",sep="")
  file_name2 <- paste("E:/circidian_algorithm/result/mouse/4h/biocycle/",tissuename[j],".tsv",sep="")
  data1 <- read.table(file_name1, header=TRUE)
  data2 <- read.table(file_name2, header=TRUE)
  gene1 <- data1[data1$Q_VALUE<0.05,1]
  gene2 <- data2[data2$Q_VALUE<0.05,1]
  temp <- c(length(gene1), length(gene2), length(intersect(gene1, gene2)),(length(intersect(gene1, gene2)) / (length(gene1) + length(gene2) - length(intersect(gene1, gene2)))),(length(intersect(gene1, gene2))/length(gene1) ))
  biocycle_table <- rbind(biocycle_table, temp)
}
write.csv(biocycle_table, "E:/circidian_algorithm/result/mouse/biocycle_table.csv")

#rain
tissuename <- c("adr","AOR","BFAT","bst","CER","HEA","HYP","KID","LIV","LUN","MUS","WAT")
rain_table <- data.frame(0,0,0,0,0)
for(j in 1:length(tissuename))
{
  file_name1 <- paste("E:/circidian_algorithm/result/mouse/2h/rain/",tissuename[j],".csv",sep="")
  file_name2 <- paste("E:/circidian_algorithm/result/mouse/4h/rain/",tissuename[j],".csv",sep="")
  data1 <- read.csv(file_name1, header=TRUE)
  data2 <- read.csv(file_name2, header=TRUE)
  gene1 <- data1[data1$bhq<0.05,2]
  gene2 <- data2[data2$bhq<0.05,2]
  temp <- c(length(gene1), length(gene2), length(intersect(gene1, gene2)),(length(intersect(gene1, gene2)) / (length(gene1) + length(gene2) - length(intersect(gene1, gene2)))),(length(intersect(gene1, gene2))/length(gene1) ))
  rain_table <- rbind(rain_table, temp)
}
write.csv(rain_table, "E:/circidian_algorithm/result/mouse/rain_table.csv")



#eJTK
tissuename <- c("adr","AOR","BFAT","bst","CER","HEA","HYP","KID","LIV","LUN","MUS","WAT")
ejtk_table <- data.frame(0,0,0,0,0)
for(j in 1:length(tissuename))
{
  file_name1 <- paste("E:/circidian_algorithm/result/mouse/2h/eJTK/",tissuename[j],"_cos24_ph00-22_by2_a02-22_by2_OTHERTEXT_jtkout_GammaP.txt",sep="")
  file_name2 <- paste("E:/circidian_algorithm/result/mouse/4h/eJTK/",tissuename[j],"_cos24_ph00-20_by4_a04-20_by4_OTHERTEXT_jtkout_GammaP.txt",sep="")
  data1 <- read.table(file_name1, header=TRUE)
  data2 <- read.table(file_name2, header=TRUE)
  gene1 <- data1[data1$GammaBH<0.05,1]
  gene2 <- data2[data2$GammaBH<0.05,1]
  temp <- c(length(gene1), length(gene2), length(intersect(gene1, gene2)),(length(intersect(gene1, gene2)) / (length(gene1) + length(gene2) - length(intersect(gene1, gene2)))),(length(intersect(gene1, gene2))/length(gene1) ))
  ejtk_table <- rbind(ejtk_table, temp)
}
write.csv(ejtk_table, "E:/circidian_algorithm/result/mouse/eJTK_table.csv")


#cicrwave
tissuename <- c("adr","AOR","BFAT","bst","CER","HEA","HYP","KID","LIV","LUN","MUS","WAT")
ejtk_table <- data.frame(0,0,0,0,0)
for(j in 1:length(tissuename))
{
  file_name1 <- paste("E:/circidian_algorithm/result/mouse/2h/circwave/",tissuename[j],".csv",sep="")
  file_name2 <- paste("E:/circidian_algorithm/result/mouse/4h/circwave/",tissuename[j],".csv",sep="")
  data1 <- read.csv(file_name1)
  data2 <- read.csv(file_name2)
  fdr1 <- p.adjust(data1$p1, method = "fdr")
  fdr2 <- p.adjust(data2$p1, method = "fdr")
  gene1 <- data1[fdr1<0.05,1]
  gene2 <- data2[fdr2<0.05,1]
  temp <- c(length(gene1), length(gene2), length(intersect(gene1, gene2)),(length(intersect(gene1, gene2)) / (length(gene1) + length(gene2) - length(intersect(gene1, gene2)))),(length(intersect(gene1, gene2))/length(gene1) ))
  ejtk_table <- rbind(ejtk_table, temp)
}
write.csv(ejtk_table, "E:/circidian_algorithm/result/mouse/circwave_table.csv")

#gtest
tissuename <- c("adr","AOR","BFAT","bst","CER","HEA","HYP","KID","LIV","LUN","MUS","WAT")
ejtk_table <- data.frame(0,0,0,0,0)
for(j in 1:length(tissuename))
{
  file_name1 <- paste("E:/circidian_algorithm/result/mouse/2h/gtest/",tissuename[j],".csv",sep="")
  file_name2 <- paste("E:/circidian_algorithm/result/mouse/4h/gtest/",tissuename[j],".csv",sep="")
  data1 <- read.csv(file_name1)
  data2 <- read.csv(file_name2)
  gene1 <- data1[data1$fdr.out.qval<0.05,2]
  gene2 <- data2[data2$fdr.out.qval<0.05,2]
  temp <- c(length(gene1), length(gene2), length(intersect(gene1, gene2)),(length(intersect(gene1, gene2)) / (length(gene1) + length(gene2) - length(intersect(gene1, gene2)))),(length(intersect(gene1, gene2))/length(gene1) ))
  ejtk_table <- rbind(ejtk_table, temp)
}
write.csv(ejtk_table, "E:/circidian_algorithm/result/mouse/gtest_table.csv")


#meta2d & rain
tissuename <- c("adr","AOR","BFAT","bst","CER","HEA","HYP","KID","LIV","LUN","MUS","WAT")
meta2d_table <- data.frame(0,0,0,0,0)
for(j in 1:length(tissuename))
{
  file_name1 <- paste("E:/circidian_algorithm/result/mouse/4h/JTK/meta2d_",tissuename[j],".csv",sep="")
  file_name2 <- paste("E:/circidian_algorithm/result/mouse/4h/rain/",tissuename[j],".csv",sep="")
  data1 <- read.csv(file_name1)
  data2 <- read.csv(file_name2)
  
  gene1 <- data1[data1$meta2d_BH.Q < 0.05, 1]
  gene2 <- data2[data2$bhq < 0.05, 2]
  temp <- c(length(gene1), length(gene2), length(intersect(gene1, gene2)),(length(intersect(gene1, gene2)) / (length(gene1) + length(gene2) - length(intersect(gene1, gene2)))),(length(intersect(gene1, gene2))/length(gene1)))

  meta2d_table <- rbind(meta2d_table, temp)
}
write.csv(meta2d_table,"E:/circidian_algorithm/result/mouse/meta2d_rain_4h.csv")



data <- read.csv("E:/circidian_algorithm/result/cyanobacteria/circwave/JTK.out.csv")
fdr <- p.adjust(data$p1, method="fdr")
sum(fdr < 0.05)
fdr[fdr<0.05]





#JTK & rain &meta2d
tissuename <- c("adr","AOR","BFAT","bst","CER","HEA","HYP","KID","LIV","LUN","MUS","WAT")
jtk_rain_table <- data.frame(0,0,0,0,0)
for(j in 1:length(tissuename))
{
  file_name1 <- paste("E:/circidian_algorithm/result/mouse/2h/JTK/meta2d_",tissuename[j],".csv",sep="")
  file_name2 <- paste("E:/circidian_algorithm/result/mouse/2h/rain/",tissuename[j],".csv",sep="")
  data1 <- read.csv(file_name1)
  data2 <- read.csv(file_name2)
  
  gene1 <- data1[data1$JTK_BH.Q < 0.05, 1]
  gene2 <- data2[data2$bhq < 0.05, 2]
  temp <- c(length(gene1), length(gene2), length(intersect(gene1, gene2)),(length(intersect(gene1, gene2)) / (length(gene1) + length(gene2) - length(intersect(gene1, gene2)))),(length(intersect(gene1, gene2))/length(gene1)))
  
  jtk_rain_table <- rbind(jtk_rain_table, temp)
}
write.csv(jtk_rain_table,"E:/circidian_algorithm/result/mouse/jtk_rain.csv")

tissuename <- c("adr","AOR","BFAT","bst","CER","HEA","HYP","KID","LIV","LUN","MUS","WAT")
jtk_rain4_table <- data.frame(0,0,0,0,0)
for(j in 1:length(tissuename))
{
  file_name1 <- paste("E:/circidian_algorithm/result/mouse/2h/JTK/meta2d_",tissuename[j],".csv",sep="")
  file_name2 <- paste("E:/circidian_algorithm/result/mouse/4h/rain/",tissuename[j],".csv",sep="")
  data1 <- read.csv(file_name1)
  data2 <- read.csv(file_name2)
  
  gene1 <- data1[data1$JTK_BH.Q < 0.05, 1]
  gene2 <- data2[data2$bhq < 0.05, 2]
  temp <- c(length(gene1), length(gene2), length(intersect(gene1, gene2)),(length(intersect(gene1, gene2)) / (length(gene1) + length(gene2) - length(intersect(gene1, gene2)))),(length(intersect(gene1, gene2))/length(gene1)))
  
  jtk_rain4_table <- rbind(jtk_rain4_table, temp)
}
write.csv(jtk_rain4_table,"E:/circidian_algorithm/result/mouse/jtk_rain4.csv")

tissuename <- c("adr","AOR","BFAT","bst","CER","HEA","HYP","KID","LIV","LUN","MUS","WAT")
jtk_rainandmeta2d_table <- data.frame(0,0,0,0,0)
for(j in 1:length(tissuename))
{
  file_name1 <- paste("E:/circidian_algorithm/result/mouse/2h/JTK/meta2d_",tissuename[j],".csv",sep="")
  file_name2 <- paste("E:/circidian_algorithm/result/mouse/4h/rain/",tissuename[j],".csv",sep="")
  file_name3 <- paste("E:/circidian_algorithm/result/mouse/4h/JTK/meta2d_",tissuename[j],".csv",sep="")
  data1 <- read.csv(file_name1)
  data2 <- read.csv(file_name2)
  data3 <- read.csv(file_name3)
  
  gene1 <- data1[data1$JTK_BH.Q < 0.05, 1]
  gene2 <- data2[data2$bhq < 0.05, 2]
  gene3 <- data3[data3$meta2d_BH.Q < 0.05, 1]
  gene2 <- union(gene2, gene3)
  temp <- c(length(gene1), length(gene2), length(intersect(gene1, gene2)),(length(intersect(gene1, gene2)) / (length(gene1) + length(gene2) - length(intersect(gene1, gene2)))),(length(intersect(gene1, gene2))/length(gene1)))
  
  jtk_rainandmeta2d_table <- rbind(jtk_rainandmeta2d_table, temp)
}
write.csv(jtk_rainandmeta2d_table,"E:/circidian_algorithm/result/mouse/jtk_rainandmeta2d.csv")


tissuename <- c("adr","AOR","BFAT","bst","CER","HEA","HYP","KID","LIV","LUN","MUS","WAT")
jtk_rainandmeta2d_table <- data.frame(0,0,0,0,0)
for(j in 1:length(tissuename))
{
  file_name1 <- paste("E:/circidian_algorithm/result/mouse/2h/JTK/meta2d_",tissuename[j],".csv",sep="")
  file_name2 <- paste("E:/circidian_algorithm/result/mouse/2h/JTK/meta2d_",tissuename[j],".csv",sep="")
  data1 <- read.csv(file_name1)
  data2 <- read.csv(file_name2)
  data3 <- read.csv(file_name3)
  
  gene1 <- data1[data1$JTK_BH.Q < 0.05, 1]
  gene2 <- data2[data2$meta2d_BH.Q < 0.05, 1]
  temp <- c(length(gene1), length(gene2), length(intersect(gene1, gene2)),(length(intersect(gene1, gene2)) / (length(gene1) + length(gene2) - length(intersect(gene1, gene2)))),(length(intersect(gene1, gene2))/length(gene1)))
  
  jtk_rainandmeta2d_table <- rbind(jtk_rainandmeta2d_table, temp)
}
write.csv(jtk_rainandmeta2d_table,"E:/circidian_algorithm/result/mouse/jtk_meta2d2h.csv")


tissuename <- c("adr","AOR","BFAT","bst","CER","HEA","HYP","KID","LIV","LUN","MUS","WAT")
jtk_rainandmeta2d_table <- data.frame(0,0,0,0,0)
for(j in 1:length(tissuename))
{
  file_name1 <- paste("E:/circidian_algorithm/result/mouse/2h/JTK/meta2d_",tissuename[j],".csv",sep="")
  file_name2 <- paste("E:/circidian_algorithm/result/mouse/4h/JTK/meta2d_",tissuename[j],".csv",sep="")
  data1 <- read.csv(file_name1)
  data2 <- read.csv(file_name2)
  data3 <- read.csv(file_name3)
  
  gene1 <- data1[data1$JTK_BH.Q < 0.05, 1]
  gene2 <- data2[data2$meta2d_BH.Q < 0.05, 1]
  temp <- c(length(gene1), length(gene2), length(intersect(gene1, gene2)),(length(intersect(gene1, gene2)) / (length(gene1) + length(gene2) - length(intersect(gene1, gene2)))),(length(intersect(gene1, gene2))/length(gene1)))
  
  jtk_rainandmeta2d_table <- rbind(jtk_rainandmeta2d_table, temp)
}
write.csv(jtk_rainandmeta2d_table,"E:/circidian_algorithm/result/mouse/jtk_meta2d4h.csv")






###############################intersect and union
tissuename <- c("adr","AOR","BFAT","bst","CER","HEA","HYP","KID","LIV","LUN","MUS","WAT")
jtk_rainandmeta2d_table <- data.frame(0,0,0,0,0)
for(j in 1:length(tissuename))
{
  file_name1 <- paste("E:/circidian_algorithm/result/mouse/2h/JTK/meta2d_",tissuename[j],".csv",sep="")
  file_name2 <- paste("E:/circidian_algorithm/result/mouse/4h/JTK/meta2d_",tissuename[j],".csv",sep="")
  file_name3 <- paste("E:/circidian_algorithm/result/mouse/4h/rain/",tissuename[j],".csv",sep="")
  file_name4 <- paste("E:/circidian_algorithm/result/mouse/4h/biocycle/",tissuename[j],".tsv",sep="")
  file_name5 <- paste("E:/circidian_algorithm/result/mouse/4h/circwave/",tissuename[j],".csv",sep="")
  data1 <- read.csv(file_name1)
  data2 <- read.csv(file_name2)
  data3 <- read.csv(file_name3)
  data4 <- read.table(file_name4,header=TRUE)
  data5 <- read.csv(file_name5)
  
  fdr <- p.adjust(data4$P_VALUE, method="fdr")
  fdr2 <- p.adjust(data5$p1, method="fdr")
  gene1 <- data1[data1$JTK_BH.Q < 0.05, 1]
  gene2 <- data2[data2$meta2d_BH.Q < 0.05, 1]
  gene3 <- data3[data3$bhq < 0.05, 2]
  gene4 <- data4[fdr < 0.05, 1]
  gene5 <- data5[fdr2 < 0.05, 1]
  inter1 <- intersect(gene2, gene3)
  inter2 <- intersect(gene2, gene4)
  inter3 <- intersect(gene2, gene5)
  inter4 <- intersect(gene3, gene4)
  inter5 <- intersect(gene3, gene5)
  inter6 <- intersect(gene4, gene5)
  union_gene <- unique(inter1,inter2,inter3,inter4,inter5,inter6)
  length(intersect(gene1,union_gene))
  gene2 <- intersect(gene1,union_gene)
  
  temp <- c(length(gene1), length(gene2), length(intersect(gene1, gene2)),(length(intersect(gene1, gene2)) / (length(gene1) + length(gene2) - length(intersect(gene1, gene2)))),(length(intersect(gene1, gene2))/length(gene1)))
  
  jtk_rainandmeta2d_table <- rbind(jtk_rainandmeta2d_table, temp)
}
write.csv(jtk_rainandmeta2d_table,"E:/circidian_algorithm/result/mouse/union.csv")




tissuename <- c("adr","AOR","BFAT","bst","CER","HEA","HYP","KID","LIV","LUN","MUS","WAT")
jtk_rainandmeta2d_table <- data.frame(0,0,0,0,0)
for(j in 1:length(tissuename))
{
  file_name1 <- paste("E:/circidian_algorithm/result/mouse/2h/JTK/meta2d_",tissuename[j],".csv",sep="")
  file_name2 <- paste("E:/circidian_algorithm/result/mouse/4h/rain/",tissuename[j],".csv",sep="")
  file_name3 <- paste("E:/circidian_algorithm/result/mouse/4h/JTK/meta2d_",tissuename[j],".csv",sep="")
  data1 <- read.csv(file_name1)
  data2 <- read.csv(file_name2)
  data3 <- read.csv(file_name3)
  
  gene1 <- data1[data1$ARS_BH.Q < 0.05, 1]
  gene2 <- data2[data2$bhq < 0.05, 2]
  gene3 <- data3[data3$meta2d_BH.Q < 0.05, 1]
  gene2 <- union(gene2, gene3)
  temp <- c(length(gene1), length(gene2), length(intersect(gene1, gene2)),(length(intersect(gene1, gene2)) / (length(gene1) + length(gene2) - length(intersect(gene1, gene2)))),(length(intersect(gene1, gene2))/length(gene1)))
  
  jtk_rainandmeta2d_table <- rbind(jtk_rainandmeta2d_table, temp)
}
write.csv(jtk_rainandmeta2d_table,"E:/circidian_algorithm/result/mouse/arser_rainandmeta2d.csv")









tissuename <- c("adr","AOR","BFAT","bst","CER","HEA","HYP","KID","LIV","LUN","MUS","WAT")
jtk_rainandmeta2d_table <- data.frame(0,0,0,0,0)
for(j in 1:length(tissuename))
{
  file_name1 <- paste("E:/circidian_algorithm/result/mouse/2h/JTK/meta2d_",tissuename[j],".csv",sep="")
  file_name2 <- paste("E:/circidian_algorithm/result/mouse/4h/JTK/meta2d_",tissuename[j],".csv",sep="")
  data1 <- read.csv(file_name1)
  data2 <- read.csv(file_name2)
  
  gene1 <- data1[data1$ARS_BH.Q < 0.05, 1]
  gene2 <- data2[data2$meta2d_BH.Q < 0.05, 1]
  temp <- c(length(gene1), length(gene2), length(intersect(gene1, gene2)),(length(intersect(gene1, gene2)) / (length(gene1) + length(gene2) - length(intersect(gene1, gene2)))),(length(intersect(gene1, gene2))/length(gene1)))
  
  jtk_rainandmeta2d_table <- rbind(jtk_rainandmeta2d_table, temp)
}
write.csv(jtk_rainandmeta2d_table,"E:/circidian_algorithm/result/mouse/arser_meta2d4h.csv")

tissuename <- c("adr","AOR","BFAT","bst","CER","HEA","HYP","KID","LIV","LUN","MUS","WAT")
jtk_rain4_table <- data.frame(0,0,0,0,0)
for(j in 1:length(tissuename))
{
  file_name1 <- paste("E:/circidian_algorithm/result/mouse/2h/JTK/meta2d_",tissuename[j],".csv",sep="")
  file_name2 <- paste("E:/circidian_algorithm/result/mouse/4h/rain/",tissuename[j],".csv",sep="")
  data1 <- read.csv(file_name1)
  data2 <- read.csv(file_name2)
  
  gene1 <- data1[data1$ARS_BH.Q < 0.05, 1]
  gene2 <- data2[data2$bhq < 0.05, 2]
  temp <- c(length(gene1), length(gene2), length(intersect(gene1, gene2)),(length(intersect(gene1, gene2)) / (length(gene1) + length(gene2) - length(intersect(gene1, gene2)))),(length(intersect(gene1, gene2))/length(gene1)))
  
  jtk_rain4_table <- rbind(jtk_rain4_table, temp)
}
write.csv(jtk_rain4_table,"E:/circidian_algorithm/result/mouse/arser_rain4.csv")




#########find best 2h
tissuename <- c("adr","AOR","BFAT","bst","CER","HEA","HYP","KID","LIV","LUN","MUS","WAT")
jtk_rain4_table <- data.frame(0,0,0,0,0)
for(j in 1:length(tissuename))
{
  file_name1 <- paste("E:/circidian_algorithm/result/mouse/2h/JTK/meta2d_",tissuename[j],".csv",sep="")
  file_name2 <- paste("E:/circidian_algorithm/result/mouse/2h/rain/",tissuename[j],".csv",sep="")
  file_name3 <- paste("E:/circidian_algorithm/result/mouse/4h/JTK/meta2d_",tissuename[j],".csv",sep="")
  file_name4 <- paste("E:/circidian_algorithm/result/mouse/4h/rain/",tissuename[j],".csv",sep="")
  data1 <- read.csv(file_name1)
  data2 <- read.csv(file_name2)
  data3 <- read.csv(file_name3)
  data4 <- read.csv(file_name4)
  
  gene1 <- data1[data1$meta2d_BH.Q < 0.05, 1]
  gene2 <- data2[data2$bhq < 0.05, 2]
  gene3 <- data1[data3$meta2d_BH.Q < 0.05, 1]
  gene4 <- data2[data4$bhq < 0.05, 2]
  gene1 <- union(gene1, gene2)
  gene2 <- union(gene3, gene4)
  temp <- c(length(gene1), length(gene2), length(intersect(gene1, gene2)),(length(intersect(gene1, gene2)) / (length(gene1) + length(gene2) - length(intersect(gene1, gene2)))),(length(intersect(gene1, gene2))/length(gene1)))
  
  jtk_rain4_table <- rbind(jtk_rain4_table, temp)
}
write.csv(jtk_rain4_table,"E:/circidian_algorithm/result/mouse/meta2drain2h_meta2drain4h.csv")


tissuename <- c("adr","AOR","BFAT","bst","CER","HEA","HYP","KID","LIV","LUN","MUS","WAT")
jtk_rain4_table <- data.frame(0,0,0,0,0)
for(j in 1:length(tissuename))
{
  file_name1 <- paste("E:/circidian_algorithm/result/mouse/2h/JTK/meta2d_",tissuename[j],".csv",sep="")
  file_name2 <- paste("E:/circidian_algorithm/result/mouse/2h/rain/",tissuename[j],".csv",sep="")
  file_name3 <- paste("E:/circidian_algorithm/result/mouse/2h/JTK/meta2d_",tissuename[j],".csv",sep="")
  file_name4 <- paste("E:/circidian_algorithm/result/mouse/4h/JTK/meta2d_",tissuename[j],".csv",sep="")
  file_name5 <- paste("E:/circidian_algorithm/result/mouse/4h/rain/",tissuename[j],".csv",sep="")
  file_name6 <- paste("E:/circidian_algorithm/result/mouse/4h/JTK/meta2d_",tissuename[j],".csv",sep="")
  data1 <- read.csv(file_name1)
  data2 <- read.csv(file_name2)
  data3 <- read.csv(file_name3)
  data4 <- read.csv(file_name4)
  data5 <- read.csv(file_name5)
  data6 <- read.csv(file_name6)
  
  gene1 <- data1[data1$JTK_BH.Q < 0.05, 1]
  gene2 <- data2[data2$bhq < 0.05, 2]
  gene3 <- data3[data3$meta2d_BH.Q < 0.05, 1]
  gene4 <- data4[data4$JTK_BH.Q < 0.05, 1]
  gene5 <- data5[data5$bhq < 0.05, 2]
  gene6 <- data6[data6$meta2d_BH.Q < 0.05, 1]
  gene11 <- intersect(gene1, gene2)
  gene12 <- intersect(gene1, gene3)
  gene13 <- intersect(gene2, gene3)
  gene21 <- intersect(gene4, gene5)
  gene22 <- intersect(gene4, gene6)
  gene23 <- intersect(gene5, gene6)
  gene1 <- unique(c(gene11,gene12,gene13))
  gene2 <- unique(c(gene21,gene22,gene23))
  temp <- c(length(gene1), length(gene2), length(intersect(gene1, gene2)),(length(intersect(gene1, gene2)) / (length(gene1) + length(gene2) - length(intersect(gene1, gene2)))),(length(intersect(gene1, gene2))/length(gene1)))
  
  jtk_rain4_table <- rbind(jtk_rain4_table, temp)
}
write.csv(jtk_rain4_table,"E:/circidian_algorithm/result/mouse/3union2h_3union4h.csv")



tissuename <- c("adr","AOR","BFAT","bst","CER","HEA","HYP","KID","LIV","LUN","MUS","WAT")
jtk_rain4_table <- data.frame(0,0,0,0,0)
for(j in 1:length(tissuename))
{
  file_name1 <- paste("E:/circidian_algorithm/result/mouse/2h/JTK/meta2d_",tissuename[j],".csv",sep="")
  file_name2 <- paste("E:/circidian_algorithm/result/mouse/2h/rain/",tissuename[j],".csv",sep="")
  file_name3 <- paste("E:/circidian_algorithm/result/mouse/4h/JTK/meta2d_",tissuename[j],".csv",sep="")
  file_name4 <- paste("E:/circidian_algorithm/result/mouse/4h/rain/",tissuename[j],".csv",sep="")
  data1 <- read.csv(file_name1)
  data2 <- read.csv(file_name2)
  data3 <- read.csv(file_name3)
  data4 <- read.csv(file_name4)
  
  gene1 <- data1[data1$JTK_BH.Q < 0.05, 1]
  gene2 <- data2[data2$bhq < 0.05, 2]
  gene3 <- data1[data3$JTK_BH.Q < 0.05, 1]
  gene4 <- data2[data4$bhq < 0.05, 2]
  gene1 <- union(gene1, gene2)
  gene2 <- union(gene3, gene4)
  temp <- c(length(gene1), length(gene2), length(intersect(gene1, gene2)),(length(intersect(gene1, gene2)) / (length(gene1) + length(gene2) - length(intersect(gene1, gene2)))),(length(intersect(gene1, gene2))/length(gene1)))
  
  jtk_rain4_table <- rbind(jtk_rain4_table, temp)
}
write.csv(jtk_rain4_table,"E:/circidian_algorithm/result/mouse/JTKrain2h_JTKrain4h.csv")


tissuename <- c("adr","AOR","BFAT","bst","CER","HEA","HYP","KID","LIV","LUN","MUS","WAT")
jtk_rain4_table <- data.frame(0,0,0,0,0)
for(j in 1:length(tissuename))
{
  file_name1 <- paste("E:/circidian_algorithm/result/mouse/2h/JTK/meta2d_",tissuename[j],".csv",sep="")
  file_name2 <- paste("E:/circidian_algorithm/result/mouse/2h/JTK/meta2d_",tissuename[j],".csv",sep="")
  file_name3 <- paste("E:/circidian_algorithm/result/mouse/4h/JTK/meta2d_",tissuename[j],".csv",sep="")
  file_name4 <- paste("E:/circidian_algorithm/result/mouse/4h/JTK/meta2d_",tissuename[j],".csv",sep="")
  data1 <- read.csv(file_name1)
  data2 <- read.csv(file_name2)
  data3 <- read.csv(file_name3)
  data4 <- read.csv(file_name4)
  
  gene1 <- data1[data1$JTK_BH.Q < 0.05, 1]
  gene2 <- data2[data2$bhq < 0.05, 2]
  gene3 <- data1[data3$JTK_BH.Q < 0.05, 1]
  gene4 <- data2[data4$bhq < 0.05, 2]
  gene1 <- union(gene1, gene2)
  gene2 <- union(gene3, gene4)
  temp <- c(length(gene1), length(gene2), length(intersect(gene1, gene2)),(length(intersect(gene1, gene2)) / (length(gene1) + length(gene2) - length(intersect(gene1, gene2)))),(length(intersect(gene1, gene2))/length(gene1)))
  
  jtk_rain4_table <- rbind(jtk_rain4_table, temp)
}
write.csv(jtk_rain4_table,"E:/circidian_algorithm/result/mouse/meta2dJTK_meta2dJTK.csv")



tissuename <- c("adr","AOR","BFAT","bst","CER","HEA","HYP","KID","LIV","LUN","MUS","WAT")
jtk_rainandmeta2d_table <- data.frame(0,0,0,0,0)
for(j in 1:length(tissuename))
{
  file_name1 <- paste("E:/circidian_algorithm/result/mouse/2h/rain/",tissuename[j],".csv",sep="")
  file_name2 <- paste("E:/circidian_algorithm/result/mouse/2h/biocycle/",tissuename[j],".tsv",sep="")
  file_name3 <- paste("E:/circidian_algorithm/result/mouse/4h/rain/",tissuename[j],".csv",sep="")
  file_name4 <- paste("E:/circidian_algorithm/result/mouse/4h/biocycle/",tissuename[j],".tsv",sep="")
  
  data1 <- read.csv(file_name1)
  data2 <- read.table(file_name2,header=TRUE)
  data3 <- read.csv(file_name3)
  data4 <- read.table(file_name4,header=TRUE)
  
  fdr <- p.adjust(data2$P_VALUE, method="fdr")
  fdr2 <- p.adjust(data4$P_VALUE, method="fdr")
  gene1 <- data1[data1$bhq < 0.05, 2]
  gene2 <- data2[fdr < 0.05, 1]
  gene3 <- data3[data3$bhq < 0.05, 2]
  gene4 <- data4[fdr2 < 0.05, 1]
  
  gene1 <- union(gene1,gene2)
  gene2 <- union(gene3,gene4)
  
  temp <- c(length(gene1), length(gene2), length(intersect(gene1, gene2)),(length(intersect(gene1, gene2)) / (length(gene1) + length(gene2) - length(intersect(gene1, gene2)))),(length(intersect(gene1, gene2))/length(gene1)))
  
  jtk_rainandmeta2d_table <- rbind(jtk_rainandmeta2d_table, temp)
}
write.csv(jtk_rainandmeta2d_table,"E:/circidian_algorithm/result/mouse/rainbiocycle_rainbiocycle.csv")


tissuename <- c("adr","AOR","BFAT","bst","CER","HEA","HYP","KID","LIV","LUN","MUS","WAT")
jtk_rainandmeta2d_table <- data.frame(0,0,0,0,0)
for(j in 1:length(tissuename))
{
  file_name1 <- paste("E:/circidian_algorithm/result/mouse/2h/JTK/meta2d_",tissuename[j],".csv",sep="")
  file_name2 <- paste("E:/circidian_algorithm/result/mouse/2h/biocycle/",tissuename[j],".tsv",sep="")
  file_name3 <- paste("E:/circidian_algorithm/result/mouse/4h/JTK/meta2d_",tissuename[j],".csv",sep="")
  file_name4 <- paste("E:/circidian_algorithm/result/mouse/4h/biocycle/",tissuename[j],".tsv",sep="")
  
  data1 <- read.csv(file_name1)
  data2 <- read.table(file_name2,header=TRUE)
  data3 <- read.csv(file_name3)
  data4 <- read.table(file_name4,header=TRUE)
  
  fdr <- p.adjust(data2$P_VALUE, method="fdr")
  fdr2 <- p.adjust(data4$P_VALUE, method="fdr")
  gene1 <- data1[data1$JTK_BH.Q < 0.05, 1]
  gene2 <- data2[fdr < 0.05, 1]
  gene3 <- data3[data3$JTK_BH.Q < 0.05, 1]
  gene4 <- data4[fdr2 < 0.05, 1]
  
  gene1 <- union(gene1,gene2)
  gene2 <- union(gene3,gene4)
  
  temp <- c(length(gene1), length(gene2), length(intersect(gene1, gene2)),(length(intersect(gene1, gene2)) / (length(gene1) + length(gene2) - length(intersect(gene1, gene2)))),(length(intersect(gene1, gene2))/length(gene1)))
  
  jtk_rainandmeta2d_table <- rbind(jtk_rainandmeta2d_table, temp)
}
write.csv(jtk_rainandmeta2d_table,"E:/circidian_algorithm/result/mouse/JTKbiocycle_JTKbiocycle.csv")


tissuename <- c("adr","AOR","BFAT","bst","CER","HEA","HYP","KID","LIV","LUN","MUS","WAT")
jtk_rainandmeta2d_table <- data.frame(0,0,0,0,0)
for(j in 1:length(tissuename))
{
  file_name1 <- paste("E:/circidian_algorithm/result/mouse/2h/JTK/meta2d_",tissuename[j],".csv",sep="")
  file_name2 <- paste("E:/circidian_algorithm/result/mouse/2h/biocycle/",tissuename[j],".tsv",sep="")
  file_name3 <- paste("E:/circidian_algorithm/result/mouse/4h/JTK/meta2d_",tissuename[j],".csv",sep="")
  file_name4 <- paste("E:/circidian_algorithm/result/mouse/4h/biocycle/",tissuename[j],".tsv",sep="")
  
  data1 <- read.csv(file_name1)
  data2 <- read.table(file_name2,header=TRUE)
  data3 <- read.csv(file_name3)
  data4 <- read.table(file_name4,header=TRUE)
  
  fdr <- p.adjust(data2$P_VALUE, method="fdr")
  fdr2 <- p.adjust(data4$P_VALUE, method="fdr")
  gene1 <- data1[data1$meta2d_BH.Q < 0.05, 1]
  gene2 <- data2[fdr < 0.05, 1]
  gene3 <- data3[data3$meta2d_BH.Q < 0.05, 1]
  gene4 <- data4[fdr2 < 0.05, 1]
  
  gene1 <- union(gene1,gene2)
  gene2 <- union(gene3,gene4)
  
  temp <- c(length(gene1), length(gene2), length(intersect(gene1, gene2)),(length(intersect(gene1, gene2)) / (length(gene1) + length(gene2) - length(intersect(gene1, gene2)))),(length(intersect(gene1, gene2))/length(gene1)))
  
  jtk_rainandmeta2d_table <- rbind(jtk_rainandmeta2d_table, temp)
}
write.csv(jtk_rainandmeta2d_table,"E:/circidian_algorithm/result/mouse/meta2dbiocycle_meta2dbiocycle.csv")




tissuename <- c("adr","AOR","BFAT","bst","CER","HEA","HYP","KID","LIV","LUN","MUS","WAT")
jtk_rain4_table <- data.frame(0,0,0,0,0)
for(j in 1:length(tissuename))
{
  file_name1 <- paste("E:/circidian_algorithm/result/mouse/2h/biocycle/",tissuename[j],".tsv",sep="")
  file_name2 <- paste("E:/circidian_algorithm/result/mouse/2h/rain/",tissuename[j],".csv",sep="")
  file_name3 <- paste("E:/circidian_algorithm/result/mouse/2h/JTK/meta2d_",tissuename[j],".csv",sep="")
  file_name4 <- paste("E:/circidian_algorithm/result/mouse/4h/biocycle/",tissuename[j],".tsv",sep="")
  file_name5 <- paste("E:/circidian_algorithm/result/mouse/4h/rain/",tissuename[j],".csv",sep="")
  file_name6 <- paste("E:/circidian_algorithm/result/mouse/4h/JTK/meta2d_",tissuename[j],".csv",sep="")
  data1 <- read.table(file_name1,header=TRUE)
  data2 <- read.csv(file_name2)
  data3 <- read.csv(file_name3)
  data4 <- read.table(file_name4,header=TRUE)
  data5 <- read.csv(file_name5)
  data6 <- read.csv(file_name6)
  data6$meta2d_BH.Q
  
  fdr <- p.adjust(data1$P_VALUE, method="fdr")
  fdr2 <- p.adjust(data4$P_VALUE, method="fdr")
  gene1 <- data1[fdr < 0.05, 1]
  gene2 <- data2[data2$bhq < 0.05, 2]
  gene3 <- data3[data3$meta2d_BH.Q < 0.05, 1]
  gene4 <- data4[fdr2 < 0.05, 1]
  gene5 <- data5[data5$bhq < 0.05, 2]
  gene6 <- data6[data6$meta2d_BH.Q < 0.05, 1]
  gene11 <- intersect(gene1, gene2)
  gene12 <- intersect(gene1, gene3)
  gene13 <- intersect(gene2, gene3)
  gene21 <- intersect(gene4, gene5)
  gene22 <- intersect(gene4, gene6)
  gene23 <- intersect(gene5, gene6)
  gene1 <- unique(c(gene11,gene12,gene13))
  gene2 <- unique(c(gene21,gene22,gene23))
  temp <- c(length(gene1), length(gene2), length(intersect(gene1, gene2)),(length(intersect(gene1, gene2)) / (length(gene1) + length(gene2) - length(intersect(gene1, gene2)))),(length(intersect(gene1, gene2))/length(gene1)))
  
  jtk_rain4_table <- rbind(jtk_rain4_table, temp)
}
write.csv(jtk_rain4_table,"E:/circidian_algorithm/result/mouse/meta2dbiocyclerain2h_meta2dbiocyclerain4h.csv")




tissuename <- c("adr","AOR","BFAT","bst","CER","HEA","HYP","KID","LIV","LUN","MUS","WAT")
jtk_rainandmeta2d_table <- data.frame(0,0,0,0,0)
for(j in 1:length(tissuename))
{
  file_name1 <- paste("E:/circidian_algorithm/result/mouse/2h/JTK/meta2d_",tissuename[j],".csv",sep="")
  file_name2 <- paste("E:/circidian_algorithm/result/mouse/2h/biocycle/",tissuename[j],".tsv",sep="")
  file_name3 <- paste("E:/circidian_algorithm/result/mouse/4h/JTK/meta2d_",tissuename[j],".csv",sep="")
  file_name4 <- paste("E:/circidian_algorithm/result/mouse/4h/biocycle/",tissuename[j],".tsv",sep="")
  
  data1 <- read.csv(file_name1)
  data2 <- read.table(file_name2,header=TRUE)
  data3 <- read.csv(file_name3)
  data4 <- read.table(file_name4,header=TRUE)
  
  fdr <- p.adjust(data2$P_VALUE, method="fdr")
  fdr2 <- p.adjust(data4$P_VALUE, method="fdr")
  gene1 <- data1[data1$ARS_BH.Q < 0.05, 1]
  gene2 <- data2[fdr < 0.05, 1]
  gene3 <- data3[data3$ARS_BH.Q < 0.05, 1]
  gene4 <- data4[fdr2 < 0.05, 1]
  
  gene1 <- union(gene1,gene2)
  gene2 <- union(gene3,gene4)
  
  temp <- c(length(gene1), length(gene2), length(intersect(gene1, gene2)),(length(intersect(gene1, gene2)) / (length(gene1) + length(gene2) - length(intersect(gene1, gene2)))),(length(intersect(gene1, gene2))/length(gene1)))
  
  jtk_rainandmeta2d_table <- rbind(jtk_rainandmeta2d_table, temp)
}
write.csv(jtk_rainandmeta2d_table,"E:/circidian_algorithm/result/mouse/arserbiocycle_arserbiocycle.csv")


tissuename <- c("adr","AOR","BFAT","bst","CER","HEA","HYP","KID","LIV","LUN","MUS","WAT")
jtk_rainandmeta2d_table <- data.frame(0,0,0,0,0)
for(j in 1:length(tissuename))
{
  file_name1 <- paste("E:/circidian_algorithm/result/mouse/2h/JTK/meta2d_",tissuename[j],".csv",sep="")
  file_name2 <- paste("E:/circidian_algorithm/result/mouse/2h/JTK/meta2d_",tissuename[j],".csv",sep="")
  file_name3 <- paste("E:/circidian_algorithm/result/mouse/4h/JTK/meta2d_",tissuename[j],".csv",sep="")
  file_name4 <- paste("E:/circidian_algorithm/result/mouse/4h/JTK/meta2d_",tissuename[j],".csv",sep="")
  
  data1 <- read.csv(file_name1)
  data2 <- read.csv(file_name2,header=TRUE)
  data3 <- read.csv(file_name3)
  data4 <- read.csv(file_name4,header=TRUE)

  gene1 <- data1[data1$ARS_BH.Q < 0.05, 1]
  gene2 <- data2[data2$JTK_BH.Q < 0.05, 1]
  gene3 <- data3[data3$ARS_BH.Q < 0.05, 1]
  gene4 <- data4[data4$JTK_BH.Q < 0.05, 1]
  
  gene1 <- union(gene1,gene2)
  gene2 <- union(gene3,gene4)
  
  temp <- c(length(gene1), length(gene2), length(intersect(gene1, gene2)),(length(intersect(gene1, gene2)) / (length(gene1) + length(gene2) - length(intersect(gene1, gene2)))),(length(intersect(gene1, gene2))/length(gene1)))
  
  jtk_rainandmeta2d_table <- rbind(jtk_rainandmeta2d_table, temp)
}
write.csv(jtk_rainandmeta2d_table,"E:/circidian_algorithm/result/mouse/jtkarser_jtkarser.csv")



tissuename <- c("adr","AOR","BFAT","bst","CER","HEA","HYP","KID","LIV","LUN","MUS","WAT")
jtk_rainandmeta2d_table <- data.frame(0,0,0,0,0)
for(j in 1:length(tissuename))
{
  file_name1 <- paste("E:/circidian_algorithm/result/mouse/2h/JTK/meta2d_",tissuename[j],".csv",sep="")
  file_name2 <- paste("E:/circidian_algorithm/result/mouse/2h/JTK/meta2d_",tissuename[j],".csv",sep="")
  file_name3 <- paste("E:/circidian_algorithm/result/mouse/4h/JTK/meta2d_",tissuename[j],".csv",sep="")
  file_name4 <- paste("E:/circidian_algorithm/result/mouse/4h/JTK/meta2d_",tissuename[j],".csv",sep="")
  
  data1 <- read.csv(file_name1)
  data2 <- read.csv(file_name2,header=TRUE)
  data3 <- read.csv(file_name3)
  data4 <- read.csv(file_name4,header=TRUE)
  
  gene1 <- data1[data1$ARS_BH.Q < 0.05, 1]
  gene2 <- data2[data2$meta2d_BH.Q < 0.05, 1]
  gene3 <- data3[data3$ARS_BH.Q < 0.05, 1]
  gene4 <- data4[data4$meta2d_BH.Q < 0.05, 1]
  
  gene1 <- union(gene1,gene2)
  gene2 <- union(gene3,gene4)
  
  temp <- c(length(gene1), length(gene2), length(intersect(gene1, gene2)),(length(intersect(gene1, gene2)) / (length(gene1) + length(gene2) - length(intersect(gene1, gene2)))),(length(intersect(gene1, gene2))/length(gene1)))
  
  jtk_rainandmeta2d_table <- rbind(jtk_rainandmeta2d_table, temp)
}
write.csv(jtk_rainandmeta2d_table,"E:/circidian_algorithm/result/mouse/meta2darser_meta2darser.csv")


tissuename <- c("adr","AOR","BFAT","bst","CER","HEA","HYP","KID","LIV","LUN","MUS","WAT")
jtk_rainandmeta2d_table <- data.frame(0,0,0,0,0)
for(j in 1:length(tissuename))
{
  file_name1 <- paste("E:/circidian_algorithm/result/mouse/2h/JTK/meta2d_",tissuename[j],".csv",sep="")
  file_name2 <- paste("E:/circidian_algorithm/result/mouse/2h/rain/",tissuename[j],".csv",sep="")
  file_name3 <- paste("E:/circidian_algorithm/result/mouse/4h/JTK/meta2d_",tissuename[j],".csv",sep="")
  file_name4 <- paste("E:/circidian_algorithm/result/mouse/4h/rain/",tissuename[j],".csv",sep="")
  
  data1 <- read.csv(file_name1)
  data2 <- read.csv(file_name2,header=TRUE)
  data3 <- read.csv(file_name3)
  data4 <- read.csv(file_name4,header=TRUE)
  
  gene1 <- data1[data1$ARS_BH.Q < 0.05, 1]
  gene2 <- data2[data2$bhq < 0.05, 2]
  gene3 <- data3[data3$ARS_BH.Q < 0.05, 1]
  gene4 <- data4[data4$bhq < 0.05, 2]
  
  gene1 <- union(gene1,gene2)
  gene2 <- union(gene3,gene4)
  
  temp <- c(length(gene1), length(gene2), length(intersect(gene1, gene2)),(length(intersect(gene1, gene2)) / (length(gene1) + length(gene2) - length(intersect(gene1, gene2)))),(length(intersect(gene1, gene2))/length(gene1)))
  
  jtk_rainandmeta2d_table <- rbind(jtk_rainandmeta2d_table, temp)
}
write.csv(jtk_rainandmeta2d_table,"E:/circidian_algorithm/result/mouse/raindarser_rainarser.csv")



tissuename <- c("adr","AOR","BFAT","bst","CER","HEA","HYP","KID","LIV","LUN","MUS","WAT")
jtk_rainandmeta2d_table <- data.frame(0,0,0,0,0)
for(j in 1:length(tissuename))
{
  file_name1 <- paste("E:/circidian_algorithm/result/archaea/meta2d/meta2d_hal_4.csv",sep="")
  file_name2 <- paste("E:/circidian_algorithm/result/archaea/rain/",'hal4',".csv",sep="")
  
  data1 <- read.csv(file_name1)
  data2 <- read.csv(file_name2,header=TRUE)

  
  gene1 <- data1[data1$meta2d_BH.Q < 0.05, 1]
  gene2 <- data2[data2$bhq < 0.05, 2]
  
  
  temp <- c(length(unique(gene1)), length(unique(gene2)), length(intersect(gene1, gene2)),length(union(gene1,gene2)))
  
  jtk_rainandmeta2d_table <- rbind(jtk_rainandmeta2d_table, temp)
}
write.csv(jtk_rainandmeta2d_table,"E:/circidian_algorithm/result/mouse/raindarser_rainarser.csv")




tissuename <- c("adr","AOR","BFAT","bst","CER","HEA","HYP","KID","LIV","LUN","MUS","WAT")
jtk_rainandmeta2d_table <- data.frame(0,0,0,0,0)
for(j in 1:length(tissuename))
{
  file_name1 <- paste("E:/circidian_algorithm/result/mouse/2h/JTK/meta2d_",tissuename[j],".csv",sep="")
  file_name2 <- paste("E:/circidian_algorithm/result/mouse/4h/biocycle/",tissuename[j],".tsv",sep="")
  file_name3 <- paste("E:/circidian_algorithm/result/mouse/4h/JTK/meta2d_",tissuename[j],".csv",sep="")
  data1 <- read.csv(file_name1)
  data2 <- read.table(file_name2,header=TRUE)
  data3 <- read.csv(file_name3)
  
  gene1 <- data1[data1$JTK_BH.Q < 0.05, 1]
  gene2 <- data2[data2$Q_VALUE < 0.05, 1]
  gene3 <- data3[data3$meta2d_BH.Q < 0.05, 1]
  gene2 <- union(gene2, gene3)
  temp <- c(length(gene1), length(gene2), length(intersect(gene1, gene2)),(length(intersect(gene1, gene2)) / (length(gene1) + length(gene2) - length(intersect(gene1, gene2)))),(length(intersect(gene1, gene2))/length(gene1)))
  
  jtk_rainandmeta2d_table <- rbind(jtk_rainandmeta2d_table, temp)
}
write.csv(jtk_rainandmeta2d_table,"E:/circidian_algorithm/result/mouse/meta2dbiocycleJTK2h_meta2dbiocycleJTK2h.csv")



tissuename <- c("adr","AOR","BFAT","bst","CER","HEA","HYP","KID","LIV","LUN","MUS","WAT")
jtk_rainandmeta2d_table <- data.frame(0,0,0,0,0)
for(j in 1:length(tissuename))
{
  file_name1 <- paste("E:/circidian_algorithm/result/mouse/2h/JTK/meta2d_",tissuename[j],".csv",sep="")
  file_name2 <- paste("E:/circidian_algorithm/result/mouse/4h/JTK/meta2d_",tissuename[j],".csv",sep="")
  file_name3 <- paste("E:/circidian_algorithm/result/mouse/4h/JTK/meta2d_",tissuename[j],".csv",sep="")
  data1 <- read.csv(file_name1)
  data2 <- read.csv(file_name2,header=TRUE)
  data3 <- read.csv(file_name3)
  JTK_BH.Q
  gene1 <- data1[data1$JTK_BH.Q < 0.05, 1]
  gene2 <- data2[data2$ARS_BH.Q < 0.05, 1]
  gene3 <- data3[data3$meta2d_BH.Q < 0.05, 1]
  gene2 <- union(gene2, gene3)
  temp <- c(length(gene1), length(gene2), length(intersect(gene1, gene2)),(length(intersect(gene1, gene2)) / (length(gene1) + length(gene2) - length(intersect(gene1, gene2)))),(length(intersect(gene1, gene2))/length(gene1)))
  
  jtk_rainandmeta2d_table <- rbind(jtk_rainandmeta2d_table, temp)
}
write.csv(jtk_rainandmeta2d_table,"E:/circidian_algorithm/result/mouse/meta2dbiocycleJTK2h_meta2dbiocycleJTK2h.csv")












###########################################mouse another 4h data
tissuename <- c("adr","AOR","BFAT","bst","CER","HEA","HYP","KID","LIV","LUN","MUS","WAT")
for(j in 1:length(tissuename))
{
  file_name <- paste("E:/GO_layers/final_photo/exp_amp_data/mouse_2h/",tissuename[j],".csv",sep="")
  data <- read.csv(file_name,check.names = FALSE)
  new_data <- data[,-c(2,4,6,8,10,12,14,16,18,20,22,24)]
  out_name <- paste("E:/GO_layers/final_photo/exp_amp_data/mouse_4h_2/",tissuename[j],".csv",sep="")
  write.csv(new_data, out_name, row.names = FALSE)
}
  
library(MetaCycle)
tissuename <- c("adr","AOR","BFAT","bst","CER","HEA","HYP","KID","LIV","LUN","MUS","WAT")
for(j in 1:length(tissuename))
{
  file_name <- paste("E:/circidian_algorithm/data/mouse/mouse_4h_2/",tissuename[j],".csv",sep="")
  data <- read.csv(file_name,check.names = FALSE)
  meta2d(infile=file_name, filestyle="csv", outdir="E:/circidian_algorithm/result/mouse/4h_2",timepoints="Line1",outRawData=TRUE,minper=20,maxper=28, cycMethod = c("ARS", "JTK", "LS"))
}

######################################################RAIN
library(rain)

set.seed(123)

tissuename <- c("adr","AOR","BFAT","bst","CER","HEA","HYP","KID","LIV","LUN","MUS","WAT")
for(j in 1:length(tissuename))
{
  file_name <- paste("E:/circidian_algorithm/data/mouse/mouse_4h_2/",tissuename[j],".csv",sep="")
  data <- read.csv(file_name)
  input_data <- t(data[,2:13])
  rainresult <- rain(input_data, period=24, deltat=4, peak.border=c(0.3,0.7),verbose=FALSE) 
  #如果有重复nr.series=n
  #results <- rain(t(menetRNASeqMouseLiver), deltat=4, period=24, nr.series=2,peak.border=c(0.3, 0.7), verbose=FALSE)
  pval <- rainresult$pVal
  bhq <- p.adjust(pval, method="fdr")
  temp <- data.frame(data[,1], pval, bhq, rainresult$phase)
  file_out_name <- paste("E:/circidian_algorithm/result/mouse/4h_2/rain/",tissuename[j],".csv",sep="")
  print(tissuename[j])
  print(sum(pval < 0.05))
  print(sum(bhq < 0.05))
  write.csv(temp, file_out_name)
}



#######################################################FIsher G test
library(GeneCycle)


tissuename <- c("adr","AOR","BFAT","bst","CER","HEA","HYP","KID","LIV","LUN","MUS","WAT")
for(j in 1:length(tissuename))
{
  file_name <- paste("E:/circidian_algorithm/data/mouse/mouse_4h_2/",tissuename[j],".csv",sep="")
  data <- read.csv(file_name)
  input_data <- t(data[,2:13])
  # p-values from Fisher's g test
  pval.caulobacter <- fisher.g.test(input_data)
  
  fdr.out <- fdrtool(pval.caulobacter, statistic="pvalue")
  # how many significant?
  temp <- data.frame(data[,1], fdr.out$pval, fdr.out$qval)
  print(tissuename[j])
  print(sum(fdr.out$pval < 0.05))
  print(sum(fdr.out$qval < 0.05))
  file_out_name <- paste("E:/circidian_algorithm/result/mouse/4h_2/gtest/",tissuename[j],".csv",sep="")
  write.csv(temp, file_out_name)
  # sum(fdr.out$qval < 0.05) # tail area-based Fdr
  # sum(fdr.out$lfdr < 0.2) # density-based local fdr
  # sum(fdr.out$lfdr < 0.05) # density-based local fdr
}




















#meta2d
tissuename <- c("adr","AOR","BFAT","bst","CER","HEA","HYP","KID","LIV","LUN","MUS","WAT")
jtk_table <- data.frame(0,0,0,0,0)
ls_table <- data.frame(0,0,0,0,0)
arser_table <- data.frame(0,0,0,0,0)
meta2d_table <- data.frame(0,0,0,0,0)
for(j in 1:length(tissuename))
{
  file_name1 <- paste("E:/circidian_algorithm/result/mouse/2h/JTK/meta2d_",tissuename[j],".csv",sep="")
  file_name2 <- paste("E:/circidian_algorithm/result/mouse/4h_2/JTK/meta2d_",tissuename[j],".csv",sep="")
  data1 <- read.csv(file_name1)
  data2 <- read.csv(file_name2)
  
  gene1 <- data1[data1$JTK_BH.Q < 0.05, 1]
  gene2 <- data2[data2$JTK_BH.Q < 0.05, 1]
  temp <- c(length(gene1), length(gene2), length(intersect(gene1, gene2)),(length(intersect(gene1, gene2)) / (length(gene1) + length(gene2) - length(intersect(gene1, gene2)))),(length(intersect(gene1, gene2))/length(gene1)))
  jtk_table <- rbind(jtk_table, temp)
  
  gene1 <- data1[data1$LS_BH.Q < 0.05, 1]
  gene2 <- data2[data2$LS_BH.Q < 0.05, 1]
  temp <- c(length(gene1), length(gene2), length(intersect(gene1, gene2)),(length(intersect(gene1, gene2)) / (length(gene1) + length(gene2) - length(intersect(gene1, gene2)))),(length(intersect(gene1, gene2))/length(gene1) ))
  ls_table <- rbind(ls_table, temp)
  
  gene1 <- data1[data1$ARS_BH.Q < 0.05, 1]
  gene2 <- data2[data2$ARS_BH.Q < 0.05, 1]
  temp <- c(length(gene1), length(gene2), length(intersect(gene1, gene2)),(length(intersect(gene1, gene2)) / (length(gene1) + length(gene2) - length(intersect(gene1, gene2)))),(length(intersect(gene1, gene2))/length(gene1) ))
  arser_table <- rbind(arser_table, temp)
  
  gene1 <- data1[data1$meta2d_BH.Q < 0.05, 1]
  gene2 <- data2[data2$meta2d_BH.Q < 0.05, 1]
  temp <- c(length(gene1), length(gene2), length(intersect(gene1, gene2)),(length(intersect(gene1, gene2)) / (length(gene1) + length(gene2) - length(intersect(gene1, gene2)))),(length(intersect(gene1, gene2))/length(gene1) ))
  meta2d_table <- rbind(meta2d_table, temp)
}
write.csv(jtk_table,"E:/circidian_algorithm/result/mouse/4h_2/result/jtk_table.csv")
write.csv(ls_table,"E:/circidian_algorithm/result/mouse/4h_2/result/ls_table.csv")
write.csv(arser_table,"E:/circidian_algorithm/result/mouse/4h_2/result/arser_table.csv")
write.csv(meta2d_table,"E:/circidian_algorithm/result/mouse/4h_2/result/meta2d_table.csv")



#biocycle
tissuename <- c("adr","AOR","BFAT","bst","CER","HEA","HYP","KID","LIV","LUN","MUS","WAT")
biocycle_table <- data.frame(0,0,0,0,0)
for(j in 1:length(tissuename))
{
  file_name1 <- paste("E:/circidian_algorithm/result/mouse/2h/biocycle/",tissuename[j],".tsv",sep="")
  file_name2 <- paste("E:/circidian_algorithm/result/mouse/4h/biocycle/",tissuename[j],".tsv",sep="")
  data1 <- read.table(file_name1, header=TRUE)
  data2 <- read.table(file_name2, header=TRUE)
  gene1 <- data1[data1$Q_VALUE<0.05,1]
  gene2 <- data2[data2$Q_VALUE<0.05,1]
  temp <- c(length(gene1), length(gene2), length(intersect(gene1, gene2)),(length(intersect(gene1, gene2)) / (length(gene1) + length(gene2) - length(intersect(gene1, gene2)))),(length(intersect(gene1, gene2))/length(gene1) ))
  biocycle_table <- rbind(biocycle_table, temp)
}
write.csv(biocycle_table, "E:/circidian_algorithm/result/mouse/biocycle_table.csv")

#rain
tissuename <- c("adr","AOR","BFAT","bst","CER","HEA","HYP","KID","LIV","LUN","MUS","WAT")
rain_table <- data.frame(0,0,0,0,0)
for(j in 1:length(tissuename))
{
  file_name1 <- paste("E:/circidian_algorithm/result/mouse/2h/rain/",tissuename[j],".csv",sep="")
  file_name2 <- paste("E:/circidian_algorithm/result/mouse/4h_2/rain/",tissuename[j],".csv",sep="")
  data1 <- read.csv(file_name1, header=TRUE)
  data2 <- read.csv(file_name2, header=TRUE)
  gene1 <- data1[data1$bhq<0.05,2]
  gene2 <- data2[data2$bhq<0.05,2]
  temp <- c(length(gene1), length(gene2), length(intersect(gene1, gene2)),(length(intersect(gene1, gene2)) / (length(gene1) + length(gene2) - length(intersect(gene1, gene2)))),(length(intersect(gene1, gene2))/length(gene1) ))
  rain_table <- rbind(rain_table, temp)
}
write.csv(rain_table, "E:/circidian_algorithm/result/mouse/4h_2/table/rain_table.csv")



#eJTK
tissuename <- c("adr","AOR","BFAT","bst","CER","HEA","HYP","KID","LIV","LUN","MUS","WAT")
ejtk_table <- data.frame(0,0,0,0,0)
for(j in 1:length(tissuename))
{
  file_name1 <- paste("E:/circidian_algorithm/result/mouse/2h/eJTK/",tissuename[j],"_cos24_ph00-22_by2_a02-22_by2_OTHERTEXT_jtkout_GammaP.txt",sep="")
  file_name2 <- paste("E:/circidian_algorithm/result/mouse/4h_2/eJTK/",tissuename[j],"_cos24_ph00-20_by4_a04-20_by4_OTHERTEXT_jtkout_GammaP.txt",sep="")
  data1 <- read.table(file_name1, header=TRUE)
  data2 <- read.table(file_name2, header=TRUE)
  gene1 <- data1[data1$GammaBH<0.05,1]
  gene2 <- data2[data2$GammaBH<0.05,1]
  temp <- c(length(gene1), length(gene2), length(intersect(gene1, gene2)),(length(intersect(gene1, gene2)) / (length(gene1) + length(gene2) - length(intersect(gene1, gene2)))),(length(intersect(gene1, gene2))/length(gene1) ))
  ejtk_table <- rbind(ejtk_table, temp)
}
write.csv(ejtk_table, "E:/circidian_algorithm/result/mouse/4h_2/table/eJTK_table.csv")



#gtest
tissuename <- c("adr","AOR","BFAT","bst","CER","HEA","HYP","KID","LIV","LUN","MUS","WAT")
ejtk_table <- data.frame(0,0,0,0,0)
for(j in 1:length(tissuename))
{
  file_name1 <- paste("E:/circidian_algorithm/result/mouse/2h/gtest/",tissuename[j],".csv",sep="")
  file_name2 <- paste("E:/circidian_algorithm/result/mouse/4h_2/gtest/",tissuename[j],".csv",sep="")
  data1 <- read.csv(file_name1)
  data2 <- read.csv(file_name2)
  gene1 <- data1[data1$fdr.out.qval<0.05,2]
  gene2 <- data2[data2$fdr.out.qval<0.05,2]
  temp <- c(length(gene1), length(gene2), length(intersect(gene1, gene2)),(length(intersect(gene1, gene2)) / (length(gene1) + length(gene2) - length(intersect(gene1, gene2)))),(length(intersect(gene1, gene2))/length(gene1) ))
  ejtk_table <- rbind(ejtk_table, temp)
}
write.csv(ejtk_table, "E:/circidian_algorithm/result/mouse/4h_2/table/gtest_table.csv")


#meta2d & rain
tissuename <- c("adr","AOR","BFAT","bst","CER","HEA","HYP","KID","LIV","LUN","MUS","WAT")
meta2d_table <- data.frame(0,0,0,0,0)
for(j in 1:length(tissuename))
{
  file_name1 <- paste("E:/circidian_algorithm/result/mouse/4h_2/JTK/meta2d_",tissuename[j],".csv",sep="")
  file_name2 <- paste("E:/circidian_algorithm/result/mouse/4h_2/rain/",tissuename[j],".csv",sep="")
  data1 <- read.csv(file_name1)
  data2 <- read.csv(file_name2)
  
  gene1 <- data1[data1$meta2d_BH.Q < 0.05, 1]
  gene2 <- data2[data2$bhq < 0.05, 2]
  temp <- c(length(gene1), length(gene2), length(intersect(gene1, gene2)),(length(intersect(gene1, gene2)) / (length(gene1) + length(gene2) - length(intersect(gene1, gene2)))),(length(intersect(gene1, gene2))/length(gene1)))
  
  meta2d_table <- rbind(meta2d_table, temp)
}
write.csv(meta2d_table,"E:/circidian_algorithm/result/mouse/4h_2/table/meta2d_rain_4h.csv")





#JTK & rain &meta2d

tissuename <- c("adr","AOR","BFAT","bst","CER","HEA","HYP","KID","LIV","LUN","MUS","WAT")
jtk_rain4_table <- data.frame(0,0,0,0,0)
for(j in 1:length(tissuename))
{
  file_name1 <- paste("E:/circidian_algorithm/result/mouse/2h/JTK/meta2d_",tissuename[j],".csv",sep="")
  file_name2 <- paste("E:/circidian_algorithm/result/mouse/4h_2/rain/",tissuename[j],".csv",sep="")
  data1 <- read.csv(file_name1)
  data2 <- read.csv(file_name2)
  
  gene1 <- data1[data1$JTK_BH.Q < 0.05, 1]
  gene2 <- data2[data2$bhq < 0.05, 2]
  temp <- c(length(gene1), length(gene2), length(intersect(gene1, gene2)),(length(intersect(gene1, gene2)) / (length(gene1) + length(gene2) - length(intersect(gene1, gene2)))),(length(intersect(gene1, gene2))/length(gene1)))
  
  jtk_rain4_table <- rbind(jtk_rain4_table, temp)
}
write.csv(jtk_rain4_table,"E:/circidian_algorithm/result/mouse/4h_2/table/jtk_rain4.csv")

tissuename <- c("adr","AOR","BFAT","bst","CER","HEA","HYP","KID","LIV","LUN","MUS","WAT")
jtk_rainandmeta2d_table <- data.frame(0,0,0,0,0)
for(j in 1:length(tissuename))
{
  file_name1 <- paste("E:/circidian_algorithm/result/mouse/2h/JTK/meta2d_",tissuename[j],".csv",sep="")
  file_name2 <- paste("E:/circidian_algorithm/result/mouse/4h_2/rain/",tissuename[j],".csv",sep="")
  file_name3 <- paste("E:/circidian_algorithm/result/mouse/4h_2/JTK/meta2d_",tissuename[j],".csv",sep="")
  data1 <- read.csv(file_name1)
  data2 <- read.csv(file_name2)
  data3 <- read.csv(file_name3)
  
  gene1 <- data1[data1$JTK_BH.Q < 0.05, 1]
  gene2 <- data2[data2$bhq < 0.05, 2]
  gene3 <- data3[data3$meta2d_BH.Q < 0.05, 1]
  gene2 <- union(gene2, gene3)
  temp <- c(length(gene1), length(gene2), length(intersect(gene1, gene2)),(length(intersect(gene1, gene2)) / (length(gene1) + length(gene2) - length(intersect(gene1, gene2)))),(length(intersect(gene1, gene2))/length(gene1)))
  
  jtk_rainandmeta2d_table <- rbind(jtk_rainandmeta2d_table, temp)
}
write.csv(jtk_rainandmeta2d_table,"E:/circidian_algorithm/result/mouse/4h_2/table/jtk_rainandmeta2d.csv")




tissuename <- c("adr","AOR","BFAT","bst","CER","HEA","HYP","KID","LIV","LUN","MUS","WAT")
jtk_rainandmeta2d_table <- data.frame(0,0,0,0,0)
for(j in 1:length(tissuename))
{
  file_name1 <- paste("E:/circidian_algorithm/result/mouse/2h/JTK/meta2d_",tissuename[j],".csv",sep="")
  file_name2 <- paste("E:/circidian_algorithm/result/mouse/4h_2/JTK/meta2d_",tissuename[j],".csv",sep="")
  data1 <- read.csv(file_name1)
  data2 <- read.csv(file_name2)
  
  gene1 <- data1[data1$JTK_BH.Q < 0.05, 1]
  gene2 <- data2[data2$meta2d_BH.Q < 0.05, 1]
  temp <- c(length(gene1), length(gene2), length(intersect(gene1, gene2)),(length(intersect(gene1, gene2)) / (length(gene1) + length(gene2) - length(intersect(gene1, gene2)))),(length(intersect(gene1, gene2))/length(gene1)))
  
  jtk_rainandmeta2d_table <- rbind(jtk_rainandmeta2d_table, temp)
}
write.csv(jtk_rainandmeta2d_table,"E:/circidian_algorithm/result/mouse/4h_2/table/jtk_meta2d4h.csv")






###############################intersect and union
tissuename <- c("adr","AOR","BFAT","bst","CER","HEA","HYP","KID","LIV","LUN","MUS","WAT")
jtk_rainandmeta2d_table <- data.frame(0,0,0,0,0)
for(j in 1:length(tissuename))
{
  file_name1 <- paste("E:/circidian_algorithm/result/mouse/2h/JTK/meta2d_",tissuename[j],".csv",sep="")
  file_name2 <- paste("E:/circidian_algorithm/result/mouse/4h_2/JTK/meta2d_",tissuename[j],".csv",sep="")
  file_name3 <- paste("E:/circidian_algorithm/result/mouse/4h_2/rain/",tissuename[j],".csv",sep="")
  file_name4 <- paste("E:/circidian_algorithm/result/mouse/4h_2/biocycle/",tissuename[j],".tsv",sep="")
  #file_name5 <- paste("E:/circidian_algorithm/result/mouse/4h_2/eJTK/",tissuename[j],"_cos24_ph00-20_by4_a04-20_by4_OTHERTEXT_jtkout_GammaP.txt",sep="")
  data1 <- read.csv(file_name1)
  data2 <- read.csv(file_name2)
  data3 <- read.csv(file_name3)
  data4 <- read.table(file_name4,header=TRUE)
  #data5 <- read.table(file_name5,header=TRUE)
  
  fdr <- p.adjust(data4$P_VALUE, method="fdr")
  
  gene1 <- data1[data1$JTK_BH.Q < 0.05, 1]
  gene2 <- data2[data2$meta2d_BH.Q < 0.05, 1]
  gene3 <- data3[data3$bhq < 0.05, 2]
  gene4 <- data4[fdr < 0.05, 1]
  #gene5 <- data5[data5$GammaBH < 0.05, 1]
  inter1 <- intersect(gene2, gene3)
  inter2 <- intersect(gene2, gene4)
 # inter3 <- intersect(gene2, gene5)
  inter4 <- intersect(gene3, gene4)
  # inter5 <- intersect(gene3, gene5)
  # inter6 <- intersect(gene4, gene5)
  if(length(inter1)==0)
  {
    inter1=0
  }
  if(length(inter2)==0)
  {
    inter2=0
  }
  if(length(inter4)==0)
  {
    inter4=0
  }
  union_gene <- unique(inter1,inter2,inter4)
  length(intersect(gene1,union_gene))

  gene2 <- intersect(gene1,union_gene)

  temp <- c(length(gene1), length(gene2), length(intersect(gene1, gene2)),(length(intersect(gene1, gene2)) / (length(gene1) + length(gene2) - length(intersect(gene1, gene2)))),(length(intersect(gene1, gene2))/length(gene1)))
  
  jtk_rainandmeta2d_table <- rbind(jtk_rainandmeta2d_table, temp)
}
write.csv(jtk_rainandmeta2d_table,"E:/circidian_algorithm/result/mouse/4h_2/table/union.csv")




tissuename <- c("adr","AOR","BFAT","bst","CER","HEA","HYP","KID","LIV","LUN","MUS","WAT")
jtk_rainandmeta2d_table <- data.frame(0,0,0,0,0)
for(j in 1:length(tissuename))
{
  file_name1 <- paste("E:/circidian_algorithm/result/mouse/2h/JTK/meta2d_",tissuename[j],".csv",sep="")
  file_name2 <- paste("E:/circidian_algorithm/result/mouse/4h_2/rain/",tissuename[j],".csv",sep="")
  file_name3 <- paste("E:/circidian_algorithm/result/mouse/4h_2/JTK/meta2d_",tissuename[j],".csv",sep="")
  data1 <- read.csv(file_name1)
  data2 <- read.csv(file_name2)
  data3 <- read.csv(file_name3)
  
  gene1 <- data1[data1$ARS_BH.Q < 0.05, 1]
  gene2 <- data2[data2$bhq < 0.05, 2]
  gene3 <- data3[data3$meta2d_BH.Q < 0.05, 1]
  gene2 <- union(gene2, gene3)
  temp <- c(length(gene1), length(gene2), length(intersect(gene1, gene2)),(length(intersect(gene1, gene2)) / (length(gene1) + length(gene2) - length(intersect(gene1, gene2)))),(length(intersect(gene1, gene2))/length(gene1)))
  
  jtk_rainandmeta2d_table <- rbind(jtk_rainandmeta2d_table, temp)
}
write.csv(jtk_rainandmeta2d_table,"E:/circidian_algorithm/result/mouse/4h_2/table/arser_rainandmeta2d.csv")









tissuename <- c("adr","AOR","BFAT","bst","CER","HEA","HYP","KID","LIV","LUN","MUS","WAT")
jtk_rainandmeta2d_table <- data.frame(0,0,0,0,0)
for(j in 1:length(tissuename))
{
  file_name1 <- paste("E:/circidian_algorithm/result/mouse/2h/JTK/meta2d_",tissuename[j],".csv",sep="")
  file_name2 <- paste("E:/circidian_algorithm/result/mouse/4h_2/JTK/meta2d_",tissuename[j],".csv",sep="")
  data1 <- read.csv(file_name1)
  data2 <- read.csv(file_name2)
  
  gene1 <- data1[data1$ARS_BH.Q < 0.05, 1]
  gene2 <- data2[data2$meta2d_BH.Q < 0.05, 1]
  temp <- c(length(gene1), length(gene2), length(intersect(gene1, gene2)),(length(intersect(gene1, gene2)) / (length(gene1) + length(gene2) - length(intersect(gene1, gene2)))),(length(intersect(gene1, gene2))/length(gene1)))
  
  jtk_rainandmeta2d_table <- rbind(jtk_rainandmeta2d_table, temp)
}
write.csv(jtk_rainandmeta2d_table,"E:/circidian_algorithm/result/mouse/4h_2/table/arser_meta2d4h.csv")

tissuename <- c("adr","AOR","BFAT","bst","CER","HEA","HYP","KID","LIV","LUN","MUS","WAT")
jtk_rain4_table <- data.frame(0,0,0,0,0)
for(j in 1:length(tissuename))
{
  file_name1 <- paste("E:/circidian_algorithm/result/mouse/2h/JTK/meta2d_",tissuename[j],".csv",sep="")
  file_name2 <- paste("E:/circidian_algorithm/result/mouse/4h_2/rain/",tissuename[j],".csv",sep="")
  data1 <- read.csv(file_name1)
  data2 <- read.csv(file_name2)
  
  gene1 <- data1[data1$ARS_BH.Q < 0.05, 1]
  gene2 <- data2[data2$bhq < 0.05, 2]
  temp <- c(length(gene1), length(gene2), length(intersect(gene1, gene2)),(length(intersect(gene1, gene2)) / (length(gene1) + length(gene2) - length(intersect(gene1, gene2)))),(length(intersect(gene1, gene2))/length(gene1)))
  
  jtk_rain4_table <- rbind(jtk_rain4_table, temp)
}
write.csv(jtk_rain4_table,"E:/circidian_algorithm/result/mouse/4h_2/table/arser_rain4.csv")











#meta2d
tissuename <- c("adr","AOR","BFAT","bst","CER","HEA","HYP","KID","LIV","LUN","MUS","WAT")
jtk_table <- data.frame(0,0,0,0,0)
ls_table <- data.frame(0,0,0,0,0)
arser_table <- data.frame(0,0,0,0,0)
meta2d_table <- data.frame(0,0,0,0,0)
for(j in 1:length(tissuename))
{
  file_name1 <- paste("E:/circidian_algorithm/result/mouse/4h/JTK/meta2d_",tissuename[j],".csv",sep="")
  file_name2 <- paste("E:/circidian_algorithm/result/mouse/4h_2/JTK/meta2d_",tissuename[j],".csv",sep="")
  data1 <- read.csv(file_name1)
  data2 <- read.csv(file_name2)
  
  gene1 <- data1[data1$JTK_BH.Q < 0.05, 1]
  gene2 <- data2[data2$JTK_BH.Q < 0.05, 1]
  temp <- c(length(gene1), length(gene2), length(intersect(gene1, gene2)),(length(intersect(gene1, gene2)) / (length(gene1) + length(gene2) - length(intersect(gene1, gene2)))),(length(intersect(gene1, gene2))/length(gene1)))
  jtk_table <- rbind(jtk_table, temp)
  
  gene1 <- data1[data1$LS_BH.Q < 0.05, 1]
  gene2 <- data2[data2$LS_BH.Q < 0.05, 1]
  temp <- c(length(gene1), length(gene2), length(intersect(gene1, gene2)),(length(intersect(gene1, gene2)) / (length(gene1) + length(gene2) - length(intersect(gene1, gene2)))),(length(intersect(gene1, gene2))/length(gene1) ))
  ls_table <- rbind(ls_table, temp)
  
  gene1 <- data1[data1$ARS_BH.Q < 0.05, 1]
  gene2 <- data2[data2$ARS_BH.Q < 0.05, 1]
  temp <- c(length(gene1), length(gene2), length(intersect(gene1, gene2)),(length(intersect(gene1, gene2)) / (length(gene1) + length(gene2) - length(intersect(gene1, gene2)))),(length(intersect(gene1, gene2))/length(gene1) ))
  arser_table <- rbind(arser_table, temp)
  
  gene1 <- data1[data1$meta2d_BH.Q < 0.05, 1]
  gene2 <- data2[data2$meta2d_BH.Q < 0.05, 1]
  temp <- c(length(gene1), length(gene2), length(intersect(gene1, gene2)),(length(intersect(gene1, gene2)) / (length(gene1) + length(gene2) - length(intersect(gene1, gene2)))),(length(intersect(gene1, gene2))/length(gene1) ))
  meta2d_table <- rbind(meta2d_table, temp)
}
write.csv(jtk_table,"E:/circidian_algorithm/result/mouse/4h_2/result/4h_jtk_table.csv")
write.csv(ls_table,"E:/circidian_algorithm/result/mouse/4h_2/result/4h_ls_table.csv")
write.csv(arser_table,"E:/circidian_algorithm/result/mouse/4h_2/result/4h_arser_table.csv")
write.csv(meta2d_table,"E:/circidian_algorithm/result/mouse/4h_2/result/4h_meta2d_table.csv")



#biocycle
tissuename <- c("adr","AOR","BFAT","bst","CER","HEA","HYP","KID","LIV","LUN","MUS","WAT")
biocycle_table <- data.frame(0,0,0,0,0)
for(j in 1:length(tissuename))
{
  file_name1 <- paste("E:/circidian_algorithm/result/mouse/4h/biocycle/",tissuename[j],".tsv",sep="")
  file_name2 <- paste("E:/circidian_algorithm/result/mouse/4h_2/biocycle/",tissuename[j],".tsv",sep="")
  data1 <- read.table(file_name1, header=TRUE)
  data2 <- read.table(file_name2, header=TRUE)
  gene1 <- data1[data1$Q_VALUE<0.05,1]
  gene2 <- data2[data2$Q_VALUE<0.05,1]
  temp <- c(length(gene1), length(gene2), length(intersect(gene1, gene2)),(length(intersect(gene1, gene2)) / (length(gene1) + length(gene2) - length(intersect(gene1, gene2)))),(length(intersect(gene1, gene2))/length(gene1) ))
  biocycle_table <- rbind(biocycle_table, temp)
}
write.csv(biocycle_table, "E:/circidian_algorithm/result/mouse/4h_biocycle_table.csv")

#rain
tissuename <- c("adr","AOR","BFAT","bst","CER","HEA","HYP","KID","LIV","LUN","MUS","WAT")
rain_table <- data.frame(0,0,0,0,0)
for(j in 1:length(tissuename))
{
  file_name1 <- paste("E:/circidian_algorithm/result/mouse/4h/rain/",tissuename[j],".csv",sep="")
  file_name2 <- paste("E:/circidian_algorithm/result/mouse/4h_2/rain/",tissuename[j],".csv",sep="")
  data1 <- read.csv(file_name1, header=TRUE)
  data2 <- read.csv(file_name2, header=TRUE)
  gene1 <- data1[data1$bhq<0.05,2]
  gene2 <- data2[data2$bhq<0.05,2]
  temp <- c(length(gene1), length(gene2), length(intersect(gene1, gene2)),(length(intersect(gene1, gene2)) / (length(gene1) + length(gene2) - length(intersect(gene1, gene2)))),(length(intersect(gene1, gene2))/length(gene1) ))
  rain_table <- rbind(rain_table, temp)
}
write.csv(rain_table, "E:/circidian_algorithm/result/mouse/4h_2/table/4h_rain_table.csv")



#eJTK
tissuename <- c("adr","AOR","BFAT","bst","CER","HEA","HYP","KID","LIV","LUN","MUS","WAT")
ejtk_table <- data.frame(0,0,0,0,0)
for(j in 1:length(tissuename))
{
  file_name1 <- paste("E:/circidian_algorithm/result/mouse/4h/eJTK/",tissuename[j],"_cos24_ph00-22_by2_a02-22_by2_OTHERTEXT_jtkout_GammaP.txt",sep="")
  file_name2 <- paste("E:/circidian_algorithm/result/mouse/4h_2/eJTK/",tissuename[j],"_cos24_ph00-20_by4_a04-20_by4_OTHERTEXT_jtkout_GammaP.txt",sep="")
  data1 <- read.table(file_name1, header=TRUE)
  data2 <- read.table(file_name2, header=TRUE)
  gene1 <- data1[data1$GammaBH<0.05,1]
  gene2 <- data2[data2$GammaBH<0.05,1]
  temp <- c(length(gene1), length(gene2), length(intersect(gene1, gene2)),(length(intersect(gene1, gene2)) / (length(gene1) + length(gene2) - length(intersect(gene1, gene2)))),(length(intersect(gene1, gene2))/length(gene1) ))
  ejtk_table <- rbind(ejtk_table, temp)
}
write.csv(ejtk_table, "E:/circidian_algorithm/result/mouse/4h_2/table/4h_eJTK_table.csv")



#gtest
tissuename <- c("adr","AOR","BFAT","bst","CER","HEA","HYP","KID","LIV","LUN","MUS","WAT")
ejtk_table <- data.frame(0,0,0,0,0)
for(j in 1:length(tissuename))
{
  file_name1 <- paste("E:/circidian_algorithm/result/mouse/4h/gtest/",tissuename[j],".csv",sep="")
  file_name2 <- paste("E:/circidian_algorithm/result/mouse/4h_2/gtest/",tissuename[j],".csv",sep="")
  data1 <- read.csv(file_name1)
  data2 <- read.csv(file_name2)
  gene1 <- data1[data1$fdr.out.qval<0.05,2]
  gene2 <- data2[data2$fdr.out.qval<0.05,2]
  temp <- c(length(gene1), length(gene2), length(intersect(gene1, gene2)),(length(intersect(gene1, gene2)) / (length(gene1) + length(gene2) - length(intersect(gene1, gene2)))),(length(intersect(gene1, gene2))/length(gene1) ))
  ejtk_table <- rbind(ejtk_table, temp)
}
write.csv(ejtk_table, "E:/circidian_algorithm/result/mouse/4h_2/table/4h_gtest_table.csv")


#meta2d & rain
tissuename <- c("adr","AOR","BFAT","bst","CER","HEA","HYP","KID","LIV","LUN","MUS","WAT")
meta2d_table <- data.frame(0,0,0,0,0)
for(j in 1:length(tissuename))
{
  file_name1 <- paste("E:/circidian_algorithm/result/mouse/4h_2/JTK/meta2d_",tissuename[j],".csv",sep="")
  file_name2 <- paste("E:/circidian_algorithm/result/mouse/4h_2/rain/",tissuename[j],".csv",sep="")
  data1 <- read.csv(file_name1)
  data2 <- read.csv(file_name2)
  
  gene1 <- data1[data1$meta2d_BH.Q < 0.05, 1]
  gene2 <- data2[data2$bhq < 0.05, 2]
  temp <- c(length(gene1), length(gene2), length(intersect(gene1, gene2)),(length(intersect(gene1, gene2)) / (length(gene1) + length(gene2) - length(intersect(gene1, gene2)))),(length(intersect(gene1, gene2))/length(gene1)))
  
  meta2d_table <- rbind(meta2d_table, temp)
}
write.csv(meta2d_table,"E:/circidian_algorithm/result/mouse/4h_2/table/meta2d_rain_4h.csv")











######################cosine
library(lme4)



######################sin
set.seed(1)  
data <- read.csv("E:/circidian_algorithm/data/mouse/mouse_4h/adr.csv",check.names = FALSE)
t <- colnames(data)
t <- as.numeric(t[-1])
data <- as.matrix(data)
num <- 0
for(i in 1:nrow(data))
{
  y <- data[1,]
  y <- as.vector(y[-1])
  res <- runif(length(t),min=-50,max=50)
  y <- y+ res
  m <- nls(y ~ a*sin(b*t+c) + mean(y),start=list(a=1, b=24,c=1))
# 计算模型的拟合优度
  if(cor(y, predict(m)) > 0.75)
  {
    num <- num + 1
  }
}


x1 <- seq(0, 50, 1)
x2 <- seq(0, 100, 2)
x <- data.frame(x1,x2)
f=function(x1, x2, a, b) {a+x1+x2^b};
result=nls(x$y~f(x$x1, x$x2, a, b), data=x, start=list(a=1, b=2));
plot(x, y) 
lines(x, predict(m), lty = 2, col = "red", lwd = 3)







############
library(nls2)
set.seed(1)  
data <- read.csv("E:/circidian_algorithm/data/mouse/mouse_4h/adr.csv",check.names = FALSE)
t <- colnames(data)
t <- as.numeric(t[-1])
data <- as.matrix(data)
num <- 0
for(i in 1:nrow(data))
{
  y <- data[1,]
  y <- as.vector(y[-1])
  res <- runif(length(t),min=-50,max=50)
  y <- y+ res
  m <- nls(y ~ a*sin(b*t+c) + mean(y),start=list(a=mean(y), b=24,c=0.5))
  # 计算模型的拟合优度
  if(cor(y, predict(m)) > 0.75)
  {
    num <- num + 1
  }
}





library(nls2)
set.seed(1)
library(lme4)
library(nlme)
tissuename <- c("adr","AOR","BFAT","bst","CER","HEA","HYP","KID","LIV","LUN","MUS","WAT")

for(j in 1:length(tissuename))
{
  file_name1 <- paste("E:/circidian_algorithm/data/mouse/mouse_4h/",tissuename[j],".csv",sep="")
  pvalue_list <- c(rep(0,35556))
  cor_list <- c(rep(0,35556))
  data <- read.csv(file_name1,check.names = FALSE)
  
  
  t <- c(18,22,26,30,34,38,42,46,50,54,58,62)
  data <- as.matrix(data[,-1])
  data <- apply(data,2,as.numeric)
  data_new <- scale(data)
  data <- data_new
  num <- 0
  
  for(i in 1:nrow(data))
  {
    y <- data[i,]
    new_data <- data.frame(y,t)
    #m <- nlme(y ~ a*sin(b*t+c) + mean(y), data=new_data, fixed = a + b + c ~ 1, random = Asym ~ 1, start=c(a=mean(y),b=24,c=6))
    
    st1 <- expand.grid(a = 2*mean(y),b = 0.26, c= -5)
    m <- nls2((y ~ a*sin(b*t+c) + mean(y)), start = st1, algorithm = "brute-force")
  
    # 计算模型的拟合优度
    cor_result <- cor.test(y, predict(m))
    rm(m)
    pvalue <- cor_result$p.value
    cor_result <- cor_result$estimate
    pvalue_list[i] <- pvalue
    cor_list[i] <- cor_result
    if(cor_result > 0.75)
    {
      num <- num + 1
    }
  }
  padj <- p.adjust(pvalue_list)
  new_dataframe <- data.frame(cor_list,pvalue_list, padj)
  new_dataframe <- new_dataframe[new_dataframe$padj<0.05,]
  
  
  print(tissuename[j])
  print(length(new_dataframe[new_dataframe$cor_list>0.75,1]))
}
# fm1 <- nlme(height ~ SSasymp(age, Asym, R0, lrc),
#             data = Loblolly,
#             fixed = Asym + R0 + lrc ~ 1,
#             random = Asym ~ 1,
#             start = c(Asym = 103, R0 = -8.5, lrc = -3.3))




















####################################################meta2d + rain
#mouse
anno <- read.csv("E:/circidian_algorithm/data/annotation/mouse/mouse.csv")
tissuename <- c("adr","AOR","BFAT","bst","CER","HEA","HYP","KID","LIV","LUN","MUS","WAT")
for(j in 1:length(tissuename))
{
  print(tissuename[j])
  file_name1 <- paste("E:/circidian_algorithm/result/mouse/4h/JTK/meta2d_",tissuename[j],".csv",sep="")
  file_name2 <- paste("E:/circidian_algorithm/result/mouse/4h/rain/",tissuename[j],".csv",sep="")
  file_name3 <- paste("E:/GO_layers/final_photo/exp_amp_data/mouse_2h_result/",tissuename[j],".csv",sep="")
  data1 <- read.csv(file_name1)
  data2 <- read.csv(file_name2)
  data3 <- read.csv(file_name3)
  gene1 <- data1[data1$meta2d_BH.Q < 0.05, 1]
  print("meta2d:")
  gene1 <- anno[ anno$probe%in%  gene1, 2]
  print(length(unique(gene1)))
  gene2 <- data2[data2$bhq < 0.05, 2]
  gene2 <- anno[anno$probe%in%  gene2, 2]
  print("rain:")
  print(length(unique(gene2)))
  gene3 <- data3[data3$BH.Q<0.05,1]
  gene3 <- anno[anno$probe%in%  gene3, 2]
  print("jtk2h:")
  print(length(unique(gene3)))
  gene <- union(gene1, gene2)
  print("union:")
  print(length(union(gene1,gene2)))
  print("intersect:")
  print(length(intersect(gene,gene3)))
  file_name <- paste("E:/circidian_algorithm/circidian_gene/mouse/",tissuename[j],".csv",sep="")
  write.csv(gene, file_name)
}





#fly

meta2d_data <- read.csv("E:/circidian_algorithm/result/fly/meta2d/meta2d_JTK.flyold.csv")
meta2d_gene <- meta2d_data[meta2d_data$meta2d_BH.Q<0.05,]$CycID
rain_data <- read.csv("E:/circidian_algorithm/result/fly/rain/old.csv")
rain_gene <- rain_data[rain_data$bhq<0.05,]$data...1.
gene <- union(meta2d_gene, rain_gene)
jtk_data <- read.csv("E:/GO_layers/final_photo/exp_amp_data/fly/JTK.flyold.csv")
jtk_gene <- jtk_data[jtk_data$ADJ.P<0.05,1]
print(length(unique(jtk_gene)))
print(length(intersect(gene,jtk_gene)))
write.csv(gene, "E:/circidian_algorithm/circidian_gene/fly/old.csv")

meta2d_data <- read.csv("E:/circidian_algorithm/result/fly/meta2d/meta2d_JTK.flyyoung.csv")
meta2d_gene <- meta2d_data[meta2d_data$meta2d_BH.Q<0.05,]$CycID
rain_data <- read.csv("E:/circidian_algorithm/result/fly/rain/young.csv")
rain_gene <- rain_data[rain_data$bhq<0.05,]$data...1.
gene <- union(meta2d_gene, rain_gene)
jtk_data <- read.csv("E:/GO_layers/final_photo/exp_amp_data/fly/JTK.flyyoung.csv")
jtk_gene <- jtk_data[jtk_data$ADJ.P<0.05,1]
print(length(unique(jtk_gene)))
print(length(intersect(gene,jtk_gene)))
write.csv(gene, "E:/circidian_algorithm/circidian_gene/fly/young.csv")

#yeast
meta2d_data <- read.csv("E:/circidian_algorithm/result/yeast/meta2d/meta2d_JTK.Sample2.csv")
meta2d_gene <- meta2d_data[meta2d_data$meta2d_BH.Q<0.05,]$CycID
rain_data <- read.csv("E:/circidian_algorithm/result/yeast/rain/sample2.csv")
rain_gene <- rain_data[rain_data$bhq<0.05,]$data...1.
gene <- union(meta2d_gene, rain_gene)
jtk_data <- read.csv("E:/GO_layers/final_photo/exp_amp_data/yeast/JTK.Sample2.csv")
jtk_gene <- jtk_data[jtk_data$ADJ.P<0.05,1]
print(length(unique(jtk_gene)))
print(length(intersect(gene,jtk_gene)))
write.csv(gene, "E:/circidian_algorithm/circidian_gene/yeast/sample2.csv")


meta2d_data <- read.csv("E:/circidian_algorithm/result/yeast/meta2d/meta2d_JTK.Sample6.csv")
meta2d_gene <- meta2d_data[meta2d_data$meta2d_BH.Q<0.05,]$CycID
rain_data <- read.csv("E:/circidian_algorithm/result/yeast/rain/sample6.csv")
rain_gene <- rain_data[rain_data$bhq<0.05,]$data...1.
jtk_data <- read.csv("E:/GO_layers/final_photo/exp_amp_data/yeast/JTK.Sample6.csv")
jtk_gene <- jtk_data[jtk_data$ADJ.P<0.05,1]
gene <- union(meta2d_gene, rain_gene)
print(length(unique(jtk_gene)))
print(length(intersect(gene,jtk_gene)))
write.csv(gene, "E:/circidian_algorithm/circidian_gene/yeast/sample6.csv")


#arabidopsis
meta2d_data <- read.csv("E:/circidian_algorithm/result/arabidopsis/meta2d/meta2d_flower.csv")
meta2d_gene <- meta2d_data[meta2d_data$meta2d_BH.Q<0.05,]$CycID
rain_data <- read.csv("E:/circidian_algorithm/result/arabidopsis/rain/flower.csv")
rain_gene <- rain_data[rain_data$bhq<0.05,]$data...1.
jtk_data <- read.csv("E:/GO_layers/final_photo/exp_amp_data/arabidopsis/flower.csv")
jtk_gene <- jtk_data[jtk_data$ADJ.P<0.05,1]
gene <- union(meta2d_gene, rain_gene)
print(length(unique(jtk_gene)))
print(length(intersect(gene,jtk_gene)))
write.csv(gene, "E:/circidian_algorithm/circidian_gene/arabidopsis/flower.csv")



meta2d_data <- read.csv("E:/circidian_algorithm/result/arabidopsis/meta2d/meta2d_JTK.Agse5612.csv")
meta2d_gene <- meta2d_data[meta2d_data$meta2d_BH.Q<0.05,]$CycID
rain_data <- read.csv("E:/circidian_algorithm/result/arabidopsis/rain/seed.csv")
rain_gene <- rain_data[rain_data$bhq<0.05,]$data...1.
jtk_data <- read.csv("E:/GO_layers/final_photo/exp_amp_data/arabidopsis/JTK.Agse5612.csv")
jtk_gene <- jtk_data[jtk_data$ADJ.P<0.05,1]
gene <- union(meta2d_gene, rain_gene)
print(length(unique(jtk_gene)))
print(length(intersect(gene,jtk_gene)))
write.csv(gene, "E:/circidian_algorithm/circidian_gene/arabidopsis/seed.csv")



#neurospora
meta2d_data <- read.csv("E:/circidian_algorithm/result/neurospora/meta2d/meta2d_JTKresult_RNA.csv")
meta2d_gene <- meta2d_data[meta2d_data$meta2d_BH.Q<0.05,]$CycID
rain_data <- read.csv("E:/circidian_algorithm/result/neurospora/rain/neurospora.csv")
rain_gene <- rain_data[rain_data$bhq<0.05,]$data...1.
jtk_data <- read.csv("E:/GO_layers/final_photo/exp_amp_data/neurospora/JTKresult_RNA.csv")
jtk_gene <- jtk_data[jtk_data$ADJ.P<0.05,1]
gene <- union(meta2d_gene, rain_gene)
print(length(unique(jtk_gene)))
print(length(intersect(gene,jtk_gene)))
write.csv(gene, "E:/circidian_algorithm/circidian_gene/neurospora/neurospora.csv")


#cyanobacteria
meta2d_data <- read.csv("E:/circidian_algorithm/result/cyanobacteria/meta2d/meta2d_JTK.LLREPLICATE2.csv")
meta2d_gene <- meta2d_data[meta2d_data$meta2d_BH.Q<0.05,]$CycID
rain_data <- read.csv("E:/circidian_algorithm/result/cyanobacteria/rain/cyanobacteria.csv")
rain_gene <- rain_data[rain_data$bhq<0.05,]$data...1.
jtk_data <- read.csv("E:/GO_layers/final_photo/exp_amp_data/cyanobacteria/JTK.LLREPLICATE2.csv")
jtk_gene <- jtk_data[jtk_data$ADJ.P<0.05,1]
gene <- union(meta2d_gene, rain_gene)
print(length(unique(jtk_gene)))
print(length(intersect(gene,jtk_gene)))
write.csv(gene, "E:/circidian_algorithm/circidian_gene/cyanobacteria/cyanobacteria.csv")




meta2d_data <- read.csv("E:/circidian_algorithm/result/archaea/meta2d/meta2d_hal_1.csv")
meta2d_gene <- meta2d_data[meta2d_data$meta2d_BH.Q<0.05,]$CycID
rain_data <- read.csv("E:/circidian_algorithm/result/archaea/rain/hal1.csv")
rain_gene <- rain_data[rain_data$bhq<0.05,]$data...1.
jtk_data <- read.csv("E:/circidian_algorithm/result/archaea/meta2d/meta2d_hal_1.csv")
jtk_gene <- jtk_data[jtk_data$JTK_pvalue <0.05,1]
gene <- union(meta2d_gene, rain_gene)
print(length(unique(jtk_gene)))
print(length(intersect(gene,jtk_gene)))
write.csv(gene, "E:/circidian_algorithm/circidian_gene/archaea/hal1.csv") 

meta2d_data <- read.csv("E:/circidian_algorithm/result/archaea/meta2d/meta2d_hal_2.csv")
meta2d_gene <- meta2d_data[meta2d_data$meta2d_BH.Q<0.05,]$CycID
rain_data <- read.csv("E:/circidian_algorithm/result/archaea/rain/hal2.csv")
rain_gene <- rain_data[rain_data$bhq<0.05,]$data...1.
jtk_data <- read.csv("E:/circidian_algorithm/result/archaea/meta2d/meta2d_hal_2.csv")
jtk_gene <- jtk_data[jtk_data$JTK_pvalue <0.05,1]
gene <- union(meta2d_gene, rain_gene)
print(length(unique(jtk_gene)))
print(length(intersect(gene,jtk_gene)))
write.csv(gene, "E:/circidian_algorithm/circidian_gene/archaea/hal2.csv") 

meta2d_data <- read.csv("E:/circidian_algorithm/result/archaea/meta2d/meta2d_hal_3.csv")
meta2d_gene <- meta2d_data[meta2d_data$meta2d_BH.Q<0.05,]$CycID
rain_data <- read.csv("E:/circidian_algorithm/result/archaea/rain/hal3.csv")
rain_gene <- rain_data[rain_data$bhq<0.05,]$data...1.
jtk_data <- read.csv("E:/circidian_algorithm/result/archaea/meta2d/meta2d_hal_3.csv")
jtk_gene <- jtk_data[jtk_data$JTK_pvalue <0.05,1]
gene <- union(meta2d_gene, rain_gene)
print(length(unique(jtk_gene)))
print(length(intersect(gene,jtk_gene)))
write.csv(gene, "E:/circidian_algorithm/circidian_gene/archaea/hal3.csv") 


meta2d_data <- read.csv("E:/circidian_algorithm/result/archaea/meta2d/meta2d_hal_4.csv")
meta2d_gene <- meta2d_data[meta2d_data$meta2d_BH.Q<0.05,]$CycID
rain_data <- read.csv("E:/circidian_algorithm/result/archaea/rain/hal4.csv")
rain_gene <- rain_data[rain_data$bhq<0.05,]$data...1.
jtk_data <- read.csv("E:/circidian_algorithm/result/archaea/meta2d/meta2d_hal_4.csv")
jtk_gene <- jtk_data[jtk_data$JTK_pvalue <0.05,1]
gene <- union(meta2d_gene, rain_gene)
print(length(unique(jtk_gene)))
print(length(intersect(gene,jtk_gene)))
write.csv(gene, "E:/circidian_algorithm/circidian_gene/archaea/hal4.csv") 




meta2d_data <- read.csv("E:/circidian_algorithm/result/human/meta2d/meta2d_data.csv")
meta2d_gene <- meta2d_data[meta2d_data$meta2d_BH.Q<0.05,]$Gene.symbol
print(length(unique(meta2d_gene)))
rain_data <- read.csv("E:/circidian_algorithm/result/human/rain/cell.csv")
rain_gene <- rain_data[rain_data$bhq<0.05,]$Gene.symbol
print(length(unique(rain_gene)))
jtk_data <- read.csv("E:/circidian_algorithm/result/human/meta2d/meta2d_data.csv")
jtk_gene <- jtk_data[jtk_data$JTK_pvalue <0.05,1]
gene <- union(meta2d_gene, rain_gene)
print(length(unique(jtk_gene)))
print(length(intersect(gene,jtk_gene)))
print(length(union(rain_gene, meta2d_gene)))
write.csv(gene, "E:/circidian_algorithm/circidian_gene/human/cell.csv") 



meta2d_data <- read.csv("E:/circidian_algorithm/result/human/meta2d/meta2d_RAW DATA.csv")
meta2d_gene <- meta2d_data[meta2d_data$meta2d_BH.Q<0.05,]$GENE_SYMBOL
print(length(unique(meta2d_gene)))
rain_data <- read.csv("E:/circidian_algorithm/result/human/rain/epistem.csv")
rain_gene <- rain_data[rain_data$bhq<0.05,]$GENE_SYMBOL
print(length(unique(rain_gene)))
jtk_data <- read.csv("E:/circidian_algorithm/result/human/meta2d/meta2d_RAW DATA.csv")
jtk_gene <- jtk_data[jtk_data$JTK_pvalue <0.05,1]
gene <- union(meta2d_gene, rain_gene)
print(length(unique(jtk_gene)))
print(length(intersect(gene,jtk_gene)))
print(length(union(rain_gene, meta2d_gene)))
write.csv(gene, "E:/circidian_algorithm/circidian_gene/human/cell.csv") 








########################################################probe
####################################################meta2d + rain
#mouse

tissuename <- c("adr","AOR","BFAT","bst","CER","HEA","HYP","KID","LIV","LUN","MUS","WAT")
for(j in 1:length(tissuename))
{
  print(tissuename[j])
  file_name1 <- paste("E:/circidian_algorithm/result/mouse/4h/JTK/meta2d_",tissuename[j],".csv",sep="")
  file_name2 <- paste("E:/circidian_algorithm/result/mouse/4h/rain/",tissuename[j],".csv",sep="")
  file_name3 <- paste("E:/GO_layers/final_photo/exp_amp_data/mouse_2h_result/",tissuename[j],".csv",sep="")
  data1 <- read.csv(file_name1)
  data2 <- read.csv(file_name2)
  data3 <- read.csv(file_name3)
  gene1 <- data1[data1$meta2d_BH.Q < 0.05, 1]
  print("meta2d:")
  
  print(length(unique(gene1)))
  gene2 <- data2[data2$bhq < 0.05, 2]
  
  print("rain:")
  print(length(unique(gene2)))
  gene3 <- data3[data3$BH.Q<0.05,1]
  
  print("jtk2h:")
  print(length(unique(gene3)))
  gene <- union(gene1, gene2)
  print("union:")
  print(length(union(gene1,gene2)))
  print("intersect:")
  print(length(intersect(gene,gene3)))
  file_name <- paste("E:/circidian_algorithm/circidian_gene/mouse/",tissuename[j],".csv",sep="")
  write.csv(gene, file_name)
}





#fly

meta2d_data <- read.csv("E:/circidian_algorithm/result/fly/meta2d/meta2d_JTK.flyold.csv")
meta2d_gene <- meta2d_data[meta2d_data$meta2d_BH.Q<0.05,]$CycID
rain_data <- read.csv("E:/circidian_algorithm/result/fly/rain/old.csv")
rain_gene <- rain_data[rain_data$bhq<0.05,]$data...1.
gene <- union(meta2d_gene, rain_gene)
jtk_data <- read.csv("E:/GO_layers/final_photo/exp_amp_data/fly/JTK.flyold.csv")
jtk_gene <- jtk_data[jtk_data$ADJ.P<0.05,1]
print(length(unique(jtk_gene)))
print(length(intersect(gene,jtk_gene)))
write.csv(gene, "E:/circidian_algorithm/circidian_gene/fly/old.csv")

meta2d_data <- read.csv("E:/circidian_algorithm/result/fly/meta2d/meta2d_JTK.flyyoung.csv")
meta2d_gene <- meta2d_data[meta2d_data$meta2d_BH.Q<0.05,]$CycID
rain_data <- read.csv("E:/circidian_algorithm/result/fly/rain/young.csv")
rain_gene <- rain_data[rain_data$bhq<0.05,]$data...1.
gene <- union(meta2d_gene, rain_gene)
jtk_data <- read.csv("E:/GO_layers/final_photo/exp_amp_data/fly/JTK.flyyoung.csv")
jtk_gene <- jtk_data[jtk_data$ADJ.P<0.05,1]
print(length(unique(jtk_gene)))
print(length(intersect(gene,jtk_gene)))
write.csv(gene, "E:/circidian_algorithm/circidian_gene/fly/young.csv")

#yeast
meta2d_data <- read.csv("E:/circidian_algorithm/result/yeast/meta2d/meta2d_JTK.Sample2.csv")
meta2d_gene <- meta2d_data[meta2d_data$meta2d_BH.Q<0.05,]$CycID
rain_data <- read.csv("E:/circidian_algorithm/result/yeast/rain/sample2.csv")
rain_gene <- rain_data[rain_data$bhq<0.05,]$data...1.
gene <- union(meta2d_gene, rain_gene)
jtk_data <- read.csv("E:/GO_layers/final_photo/exp_amp_data/yeast/JTK.Sample2.csv")
jtk_gene <- jtk_data[jtk_data$ADJ.P<0.05,1]
print(length(unique(jtk_gene)))
print(length(intersect(gene,jtk_gene)))
write.csv(gene, "E:/circidian_algorithm/circidian_gene/yeast/sample2.csv")


meta2d_data <- read.csv("E:/circidian_algorithm/result/yeast/meta2d/meta2d_JTK.Sample6.csv")
meta2d_gene <- meta2d_data[meta2d_data$meta2d_BH.Q<0.05,]$CycID
rain_data <- read.csv("E:/circidian_algorithm/result/yeast/rain/sample6.csv")
rain_gene <- rain_data[rain_data$bhq<0.05,]$data...1.
jtk_data <- read.csv("E:/GO_layers/final_photo/exp_amp_data/yeast/JTK.Sample6.csv")
jtk_gene <- jtk_data[jtk_data$ADJ.P<0.05,1]
gene <- union(meta2d_gene, rain_gene)
print(length(unique(jtk_gene)))
print(length(intersect(gene,jtk_gene)))
write.csv(gene, "E:/circidian_algorithm/circidian_gene/yeast/sample6.csv")


#arabidopsis
meta2d_data <- read.csv("E:/circidian_algorithm/result/arabidopsis/meta2d/meta2d_flower.csv")
meta2d_gene <- meta2d_data[meta2d_data$meta2d_BH.Q<0.05,]$CycID
rain_data <- read.csv("E:/circidian_algorithm/result/arabidopsis/rain/flower.csv")
rain_gene <- rain_data[rain_data$bhq<0.05,]$data...1.
jtk_data <- read.csv("E:/GO_layers/final_photo/exp_amp_data/arabidopsis/flower.csv")
jtk_gene <- jtk_data[jtk_data$ADJ.P<0.05,1]
gene <- union(meta2d_gene, rain_gene)
print(length(unique(jtk_gene)))
print(length(intersect(gene,jtk_gene)))
write.csv(gene, "E:/circidian_algorithm/circidian_gene/arabidopsis/flower.csv")



meta2d_data <- read.csv("E:/circidian_algorithm/result/arabidopsis/meta2d/meta2d_JTK.Agse5612.csv")
meta2d_gene <- meta2d_data[meta2d_data$meta2d_BH.Q<0.05,]$CycID
rain_data <- read.csv("E:/circidian_algorithm/result/arabidopsis/rain/seed.csv")
rain_gene <- rain_data[rain_data$bhq<0.05,]$data...1.
jtk_data <- read.csv("E:/GO_layers/final_photo/exp_amp_data/arabidopsis/JTK.Agse5612.csv")
jtk_gene <- jtk_data[jtk_data$ADJ.P<0.05,1]
gene <- union(meta2d_gene, rain_gene)
print(length(unique(jtk_gene)))
print(length(intersect(gene,jtk_gene)))
write.csv(gene, "E:/circidian_algorithm/circidian_gene/arabidopsis/seed.csv")



#neurospora
meta2d_data <- read.csv("E:/circidian_algorithm/result/neurospora/meta2d/meta2d_JTKresult_RNA.csv")
meta2d_gene <- meta2d_data[meta2d_data$meta2d_BH.Q<0.05,]$CycID
rain_data <- read.csv("E:/circidian_algorithm/result/neurospora/rain/neurospora.csv")
rain_gene <- rain_data[rain_data$bhq<0.05,]$data...1.
jtk_data <- read.csv("E:/GO_layers/final_photo/exp_amp_data/neurospora/JTKresult_RNA.csv")
jtk_gene <- jtk_data[jtk_data$ADJ.P<0.05,1]
gene <- union(meta2d_gene, rain_gene)
print(length(unique(jtk_gene)))
print(length(intersect(gene,jtk_gene)))
write.csv(gene, "E:/circidian_algorithm/circidian_gene/neurospora/neurospora.csv")


#cyanobacteria
meta2d_data <- read.csv("E:/circidian_algorithm/result/cyanobacteria/meta2d/meta2d_JTK.LLREPLICATE2.csv")
meta2d_gene <- meta2d_data[meta2d_data$meta2d_BH.Q<0.05,]$CycID
rain_data <- read.csv("E:/circidian_algorithm/result/cyanobacteria/rain/cyanobacteria.csv")
rain_gene <- rain_data[rain_data$bhq<0.05,]$data...1.
jtk_data <- read.csv("E:/GO_layers/final_photo/exp_amp_data/cyanobacteria/JTK.LLREPLICATE2.csv")
jtk_gene <- jtk_data[jtk_data$ADJ.P<0.05,1]
gene <- union(meta2d_gene, rain_gene)
print(length(unique(jtk_gene)))
print(length(intersect(gene,jtk_gene)))
write.csv(gene, "E:/circidian_algorithm/circidian_gene/cyanobacteria/cyanobacteria.csv")




meta2d_data <- read.csv("E:/circidian_algorithm/result/archaea/meta2d/meta2d_hal_1.csv")
meta2d_gene <- meta2d_data[meta2d_data$meta2d_BH.Q<0.05,]$CycID
rain_data <- read.csv("E:/circidian_algorithm/result/archaea/rain/hal1.csv")
rain_gene <- rain_data[rain_data$bhq<0.05,]$data...1.
jtk_data <- read.csv("E:/circidian_algorithm/result/archaea/meta2d/meta2d_hal_1.csv")
jtk_gene <- jtk_data[jtk_data$JTK_pvalue <0.05,1]
gene <- union(meta2d_gene, rain_gene)
print(length(unique(jtk_gene)))
print(length(intersect(gene,jtk_gene)))
write.csv(gene, "E:/circidian_algorithm/circidian_gene/archaea/hal1.csv") 

meta2d_data <- read.csv("E:/circidian_algorithm/result/archaea/meta2d/meta2d_hal_2.csv")
meta2d_gene <- meta2d_data[meta2d_data$meta2d_BH.Q<0.05,]$CycID
rain_data <- read.csv("E:/circidian_algorithm/result/archaea/rain/hal2.csv")
rain_gene <- rain_data[rain_data$bhq<0.05,]$data...1.
jtk_data <- read.csv("E:/circidian_algorithm/result/archaea/meta2d/meta2d_hal_2.csv")
jtk_gene <- jtk_data[jtk_data$JTK_pvalue <0.05,1]
gene <- union(meta2d_gene, rain_gene)
print(length(unique(jtk_gene)))
print(length(intersect(gene,jtk_gene)))
write.csv(gene, "E:/circidian_algorithm/circidian_gene/archaea/hal2.csv") 

meta2d_data <- read.csv("E:/circidian_algorithm/result/archaea/meta2d/meta2d_hal_3.csv")
meta2d_gene <- meta2d_data[meta2d_data$meta2d_BH.Q<0.05,]$CycID
rain_data <- read.csv("E:/circidian_algorithm/result/archaea/rain/hal3.csv")
rain_gene <- rain_data[rain_data$bhq<0.05,]$data...1.
jtk_data <- read.csv("E:/circidian_algorithm/result/archaea/meta2d/meta2d_hal_3.csv")
jtk_gene <- jtk_data[jtk_data$JTK_pvalue <0.05,1]
gene <- union(meta2d_gene, rain_gene)
print(length(unique(jtk_gene)))
print(length(intersect(gene,jtk_gene)))
write.csv(gene, "E:/circidian_algorithm/circidian_gene/archaea/hal3.csv") 


meta2d_data <- read.csv("E:/circidian_algorithm/result/archaea/meta2d/meta2d_hal_4.csv")
meta2d_gene <- meta2d_data[meta2d_data$meta2d_BH.Q<0.05,]$CycID
rain_data <- read.csv("E:/circidian_algorithm/result/archaea/rain/hal4.csv")
rain_gene <- rain_data[rain_data$bhq<0.05,]$data...1.
jtk_data <- read.csv("E:/circidian_algorithm/result/archaea/meta2d/meta2d_hal_4.csv")
jtk_gene <- jtk_data[jtk_data$JTK_pvalue <0.05,1]
gene <- union(meta2d_gene, rain_gene)
print(length(unique(jtk_gene)))
print(length(intersect(gene,jtk_gene)))
write.csv(gene, "E:/circidian_algorithm/circidian_gene/archaea/hal4.csv") 




meta2d_data <- read.csv("E:/circidian_algorithm/result/human/meta2d/meta2d_data.csv")
meta2d_gene <- meta2d_data[meta2d_data$meta2d_BH.Q<0.05,]$CycID
print(length(unique(meta2d_gene)))
rain_data <- read.csv("E:/circidian_algorithm/result/human/rain/cell.csv")
rain_gene <- rain_data[rain_data$bhq<0.05,]$X
print(length(unique(rain_gene)))
jtk_data <- read.csv("E:/circidian_algorithm/result/human/meta2d/meta2d_data.csv")
jtk_gene <- jtk_data[jtk_data$JTK_pvalue <0.05,1]
gene <- union(meta2d_gene, rain_gene)
print(length(unique(jtk_gene)))
print(length(intersect(gene,jtk_gene)))
print(length(union(rain_gene, meta2d_gene)))
write.csv(gene, "E:/circidian_algorithm/circidian_gene/human/cell.csv") 



meta2d_data <- read.csv("E:/circidian_algorithm/result/human/meta2d/meta2d_RAW DATA.csv")
meta2d_gene <- meta2d_data[meta2d_data$meta2d_BH.Q<0.05,]$CycID
print(length(unique(meta2d_gene)))
rain_data <- read.csv("E:/circidian_algorithm/result/human/rain/epistem.csv")
rain_gene <- rain_data[rain_data$bhq<0.05,]$X
print(length(unique(rain_gene)))
jtk_data <- read.csv("E:/circidian_algorithm/result/human/meta2d/meta2d_RAW DATA.csv")
jtk_gene <- jtk_data[jtk_data$JTK_pvalue <0.05,1]
gene <- union(meta2d_gene, rain_gene)
print(length(unique(jtk_gene)))
print(length(intersect(gene,jtk_gene)))
print(length(union(rain_gene, meta2d_gene)))
write.csv(gene, "E:/circidian_algorithm/circidian_gene/human/epistem.csv") 






##############################human
library(MetaCycle)
meta2d(infile="E:/circidian_algorithm/data/原始数据/human/GSE13949/data.csv", filestyle="csv", outdir="E:/circidian_algorithm/result/human/meta2d",timepoints="Line1", outIntegration="onlyIntegration")



library(rain)

set.seed(123)
data <- read.csv("E:/circidian_algorithm/data/原始数据/human/GSE13949/data.csv")
input_data <- t(data[,2:49])
rainresult <- rain(input_data, period=24, deltat=1, peak.border=c(0.3,0.7),verbose=FALSE) 
pval <- rainresult$pVal
bhq <- p.adjust(pval, method="fdr")
temp <- data.frame(data[,1], pval, bhq, rainresult$phase)
print(sum(pval < 0.05))
print(sum(bhq < 0.05))
write.csv(temp, "E:/circidian_algorithm/result/human/rain/cell.csv")


library(MetaCycle)
meta2d(infile="E:/circidian_algorithm/data/原始数据/human/gse50631/RAW DATA.csv", filestyle="csv", outdir="E:/circidian_algorithm/result/human/meta2d",timepoints="Line1", outIntegration="onlyIntegration")


data <- read.csv("E:/circidian_algorithm/data/原始数据/human/gse50631/RAW DATA.csv")
input_data <- t(data[,2:31])
rainresult <- rain(input_data, period=24, deltat=5, nr.series=3,peak.border=c(0.3,0.7),verbose=FALSE) 
pval <- rainresult$pVal
bhq <- p.adjust(pval, method="fdr")
temp <- data.frame(data[,1], pval, bhq, rainresult$phase)
print(sum(pval < 0.05))
print(sum(bhq < 0.05))
write.csv(temp, "E:/circidian_algorithm/result/human/rain/epistem.csv")


########################neurospora

library(MetaCycle)
meta2d(infile="E:/circidian_algorithm/data/原始数据/Neurospora/gse113845/DATA.csv", filestyle="csv", outdir="E:/circidian_algorithm/result/neurospora/meta2d",timepoints="Line1", minper=1,maxper=28, outIntegration="onlyIntegration")


data <- read.csv("E:/circidian_algorithm/data/原始数据/Neurospora/gse113845/DATA.csv")
input_data <- t(data[,2:12])
rainresult <- rain(input_data, period=24, deltat=2,peak.border=c(0.3,0.7),verbose=FALSE) 
pval <- rainresult$pVal
bhq <- p.adjust(pval, method="fdr")
temp <- data.frame(data[,1], pval, bhq, rainresult$phase)
print(sum(pval < 0.05))
print(sum(bhq < 0.05))
write.csv(temp, "E:/circidian_algorithm/result/neurospora/rain/2.csv")






#####################################################################################mouse abosulte expression change
library(ggplot2)
tissuename <- c("adr","AOR","BFAT","bst","CER","HEA","HYP","KID","LIV","LUN","MUS","WAT")
for(j in 1:12)
{
  file_name <- paste("E:/circidian_algorithm/result/mouse/4h/JTK/meta2d_",tissuename[j],".csv",sep="")
  print(tissuename[j])
  mouse_data <- read.csv(file_name)
  name <- mouse_data[mouse_data$meta2d_BH.Q<0.05,1]
  line1 = abs(mouse_data[,25] - mouse_data[,24])
  line2 = abs(mouse_data[,26] - mouse_data[,25])
  line3 = abs(mouse_data[,27] - mouse_data[,26])
  line4 = abs(mouse_data[,28] - mouse_data[,27])
  line5 = abs(mouse_data[,29] - mouse_data[,28])
  line6 = abs(mouse_data[,30] - mouse_data[,29])
  line7 = abs(mouse_data[,31] - mouse_data[,30])
  line8 = abs(mouse_data[,32] - mouse_data[,31])
  line9 = abs(mouse_data[,33] - mouse_data[,32])
  line10 = abs(mouse_data[,34] - mouse_data[,33])
  line11 = abs(mouse_data[,35] - mouse_data[,34])
  line12 = abs(mouse_data[,24] - mouse_data[,35])
  abs_list <- c()
  for(i in 1:nrow(mouse_data))
  {
    absk = max(c(line1[i],line2[i],line3[i],line4[i],line5[i],line6[i],line7[i],line8[i],line9[i],line10[i],line11[i],line12[i]))
    abs_list[i] <- absk
  }
  abs_name <- mouse_data[order(abs_list,decreasing = TRUE),1]
  abs_name <- abs_name[1:floor(length(mouse_data[,1])/1.25)]
  print(length(name))
  print(length(intersect(name,abs_name)))
}

tissuename <- c("adr","AOR","BFAT","bst","CER","HEA","HYP","KID","LIV","LUN","MUS","WAT")
for(j in 1:12)
{
file_name <- paste("E:/circidian_algorithm/result/mouse/4h/JTK/meta2d_",tissuename[j],".csv",sep="")
print(tissuename[j])
mouse_data <- read.csv(file_name)
name <- mouse_data[mouse_data$meta2d_BH.Q<0.05,1]
line1 = abs(mouse_data[,25] - mouse_data[,24])
line2 = abs(mouse_data[,26] - mouse_data[,25])
line3 = abs(mouse_data[,27] - mouse_data[,26])
line4 = abs(mouse_data[,28] - mouse_data[,27])
line5 = abs(mouse_data[,29] - mouse_data[,28])
line6 = abs(mouse_data[,30] - mouse_data[,29])
line7 = abs(mouse_data[,31] - mouse_data[,30])
line8 = abs(mouse_data[,32] - mouse_data[,31])
line9 = abs(mouse_data[,33] - mouse_data[,32])
line10 = abs(mouse_data[,34] - mouse_data[,33])
line11 = abs(mouse_data[,35] - mouse_data[,34])
line12 = abs(mouse_data[,24] - mouse_data[,35])
abs_list <- c()
k1_list <- c()
for(i in 1:nrow(mouse_data))
{
  absk = mean(c(line1[i],line2[i],line3[i],line4[i],line5[i],line6[i],line7[i],line8[i],line9[i],line10[i],line11[i],line12[i]))
  abs_list[i] <- absk
  k1_list[i] <- mean(as.numeric(mouse_data[i,24:35]))
}
plot_data <- data.frame(mouse_data[,1], abs_list, k1_list, 0)
plot_data[plot_data$mouse_data...1. %in% name,]$X0 <- 1
plot_data$abs_list <- log(plot_data$abs_list)
plot_data$k1_list <- log(plot_data$k1_list)
#ggplot(plot_data, aes(x=k1_list,y=abs_list)) + geom_point(aes(colour=factor(plot_data$X0)))+scale_color_manual(values =c('grey','red'))

}





file_name <- paste("E:/circidian_algorithm/result/mouse/4h/JTK/meta2d_",tissuename[j],".csv",sep="")
print(tissuename[j])
mouse_data <- read.csv(file_name)
name <- mouse_data[mouse_data$meta2d_BH.Q<0.05,1]
line1 = abs(mouse_data[,25] - mouse_data[,24])
line2 = abs(mouse_data[,26] - mouse_data[,25])
line3 = abs(mouse_data[,27] - mouse_data[,26])
line4 = abs(mouse_data[,28] - mouse_data[,27])
line5 = abs(mouse_data[,29] - mouse_data[,28])
line6 = abs(mouse_data[,30] - mouse_data[,29])
line7 = abs(mouse_data[,31] - mouse_data[,30])
line8 = abs(mouse_data[,32] - mouse_data[,31])
line9 = abs(mouse_data[,33] - mouse_data[,32])
line10 = abs(mouse_data[,34] - mouse_data[,33])
line11 = abs(mouse_data[,35] - mouse_data[,34])
line12 = abs(mouse_data[,24] - mouse_data[,35])
abs_list <- c()
k1_list <- c()
plot_frame <- data.frame(0,0,0)
for(i in 1:nrow(mouse_data))
{
  temp <- data.frame(0,0,0)
  temp[1,] <- c(mouse_data[i,1],mouse_data[i,24],line1[i])
  temp <- data.frame(nrow=11,ncol=3)
  temp <- data.frame(c(mouse_data[i,1],mouse_data[i,25],line2[i]), c(mouse_data[i,1],mouse_data[i,26],line3[i]),c(mouse_data[i,1],mouse_data[i,27],line4[i]),c(mouse_data[i,1],mouse_data[i,28],line5[i]),
                     c(mouse_data[i,1],mouse_data[i,29],line6[i]),c(mouse_data[i,1],mouse_data[i,30],line7[i]),c(mouse_data[i,1],mouse_data[i,31],line8[i]),c(mouse_data[i,1],mouse_data[i,32],line9[i]),
                     c(mouse_data[i,1],mouse_data[i,33],line10[i]),c(mouse_data[i,1],mouse_data[i,34],line11[i]),c(mouse_data[i,1],mouse_data[i,35],line12[i]),nrow=11,ncol=3)
  temp <- rbind(temp, c(mouse_data[i,1],mouse_data[i,25],line2[i]))
  temp <- rbind(temp, c(mouse_data[i,1],mouse_data[i,26],line3[i]))
  temp <- rbind(temp, c(mouse_data[i,1],mouse_data[i,27],line4[i]))
  temp <- rbind(temp, c(mouse_data[i,1],mouse_data[i,28],line5[i]))
  temp <- rbind(temp, c(mouse_data[i,1],mouse_data[i,29],line6[i]))
  temp <- rbind(temp, c(mouse_data[i,1],mouse_data[i,30],line7[i]))
  temp <- rbind(temp, c(mouse_data[i,1],mouse_data[i,31],line8[i]))
  temp <- rbind(temp, c(mouse_data[i,1],mouse_data[i,32],line9[i]))
  temp <- rbind(temp, c(mouse_data[i,1],mouse_data[i,33],line10[i]))
  temp <- rbind(temp, c(mouse_data[i,1],mouse_data[i,34],line11[i]))
  temp <- rbind(temp, c(mouse_data[i,1],mouse_data[i,35],line12[i]))
  absk = mean(c(line1[i],line2[i],line3[i],line4[i],line5[i],line6[i],line7[i],line8[i],line9[i],line10[i],line11[i],line12[i]))
  abs_list[i] <- absk
  k1_list[i] <- mean(as.numeric(mouse_data[i,24:35]))
  plot_frame <- rbind(plot_frame, temp)
}
plot_frame <- plot_frame[-1,]
plot_data <- cbind(plot_frame, 0)
plot_data[plot_data$X0 %in% name,]$`0` <- 1
plot_data$X0.1 <- log(plot_data$X0.1)
plot_data$X0.2 <- log(plot_data$X0.2)
ggplot(plot_data, aes(x=X0.1,y=X0.2)) + geom_point(aes(colour=factor(plot_data$`0`)))+scale_color_manual(values =c('grey','red'))
cor(plot_data$k1_list,plot_data$abs_list)


















################################fly
library(ggplot2)


file_name <- paste("E:/circidian_algorithm/result/fly/meta2d/meta2d_","JTK.flyold.csv",sep="")
mouse_data <- read.csv(file_name)
mouse_data <- na.omit(mouse_data)
name <- mouse_data[mouse_data$meta2d_BH.Q<0.05,1]
line1 = abs(mouse_data[,25] - mouse_data[,24])
line2 = abs(mouse_data[,26] - mouse_data[,25])
line3 = abs(mouse_data[,27] - mouse_data[,26])
line4 = abs(mouse_data[,28] - mouse_data[,27])
line5 = abs(mouse_data[,29] - mouse_data[,28])
line6 = abs(mouse_data[,30] - mouse_data[,29])
line7 = abs(mouse_data[,31] - mouse_data[,30])
line8 = abs(mouse_data[,32] - mouse_data[,31])
line9 = abs(mouse_data[,33] - mouse_data[,32])
line10 = abs(mouse_data[,34] - mouse_data[,33])
line11 = abs(mouse_data[,35] - mouse_data[,34])
line12 = abs(mouse_data[,24] - mouse_data[,35])
abs_list <- c()
k1_list <- c()
for(i in 1:nrow(mouse_data))
{
  absk = mean(c(line1[i],line2[i],line3[i],line4[i],line5[i],line6[i],line7[i],line8[i],line9[i],line10[i],line11[i],line12[i]))
  abs_list[i] <- absk
  k1_list[i] <- mean(as.numeric(mouse_data[i,24:35]))
}
plot_data <- data.frame(mouse_data[,1], abs_list, k1_list, 0)
plot_data[plot_data$mouse_data...1. %in% name,]$X0 <- 1
plot_data$abs_list <- log(plot_data$abs_list)
plot_data$k1_list <- log(plot_data$k1_list)
ggplot(plot_data, aes(x=k1_list,y=abs_list)) + geom_point(aes(colour=factor(plot_data$X0)))+scale_color_manual(values =c('grey','red'))
t.test(plot_data[plot_data$X0==1,2], plot_data[plot_data$X0==0,2])
cor.test(abs_list, k1_list)

  
file_name <- paste("E:/circidian_algorithm/result/fly/meta2d/meta2d_","JTK.flyyoung.csv",sep="")
mouse_data <- read.csv(file_name)
mouse_data <- na.omit(mouse_data)
name <- mouse_data[mouse_data$meta2d_BH.Q<0.05,1]
line1 = abs(mouse_data[,25] - mouse_data[,24])
line2 = abs(mouse_data[,26] - mouse_data[,25])
line3 = abs(mouse_data[,27] - mouse_data[,26])
line4 = abs(mouse_data[,28] - mouse_data[,27])
line5 = abs(mouse_data[,29] - mouse_data[,28])
line6 = abs(mouse_data[,30] - mouse_data[,29])
line7 = abs(mouse_data[,31] - mouse_data[,30])
line8 = abs(mouse_data[,32] - mouse_data[,31])
line9 = abs(mouse_data[,33] - mouse_data[,32])
line10 = abs(mouse_data[,34] - mouse_data[,33])
line11 = abs(mouse_data[,35] - mouse_data[,34])
line12 = abs(mouse_data[,24] - mouse_data[,35])
abs_list <- c()
k1_list <- c()
for(i in 1:nrow(mouse_data))
{
  absk = mean(c(line1[i],line2[i],line3[i],line4[i],line5[i],line6[i],line7[i],line8[i],line9[i],line10[i],line11[i],line12[i]))
  abs_list[i] <- absk
  k1_list[i] <- mean(as.numeric(mouse_data[i,24:35]))
}
plot_data <- data.frame(mouse_data[,1], abs_list, k1_list, 0)
plot_data[plot_data$mouse_data...1. %in% name,]$X0 <- 1
plot_data$abs_list <- log(plot_data$abs_list)
plot_data$k1_list <- log(plot_data$k1_list)
ggplot(plot_data, aes(x=k1_list,y=abs_list)) + geom_point(aes(colour=factor(plot_data$X0)))+scale_color_manual(values =c('grey','red'))
t.test(plot_data[plot_data$X0==1,2], plot_data[plot_data$X0==0,2])
cor.test(abs_list, k1_list)





file_name <- paste("E:/circidian_algorithm/result/yeast/meta2d/meta2d_","JTK.Sample2.csv",sep="")
mouse_data <- read.csv(file_name)
mouse_data <- na.omit(mouse_data)
name <- mouse_data[mouse_data$meta2d_BH.Q<0.05,1]
line1 = abs(mouse_data[,25] - mouse_data[,24])
line2 = abs(mouse_data[,26] - mouse_data[,25])
line3 = abs(mouse_data[,27] - mouse_data[,26])
line4 = abs(mouse_data[,28] - mouse_data[,27])
line5 = abs(mouse_data[,29] - mouse_data[,28])
line6 = abs(mouse_data[,30] - mouse_data[,29])
line7 = abs(mouse_data[,31] - mouse_data[,30])
line8 = abs(mouse_data[,32] - mouse_data[,31])
line9 = abs(mouse_data[,33] - mouse_data[,32])
line10 = abs(mouse_data[,34] - mouse_data[,33])
line11 = abs(mouse_data[,35] - mouse_data[,34])
line12 = abs(mouse_data[,36] - mouse_data[,35])
line13 = abs(mouse_data[,37] - mouse_data[,36])
line14 = abs(mouse_data[,38] - mouse_data[,37])
line15 = abs(mouse_data[,39] - mouse_data[,38])
line16 = abs(mouse_data[,40] - mouse_data[,39])
line17 = abs(mouse_data[,41] - mouse_data[,40])
line18 = abs(mouse_data[,42] - mouse_data[,41])
line19 = abs(mouse_data[,43] - mouse_data[,42])
line20 = abs(mouse_data[,24] - mouse_data[,43])
abs_list <- c()
k1_list <- c()
for(i in 1:nrow(mouse_data))
{
  absk = mean(c(line1[i],line2[i],line3[i],line4[i],line5[i],line6[i],line7[i],line8[i],line9[i],line10[i],line11[i],line12[i],line13[i],line14[i],line15[i],line16[i],line17[i],line18[i],line19[i],line20[i]))
  abs_list[i] <- absk
  k1_list[i] <- mean(as.numeric(mouse_data[i,24:35]))
}
plot_data <- data.frame(mouse_data[,1], abs_list, k1_list, 0)
plot_data[plot_data$mouse_data...1. %in% name,]$X0 <- 1
plot_data$abs_list <- log(plot_data$abs_list)
plot_data$k1_list <- log(plot_data$k1_list)
ggplot(plot_data, aes(x=k1_list,y=abs_list)) + geom_point(aes(colour=factor(plot_data$X0)))+scale_color_manual(values =c('grey','red'))
t.test(plot_data[plot_data$X0==1,2], plot_data[plot_data$X0==0,2])
cor.test(abs_list, k1_list)


file_name <- paste("E:/circidian_algorithm/result/yeast/meta2d/meta2d_","JTK.Sample6.csv",sep="")
mouse_data <- read.csv(file_name)
mouse_data <- na.omit(mouse_data)
name <- mouse_data[mouse_data$meta2d_BH.Q<0.05,1]
line1 = abs(mouse_data[,25] - mouse_data[,24])
line2 = abs(mouse_data[,26] - mouse_data[,25])
line3 = abs(mouse_data[,27] - mouse_data[,26])
line4 = abs(mouse_data[,28] - mouse_data[,27])
line5 = abs(mouse_data[,29] - mouse_data[,28])
line6 = abs(mouse_data[,30] - mouse_data[,29])
line7 = abs(mouse_data[,31] - mouse_data[,30])
line8 = abs(mouse_data[,32] - mouse_data[,31])
line9 = abs(mouse_data[,33] - mouse_data[,32])
line10 = abs(mouse_data[,34] - mouse_data[,33])
line11 = abs(mouse_data[,35] - mouse_data[,34])
line12 = abs(mouse_data[,36] - mouse_data[,35])
line13 = abs(mouse_data[,37] - mouse_data[,36])
line14 = abs(mouse_data[,38] - mouse_data[,37])
line15 = abs(mouse_data[,39] - mouse_data[,38])
line16 = abs(mouse_data[,40] - mouse_data[,39])
line17 = abs(mouse_data[,41] - mouse_data[,40])
line18 = abs(mouse_data[,42] - mouse_data[,41])
line19 = abs(mouse_data[,43] - mouse_data[,42])
line20 = abs(mouse_data[,24] - mouse_data[,43])
abs_list <- c()
k1_list <- c()
for(i in 1:nrow(mouse_data))
{
  absk = mean(c(line1[i],line2[i],line3[i],line4[i],line5[i],line6[i],line7[i],line8[i],line9[i],line10[i],line11[i],line12[i],line13[i],line14[i],line15[i],line16[i],line17[i],line18[i],line19[i],line20[i]))
  abs_list[i] <- absk
  k1_list[i] <- mean(as.numeric(mouse_data[i,24:35]))
}
plot_data <- data.frame(mouse_data[,1], abs_list, k1_list, 0)
plot_data[plot_data$mouse_data...1. %in% name,]$X0 <- 1
plot_data$abs_list <- log(plot_data$abs_list)
plot_data$k1_list <- log(plot_data$k1_list)
ggplot(plot_data, aes(x=k1_list,y=abs_list)) + geom_point(aes(colour=factor(plot_data$X0)))+scale_color_manual(values =c('grey','red'))
t.test(plot_data[plot_data$X0==1,2], plot_data[plot_data$X0==0,2])
cor.test(abs_list, k1_list)







file_name <- paste("E:/circidian_algorithm/result/arabidopsis/meta2d/meta2d_","JTK.Agse5612.csv",sep="")
mouse_data <- read.csv(file_name)
mouse_data <- na.omit(mouse_data)
name <- mouse_data[mouse_data$meta2d_BH.Q<0.05,1]
line1 = abs(mouse_data[,25] - mouse_data[,24])
line2 = abs(mouse_data[,26] - mouse_data[,25])
line3 = abs(mouse_data[,27] - mouse_data[,26])
line4 = abs(mouse_data[,28] - mouse_data[,27])
line5 = abs(mouse_data[,29] - mouse_data[,28])
line6 = abs(mouse_data[,30] - mouse_data[,29])
line7 = abs(mouse_data[,31] - mouse_data[,30])
line8 = abs(mouse_data[,32] - mouse_data[,31])
line9 = abs(mouse_data[,33] - mouse_data[,32])
line10 = abs(mouse_data[,34] - mouse_data[,33])
line11 = abs(mouse_data[,35] - mouse_data[,34])
line12 = abs(mouse_data[,36] - mouse_data[,35])
line13 = abs(mouse_data[,24] - mouse_data[,36])
abs_list <- c()
k1_list <- c()
for(i in 1:nrow(mouse_data))
{
  absk = mean(c(line1[i],line2[i],line3[i],line4[i],line5[i],line6[i],line7[i],line8[i],line9[i],line10[i],line11[i],line12[i],line13[i]))
  abs_list[i] <- absk
  k1_list[i] <- mean(as.numeric(mouse_data[i,24:35]))
}
plot_data <- data.frame(mouse_data[,1], abs_list, k1_list, 0)
plot_data[plot_data$mouse_data...1. %in% name,]$X0 <- 1
plot_data$abs_list <- log(plot_data$abs_list)
plot_data$k1_list <- log(plot_data$k1_list)
ggplot(plot_data, aes(x=k1_list,y=abs_list)) + geom_point(aes(colour=factor(plot_data$X0)))+scale_color_manual(values =c('grey','red'))
t.test(plot_data[plot_data$X0==1,2], plot_data[plot_data$X0==0,2])
cor.test(abs_list, k1_list)



file_name <- paste("E:/circidian_algorithm/result/neurospora/meta2d/meta2d_","JTKresult_RNA.csv",sep="")
mouse_data <- read.csv(file_name)
mouse_data <- na.omit(mouse_data)
name <- mouse_data[mouse_data$meta2d_BH.Q<0.05,1]
line1 = abs(mouse_data[,25] - mouse_data[,24])
line2 = abs(mouse_data[,26] - mouse_data[,25])
line3 = abs(mouse_data[,27] - mouse_data[,26])
line4 = abs(mouse_data[,28] - mouse_data[,27])
line5 = abs(mouse_data[,29] - mouse_data[,28])
line6 = abs(mouse_data[,30] - mouse_data[,29])
line7 = abs(mouse_data[,31] - mouse_data[,30])
line8 = abs(mouse_data[,32] - mouse_data[,31])
line9 = abs(mouse_data[,33] - mouse_data[,32])
line10 = abs(mouse_data[,34] - mouse_data[,33])
line11 = abs(mouse_data[,35] - mouse_data[,34])
line12 = abs(mouse_data[,24] - mouse_data[,35])
abs_list <- c()
k1_list <- c()
for(i in 1:nrow(mouse_data))
{
  absk = mean(c(line1[i],line2[i],line3[i],line4[i],line5[i],line6[i],line7[i],line8[i],line9[i],line10[i],line11[i],line12[i]))
  abs_list[i] <- absk
  k1_list[i] <- mean(as.numeric(mouse_data[i,24:35]))
}
plot_data <- data.frame(mouse_data[,1], abs_list, k1_list, 0)
plot_data[plot_data$mouse_data...1. %in% name,]$X0 <- 1
plot_data$abs_list <- log(plot_data$abs_list)
plot_data$k1_list <- log(plot_data$k1_list)
ggplot(plot_data, aes(x=k1_list,y=abs_list)) + geom_point(aes(colour=factor(plot_data$X0)))+scale_color_manual(values =c('grey','red'))
t.test(plot_data[plot_data$X0==1,2], plot_data[plot_data$X0==0,2])
cor.test(abs_list, k1_list)








file_name <- paste("E:/circidian_algorithm/result/archaea/meta2d/meta2d_","hal_1.csv",sep="")
mouse_data <- read.csv(file_name)
mouse_data <- na.omit(mouse_data)
name <- mouse_data[mouse_data$meta2d_BH.Q<0.05,1]
line1 = abs(mouse_data[,25] - mouse_data[,24])
line2 = abs(mouse_data[,26] - mouse_data[,25])
line3 = abs(mouse_data[,27] - mouse_data[,26])
line4 = abs(mouse_data[,28] - mouse_data[,27])
line5 = abs(mouse_data[,29] - mouse_data[,28])
line6 = abs(mouse_data[,30] - mouse_data[,29])
line7 = abs(mouse_data[,31] - mouse_data[,30])
line8 = abs(mouse_data[,32] - mouse_data[,31])
line9 = abs(mouse_data[,33] - mouse_data[,32])
line10 = abs(mouse_data[,34] - mouse_data[,33])
line11 = abs(mouse_data[,35] - mouse_data[,34])
line12 = abs(mouse_data[,36] - mouse_data[,35])
line13 = abs(mouse_data[,37] - mouse_data[,36])
line14 = abs(mouse_data[,38] - mouse_data[,37])
line15 = abs(mouse_data[,39] - mouse_data[,38])
line16 = abs(mouse_data[,40] - mouse_data[,39])
line17 = abs(mouse_data[,41] - mouse_data[,40])
line18 = abs(mouse_data[,42] - mouse_data[,41])
line19 = abs(mouse_data[,24] - mouse_data[,42])
abs_list <- c()
k1_list <- c()
for(i in 1:nrow(mouse_data))
{
  absk = mean(c(line1[i],line2[i],line3[i],line4[i],line5[i],line6[i],line7[i],line8[i],line9[i],line10[i],line11[i],line12[i],line13[i],line14[i],line15[i],line16[i],line17[i],line18[i],line19[i]))
  abs_list[i] <- absk
  k1_list[i] <- mean(as.numeric(mouse_data[i,24:35]))
}
plot_data <- data.frame(mouse_data[,1], abs_list, k1_list, 0)
plot_data[plot_data$mouse_data...1. %in% name,]$X0 <- 1
plot_data$abs_list <- log(plot_data$abs_list)
plot_data$k1_list <- log(plot_data$k1_list)
ggplot(plot_data, aes(x=k1_list,y=abs_list)) + geom_point(aes(colour=factor(plot_data$X0)))+scale_color_manual(values =c('grey','red'))
t.test(plot_data[plot_data$X0==1,2], plot_data[plot_data$X0==0,2])
cor.test(abs_list, k1_list)



file_name <- paste("E:/circidian_algorithm/result/archaea/meta2d/meta2d_","hal_2.csv",sep="")
mouse_data <- read.csv(file_name)
mouse_data <- na.omit(mouse_data)
name <- mouse_data[mouse_data$meta2d_BH.Q<0.05,1]
line1 = abs(mouse_data[,25] - mouse_data[,24])
line2 = abs(mouse_data[,26] - mouse_data[,25])
line3 = abs(mouse_data[,27] - mouse_data[,26])
line4 = abs(mouse_data[,28] - mouse_data[,27])
line5 = abs(mouse_data[,29] - mouse_data[,28])
line6 = abs(mouse_data[,30] - mouse_data[,29])
line7 = abs(mouse_data[,31] - mouse_data[,30])
line8 = abs(mouse_data[,32] - mouse_data[,31])
line9 = abs(mouse_data[,33] - mouse_data[,32])
line10 = abs(mouse_data[,34] - mouse_data[,33])
line11 = abs(mouse_data[,35] - mouse_data[,34])
line12 = abs(mouse_data[,36] - mouse_data[,35])
line13 = abs(mouse_data[,37] - mouse_data[,36])
line14 = abs(mouse_data[,38] - mouse_data[,37])
line15 = abs(mouse_data[,39] - mouse_data[,38])
line16 = abs(mouse_data[,40] - mouse_data[,39])
line17 = abs(mouse_data[,41] - mouse_data[,40])
line18 = abs(mouse_data[,42] - mouse_data[,41])
line19 = abs(mouse_data[,24] - mouse_data[,42])
abs_list <- c()
k1_list <- c()
for(i in 1:nrow(mouse_data))
{
  absk = mean(c(line1[i],line2[i],line3[i],line4[i],line5[i],line6[i],line7[i],line8[i],line9[i],line10[i],line11[i],line12[i],line13[i],line14[i],line15[i],line16[i],line17[i],line18[i],line19[i]))
  abs_list[i] <- absk
  k1_list[i] <- mean(as.numeric(mouse_data[i,24:35]))
}
plot_data <- data.frame(mouse_data[,1], abs_list, k1_list, 0)
plot_data[plot_data$mouse_data...1. %in% name,]$X0 <- 1
plot_data$abs_list <- log(plot_data$abs_list)
plot_data$k1_list <- log(plot_data$k1_list)
ggplot(plot_data, aes(x=k1_list,y=abs_list)) + geom_point(aes(colour=factor(plot_data$X0)))+scale_color_manual(values =c('grey','red'))
t.test(plot_data[plot_data$X0==1,2], plot_data[plot_data$X0==0,2])
cor.test(abs_list, k1_list)



file_name <- paste("E:/circidian_algorithm/result/archaea/meta2d/meta2d_","hal_4.csv",sep="")
mouse_data <- read.csv(file_name)
mouse_data <- na.omit(mouse_data)
name <- mouse_data[mouse_data$meta2d_BH.Q<0.05,1]
line1 = abs(mouse_data[,25] - mouse_data[,24])
line2 = abs(mouse_data[,26] - mouse_data[,25])
line3 = abs(mouse_data[,27] - mouse_data[,26])
line4 = abs(mouse_data[,28] - mouse_data[,27])
line5 = abs(mouse_data[,29] - mouse_data[,28])
line6 = abs(mouse_data[,30] - mouse_data[,29])
line7 = abs(mouse_data[,31] - mouse_data[,30])
line8 = abs(mouse_data[,32] - mouse_data[,31])
line9 = abs(mouse_data[,33] - mouse_data[,32])
line10 = abs(mouse_data[,34] - mouse_data[,33])
line11 = abs(mouse_data[,35] - mouse_data[,34])
line12 = abs(mouse_data[,36] - mouse_data[,35])
line13 = abs(mouse_data[,37] - mouse_data[,36])
line14 = abs(mouse_data[,38] - mouse_data[,37])
line15 = abs(mouse_data[,39] - mouse_data[,38])
line16 = abs(mouse_data[,40] - mouse_data[,39])
line17 = abs(mouse_data[,41] - mouse_data[,40])
line18 = abs(mouse_data[,42] - mouse_data[,41])
line19 = abs(mouse_data[,24] - mouse_data[,42])
abs_list <- c()
k1_list <- c()
for(i in 1:nrow(mouse_data))
{
  absk = mean(c(line1[i],line2[i],line3[i],line4[i],line5[i],line6[i],line7[i],line8[i],line9[i],line10[i],line11[i],line12[i],line13[i],line14[i],line15[i],line16[i],line17[i],line18[i],line19[i]))
  abs_list[i] <- absk
  k1_list[i] <- mean(as.numeric(mouse_data[i,24:35]))
}
plot_data <- data.frame(mouse_data[,1], abs_list, k1_list, 0)
plot_data[plot_data$mouse_data...1. %in% name,]$X0 <- 1
plot_data$abs_list <- log(plot_data$abs_list)
plot_data$k1_list <- log(plot_data$k1_list)
ggplot(plot_data, aes(x=k1_list,y=abs_list)) + geom_point(aes(colour=factor(plot_data$X0)))+scale_color_manual(values =c('grey','red'))
t.test(plot_data[plot_data$X0==1,2], plot_data[plot_data$X0==0,2])
cor.test(abs_list, k1_list)
