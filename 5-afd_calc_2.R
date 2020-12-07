#####################################################
# Lake Michigan vs. Lake Champlain - calculate Z(AFD)
#####################################################

afdlmlchwe<-read.csv(file="Z:/gatksep/rnavdna/pmarinus1/dp5af0025VF/popmuscle-wo408/afdlmlchwe.csv")
avgafd_lmlc_hwe<-colMeans(as.matrix(afdlmlchwe$afd))
sdafd_lmlc_hwe<-sd(afdlmlchwe$afd)
afdlmlchwe$Zafd<-(afdlmlchwe$afd-avgafd_lmlc_hwe)/sdafd_lmlc_hwe
write.csv(afdlmlchwe,"Z:/gatksep/rnavdna/pmarinus1/dp5af0025VF/popmuscle-wo408/afdlmlchwe_zfst.csv")

########################################################
# Lake Michigan vs. Connecticut River - calculate Z(AFD)
########################################################

afdlmcthwe<-read.csv(file="Z:/gatksep/rnavdna/pmarinus1/dp5af0025VF/popmuscle-wo408/afdlmcthwe.csv")
avgafd_lmct_hwe<-colMeans(as.matrix(afdlmcthwe$afd))
sdafd_lmct_hwe<-sd(afdlmcthwe$afd)
afdlmcthwe$Zafd<-(afdlmcthwe$afd-avgafd_lmct_hwe)/sdafd_lmct_hwe
write.csv(afdlmcthwe,"Z:/gatksep/rnavdna/pmarinus1/dp5af0025VF/popmuscle-wo408/afdlmcthwe_zfst.csv")

#########################################################
# Lake Champlain vs. Connecticut River - calculate Z(AFD)
#########################################################

afdlccthwe<-read.csv(file="Z:/gatksep/rnavdna/pmarinus1/dp5af0025VF/popmuscle-wo408/afdlccthwe.csv")
avgafd_lcct_hwe<-colMeans(as.matrix(afdlccthwe$afd))
sdafd_lcct_hwe<-sd(afdlccthwe$afd)
afdlccthwe$Zafd<-(afdlccthwe$afd-avgafd_lcct_hwe)/sdafd_lcct_hwe
write.csv(afdlccthwe,"Z:/gatksep/rnavdna/pmarinus1/dp5af0025VF/popmuscle-wo408/afdlccthwe_zfst.csv")
