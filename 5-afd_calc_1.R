###################################
## Lake Michigan vs. Lake Champlain
###################################

flmlc<-read.csv("Z:/gatksep/rnavdna/pmarinus1/dp5af0025VF/popmuscle-wo408/flmlc.csv")
flmlc$afd<-(abs(flmlc$flm.freq1-flmlc$flc.freq1)+abs(flmlc$flm.freq2-flmlc$flc.freq2)+abs(flmlc$flm.freq3-flmlc$flc.freq3)+abs(flmlc$flm.freq4-flmlc$flc.freq4))/2
write.csv(flmlc,"Z:/gatksep/rnavdna/pmarinus1/dp5af0025VF/popmuscle-wo408/afdlmlc.csv")

fstlmlchwe<-read.csv("Z:/gatksep/rnavdna/pmarinus1/dp5af0025VF/popmuscle-wo408/fstlmlchwe.csv")
hwelmlc<-data.frame(fstlmlchwe$locus)
afdlmlchwe<-merge(flmlc,hwelmlc,by.x="flm.locus",by.y="fstlmlchwe.locus")
write.csv(afdlmlchwe,"Z:/gatksep/rnavdna/pmarinus1/dp5af0025VF/popmuscle-wo408/afdlmlchwe.csv")

######################################
## Lake Michigan vs. Connecticut River
######################################

flmct<-read.csv("Z:/gatksep/rnavdna/pmarinus1/dp5af0025VF/popmuscle-wo408/flmct.csv")
flmct$afd<-(abs(flmct$flm.freq1-flmct$fct.freq1)+abs(flmct$flm.freq2-flmct$fct.freq2)+abs(flmct$flm.freq3-flmct$fct.freq3)+abs(flmct$flm.freq4-flmct$fct.freq4))/2
write.csv(flmct,"Z:/gatksep/rnavdna/pmarinus1/dp5af0025VF/popmuscle-wo408/afdlmct.csv")

fstlmcthwe<-read.csv("Z:/gatksep/rnavdna/pmarinus1/dp5af0025VF/popmuscle-wo408/fstlmcthwe.csv")
hwelmct<-data.frame(fstlmcthwe$locus)
afdlmcthwe<-merge(flmct,hwelmct,by.x="flm.locus",by.y="fstlmcthwe.locus")
write.csv(afdlmcthwe,"Z:/gatksep/rnavdna/pmarinus1/dp5af0025VF/popmuscle-wo408/afdlmcthwe.csv")

#######################################
## Lake Champlain vs. Connecticut River
#######################################

flcct<-read.csv("Z:/gatksep/rnavdna/pmarinus1/dp5af0025VF/popmuscle-wo408/flcct.csv")
flcct$afd<-(abs(flcct$flc.freq1-flcct$fct.freq1)+abs(flcct$flc.freq2-flcct$fct.freq2)+abs(flcct$flc.freq3-flcct$fct.freq3)+abs(flcct$flc.freq4-flcct$fct.freq4))/2
write.csv(flcct,"Z:/gatksep/rnavdna/pmarinus1/dp5af0025VF/popmuscle-wo408/afdlcct.csv")

fstlccthwe<-read.csv("Z:/gatksep/rnavdna/pmarinus1/dp5af0025VF/popmuscle-wo408/fstlccthwe.csv")
hwelcct<-data.frame(fstlccthwe$locus)
afdlccthwe<-merge(flcct,hwelcct,by.x="flc.locus",by.y="fstlccthwe.locus")
write.csv(afdlccthwe,"Z:/gatksep/rnavdna/pmarinus1/dp5af0025VF/popmuscle-wo408/afdlccthwe.csv")
