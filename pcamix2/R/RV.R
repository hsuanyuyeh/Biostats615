# @name RV
# @title  RV coefficient
# @description Computes the  RV coefficient between matrices.
# @param liste.mat a list of k matrices.
# @return  a (k,k) matrix containing the RV coefficients.
# @examples 
# V0<-c("a","b","a","a","b")
# V01<-c("c","d","e","c","e")
# V1<-c(5,4,2,3,6)
# V2<-c(8,15,4,6,5)
# V3<-c(4,12,5,8,7)
# V4<-c("vert","vert","jaune","rouge","jaune")
# V5<-c("grand","moyen","moyen","petit","grand")
# G1<-data.frame(V0,V01,V1)
# G2<-data.frame(V2,V3)
# G3<-data.frame(V4,V5)
# liste.mat<-list(G1,G2,G3)
# RV(liste.mat)
# @keywords internal 
RV<-function(liste.mat){
  Lg<-Lg(liste.mat)
  RV<-matrix(NA,ncol=ncol(Lg),nrow=nrow(Lg))
  for (i in 1: nrow(Lg)){
    for (j in 1:ncol(Lg)){
      RV[i,j]<-Lg[i,j]/sqrt(Lg[i,i]*Lg[j,j])
    }
  }
  rownames(RV)<-colnames(RV)<-rownames(Lg)
  return(RV)
}