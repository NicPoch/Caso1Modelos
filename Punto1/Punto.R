#-------------Autores------------------------------------------
#Juan Fernando Castañeda
#Nicolás Poch
#------------Paquetes------------------------------------------
library(markovchain)
library(expm)
#-----------Definir Variables----------------------------------
gamma<-5
p<-0.5
mu<-5
theta<-5
k<-3
m<-3
r<-0.5
N<-5
#----------------Definir Estados------------------------
S_p<-0:N
S_s<-0:N
S_e<-c()
for(x in S_p)
{
  for(y in S_s)
  {
    if((x+y)<=N)
    {
      S_e<-c(S_e,paste(x,y,sep=","))
    }
  }
}
#----------------------Generar Matriz Q-----------------
matrizQ<-matrix(0,ncol=length(S_e),nrow=length(S_e))
dimnames(matrizQ)<-list(S_e,S_e)
for(n in S_e)
{
  i<-as.numeric(strsplit(n,split = ",")[[1]][1])
  j<-as.numeric(strsplit(n,split = ",")[[1]][2])
  for(n_1 in S_e)
  {
    iP<-as.numeric(strsplit(n_1,split = ",")[[1]][1])
    jP<-as.numeric(strsplit(n_1,split = ",")[[1]][2])
    if((iP==i+1)&(j==jP))
    {
      matrizQ[n,n_1]<-p*(1/gamma)
    }
    else if((iP==i-1)&(jP==j+1))
    {
      matrizQ[n,n_1]<-mu*min(c(i,k))*theta
    }
    else if((iP==i+1)&(jP==j-1))
    {
      matrizQ[n,n_1]<-r*theta*min(c(m,j))
    }
    else if((iP==i)&(jP==j-1))
    {
      matrizQ[n,n_1]<-(1-r)*theta*min(c(m,j))
    }
    else
    {
      matrizQ[n,n_1]<-0
    }
  }
}
for(S in S_e)
{
  matrizQ[S,S]<-(-1)*sum(matrizQ[S,])
}
View(matrizQ)
#---------------Generar CMTC-------------------------------
CMTC<-new("ctmc",states=as.character(S_e),generator=matrizQ)