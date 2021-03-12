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
#--------------Punto B-------------------------------------
a<-2
b<-2
alpha<-matrix(0,ncol=length(S_e),nrow=1)
colnames(alpha)<-S_e
index<-paste(a,b,sep=",")
alpha[1,index]<-1
t<-5
pi_t<-alpha%*%expm((matrizQ*t))
cantidades<-matrix(0,nrow=length(S_e),ncol=1)
row.names(cantidades)<-S_e
for(S in S_e)
{
  cantidades[S,1]<-as.numeric(strsplit(S,split=",")[[1]][1])+as.numeric(strsplit(S,split=",")[[1]][2])
}
ET<-(pi_t%*%cantidades)[1,1]
print(ET)
#-----------------Punto C---------------------------------
gamma<-5
p<-0.35
mu<-2
theta<-1.5
k<-2
m<-2
r<-0.2
N<-8
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