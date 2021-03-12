#-------------Autores-----------------------------------
#Juan Fernando castañeda
#Nicolás Poch
#-----------Paquetes------------------------------------
library(markovchain)
library(fitdistrplus)
library(expm)
#---------------Variables-------------------------------
lambda<-14
M<-20
h<-10
#-----------------Generar Matriz P----------------------
S_x<-0:M
matrizP<-matrix(0,ncol=length(S_x),nrow=length(S_x))
for(i in S_x)
{
  for(j in S_x)
  {
    if((i<=h) & (j<M))
    {
      matrizP[i+1,j+1]<-dpois(j,lambda)
    }
    else if((i<=h) & (j==M))
    {
      matrizP[i+1,j+1]<-ppois(j-1,lambda,lower.tail = F)
    }
    else if((i>h) & (j<M) & (j>=(i-h)))
    {
      matrizP[i+1,j+1]<-dpois(j-(i-h),lambda)
    }
    else if((i>h) & (j==M) & (j>=(i-h)))
    {
      matrizP[i+1,j+1]<-ppois(j-(i-h)-1,lambda,lower.tail = F)
    }
    else
    {
      matrizP[i+1,j+1]<-0
    }
  }
}
dimnames(matrizP)<-list(S_x,S_x)
View(matrizP)
#---------------------Generar CMTD---------------------------------------
CMTD<-new("markovchain",states=as.character(S_x),transitionMatrix=matrizP)
#Valor esperado de lotes
p_estables<-steadyStates(CMTD)
EL<-(p_estables%*%matrix(S_x,ncol = 1,nrow=length(S_x)))[1,1]
#Valor esperado de lotes que se lleva el camión
print("El camión recoje:")
#hay menos que la capacidad máxima del camión
if(EL<h)
{
  print(EL)
}
#hay más que la capacidad máxima del camión
if(EL>=h)
{
  print(h)
}
#-----------------Perdida de imagen---------------------------------------
beta=1
alpha=rep(0,M+1)
alpha[h+1]<-1
alpha<-matrix(alpha,ncol=M+1,nrow=1)
colnames(alpha)<-S_x
puntosPerdidos<-0
for(n in 0:4)
{
  pi_n<-alpha%*%(matrizP%^%n)
  EL_n<-(pi_n%*%matrix(S_x,ncol = 1,nrow=length(S_x)))[1,1]
  #si el valor esperado al final de la n-esima semana es mayor a h se suma un punto de perdida
  if(EL_n>h)
  {
    puntosPerdidos<-puntosPerdidos+beta
  }
}
print("Puntos perdidos en 4 semanas:")
print(puntosPerdidos)