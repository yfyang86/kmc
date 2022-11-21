library(survival)

LL= 50
beta0 <- 3.218
beta1 <- -0.0145

stanford5 <- stanford2[!is.na(stanford2$t5), ]

beta.grid <- function(x0,range,n0,type="sq",u=5){
	n0 = as.double(n0)
	if (type=="sq"){
		o1 <- c(
		-range*(u*(n0:1)^2)/(u*n0^2),0,
		range*(u*(1:n0)^2)/(u*n0^2)
		)
	}else{
	if (type=='sqrt'){
		o1 <- c(
		-range*(u*sqrt(n0:1))/(u*sqrt(n0)),0,
		range*(u*sqrt(1:n0))/(u*sqrt(n0)))
		}else{
		o1=c(
		-range*(n0:1)/n0,
		0,
		range*(1:n0)/n0
		)
		}  
	}
	return(
		x0+o1
		);
}

beta.0 <- beta.grid(beta0, 0.05, LL,"l")
beta.1 <- beta.grid(beta1,.003,LL,"l")

set.seed(1234)

y=log10(stanford5$time)+runif(152)/1000

d <- stanford5$status

oy = order(y,-d)
d=d[oy]
y=y[oy]
x=cbind(1,stanford5$age)[oy,]

ZZ=matrix(0,2*LL+1,2*LL+1)

library(kmc)
tic=0
for(jj in 1:(2*LL+1)){
for(ii in 1:(2*LL+1)){
  beta=c(beta.0[ii],beta.1[jj])
  ZZ[jj,ii]=kmc.bjtest(y,d,x=x,beta=beta,init.st="naive")$"-2LLR"
}
}
ZZ2<-ZZ
ZZ[ZZ<0]=NA ## when KMC.BJTEST fails to converge, it'll return a negative value.
ZZ[ZZ>0.5]=NA

range(ZZ,finite=T) -> zlim
floor.d<-function(x,n=4){floor(x*10^n)/(10^n)}

#postscript("fig2_1.eps",width=7,height=7)
contour(
  y=beta.0,
  x=beta.1,
  ZZ,
  zlim=c(0,1),
  levels=unique(floor.d(
		beta.grid(x0=mean(zlim),range=diff(zlim)/2, n0=15,type="sqrt",u=10),
		4)),
  ylab="Intercept",
  xlab=expression(beta[Age])
  ) 

#dev.off()
