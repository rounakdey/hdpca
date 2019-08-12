

#n=No. of Samples, p=No. of features, m=No. of distant spikes

#d-estimator for the spikes
d.eval<-function(samp.eval,m,p,n)	
{
	gamma=p/n
	spike.est<-rep(0,m)
	for(i in 1:m)
	{
		temp=0
		for(j in (m+1):length(samp.eval))
			temp=temp+samp.eval[j]/(samp.eval[i]-samp.eval[j])
		spike.est[i]<-samp.eval[i]/(1+(gamma/(p-m))*temp)
	}
	return(spike.est)
}

#d-estimator for the angle between eigenvectors
d.angle<-function(spike.est,samp.eval,m,p,n)
{
	gamma=p/n
	angle.est<-rep(0,m)
	for(i in 1:m)
	{
		temp=0
		for(j in (m+1):length(samp.eval))
			temp=temp+samp.eval[j]/(samp.eval[i]-samp.eval[j])^2
		angle.est[i]<-sqrt(1/(1+(gamma*spike.est[i]/(p-m))*temp))
	}
	return(angle.est)
}

#psi function minus the sample eigenvalue
psi<-function(lambda,samp.eval,pop.nonsp,m,p,n)
{
	gamma=p/n
	temp<-0
	for(i in (m+1):length(which(pop.nonsp>0)))
	{
		temp<-temp+pop.nonsp[i]/(lambda-pop.nonsp[i])
	}
	temp<-1+temp*gamma/(p-m)
	return(lambda*temp-samp.eval)
}

#First derivative of the psi function
psiprime<-function(lambda,pop.nonsp,m,p,n)
{
	gamma=p/n
	temp<-0
	for(i in (m+1):length(which(pop.nonsp>0)))
	{
		temp<-temp+(pop.nonsp[i]/(lambda-pop.nonsp[i]))^2
	}
	return(1-temp*gamma/(p-m))
}


#lambda estimator for the spikes
l.eval<-function(samp.eval,pop.nonsp,m,p,n)
{
	spike.est<-rep(-100,m)	#A negative value in the final output would mean this is not a distant spike
	alpha<-uniroot(psiprime,interval=c(pop.nonsp[1]+1/p,samp.eval[1]-1/p),pop.nonsp,m,p,n)$root	#Cutoff for distant spikes
	for(i in 1:m)
	{
		psi.alpha<-psi(alpha,samp.eval[i],pop.nonsp,m,p,n)	#This is actually psi(alpha)-d
		if(psi.alpha<0)
			spike.est[i]<-uniroot(psi,interval=c(alpha,samp.eval[1]-1/p),samp.eval[i],pop.nonsp,m,p,n)$root
	}
	return(spike.est)
}

#lambda estimator for the angle between eigenvectors
l.angle<-function(spike.est,samp.eval,pop.nonsp,m,p,n)
{
	gamma=p/n
	angle.est<-rep(0,m)
	for(i in 1:m)
	{
		r2.est<-psiprime(spike.est[i],pop.nonsp,m,p,n)
		angle.est[i]<-sqrt(spike.est[i]*r2.est/samp.eval[i])
	}
	return(angle.est)
}


#OSP based estimator for the spikes
osp.eval<-function(samp.eval,m,p,n)
{
	gamma<-p/n
	a1<-0
	d<-p*samp.eval/sum(samp.eval)
	repeat{
		temp<-(d[1:m]+1-gamma)^2-4*d[1:m]
		m<-max(which(temp>=0))
		lambda<-(d[1:m]+1-gamma+sqrt((d[1:m]+1-gamma)^2-4*d[1:m]))/2
		a<-sum(lambda)+p-m
		d<-a*samp.eval/sum(samp.eval)
		if(abs(a-a1)<0.001) break
		a1<-a
	}
	nonspikes.est<-sum(samp.eval)/a
	spike.est<-lambda*nonspikes.est
	return(list(spikes=spike.est,nonspikes=nonspikes.est,n.spikes=m))
}


#OSP based estimator for the angle between eigenvectors
osp.angle<-function(spike.est,p,n)
{
	gamma=p/n
	temp1=1-gamma/(spike.est-1)^2
	temp2=1+gamma/(spike.est-1)
	return(sqrt(temp1/temp2))
}


#All functions used in Karoui's algorithm
#Inverse Steiltjes (Inverse w.r.t. z)
invsteiltjes<-function(v,z,samp.scaled,m,n)
{
	rep<-1
	repeat{
	vs=sum(1/(samp.scaled-z))/(n-m)
	vs<-c(Re(vs),Im(vs))
	v1<-Re(v)
	v2<-Im(v)
	J<-matrix(0,2,2)
	J[1,1]<-sum(((samp.scaled-Re(z))^2-Im(z)^2)/((samp.scaled-Re(z))^2+Im(z)^2)^2)/(n-m)
	J[2,2]<-J[1,1]
	J[1,2]<--sum(2*Im(z)*(samp.scaled-Re(z))/((samp.scaled-Re(z))^2+Im(z)^2)^2)/(n-m)
	J[2,1]<--J[1,2]
	znew<-c(Re(z),Im(z))-solve(J)%*%(vs-c(v1,v2))
	znew<-complex(real=znew[1],imaginary=znew[2])
	if(norm(as.matrix(vs-c(v1,v2)),type='f')<0.00001 || rep>1000) break
	z<-znew
	rep=rep+1
	}
	if(rep>1000)
	{
		return(1)
	} else {
		return(znew)
	}
}

#Produce quantiles as non-spikes from the estimated LSD
quant<-function(q,t,cuml)
{
	if (q>1 || q<0)
	{
		return(NULL)
	} else {
		for(i in 1:length(t))
		{
			if (q<=cuml[i]) break
		}
		if (q==cuml[i])
		{
			for(j in (i+1):length(t))
			{
				if(cuml[j]>cuml[i]) break
			}
			if(i==j-1)
			{
				return(t[i])
			} else {
				return(runif(1,t[i],t[j-1]))
			}
		} else {
			return(t[i])
		}
	}
}

#Loss function
lossf<-function(pop,eimax,z,v,gamma)
{
	##Loss function
	p<-length(pop)
	t<-pop[which(pop>0)]/eimax
	w<-rep(1/p,length(which(pop>0)))
	e<-rep(0,p)
	for (i in 1:length(z))
	{
		term1<-sum(w*t/(1+v[i]*t))
		e[i]=1/v[i]+z[i]-gamma*term1
	}
	L=max(abs(Re(e)),abs(Im(e))) #L is given by this
	return(L)
}

#Smoothing the LSD
popsmooth<-function(pop,t1,z,v,gamma,eimax,smooth)
{	
	p<-length(pop)
	zbak<-z
	vbak<-v
	loss1<-1
	loss2<-1
	if(smooth==TRUE){
	selz<-sample(1:length(z),200)	#Thinning - otherwise takes a lot of time
	z<-z[selz]
	v<-v[selz]
	bwinc<-round(p/length(t1),0)
	bw<-bwinc
	loss1=lossf(pop,eimax,z,v,gamma)
	pop2<-ksmooth(1:length(pop),pop,bandwidth=bw)$y
	loss2=lossf(pop2,eimax,z,v,gamma)
	rep<-1
	repeat{
		bw<-bw+bwinc
		pop2<-ksmooth(1:length(pop),pop,bandwidth=bw)$y
		loss3=lossf(pop2,eimax,z,v,gamma)	
		if(loss3>loss2*(1+10^-7) || rep>10) break
		loss2<-loss3
		rep<-rep+1
	}
	bw<-bw-bwinc
	if(loss2<loss1) pop<-ksmooth(1:length(pop),pop,bandwidth=bw)$y
	}
	loss=lossf(pop,eimax,zbak,vbak,gamma)	
	return(c(pop,loss))
}

findv2hi<-function(samp.scaled,lo,hi,n,m)
{
	mid=floor((lo+hi)/2)
	v<-complex(real=0,imaginary=mid)
	z<-complex(real=mean(samp.scaled),imaginary=1/mid)

	flag.err<-1
	try({
	z<-invsteiltjes(v,z,samp.scaled,m,n)
	flag.err<-0
	},silent=TRUE)
	if(flag.err==1)
	{
		return(findv2hi(samp.scaled,lo,mid,n,m))
	} else {
		if(mid==lo) return(mid)
		return(findv2hi(samp.scaled,mid,hi,n,m))
	}

}

karoui.nonsp<-function(samp.eval,m,p,n,smooth=TRUE)
{
gamma<-p/n
eimax=samp.eval[m+1]	#Largest sample eigenvalue

samp.scaled=samp.eval[(m+1):n]/eimax
rev=seq(0,1,length.out=50)

v2hi<-findv2hi(samp.scaled,0.0001,100,n,m)
imv=seq(0.0001,v2hi,length.out=50)
v=matrix(0,length(rev),length(imv))
z<-matrix(mean(samp.scaled)+1i,length(rev),length(imv))
for(i in 1:length(rev))
{
	for(j in 1:length(imv))
	{
		v[i,j]=complex(real=rev[i],imaginary=imv[j])
		if(i==1)
		{
			z[i,j]<-complex(real=mean(samp.scaled),imaginary=1/Im(v[i,j]))
		}
		if(i==2 && j>1) 		z[i,j]<-z[i,(j-1)]
		if(i>2) z[i,j]<-z[(i-1),j]
		flag.err<-1
		try({
		z[i,j]<-invsteiltjes(v[i,j],z[i,j],samp.scaled,m,n)
		flag.err<-0
		},silent=TRUE)
		if(flag.err==1)
		{
		try({
		z[i,j]<-complex(real=mean(samp.scaled),imaginary=1/Im(v[i,j]))
		z[i,j]<-invsteiltjes(v[i,j],z[i,j],samp.scaled,m,n)
		flag.err<-0
		},silent=TRUE)
		}
		if(flag.err==1)
		{
			z[i,j]<-NA
		}
	}
}
z<-na.omit(as.vector(z))
isImz0<-which(Im(z)==0)
if(length(isImz0)>0) z<-z[-isImz0]
if(length(z)<0.5*nrow(v)*ncol(v)) stop("The inverse stieltjes transform does not converge in Karoui's algorithm")
zadd<-complex(real=c(seq(1.001,4,length.out=600)),imaginary=c(1e-2,1e-3))
z<-c(z,zadd)
dist<-NULL
for(i in 1:length(z))	dist<-c(dist,min(abs(Re(z[i])-samp.scaled)))
infd<-which(dist<10^-8)
if(length(infd)>0)	z<-z[-infd]

#Recalculating v from z as the Steiljes transform
v<-rep(0,length(z))
for(i in 1:length(z))
{
	v[i]=sum(1/(samp.scaled-z[i]))/(n-m) #These are the Steiljes transforms
}

print("Stieltjes transformations calculated")


#Linear Programming problem
#t<-c(seq(0,0.9,length.out=100),seq(0.9,1,length.out=300))
t<-seq(0,1,length.out=200)
ai<-seq(0,1-2^(-7),2^(-7))
bi<-seq(0+2^(-7),1,2^(-7))
a=c(rep(0,(length(t)+length(ai))),1)
A1<-matrix(NA,length(z)*2,length(a))
A2<-matrix(NA,length(z)*2,length(a))
A3=c(rep(1,(length(t)+length(ai))),0)
b3=1
b1=NULL
b2=NULL
for(i in 1:length(z))
{
z1=z[i]
v2=v[i] #This is the Steiljes transform for z1

term1<-Re(v2)/(Mod(v2)^2)+Re(z1)
term2<- -Im(v2)/(Mod(v2)^2)+Im(z1)
a1<-(p/sum(n))*t/(1+t*v2)
a2<-(p/sum(n))*(1/v2-log((1+v2*bi)/(1+v2*ai))/((bi-ai)*(v2^2)))
a11<-Re(a1)
a12<-Re(a2)
a21<--Im(a1)
a22<--Im(a2)
A1[(i-1)*2+1,]=c(a11,a12,-1)
A1[(i-1)*2+2,]=c(a21,a22,-1)
A2[(i-1)*2+1,]=c(a11,a12,1)
A2[(i-1)*2+2,]=c(a21,a22,1)
b1=c(b1,term1,-term2)
b2=c(b2,term1,-term2)
if(term1<0 || term2>0) break
}
sign<-c(rep("<=",nrow(A1)),rep(">=",nrow(A2)),rep("=",1))
flag.err<-0
try({
simp<-lpSolve::lp(direction="min",a,rbind(A1,A2,A3),sign,c(b1,b2,b3))
flag.err<-1
})

if(flag.err==0)	stop("The solution for the linear programming problem is not available in Karoui's algorithm")
if(simp$status!=0)
{
	flag.err<-0
	try({
	simp=boot::simplex(a,A1,b1,A2,b2,A3,b3)
	flag.err<-2
	})
	if(flag.err==0)	stop("The solution for the linear programming problem is not available in Karoui's algorithm")
	flag.err<-ifelse(simp$solved==1,2,0)
}
if(flag.err==0)	stop("The solution for the linear programming problem is not available in Karoui's algorithm")

print("Spectrum estimated using linear programming")
oldw <- getOption("warn")
options(warn = -1)
t1<-c(t,ai)
t1<-t1[order(t1)]
if(flag.err==1)
{
	w1=simp$solution[1:length(t)]
	w2=simp$solution[(length(t)+1):(length(t)+length(ai))]
} else {
	w1=simp$soln[1:length(t)]
	w2=simp$soln[(length(t)+1):(length(t)+length(ai))]
}
cuml=rep(0,length(t1))			#Distribution function (Sample)
for(j in 1:length(t1))
{
	a<-sum(w1[which(t<=t1[j])])
	b<-sum(w2[which(bi<=t1[j])])
	if (max(which(bi<=t1[j]))==length(bi) || max(which(bi<=t1[j]))==-Inf)
	{
		d<-0
	} else {
		c<-max(which(bi<=t1[j]))+1
		d<-w2[c]*(t1[j]-ai[c])/(bi[c]-ai[c])
	}
	cuml[j]=a+b+d
}
options(warn = oldw)

#Obtain Quantiles
est<-rep(0,p-m)
for(i in 1:(p-m))
{
	est[i]<-quant((i-1)/(p-m),t1,cuml)
}
est<-est*eimax
pop=rev(est)

#Smoothing the non-spikes

pop<-popsmooth(pop,t1,z,v,gamma,eimax,smooth)
if(smooth==TRUE)	print("Spectrum smoothing completed")


pop<-c(rep(pop[1],m),pop)

return(pop)
}


select.nspike<-function(samp.eval,p,n,n.spikes.max,evals.out=FALSE,smooth=TRUE)
{
	m<-n.spikes.max
	samp.eval<-samp.eval[order(samp.eval,decreasing=TRUE)]

	if(length(samp.eval)<n-1 || length(samp.eval)>n)
	{
		stop("samp.eval must have length n or (n-1)")
	} else {
		if(length(samp.eval)==n-1)	n=n-1
		if(length(samp.eval)==n & abs(samp.eval[n]-samp.eval[n-1])<=10^-5*samp.eval[n-1])
		{
			samp.eval<-samp.eval[-n]
			n<-n-1
		}
	}
repeat
{
pop.nonsp<-karoui.nonsp(samp.eval,m,p,n,smooth)
loss<-pop.nonsp[p+1]
pop.nonsp<-pop.nonsp[-c(p+1)]
eval.l<-l.eval(samp.eval,pop.nonsp,m,p,n)

if(min(eval.l)>0)
{
	break
} else {
	m=min(which(eval.l==min(eval.l)))-1
	print(paste("Restarting algorithm with",m,"spikes"))
}
if(m==0)	stop("No distant spike was found")
}
if(evals.out)
{
	return(list(n.spikes=m,spikes=eval.l,nonspikes=pop.nonsp[(m+1):p],loss=loss))
} else {
	return(list(n.spikes=m))
}
}


hdpc_est<-function(samp.eval,p,n,method=c("d.gsp","l.gsp","osp"),n.spikes,n.spikes.max,n.spikes.out,nonspikes.out=FALSE,smooth=TRUE)
{
	method<-match.arg(method)
	
	samp.eval<-samp.eval[order(samp.eval,decreasing=TRUE)]
	if(length(samp.eval)<n-1 || length(samp.eval)>n)
	{
		stop("samp.eval must have length n or (n-1)")
	} else {
		if(length(samp.eval)==n-1)	n=n-1
		if(length(samp.eval)==n & abs(samp.eval[n]-samp.eval[n-1])<=10^-5*samp.eval[n-1])
		{
			samp.eval<-samp.eval[-n]
			n<-n-1
		}
	}

	flag.sp.lambda<-0
	if(missing(n.spikes))
	{
		if(missing(n.spikes.max))
		{
			stop("Both n.spikes and n.spikes.max cannot be missing")
		} else {
			spike.select<-select.nspike(samp.eval,p,n,n.spikes.max,evals.out=TRUE,smooth=smooth)
			n.spikes<-spike.select$n.spikes
			eval.l<-spike.select$spikes
			pop.nonsp<-spike.select$nonspikes
			loss<-spike.select$loss
			pop.nonsp<-c(rep(pop.nonsp[1],n.spikes),pop.nonsp)
			angle.l<-l.angle(eval.l,samp.eval,pop.nonsp,n.spikes,p,n)
			corr.l<-angle.l*sqrt(samp.eval[1:n.spikes]/eval.l)
			shrink.l<-eval.l/samp.eval[1:n.spikes]
			flag.sp.lambda<-1
		}
	}
	if(method=="l.gsp")
	{
		if(missing(n.spikes.out))	n.spikes.out<-n.spikes
		if(n.spikes.out>n.spikes)	stop("n.spikes.out must be smaller than n.spikes")
		if(flag.sp.lambda==0)
		{
			pop.nonsp<-karoui.nonsp(samp.eval,n.spikes,p,n,smooth=smooth)
			loss<-pop.nonsp[p+1]
			pop.nonsp<-pop.nonsp[-c(p+1)]
			eval.l<-l.eval(samp.eval,pop.nonsp,n.spikes,p,n)
			if(min(eval.l)<0)	stop("n.spikes is too large, the spectrum does not have that many distant spikes")
			angle.l<-l.angle(eval.l,samp.eval,pop.nonsp,n.spikes,p,n)
			corr.l<-angle.l*sqrt(samp.eval[1:n.spikes]/eval.l)
			shrink.l<-eval.l/samp.eval[1:n.spikes]
		}
		if(nonspikes.out)
		{
			return(list(spikes=eval.l[1:n.spikes.out],n.spikes=n.spikes,
			angles=angle.l[1:n.spikes.out],correlations=corr.l[1:n.spikes.out],
			shrinkage=shrink.l[1:n.spikes.out],loss=loss,nonspikes=pop.nonsp[(n.spikes+1):p]))
		} else {
			return(list(spikes=eval.l[1:n.spikes.out],n.spikes=n.spikes,
			angles=angle.l[1:n.spikes.out],correlations=corr.l[1:n.spikes.out],
			shrinkage=shrink.l[1:n.spikes.out],loss=loss))
		}
	} else if(method=="d.gsp") {
		if(missing(n.spikes.out))	n.spikes.out<-n.spikes
		if(n.spikes.out>n.spikes)	stop("n.spikes.out must be smaller than n.spikes")
		eval.d<-d.eval(samp.eval,n.spikes,p,n)
		angle.d<-d.angle(eval.d,samp.eval,n.spikes,p,n)
		corr.d<-angle.d*sqrt(samp.eval[1:n.spikes]/eval.d)
		shrink.d<-eval.d/samp.eval[1:n.spikes]
		return(list(spikes=eval.d[1:n.spikes.out],n.spikes=n.spikes,
		angles=angle.d[1:n.spikes.out],correlations=corr.d[1:n.spikes.out],
		shrinkage=shrink.d[1:n.spikes.out]))
	} else if(method=="osp") {
		osp.out<-osp.eval(samp.eval,n.spikes,p,n)
		eval.osp<-osp.out$spikes
		pop.nonsp<-osp.out$nonspikes
		n.spikes<-osp.out$n.spikes
		osp.scaled<-eval.osp/pop.nonsp
		if(missing(n.spikes.out))	n.spikes.out<-n.spikes
		if(n.spikes.out>n.spikes)	stop("n.spikes.out must be smaller than n.spikes")
		angle.osp<-osp.angle(osp.scaled,p,n)
		corr.osp<-angle.osp*sqrt(samp.eval[1:n.spikes]/eval.osp)
		shrink.osp<-(osp.scaled-1)/(osp.scaled+p/n-1)
		if(nonspikes.out)
		{
			return(list(spikes=eval.osp[1:n.spikes.out],n.spikes=n.spikes,
			angles=angle.osp[1:n.spikes.out],correlations=corr.osp[1:n.spikes.out],
			shrinkage=shrink.osp[1:n.spikes.out],nonspikes=pop.nonsp))
		} else {
			return(list(spikes=eval.osp[1:n.spikes.out],n.spikes=n.spikes,
			angles=angle.osp[1:n.spikes.out],correlations=corr.osp[1:n.spikes.out],
			shrinkage=shrink.osp[1:n.spikes.out]))
		}
	}
}

pc_adjust<-function(train.eval,p,n,test.scores,method=c("d.gsp","l.gsp","osp"),n.spikes,n.spikes.max,smooth=TRUE)
{
	test.scores<-as.matrix(test.scores)
	if(!missing(n.spikes))
	{
		shrinkage<-hdpc_est(train.eval,p,n,method=method,n.spikes=n.spikes,smooth=smooth)$shrinkage
	} else {
		shrinkage<-hdpc_est(train.eval,p,n,method=method,n.spikes.max=n.spikes.max,smooth=smooth)$shrinkage
	}
	m<-min(length(shrinkage),ncol(test.scores))
	for(i in 1:m)
		test.scores[,i]<-test.scores[,i]/shrinkage[i]
	print(paste("Top ",m," PC scores have been adjusted",sep=""))
	return(test.scores)
}

