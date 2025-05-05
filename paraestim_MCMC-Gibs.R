library(BEST)
source("BESTmcmc2.R")
environment(BESTmcmc2)=asNamespace('BEST')

# cnv0=read.delim("CNV.txt",sep="\t")
# exp0=read.delim("Expression.txt",sep="\t")
args=commandArgs(T)
cnv0=read.delim(args[1],sep="\t")
exp0=read.delim(args[2],sep="\t")
fnm0=sub(".txt","",args[1])
gnms=cnv0$genes
lnm=colnames(cnv0)
ps=match("genes",lnm)
cnv1=cnv0[,-ps]
exp1=exp0[,-ps]
ngn=nrow(cnv1)

fnm=paste0("BayesEst2-",fnm0,".txt")
sink(fnm)
cat("GenesENS\tExpLoss\tExpNor\tExpGain\tRateLoss\tRateGain\n")

# ll=c("less","normal","over")


sg1=0
sg2=1
sg3=-1
ps1=sg1+2
ps2=sg2+2
ps3=sg3+2

mat_rtover=c()
mat_rtless=c()
mat_mdover=c()
mat_mdless=c()
mat_mdnorm=c()

addvar<-function(vec){
	# vec1=round(vec,digit=2)
	varr=0.0001
	lnn=length(vec)
	if(lnn>1 && sum(vec==vec[1])==lnn){
		# vec[1]=vec[1]*(1+varr)
		pp=which.max(vec)
		vec[pp]=vec[pp]+varr
	}
	return(vec)
}

bg=1
nt=ngn
for(i in bg:nt){
	ct=cnv1[i,]
	et=exp1[i,]
	gt=et[ct==sg1]
	gp1=gt[!is.na(gt)]
	gt=et[ct==sg2]
	gp2=gt[!is.na(gt)]
	gt=et[ct==sg3]
	gp3=gt[!is.na(gt)]
	md1=median(gp1)
	md2=median(gp2)
	md3=median(gp3)
	if(length(gp1)<2 || length(gp2)<2 || length(gp3)<2 || sum(et)==0){
		ss=paste(c(gnms[i],rep("NA",5)),collapse="\t")
		cat(ss)
		cat("\n")
		next
	}else{
		gp1=addvar(gp1)
	}
	if(length(gp2)==0){
		diffover=NA
	}else{
		gp2=addvar(gp2)
		resgain <- BESTmcmc2(gp1, gp2,parallel=TRUE)
		diffover <- (resgain$mu2 - resgain$mu1)
	}
	if(length(gp3)==0){
		diffless=NA
	}else{
		gp3=addvar(gp3)
		resless <- BESTmcmc2(gp1, gp3,parallel=TRUE)
		diffless <- (resless$mu2 - resless$mu1)
	}

	rtgain <- mean(diffover > 0)
	rtless <- mean(diffless < 0)
	ss=paste(gnms[i],md3,md1,md2,rtless,rtgain,sep="\t")
	cat(ss)
	cat("\n")
	# mat_rtover=c(mat_rtover,rtgain)
	# mat_rtless=c(mat_rtless,rtless)
	# mat_mdover=c(mat_mdover,md2) ###,na.rm=T
	# mat_mdnorm=c(mat_mdnorm,md1)
	# mat_mdless=c(mat_mdless,md3)

}
sink()
# dt=data.frame(GenesENS=gnms[1:nt],ExpLoss=mat_mdless,ExpNor=mat_mdnorm,ExpGain=mat_mdover,RateLoss=mat_rtless,RateGain=mat_rtover)
# write.table(dt,file=paste0("BayesEst-",fnm,".txt"),quote=F,row.names=F,col.names=T,sep="\t")
# paste0(ll[ps2],"-vs-",ll[ps1],"_BayesEst.txt")


