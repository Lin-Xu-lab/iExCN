library("org.Hs.eg.db")
args=commandArgs(T)
cutf=args[2]

dat=read.table(args[1],sep="\t",header=T,as.is=T)
n=nrow(dat)

if(cutf>0)
{
	fnm=paste0("BayesEst_CancerRelatedGene_cutf",cutf,".xls")
	ps=c()
	for(i in 1:n){
		dt=dat[i,-1]
		
		if(sum(is.na(dt))==0 && dt[2]>dt[1] && dt[3]>dt[2] && dt[5]>cutf && dt[4]>cutf){
				ps=c(ps,i)
		}

	}

	mat=dat[ps,]
}else{
	fnm=paste0("BayesEst_CancerRelatedGene_all.xls")
	mat=dat
}

# write.table(mat,file=fnm,sep="\t",col.names=T,row.names=F,quote=F)
ids=mat[,1]
id1=sub("(ENSG\\d+)\\-?(\\d*)","\\1",ids)
id2=sub("(ENSG\\d+)\\-?(\\d*)","\\2",ids)
symbols <- mapIds(org.Hs.eg.db, keys = id1, keytype = "ENSEMBL", column="SYMBOL")
gnms=c()
for(k in 1:length(id2)){
	t=id2[k]
	if(t==""){
		gnms=c(gnms,symbols[k])
	}else{
		gnms=c(gnms,paste0(symbols[k],"*",id2[k]))
	}
}
mat$Symbol=gnms
write.table(mat,file=fnm,sep="\t",col.names=T,row.names=F,quote=F)

############################################################################
################### negative gene selection ################################
############################################################################
library("org.Hs.eg.db")
args=commandArgs(T)
cutf=args[2]
fl=args[1]
# fl="BayesEst2-CNV-all_0713.xls" 
cutf=0.5

dat=read.table(fl,sep="\t",header=T,as.is=T)
n=nrow(dat)
fnm=paste0("BayesEst_CancerRelatedGene_negGenelist1-",cutf,".xls")
ps=c()
for(i in 1:n){
	dt=dat[i,-1]
	
	if(sum(is.na(dt))==0 && dt[5]<cutf && dt[4]<cutf){
			ps=c(ps,i)
	}

}
mat=dat[ps,]
ids=mat[,1]
id1=sub("(ENSG\\d+)\\-?(\\d*)","\\1",ids)
id2=sub("(ENSG\\d+)\\-?(\\d*)","\\2",ids)
symbols <- mapIds(org.Hs.eg.db, keys = id1, keytype = "ENSEMBL", column="SYMBOL")
gnms=c()
for(k in 1:length(id2)){
	t=id2[k]
	if(t==""){
		gnms=c(gnms,symbols[k])
	}else{
		gnms=c(gnms,paste0(symbols[k],"*",id2[k]))
	}
}
mat$Symbol=gnms
write.table(mat,file=fnm,sep="\t",col.names=T,row.names=F,quote=F)


flps="BayesEst_CancerRelatedGene_cutf0.9.xls" 
flng="BayesEst_CancerRelatedGene_negGenelist1-0.5.xls"
datps=read.table(flps,sep="\t",header=T,as.is=T)
datng=read.table(flng,sep="\t",header=T,as.is=T)

datexp=read.table("Expression.txt",sep="\t",header=T,row.names=1)
datval=rowSums(datexp)
nmps=datps[,1]
rnm=rownames(datexp)
vv=match(nmps,rnm)
psval=datval[vv]
nps=length(psval)

nmng=datng[,1]
vv=match(nmng,rnm)
ngval=datval[vv]
nng=length(ngval)

library("LaplacesDemon")
ds=1000
for(i in 1:1000){
	smp=sample(1:nng,nps)
	y=ngval[smp]
	KL <- KLD(psval, y)
	dt=KL$mean.sum.KLD
	if(dt<ds){
		ngot=names(y)
		ds=dt
	}
}

psval2=psval[psval<1000000]
vv=(ngval>100 & ngval<1000)
ngval_t1=sample(ngval[vv],150)
vv=(ngval>1000 & ngval<10000)
ngval_t2=ngval[vv]
vv=ngval>1 & ngval<100
ngval_t3=sample(ngval[vv],100)
vv=ngval>10000 & ngval<1000000
ngval_t4=sample(ngval[vv],69)
y=c(ngval_t1,ngval_t2,ngval_t3,ngval_t4)
# y=ngval_t1
# mx=ceiling(max(c(psval,y)))
pdf("negativeGene_distribution.pdf")
y1=log10(psval)
y2=log10(y)
pv=ks.test(y1,y1)$p.value
par(mfrow=c(2,1),xpd=T,ypd=T)
d1=density(log10(psval2))
plot(d1,main="iExCN-predicted genes",
col="darkmagenta",xlim=c(0,8))
d2=density(log10(y))
plot(d2,main="Negative Gene Expression",
col="darkgreen",xlim=c(0,8))
text(7,0.7,paste0("KS test pvalue:",pv))
# df=data.frame(Expval=c(y1,y2),
# 	Group=c(rep("Positive",length(psval2)),rep("Negative",length(y))),
# 	ylab="Log Expression")

# pv=t.test(y1,y2)$p.value
# boxplot(Expval~Group, data=df)
dev.off()
write.table(names(y),file="negativeGene.txt",sep="\n",
	quote=F,row.names=F,col.names=F)





