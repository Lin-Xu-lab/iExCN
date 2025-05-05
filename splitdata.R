# cnv0=read.delim("CNV.txt",sep="\t")
# exp0=read.delim("Expression.txt",sep="\t")
args=commandArgs(T)
cnv00=read.delim("CNV.txt",sep="\t",header=T)
exp00=read.delim("Expression.txt",sep="\t",header=T)
ngn=nrow(cnv00)

gnms=cnv00[,1]
cnv0=cnv00[,-1]
exp0=exp00[,-1]

gls=as.character(gnms)
lvs=levels(gnms)
for(k in 1:length(lvs)){
	vv=gls%in%lvs[k]
	if(sum(vv)>1){
		ps=which(vv)
		for(j in 2:sum(vv))
		   gls[ps[j]]=paste0(lvs[k],"-",j)
	}
}
cnv0$genes=gls
exp0$genes=gls


dv=args[1]
bg=seq(1,ngn,dv)
t1=bg-1
ed=c(t1[-1],ngn)
# rg=cbind(bg,ed)

nn=length(bg)
for(i in 1:nn){
	nm1=paste0("CNV-",i,".txt")
	tab1=cnv0[bg[i]:ed[i],]
	write.table(tab1,file=nm1,sep="\t",quote=F,row.names=F,col.names=T)
	nm2=paste0("Expression-",i,".txt")
	tab2=exp0[bg[i]:ed[i],]
	write.table(tab2,file=nm2,sep="\t",quote=F,row.names=F,col.names=T)
	
}

