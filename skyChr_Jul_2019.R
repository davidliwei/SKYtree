compareAB.cnv <- function(sA,sB,nrt.w=1,NRT=F,homo.del=1,aneu.thr=45){
	## For two samples, sA and sB (a vector with chromosome numbers), decide the minimum steps that sB (the target) coming from sA
	Nc=length(sA)
	Nl=min(length(!is.na(sA)),length(!is.na(sB)))
    weights=c(rep(1,Nc-NRT),rep(nrt.w,NRT))
    is.double=F
    tmp=sA-sB
    tmp=tmp[1:(Nc-NRT)]
    D0=sum(abs(tmp),na.rm=T)
    tmp=2*sA-sB
    tmp=tmp[1:(Nc-NRT)]
    D1=sum(abs(tmp),na.rm=T)
    tmp=sA-sB*2
    tmp=tmp[1:(Nc-NRT)]
    D2=sum(abs(tmp),na.rm=T)
	if(D2<D1&D2<D0)return(9999999)
	if(D0>D1)is.double=T
	st=Nc-NRT+1
	vv21=which(sB[st:Nc]==2&sA[st:Nc]==1)	## NRT doubling
	vv20=which(sB[st:Nc]==2&sA[st:Nc]==0)	## NRT gain and doubling
	vv10=which(sB[st:Nc]==1&sA[st:Nc]==0)	## NRT gain
	vv12=which(sB[st:Nc]==1&sA[st:Nc]==2)	## NRT random loss
	vv02=which(sB[st:Nc]==0&sA[st:Nc]==2)	## NRT total loss
	vv01=which(sB[st:Nc]==0&sA[st:Nc]==1)	## NRT loss
	if(!is.double){
		if(sum(sA[1:(st-1)],na.rm=T)>=aneu.thr){
			## both aneuploid
			Dm=D0+length(vv21)*nrt.w+length(vv20)*nrt.w*2+length(vv10)*nrt.w+length(vv12)+length(vv02)*2+length(vv01)
		}else{
			## both euploid
			Dm=D0+length(vv21)*nrt.w+length(vv20)*nrt.w*2+length(vv10)*nrt.w+length(vv12)+length(vv02)*2+length(vv01)
		}
	}else{
		## sA euploid, sB aneuploid
#         vv24=which(sA[1:(Nc-NRT)]==2&sB[1:(Nc-NRT)]==2)
#         D1=0
#         for(i in 1:Nc){
#             aa=sA[i]
#             bb=sB[i]
#             tmp=min(c(abs(bb-2*aa),abs(bb-2*(aa-1)),abs(bb-2*(aa+1))))
#             D1=D1+tmp
#         }
		#Dm=D1+length(vv20)*nrt.w+length(vv10)*nrt.w+length(vv12)+length(vv02)*2+length(vv01)+1-length(vv24)
		Dm=D1+length(vv20)*nrt.w+length(vv10)*nrt.w+length(vv12)+length(vv02)*2+length(vv01)+1
	}
	return(Dm/Nl)
	if(D0<Dm)return(D0)
	bfD=rep(0,Nc)
	names(bfD)=names(sA)
	bfD=ceiling(sB/2)-sA
    if(!NRT)nrt.id=rep(1,Nc) else nrt.id=c(rep(1,Nc-NRT),rep(0,NRT))
    vv=which(sA==0&sB!=0)
    vv.nnrt=which(nrt.id[vv]==1)
	bfD[vv[vv.nnrt]]=homo.del	## homozygous deletion restriction
	afD=sB-2*(sA+bfD)
	return((sum((abs(bfD)+abs(afD))*w,na.rm=T)+1)/Nl)
}


computeDistMat <- function(x,nrt.w=1,weights=NULL,homo.del=100,mut.del=100,NRT=F,type='cnv',MUT=0,gender='M'){
	## x: input chromosome number matrix
	## NRT: specific to Dr. Zhu's project. Either False or integer number, the number of NRT markers added to the end
	Nc=ncol(x)
	if(!NRT){
		if(type=='cnv'){if(gender=='M')x=rbind(c(rep(2,Nc-2),1,1),x) else {if(gender=='F')x=rbind(rep(2,Nc-1),2,x) else x=rbind(rep(2,Nc),x)}}
		if(type=='baf')x=rbind(rep(1,Nc),x)
		if(type=='mut'|type=='gt')x=rbind(rep(0,Nc),x)	# add a euploid sample to be the ancestry
		if(type=='both')x=rbind(c(rep(1,Nc-MUT),rep(0,MUT)),x)
	}
	else{
		if(type=='cnv'){if(gender=='M')tmp=c(rep(2,Nc-NRT-2),1,1,rep(0,NRT)) else{ if(gender=='F')tmp=c(rep(2,Nc-NRT-1),2,rep(0,NRT)) else tmp=c(rep(2,Nc-NRT),rep(0,NRT))}}
		if(type=='baf')tmp=c(rep(1,(Nc-NRT)),rep(0,NRT))
		if(type=='mut'|type=='gt')tmp=c(rep(0,Nc))
		if(type=='both')tmp=c(rep(1,Nc-MUT),rep(0,MUT))
		x=rbind(tmp,x)
	}
	rownames(x)[1]='CONTROL'
	Ns=nrow(x)
	dist.mat=matrix(NA,Ns,Ns)
	rownames(dist.mat)=colnames(dist.mat)=rownames(x)
	for(i in 1:Ns){
		cat(i,'.',sep='')
		if(type=='cnv'|type=='baf')dist.mat[i,1:Ns]=unlist(lapply(1:Ns,function(y)compareAB.cnv(x[i,],x[y,],nrt.w=nrt.w,homo.del=homo.del,NRT=NRT)))
		if(type=='mut')dist.mat[i,1:Ns]=unlist(lapply(1:Ns,function(y)compareAB.mut(x[i,],x[y,],weights=weights,mut.del=mut.del)))
		if(type=='gt')dist.mat[i,1:Ns]=unlist(lapply(1:Ns,function(y)compareAB.gt(x[i,],x[y,],weights=weights,mut.del=mut.del)))
		if(type=='both'){
			Ncnv=1:(Nc-MUT)
			Nmut=(Nc-MUT+1):Nc
			if(is.null(weights))weights=rep(1,Nc)
			dd1=unlist(lapply(1:Ns,function(y)compareAB.cnv(x[i,Ncnv],x[y,Ncnv],weights=weights[Ncnv],homo.del=homo.del,NRT=NRT)))
			dd2=unlist(lapply(1:Ns,function(y)compareAB.mut(x[i,Nmut],x[y,Nmut],weights=weights[Nmut],mut.del=mut.del)))
			dist.mat[i,1:Ns]=(dd1*length(Ncnv)+dd2*length(Nmut))/Nc
		}
	}
	return(dist.mat)
}