source('skyChr_Jun_2019.R')

library(ape)
#tmp=read.csv('NRT-edited-Feb06.txt',header=T,sep='\t',stringsAsFactors=F)

# mouse 1, P2
m_label='mouse1P2'
tmp=read.csv('data/mouse_1_P2.txt',header=T,sep='\t',stringsAsFactors=F)

NRTmat=tmp[,2:ncol(tmp)]
rownames(NRTmat)=tmp[,1]

# mouse 1, SVZP1
m_label='mouse1SVZP1'
tmp=read.csv('data/mouse1_SVZP1.txt',header=T,sep='\t',stringsAsFactors=F)

NRTmat=tmp[,2:ncol(tmp)]
rownames(NRTmat)=tmp[,1]




# mouse 3
m_label='mouse3'
tmp=read.csv('data/Mouse 3_Chromosomes - 3.30.14.txt',header=T,sep='\t',stringsAsFactors=F)

NRTmat=tmp[,c(2:ncol(tmp))]
rownames(NRTmat)=tmp[,1]

# mouse 4
m_label='mouse4'
tmp=read.csv('data/Mouse 4-2014.05.27.txt',header=T,sep='\t',stringsAsFactors = F)

NRTmat=tmp[,c(2:ncol(tmp))]
rownames(NRTmat)=tmp[,1]

# mouse 5
m_label='mouse5'
tmp=read.csv('data/mouse5.txt',header=T,sep='\t',stringsAsFactors = F)

NRTmat=tmp[,c(2:ncol(tmp))]
rownames(NRTmat)=tmp[,1]

# mouse 6
m_label='mouse6'
tmp=read.table('data/mouse6.txt',header = T,sep='\t',stringsAsFactors = F,comment.char = '',na.strings = '?')

NRTmat=tmp[,c(2:ncol(tmp))]
rownames(NRTmat)=tmp[,1]



# start processing ####
NRTmat=as.matrix(NRTmat)
#vv=which(apply(NRTmat[,21:48],2,sum)>1)
#NRTmat=NRTmat[,c(1:20,vv+20)]

metaphase.types=as.character(tmp[,1])
names(metaphase.types)=as.character(tmp[,2])
##color.list=c('white','cyan','hotpink1','red','deeppink4','violet','deepskyblue1','deepskyblue3','dodgerblue1','dodgerblue3','darkseagreen1','darkolivegreen1','pink','darkolivegreen3','darkolivegreen4')
##names(color.list)=c('control',as.character(unique(metaphase.types)))
nrt.w=5
nrt.types=rownames(NRTmat)
  

if(m_label=='mouse1P2'){
  nrt=15
  
  mat.P2=computeDistMat(NRTmat,nrt.w,NRT=nrt,type='cnv',gender='F')
  mat0list=list(mat.P2)
  mat0list_label=c('P2')
  
  
} else if(m_label=='mouse1SVZP1'){
  nrt=15
  
  mat.SVZP1=computeDistMat(NRTmat,nrt.w,NRT=nrt,type='cnv',gender='F')
  mat0list=list(mat.SVZP1)
  mat0list_label=c('SVZP1')
  
  
}else if(m_label=='mouse2'){
  
  nrt=28
  
  mat.SVZR=computeDistMat(NRTmat[grep('SVZ\\(R\\)',nrt.types),],nrt.w,NRT=nrt,type='cnv',gender='F')
  mat.SVZL=computeDistMat(NRTmat[grep('SVZ\\(L\\)',nrt.types),],nrt.w,NRT=nrt,type='cnv',gender='F')
  mat.T1=computeDistMat(NRTmat[grep('T\\(1\\)',nrt.types),],nrt.w,NRT=nrt,type='cnv',gender='F')
  mat.T2=computeDistMat(NRTmat[grep('T\\(2\\)',nrt.types),],nrt.w,NRT=nrt,type='cnv',gender='F')
  #mat.new= computeDistMat(NRTmat,nrt.w,NRT=28,type='cnv')
  mat0=mat.SVZR
  mat0=mat.T2
  
  mat0list=list(mat.SVZR,mat.SVZL,mat.T1,mat.T2)
  mat0list_label=c('SVZR','SVZL','T1','T2')
}else if(m_label=='mouse2_separate_sample_SVZR'){
  nrt=8
  
  mat.SVZR=computeDistMat(NRTmat,nrt.w,NRT=nrt,type='cnv',gender='F')
  mat0list=list(mat.SVZR)
  mat0list_label=c('SVZR')
  
  
}else if(m_label=='mouse3'){
  
  nrt=19
  
  mat.SVZ=computeDistMat(NRTmat[grep('SVZ\\(R\\)',nrt.types),],nrt.w,NRT=nrt,type='cnv',gender='F')
  mat.T1=computeDistMat(NRTmat[grep('T\\(L\\)',nrt.types),],nrt.w,NRT=nrt,type='cnv',gender='F')
  mat.T2=computeDistMat(NRTmat[grep('T\\(R\\)',nrt.types),],nrt.w,NRT=nrt,type='cnv',gender='F')
  #mat.new= computeDistMat(NRTmat,nrt.w,NRT=28,type='cnv')
  mat0=mat.SVZ
  
  
  mat0list=list(mat.SVZ,mat.T1,mat.T2)
  mat0list_label=c('SVZ','T1','T2')
}else if(m_label=='mouse4'){
  nrt=23
  
  mat.SVZL=computeDistMat(NRTmat[grep('^SVZL',nrt.types),],nrt.w,NRT=nrt,type='cnv',gender='F')
  mat.SVZR=computeDistMat(NRTmat[grep('^SVZR',nrt.types),],nrt.w,NRT=nrt,type='cnv',gender='F')
  mat.T1=computeDistMat(NRTmat[grep('^T1',nrt.types),],nrt.w,NRT=nrt,type='cnv',gender='F')
  
  
  mat0list=list(mat.SVZL,mat.SVZR,mat.T1)
  mat0list_label=c('SVZL','SVZR','T1')
}else if(m_label=='mouse5'){
  nrt=9
  NRTmat[is.na(NRTmat)]=2
  
  mat.SVZRP1=computeDistMat(NRTmat[grep('SVZRP1',nrt.types),],nrt.w,NRT=nrt,type='cnv',gender='F')
  mat.tp3=computeDistMat(NRTmat[grep('TP3',nrt.types),],nrt.w,NRT=nrt,type='cnv',gender='F')
  
  mat0list=list(mat.SVZRP1,mat.tp3)
  mat0list_label=c('SVZRP1','TP3')
}else if(m_label=='mouse6'){
  
  nrt=3
  for(i in 1:4){
    
    NRTmat[is.na(NRTmat[,i]),i]=0
    
  }
  NRTmat[is.na(NRTmat[,5]),5]=2
  
  mat.SVZRP2=computeDistMat(NRTmat[grep('SVZ\\(R\\)',nrt.types),],nrt.w,NRT=nrt,type='cnv',gender='F')
  mat.SVZLP4=computeDistMat(NRTmat[grep('SVZ\\(L\\)',nrt.types),],nrt.w,NRT=nrt,type='cnv',gender='F')
  mat.T_P4=computeDistMat(NRTmat[grep('T_P4',nrt.types),],nrt.w,NRT=nrt,type='cnv',gender='F')
  mat.TP4=computeDistMat(NRTmat[grep('TP4',nrt.types),],nrt.w,NRT=nrt,type='cnv',gender='F')
  
  
  mat0list=list(mat.SVZRP2,mat.SVZLP4,mat.T_P4,mat.TP4)
  mat0list_label=c('SVZRP2','SVZLP4','T_P4','TP4')
}

# plot tree figures


for(s_i in 1:length(mat0list)){
  mat0=mat0list[[s_i]]
  sample_label=mat0list_label[s_i]
  pdf(file=paste('tree/',m_label,'_',sample_label,'.pdf',sep=''),width=8,height=12)
  
  for(i in 1:nrow(mat0)) for(j in 1:nrow(mat0))mat0[i,j]=min(mat0[i,j],mat0[j,i])
  #mat0=as.dist(mat0)
  tmp.tr=nj(mat0)
  tmp.tr.rooted=root(tmp.tr,1)
  
  # 
  tmp.tr.rooted$edge.length[tmp.tr.rooted$edge.length<0]=0
  par(mar=c(2,2,2,2))
  plot(tmp.tr.rooted)
  
  dev.off()
  
  ## export to MEGA
  
  ## export to MEGA
  mega_outfile=paste('tree/',m_label,'_',sample_label,'_mega.meg',sep='')
  tmp.dd=mat0
  
  write('#Mega',file=mega_outfile)
  write(paste('!Title ',m_label,'_',sample_label,';',sep=''),file=mega_outfile,append=T)
  
  write('',file=mega_outfile,append=T)
  rnames=gsub('#','-',rownames(mat0))
  rnames=gsub(' ','-',rnames)
  write(paste('#',rnames,sep=''),append=T,file=mega_outfile)
  for(i in 2:nrow(tmp.dd)){
    write(paste(tmp.dd[i,1:(i-1)],collapse='\t'),append=T,file=mega_outfile)
  }
}

# other codes 

if(F){
  
  ## rotation for SVZ-R
  tmp.edge=tmp.tr.rooted$edge
  tmp.edge[34:50,]=tmp.tr.rooted$edge[63:79,]
  tmp.edge[51:79,]=tmp.tr.rooted$edge[34:62,]
  
  tmp.tr.rooted$edge=tmp.edge
  
  # rotation for T1
  #tmp.edge=tmp.tr.rooted$edge
  #tmp.edge[2,]=tmp.tr.rooted$edge[55,]
  #tmp.edge[3:55,]=tmp.tr.rooted$edge[2:54,]
  #tmp.tr.rooted$edge=tmp.edge
  
  #tmp.tr.rooted1=rotate(tmp.tr.rooted,node=3)
  par(mar=c(1,0,0,0))
  tmp.lab=tmp.tr.rooted$tip.label
  tmp.tr.rooted$tip.label=gsub('SVZ\\(R\\)P2','',tmp.tr.rooted$tip.label)
  tmp.tr.rooted$tip.label= gsub('T\\(1\\)P1','',tmp.tr.rooted$tip.label)
  tmp.tr.rooted$tip.label= gsub('T\\(2\\)P2','',tmp.tr.rooted$tip.label)
  plot(tmp.tr.rooted,cex=0.8,tip.color=color.list[metaphase.types[tmp.lab]],font=2,direction='downwards',use.edge.length=T)
  ll=tmp.lab
  ll=ll[2:length(ll)]
  xlimits=par('xaxp')[1:2]
  ylimits=par('yaxp')[1:2]
  legend(xlimits[2]*0.05,ylimits[2]*0.9,legend=unique(metaphase.types[ll]),col=color.list[unique(metaphase.types[ll])],pch=19,pt.cex=1,cex=0.8)
  
  ## export to MEGA
  tmp.dd=mat0
  write(paste('#',gsub('#','-',rownames(mat0)),sep=''),file='tmp.txt')
  for(i in 2:nrow(tmp.dd)){
    write(paste(tmp.dd[i,1:(i-1)],collapse='\t'),append=T,file='tmp.txt')
  }

}
