



rm(list=ls())

source('skyChr_Jul_2019.R')

library(ape)



# mouse 2
m_label='mouse2'
#tmp=read.csv('../../Mouse 2 Chromosomes - 1.27.14.txt',header=T,sep='\t',stringsAsFactors=F)

if(TRUE){
  # dataset 1-3
  tmp=read.csv('data/NRT-edited-Feb06.txt',header=T,sep='\t',stringsAsFactors=F)
  #tmp=read.csv('OtherInputData/Mouse 2 Chromosomes - 1.27.14-edited.txt',header=T,sep='\t',stringsAsFactors=F)
  #tmp=read.csv('OtherInputData/Mouse 2 Chromosomes - 1.27.14.txt',header=T,sep='\t',stringsAsFactors=F)
  
  
  
  NRTmat=tmp[,3:ncol(tmp)]
  rownames(NRTmat)=tmp[,2]
  
  
  nrt=28
  
}



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


if(m_label=='mouse2'){
  
  mat.SVZR=computeDistMat(NRTmat[grep('SVZ\\(R\\)',nrt.types),],nrt.w,NRT=nrt,type='cnv',gender='M')
  mat.SVZL=computeDistMat(NRTmat[grep('SVZ\\(L\\)',nrt.types),],nrt.w,NRT=nrt,type='cnv',gender='M')
  mat.T1=computeDistMat(NRTmat[grep('T\\(1\\)',nrt.types),],nrt.w,NRT=nrt,type='cnv',gender='M')
  mat.T2=computeDistMat(NRTmat[grep('T\\(2\\)',nrt.types),],nrt.w,NRT=nrt,type='cnv',gender='M')
  #mat.new= computeDistMat(NRTmat,nrt.w,NRT=28,type='cnv')
  mat0=mat.SVZR
  mat0=mat.T2
  
  mat0list=list(mat.SVZR,mat.SVZL,mat.T1,mat.T2)
  mat0list_label=c('SVZR','SVZL','T1','T2')
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
