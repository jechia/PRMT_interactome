library(RMySQL)
library(stringr)
library(RColorBrewer)
library(ggplot2)
library(reshape2)
library(biomaRt)
library(future.apply)
library(pheatmap)
library(GeneOverlap)
library(venn)
library(Cairo)
library(plyr)

library(DBI)
library(pool)

#data breaks according to the quantile
quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

#read miso result
setwd("~/whh/miso_summary")
sample.l<-c('18R127','18R128','18R129','18R130','18R131','18R132','18R133','18R134','18R135','18R136','18R137','18R138')
sample_name.l<-c('Control_1','PRMT1_1','PRMT1_2','PRMT3_1','PRMT3_2','PRMT4_1','PRMT4_2','PRMT5_1','PRMT5_2','PRMT6_1','PRMT7_1','PRMT7_2')
event.l<-c('A3SS','A5SS','MXE','RI','SE')
sample.n<-c('PRMT1','PRMT3','PRMT4','PRMT5','PRMT6','PRMT7')
rep.l<-c(1,1,1,1,0,1)

#read miso summary file set path 
read_miso_summary<-function(sample_list,event_list,sample_name_list,replicate=T){
  all_data.l<-sapply(1:length(sample_list),function(i) {
    data<-future_sapply(1:length(event_list),function(j) {
      file_name<-paste(sample_list[i],event_list[j],sep='_')
      print(file_name)
      event_result<-read.delim(file_name,header=T)
      event.m<-matrix(unlist(apply(event_result,1,filter_miso,sample_name=sample_name_list[i],event_name=event_list[j])),ncol=4,byrow=T)
      colnames(event.m)<-c("sample_name","event_type","event_name","psi")
      return(event.m)
    })
    df<-Reduce(rbind,data)
    return(df)
  })
  all_data.df<-Reduce(rbind,all_data.l)
  all_data.df<-as.data.frame(all_data.df)
  all_data.df$psi<-as.numeric(as.character(all_data.df$psi))
  #filter replicate
  if(replicate){
    exp_design<-matrix(unlist(lapply(as.character(all_data.df$sample_name),function(x) unlist(strsplit(x,split = '_')))),ncol=2,byrow=T)
    new_data.df<-data.frame(exp_name=exp_design[,1],rep=exp_design[,2],subset(all_data.df,select = -sample_name))
    new_data.l<-lapply(levels(new_data.df$exp_name),filter_rep,data_frame=new_data.df)
    data.df<-Reduce(rbind,new_data.l)
  }else{
    data.df<-all_data.df
    colnames(data.df)<-c('exp_name','event_type','event_name','psi')
  }
  return(data.df)
}
#replicate filter
filter_rep<-function(x,data_frame){
  df=subset(data_frame,exp_name==x)
  if(length(unique(df$rep))>1){
    rep_event.l<-lapply(levels(df$rep),function(y,dt_frame) subset(dt_frame,rep==y)$event_name , dt_frame=df)
    rep_event<-Reduce(intersect,rep_event.l)
    rep_data.df<-subset(df,event_name%in%rep_event)
    Var<-aggregate( psi ~ .,data=subset(rep_data.df,select=-rep),FUN=var)
    idx<-which(Var$psi<0.05)
    re<-aggregate( psi ~ .,subset(rep_data.df[-idx,],select=-rep),mean)
  }
  else re=subset(df,select=-rep)
  return(re)
}
#differential event filter
filter_diff<-function(x,data_frame,ct,cut_off) {
  event.idx<-intersect(subset(data_frame,exp_name==x)$event_name,ct)
  diff.l<-future_sapply(event.idx,function(e,dt_frame,s_name){
    as.numeric(subset(dt_frame,exp_name==s_name&event_name==e)$psi)-as.numeric(subset(dt_frame,exp_name=='Control'&event_name==e)$psi)
  },dt_frame=data_frame,s_name=x)
  idx<-which(abs(diff.l)>cut_off)
  return(event.idx[idx])
}
#filter miso confidence interval w/o read count
filter_miso<-function(psi.l,sample_name,event_name) {
  psi<-as.vector(psi.l)
  if(as.numeric(psi[4])-as.numeric(psi[3])<=0.3)
    return(c(sample_name,event_name,psi[1],psi[2]))
  else
    return(NULL)
}

#filter read count
filter_miso<-function(psi.l,sample_name,event_name) {
  psi<-as.vector(psi.l)
  count.s<-sum(as.numeric(sapply(unlist(strsplit(as.character(psi[7]),split = ",",fixed = T)),function(x) strsplit(x,split=":")[[1]][2])))
  if((as.numeric(psi[4])-as.numeric(psi[3])<=0.3)&count.s>50)
    return(c(sample_name,event_name,psi[1],psi[2]))
  else
    return(NULL)
}

data.df<-read_miso_summary(sample_list = sample.l,event_list = event.l,sample_name_list = sample_name.l)
#find differential event list for each sample
control<-subset(data.df,exp_name=="Control")$event_name
diff_event.l<-lapply(levels(data.df$exp_name)[-1],filter_diff,data_frame=data.df,ct=control,cut_off=0.1)

inter_event<-Reduce(intersect,diff_event.l)
inter.df=subset(data.df,event_name%in%inter_event)

#venn plot
#venn(diff_event.l, ilabels = TRUE,zcolor = "style", snames="PRMT1,PRMT3,PRMT4,PRMT5,PRMT6,PRMT7",size = 25, cexil = 2, cexsn = 2.5);

#annotation
annotate_event<-function(e.m,dt.df,inter=F,connection){
  if(inter){
    e.n<-e.m
    event.n<-unlist(strsplit(e.n,split=":"))
    strand<-event.n[length(event.n)]
    event<-unlist(strsplit(event.n,split="-"))
    eve<-unlist(strsplit(event,split="|",fixed = T))
    chr<-event[1]
    s1<-eve[2]
    s2<-eve[length(eve)-1]
    if(s1>s2){
      start<-s2
      end<-s1
    }else{
      end<-s2
      start<-s1
    }
    name<-ucsc_query(chr,start,end,strand,connection)
    event.m<-subset(dt.df,event_name==e.n)
    delta_psi<-mean(subset(dt.df,event_name==e.n&exp_name%in%sample.n)$psi)-subset(dt.df,event_name==e.n&exp_name=="Control")$psi
    event_type<-as.character(event.m[1,2])
    return(c(chr,start,end,name,strand,e.n,event_type,delta_psi))
  }else{
    e.n<-as.character(e.m)
    event.n<-unlist(strsplit(e.n,split=":"))
    strand<-event.n[length(event.n)]
    event<-unlist(strsplit(event.n,split="-"))
    eve<-unlist(strsplit(event,split="|",fixed = T))
    chr<-event[1]
    s1<-eve[2]
    s2<-eve[length(eve)-1]
    if(s1>s2){
      start<-s2
      end<-s1
    }else{
      end<-s2
      start<-s1
    }
    name<-ucsc_query(chr,start,end,strand,connection)
    event.l<-subset(dt.df,event_name==e.n)
    event_type<-as.character(event.l[2]$event_type)
    return(c(chr,start,end,name,strand,e.n,event_type,event.l[4]))
  }
}
ucsc_query <- function(chr,start,end,strand,connection) {
  refGene <- dbGetQuery(connection,
                   stringr::str_interp(
                     "SELECT DISTINCT name2
      FROM refGene
      WHERE chrom = '${chr}' AND txStart <= ${end} AND txEnd >= ${start} AND strand='${strand}'"))
  if(nrow(refGene)==1){
    refGene<-refGene$name2
  }else if(nrow(refGene)>1){
    refGene<-lapply(refGene,paste,collapse=" ")$name2
  }else{
    refGene<-"NA"
  }
  return(refGene)
}
con_ucsc <- dbPool(drv = RMySQL::MySQL(), db = "hg19", user = "genome", host = "genome-mysql.soe.ucsc.edu")
lapply(1:length(diff_event.l),function(i,e.l,df,con) {
  data.raw<-future_sapply(e.l[[i]],annotate_event,inter=F,dt.df=subset(data.df,exp_name==sample.n[i]),connection=con)
  data<-matrix(unlist(data.raw),ncol=8,byrow=T)
  colnames(data)<-c('chr','start','end','gene name','strand','event name','event type','delta_psi')
  filename<-paste(sample.n[i],"all_event","csv",sep = ".")
  write.csv(data,filename,row.names = F,col.names = T)
} ,e.l=diff_event.l,df=data.df,con=con_ucsc)

data.raw<-future_sapply(inter_event,annotate_event,inter=T,dt.df=inter.df,connection=con_ucsc)
data<-matrix(unlist(data.raw),ncol=8,byrow=T)
colnames(data)<-c('chr','start','end','gene name','strand','event name','event type','delta_psi')
filename<-paste("inter_event","csv",sep = ".")
write.csv(data,filename,row.names = F,quote = F)
poolClose(con_ucsc)

total_event.l<-unique(unlist(diff_event.l))

read_raw_miso<-function(sample_list,event_list,sample_name_list,replicate=T,total_event){
  all_data.l<-lapply(1:length(sample_list),function(i) {
    data<-future_sapply(1:length(event_list),function(j) {
      file_name<-paste(sample_list[i],event_list[j],sep='_')
      event_result<-read.delim(file_name,header=T)
      event.m<-matrix(unlist(apply(event_result,1,filter_raw_miso,sample_name=sample_name_list[i],event_name=event_list[j])),ncol=4,byrow=T)
      colnames(event.m)<-c("sample_name","event_type","event_name","psi")
      event.df<-as.matrix(event.m[event.m[,3]%in%total_event,])
      return(event.df)
    })
    df<-Reduce(rbind,data)
    return(df)
  })
  all_data.df<-Reduce(rbind,all_data.l)
  all_data.df<-as.data.frame(all_data.df)
  all_data.df$psi<-as.numeric(as.character(all_data.df$psi))
  #filter replicate
  if(replicate){
    exp_design<-matrix(unlist(lapply(as.character(all_data.df$sample_name),function(x) unlist(strsplit(x,split = '_')))),ncol=2,byrow=T)
    new_data.df<-data.frame(exp_name=exp_design[,1],rep=exp_design[,2],subset(all_data.df,select = -sample_name))
    new_data.l<-lapply(levels(new_data.df$exp_name),filter_rep,data_frame=new_data.df)
    data.df<-Reduce(rbind,new_data.l)
  }else{
    data.df<-all_data.df
    colnames(data.df)<-c('exp_name','event_type','event_name','psi')
  }
  return(data.df)
}
filter_raw_miso<-function(psi.l,sample_name,event_name) {
  psi<-as.vector(psi.l)
  return(c(sample_name,event_name,psi[1],psi[2]))
}

plan(multiprocess) 
all.df<-read_raw_miso(sample_list = sample.l,event_list = event.l,sample_name_list = sample_name.l,total_event = total_event.l)
colnames(all.df)<-c("exp_name","event_type","event_name","delta_psi")
all.m<-acast(all.df, event_name~exp_name, value.var="delta_psi")

all_psi.m<-apply(all.m[,-1],2,function(x) x-all.m[,1])

check<-apply(all_psi.m,1,function(x) max(abs(x))<0.1)

all<-all_psi.m[-which(check),]

re.l<-apply(all,2,function(x) names(x[abs(x)>0.1]))
venn(re.l,zcolor = "style",size = 20, cexil = 2, cexsn = 1,borders = F)

all.sort<-apply(all,1,function(x){
  a<-ifelse(abs(x)>0.1,sign(x),0)
  s<-sign(sum(x))
  n<-sum(abs(a))
  re<-s*n+s*abs(sum(a*c(6,5,1,3,4,2))/30)
  if(abs(re)<1){
    idx<-which(abs(x)==max(abs(x)))
    s<-sign(x[idx])
    return(s*n+s*abs(sum(a*c(6,5,1,3,4,2))/30))
  }
  else
    return(re)
})

names(all.sort)<-rownames(all)
row.idx<-names(sort(all.sort,decreasing = T))

quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

mat_breaks <- quantile_breaks(all_psi.m, n = 100)
cc = colorRampPalette(rev(brewer.pal(n = 7, 
                                     name = "RdYlBu")))	

a<-sapply(floor(all.sort[row.idx]),function(x) {
  ifelse(x < 0,abs(x)-1,x)})


row_anno<-data.frame(Var1 = factor(a))
col_anno<-c("#66c2a5","#fc8d62","#8da0cb","#e78ac3","#a6d854","#ffd92f")
names(col_anno)<-1:6
pheatmap(all_psi.m[row.idx,c(1,2,5,4,6,3)],show_rownames = F,cluster_rows = F,cluster_cols = F,breaks = mat_breaks,annotation_row = row_anno,annotation_colors=list(Var1=col_anno))  

plot.df<-Reduce(rbind,lapply(1:6,function(x){
  count.df<-subset(all.df,exp_name==names(re.l)[x]&event_name%in%re.l[[x]])
  count.l<-sapply(event.l,function(x)nrow(subset(count.df,event_type==x)))
  df<-data.frame(exp_name=rep(names(re.l)[x],5),event_type=event.l,count=count.l)
  return(df)
  }))

plot.df <- ddply(plot.df,.(exp_name), transform, pos = (sum(count)-cumsum(count)+0.5*count))

ggplot(plot.df,aes(x=exp_name,y=count,fill=event_type))+geom_bar(stat="identity")+ geom_text(data=plot.df,aes(x=exp_name,y=pos,label = count), size = 4,show.legend = F)+scale_fill_manual(values=c("#8dd3c7","#ffffb3","#bebada","#80b1d3","#fb8072"))+theme_classic()
  scale_fill_manual(values=c("#fed9a6","#b3cde3","#ccebc5","#decbe4","#fbb4ae"))

read_miso_raw<-function(sample_name_list,total_event){
  all_data.l<-lapply(1:length(sample_name_list),function(i) {
    file_name<-paste(sample_name_list[i],"all_event.csv",sep='.')
    print(file_name)
    event_result<-read.csv(file_name,header=T)
    event.m<-matrix(unlist(apply(event_result,1,filter_miso_raw,sample_name=sample_name_list[i])),ncol=4,byrow=T)
    colnames(event.m)<-c("sample_name","event_type","event_name","gene_name")
    event.df<-as.matrix(event.m[event.m[,3]%in%total_event,])
    return(event.df)
  })
  all_data.df<-Reduce(rbind,all_data.l)
  all_data.df<-as.data.frame(all_data.df)
  colnames(all_data.df)<-c('exp_name','event_type','event_name','gene_name')
  return(all_data.df)
}
filter_miso_raw<-function(psi.l,sample_name) {
  psi<-as.vector(psi.l)
  return(c(sample_name,psi[7],psi[6],psi[4]))
}
miso_gene.df<-read_miso_raw(sample.n,miso.n)

miso.m<-all

miso.n<-rownames(all)
miso.df<-cbind(t(sapply(miso.n,function(x) c(as.character(unique(subset(miso_gene.df,event_name==x)$event_type)),as.character(unique(subset(miso_gene.df,event_name==x)$gene_name))))),miso.n)
colnames(miso.df)<-c("event_type","gene_name","event_name")
rownames(miso.df)<-NULL

setwd("~/whh/rMATs_diff")
read_rmats_diff<-function(sample_list,event_list,rep_list=rep.l){
  all_data.l<-sapply(1:length(sample_list),function(i) {
    data<-future_sapply(1:length(event_list),function(j) {
      file.n<-paste(sample_list[i],event_list[j],sep='_')
      file_name<-paste(file.n,"diff",sep=".")
      print(file_name)
      event_result<-read.delim(file_name,header=F)
      event.m<-matrix(unlist(apply(event_result,1,filter_rmats,sample_name=sample_list[i],event_type=event_list[j],ifrep=rep_list[i])),ncol=5,byrow=T)
      colnames(event.m)<-c("sample_name","event_type","event_name","gene_name","delta_psi")
      return(event.m)
    })
    df<-Reduce(rbind,data)
    return(df)
  })
  all_data.df<-Reduce(rbind,all_data.l)
  all_data.df<-as.data.frame(all_data.df)
  all_data.df$delta_psi<-as.numeric(as.character(all_data.df$delta_psi))
  data.df<-all_data.df
  colnames(data.df)<-c('exp_name','event_type','event_name',"gene_name",'delta_psi')
  return(data.df)
}
filter_rmats<-function(psi.l,sample_name,event_type,ifrep){
  psi<-as.vector(psi.l)
  if(ifrep)
    count.s<-sum(as.numeric(psi[13]),as.numeric(psi[14]),as.numeric(unlist(strsplit(as.character(psi[15]),split = ','))),as.numeric(unlist(strsplit(as.character(psi[16]),split = ','))))
  else
    count.s<-sum(as.numeric(psi[13]),as.numeric(psi[14]),as.numeric(psi[15]),as.numeric(psi[16]))
  if((count.s>50) & (abs(as.numeric(psi[23]))>0.1)){
    eve.n<-paste(as.character(psi[6:11]),collapse = ":")
    strand<-as.character(psi[5])
    chr<-as.character(psi[4])
    event.n<-paste(chr,eve.n,strand,sep = ":")
    event_name<-gsub(" ","",event.n)
    return(c(sample_name,event_type,event_name,psi[3],as.numeric(psi[23])))
  }
  else
    return(NULL)
}

sample.l<-c('PRMT1','PRMT3','PRMT4','PRMT5','PRMT6','PRMT7')
rep.l<-c(1,1,1,1,0,1)
event.l<-c('A3SS','A5SS','RI','SE')
plan(multiprocess) 
data.df<-read_rmats_diff(sample_list = sample.l,event_list = event.l,rep_list = rep.l)
diff_event.l<-lapply(levels(data.df$exp_name),function(x) subset(data.df,exp_name==x)$event_name)
inter_event<-Reduce(intersect,diff_event.l)
inter.df=subset(data.df,event_name%in%inter_event)

CairoJPEG(file="venn_rmats.jpeg",width=1000,height=1000)
venn(diff_event.l, ilabels = TRUE,zcolor = "style", snames="PRMT1,PRMT3,PRMT4,PRMT5,PRMT6,PRMT7",size = 25, cexil = 2, cexsn = 2.5);
dev.off()

lapply(1:length(diff_event.l),function(i,e.l,df) {
  data.raw<-future_sapply(e.l[[i]],annotate_event,inter=F,dt.df=subset(data.df,exp_name==sample.n[i]))
  data<-matrix(unlist(data.raw),ncol=8,byrow=T)
  colnames(data)<-c('chr','start','end','gene name','strand','event name','event type','delta_psi')
  filename<-paste(sample.n[i],"event","csv",sep = ".")
  write.csv(data,filename,row.names = F,col.names = T)
} ,e.l=diff_event.l,df=data.df)

#write unique event
lapply(sample.l,function(x,df){
  print(x)
  uni.eve<-setdiff(subset(df,exp_name==x)$event_name,subset(df,exp_name!=x)$event_name)
  print(length(uni.eve))
  out.data<-subset(df,event_name %in% uni.eve)
  filename<-paste(x,"event","csv",sep = ".")
  write.csv(out.data,filename,row.names = F,col.names = T)
},df=data.df)
#write all event
lapply(sample.l,function(x,df){
  print(x)
  out.data<-subset(df,exp_name==x)
  print(nrow(out.data))
  filename<-paste(x,"all_event","csv",sep = ".")
  write.csv(out.data,filename,row.names = F,col.names = T)
},df=data.df)


total_event.l<-as.character(unique(data.df$event_name))
read_rmats_raw<-function(sample_list,event_list,rep_list=rep.l,total_event){
  all_data.l<-sapply(1:length(sample_list),function(i) {
    data<-future_sapply(1:length(event_list),function(j) {
      file.n<-paste(sample_list[i],event_list[j],sep='_')
      file_name<-paste(file.n,"txt",sep=".")
      print(file_name)
      event_result<-read.delim(file_name,header=T)
      event.m<-matrix(unlist(apply(event_result,1,filter_raw_rmats,sample_name=sample_list[i],event_type=event_list[j],ifrep=rep_list[i],event.t=total_event)),ncol=3,byrow=T)
      colnames(event.m)<-c("sample_name","event_name","delta_psi")
      return(event.m)
    })
    df<-Reduce(rbind,data)
    return(df)
  })
  all_data.df<-Reduce(rbind,all_data.l)
  all_data.df<-as.data.frame(all_data.df)
  all_data.df$delta_psi<-as.numeric(as.character(all_data.df$delta_psi))
  data.df<-all_data.df
  colnames(data.df)<-c('exp_name','event_name','delta_psi')
  return(data.df)
}
filter_raw_rmats<-function(psi.l,sample_name,event_type,ifrep,event.t){
  psi<-as.vector(psi.l)
  event.n<-paste(as.character(psi[4]),paste(as.character(psi[6:11]),collapse = ":"),as.character(psi[5]),sep = ":")
  event_name<-gsub(" ","",event.n)
  if(event_name %in% event.t){
    return(c(sample_name,event_name,as.numeric(psi[23])))
  }
  else
    return(NULL)
}
plan(multiprocess)
all.df<-read_rmats_raw(sample_list = sample.l,event_list = event.l,rep_list = rep.l,total_event = total_event.l)
all.m<-acast(all.df, event_name~exp_name, value.var="delta_psi")
all.m[is.na(all.m)]<-0



rmats.m<- -all.m

rmats.n<-rownames(all.m)
rmats.df<-cbind(t(sapply(rmats.n,function(x) c(as.character(unique(subset(data.df,event_name==x)$event_type)),as.character(unique(subset(data.df,event_name==x)$gene_name))))),rmats.n)
colnames(rmats.df)<-c("event_type","gene_name","event_name")
rownames(rmats.df)<-NULL

new_rmats_name<-matrix(unlist(apply(miso.df[miso.df[,1]!="MXE",],1,function(x) miso2rmats(x[3],x[1]))),ncol=2,byrow = T)
inter.n<-new_rmats_name[new_rmats_name[,2]%in%rmats.df[,3],]
event_name.df<-rbind(miso.df,rmats.df)
merge.m<-rbind(miso.m,rmats.m[!rownames(rmats.m)%in%inter.n[,2],])

miso2rmats<-function(name,tp){
  event.n<-unlist(strsplit(name,split=":"))
  strand<-event.n[length(event.n)]
  event<-unlist(strsplit(event.n,split="-"))
  eve<-unlist(strsplit(event,split="|",fixed = T))
  if(tp=="A3SS"){
    if(strand=="-"){
      s1<-as.character(as.numeric(eve[7])-1)
      s5<-as.character(as.numeric(eve[2])-1)
      re<-paste(eve[1],s1,eve[6],s1,eve[5],s5,eve[3],strand,sep = ":")
    }
    if(strand=="+"){
      s1<-as.character(as.numeric(eve[5])-1)
      s3<-as.character(as.numeric(eve[6])-1)
      s5<-as.character(as.numeric(eve[2])-1)
      re<-paste(eve[1],s1,eve[7],s3,eve[7],s5,eve[3],strand,sep = ":")
    } 
  }
  if(tp=="A5SS"){
    if(strand=="-"){
      s1<-as.character(as.numeric(eve[3])-1)
      s3<-as.character(as.numeric(eve[4])-1)
      s5<-as.character(as.numeric(eve[6])-1)
      re<-paste(eve[1],s1,eve[2],s3,eve[2],s5,eve[7],strand,sep = ":")
    }
    if(strand=="+"){
      s1<-as.character(as.numeric(eve[2])-1)
      s5<-as.character(as.numeric(eve[6])-1)
      re<-paste(eve[1],s1,eve[4],s1,eve[3],s5,eve[7],strand,sep = ":")
    }
  }
  if(tp=="RI"){
    if(strand=="+"){
      s1<-as.character(as.numeric(eve[2])-1)
      s5<-as.character(as.numeric(eve[5])-1)
      re<-paste(eve[1],s1,eve[6],s1,eve[3],s5,eve[6],strand,sep = ":")
    }
    if(strand=="-"){
      s1<-as.character(as.numeric(eve[6])-1)
      s5<-as.character(as.numeric(eve[3])-1)
      re<-paste(eve[1],s1,eve[2],s1,eve[5],s5,eve[2],strand,sep = ":")
    }
  }
  if(tp=="SE"){
    if(strand=="+"){
      s1<-as.character(as.numeric(eve[5])-1)
      s3<-as.character(as.numeric(eve[2])-1)
      s5<-as.character(as.numeric(eve[8])-1)
      re<-paste(eve[1],s1,eve[6],s3,eve[3],s5,eve[9],strand,sep = ":")
    }
    if(strand=="-"){
      s1<-as.character(as.numeric(eve[5])-1)
      s5<-as.character(as.numeric(eve[2])-1)
      s3<-as.character(as.numeric(eve[8])-1)
      re<-paste(eve[1],s1,eve[6],s3,eve[9],s5,eve[3],strand,sep = ":")
    }
  }
  return(c(name,re))
}

re<-sapply(intersect(rmats.df[,2],miso.df[,2]),function(x){
  a=miso.df[miso.df[,2]==x,1]
  b=rmats.df[rmats.df[,2]==x,1]
  if(a%in%b){
    print(x)
    print(intersect(a,b))
    print(miso.df[miso.df[,2]==x,3])
  }
})

check<-apply(merge.m,1,function(x) max(abs(x))<0.1)


re.l<-apply(merge.m,2,function(x) names(x[abs(x)>0.1]))
venn(re.l,zcolor = "style",size = 20, cexil = 2, cexsn = 1,borders = F)

all.sort<-apply(merge.m,1,function(x){
  a<-ifelse(abs(x)>0.1,sign(x),0)
  s<-sign(sum(x))
  n<-sum(abs(a))
  re<-s*n+s*abs(sum(a*c(6,5,1,3,4,2))/30)
  if(abs(re)<1){
    idx<-which(abs(x)==max(abs(x)))
    s<-sign(x[idx])
    return(s*n+s*abs(sum(a*c(6,5,1,3,4,2))/30))
  }
  else
    return(re)
})

names(all.sort)<-rownames(merge.m)
row.idx<-names(sort(all.sort,decreasing = T))

mat_breaks <- quantile_breaks(merge.m, n = 100)

a<-sapply(floor(all.sort[row.idx]),function(x) {
  ifelse(x < 0,abs(x)-1,x)})


row_anno<-data.frame(Var1 = factor(a))
col_anno<-c("#66c2a5","#fc8d62","#8da0cb","#e78ac3","#a6d854","#ffd92f")
names(col_anno)<-1:6
pheatmap(merge.m[row.idx,c(1,2,5,4,6,3)],show_rownames = F,cluster_rows = F,cluster_cols = F,color=cc(100),breaks = mat_breaks,annotation_row = row_anno,annotation_colors=list(Var1=col_anno))  


plot.df<-Reduce(rbind,lapply(1:6,function(x){
  count.df<-event_name.df[event_name.df[,3]%in%re.l[[x]],]
  count.l<-sapply(event.l,function(x) length(which(count.df[,1]==x)))
  df<-data.frame(exp_name=rep(names(re.l)[x],5),event_type=event.l,count=count.l)
  return(df)
}))

ggplot(plot.df,aes(x=exp_name,y=count,fill=event_type))+geom_bar(stat="identity")+ geom_text(data=plot.df,aes(x=exp_name,y=pos,label = count), size = 4,show.legend = F)+scale_fill_manual(values=c("#8dd3c7","#ffffb3","#bebada","#80b1d3","#fb8072"))+theme_classic()
scale_fill_manual(values=c("#fed9a6","#b3cde3","#ccebc5","#decbe4","#fbb4ae"))