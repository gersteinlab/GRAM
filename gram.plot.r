args = commandArgs(trailingOnly=TRUE)

require("reshape2")
require("gplots")
require("grid")

thr.em = 0.45

gram.file = args[1]
funseq.file = args[2]
output.path = args[3]

setwd(output.path)

#snp.regions = read.table("SNP_regions_selected.txt", header=TRUE, stringsAsFactors=FALSE, sep="\t")
### Load and process data
## TODO: check if there is patient id
# (1) no patient id or patient number > 5: same color
# (2) patient number < 5: color by patient

#gram.file='../output/eQTL_eqtl.merged.gram.190404.nc.bed'
#funseq.file='../output/eQTL_eqtl.merged.funseq.190404.nc.bed'

smallSize=FALSE
multiRegion=FALSE

#### Load GRAM output
# Data format: chr st ed allele1 allele2 name score
dat.gram = read.table(gram.file, stringsAsFactors=FALSE, header=TRUE)
colnames(dat.gram)=c("chr", "st", "ed", "ref", "mut", "sampleid", 
                     "snpid", "ref.enhAct", "alt.enhAct", "logodds", "expr.rf", "vodds.rf", "score")


get_region_dict <- function(dat, region.len=1000000){
  # Automatically identify regions for plotting
  dat.pos = unique(dat[,c("chr","st","ed","snpid")])
  region.list = list()
  for(chr in unique(dat.pos$chr)){
    dat.pos.chr = dat.pos[dat.pos$chr==chr,]
    pos.min = min(dat.pos.chr$ed)
    pos.max = max(dat.pos.chr$ed)
    seg.len = min(1000000, pos.max-pos.min+1)
    pos.seg = seq(pos.min,pos.max,by=seg.len)
    region.chr = data.frame(st=pos.min, ed=pos.min+seg.len)
    rownames(region.chr) = paste(chr,"-",pos.min,"-",pos.max,sep="")
    region.list[[chr]] = region.chr
  }
  
  region.var = NULL
  for(i in seq(1, nrow(dat.pos))){
    chr = dat.pos[i,]$chr
    pos = dat.pos[i,]$ed
    
    region.chr = region.list[[chr]]
    for(j in seq(1, nrow(region.chr))){
      if(pos>=region.chr[j,1] && pos<=region.chr[j,2]) region.var = c(region.var, rownames(region.chr)[j])
    }
    
  }
  names(region.var) = dat.pos$snpid
  return(region.var)
}
  
region.var = get_region_dict(dat.gram)
dat.gram$region = region.var[dat.gram$snpid]

if(length(unique(region.var))>1) multiRegion=TRUE
if(length(unique(dat.gram$sample))<=5){smallSize=TRUE}

#### Load Funseq output
# Data format: chr st ed ref mut name info
dat.funseq.raw = read.table(funseq.file, stringsAsFactors=FALSE,
                            col.names=c("chr","st","ed","ref","mut","name","info"),comment.char = "$",sep="\t")

dat.funseq = dat.funseq.raw[c("chr","ed")]
dat.funseq$score = unlist(lapply(dat.funseq.raw$info, 
                                 function(x){l=strsplit(x, ";")[[1]];return(l[length(l)-2])}))

dat.gram = dat.gram[(dat.funseq$score != "."),]
dat.funseq = dat.funseq[(dat.funseq$score != "."),]

dat.funseq$score = as.numeric(dat.funseq$score)
colnames(dat.funseq) = c("chr","ed","score")


########## Plot
### Heatmap
if(!smallSize){
  pdf("figure/heatmap.pdf", height=10, width=8)
  dat.heatmap = acast(dat.gram, sampleid~snpid, value.var="score")
  min.dat.heatmap = min(dat.heatmap[!is.na(dat.heatmap)])
  max.dat.heatmap = max(dat.heatmap[!is.na(dat.heatmap)])
  col.ratio = max.dat.heatmap/min.dat.heatmap
  num.col1 = floor(32/(1 + col.ratio)) # size of color palette for "NA" (0) values
  num.col2 = 32 - num.col1
  
  dat.heatmap[is.na(dat.heatmap)] = 0
  #dat.heatmap[is.na(dat.heatmap)] = min(dat.heatmap[!is.na(dat.heatmap)])
  heatmap.2(t(dat.heatmap), col=c(colorRampPalette(c("gray80"))(num.col1),colorRampPalette(c("gray99", "red"))(num.col2)),
            trace="none", lhei = c(1.5,7), cexRow=0.5, cexCol=0.35)
  dev.off()
}


### Write data table
if(multiRegion){
  max.score = aggregate(score~sampleid+region, dat.gram, max)
  colnames(max.score)[3] = "max.score"
  rownames(max.score) = paste(max.score$sampleid,max.score$region,sep="_")
  
  variant.count = as.data.frame(table(dat.gram$sampleid, dat.gram$region))
  rownames(variant.count) = paste(variant.count$Var1,variant.count$Var2,sep="_")
  
  max.score$variant.count = variant.count[rownames(max.score),"Freq"]
  
  max.score = max.score[c("sampleid","region","variant.count", "max.score")]
  max.score$file.name = paste("score_",max.score$region,"_",max.score$sampleid,".pdf", sep="")
} else{
  max.score = aggregate(score~sampleid, dat.gram, max)
  colnames(max.score)[2] = "max.score"
  rownames(max.score) = max.score$sampleid
  
  variant.count = as.data.frame(table(dat.gram$sampleid))
  rownames(variant.count) = variant.count$Var1
  
  max.score$variant.count = variant.count[rownames(max.score),"Freq"]
  
  max.score = max.score[c("sampleid","variant.count", "max.score")]
  max.score$file.name = paste("score_",max.score$sampleid,".pdf", sep="")
}

max.score = max.score[order(max.score$max.score, decreasing=TRUE),]
write.table(max.score, "result_list.txt", quote=FALSE, row.names=FALSE)

### Track plot
plot.w = 8
plot.h = 4

region.list = unique(dat.gram$region)
sample.list = unique(dat.gram$sampleid)

plot_score_track <- function(dat.plot, xlim, ylim, main, ylab, thr="none"){
  # var.list: the list of vars to highlight
  dat.top.plot = dat.plot[c("ed","snpid","score")]
  #dat.top.plot$snpid.trimmed = unlist(lapply(dat.top.plot$snpid, function(x){l=strsplit(x, "_")[[1]];n=length(l);if(n==1){return(l);}else{
  #  return(paste(as.character(l[-c(n,(n-1),(n-2))]),collapse="_"))}}))
  
  dat.top.plot$snpid.trimmed = dat.top.plot$snpid
  
  dat.top.plot.top = dat.top.plot[order(dat.top.plot$score, decreasing=TRUE)[1:(min(3,nrow(dat.top.plot)))],]
  var.list = dat.top.plot.top$snpid[dat.top.plot.top$score>=thr.em]
  
  color.panel = c("red","purple","green","cyan") 
  color.panel.sub = color.panel[1:length(var.list)]
  names(color.panel.sub) = var.list
  
  col.list = color.panel.sub[dat.top.plot$snpid] # Coloring: special color for highlighted points
  col.list[is.na(col.list)] = "gray"
  
  pch.list = rep(1, nrow(dat.plot)) # Point type: filled for highlighted points
  pch.list[dat.top.plot$snpid %in% var.list] = 16
  
  dat.label = dat.top.plot[dat.top.plot$snpid %in% var.list,] # vars that need labelling
  
  #gg = ggplot(dat.top.plot, aes(x=ed, y=score)) + geom_point(shape=pch.list, color=col.list, size=2) +
  #  xlim(xlim[1], xlim[2])+
  #  ylim(ylim[1], ylim[2])+
  #  theme_bw() +
  #  theme(axis.line = element_line(colour = "black"),
  #        panel.grid.major = element_blank(),
  #        panel.grid.minor = element_blank(),
  #        panel.border = element_blank(),
  #        panel.background = element_blank()) 
  
  plot(dat.top.plot$ed, dat.top.plot$score, xlab=region, ylab=ylab, xlim=xlim, ylim=ylim,
       col=col.list, pch=pch.list, cex=1.5, xaxt='n',yaxt='n', axes=FALSE)
  axis(side=1, at=seq(xlim[1], xlim[2], length.out = 5))
  axis(side=2)
  if(sum(dat.top.plot$snpid %in% var.list) > 0) text(dat.label$ed, dat.label$score, dat.label$snpid.trimmed, cex=1, pos=4, offset=0.5)
  mtext(main, side=3, line=-2)
  
  if(thr != "none"){
    abline(h=thr, col="red",lty=2)
  }
}

plot_tracks <- function(dat.gram.sub, dat.funseq.sub){
  mat.gram.sub = acast(dat.gram.sub, snpid~sampleid, value.var="score")
  
  gram.mean = apply(mat.gram.sub, 1, function(x){mean(x[!is.na(x)])})
  #gram.mean = apply(mat.gram.sub, 1, function(x){max(x[!is.na(x)])})
  #gram.score = apply(mat.gram.sub, 1, function(x){max(x[!is.na(x)])})
  
  gram.top.mean = sort(gram.mean, decreasing=TRUE)[1:min(3,length(gram.mean))]
  gram.top.mean = gram.top.mean[gram.top.mean>thr.em]
  gram.top.mean.var = names(gram.top.mean)
  
  #dat.gram.top.var = dat.gram.sub[dat.gram.sub$snpid==gram.top.mean.var[1],]
  #dat.gram.top.var = dat.gram.top.var[dat.gram.top.var$score>=thr.em,]
  #sample.top.var = dat.gram.top.var[order(dat.gram.top.var$score,decreasing=TRUE),]$sampleid[1:min(5,nrow(dat.gram.top.var))]
  
  pos=dat.gram.sub$ed
  pos.left = min(pos)
  pos.right = max(pos)
  
  
  ### Calculate plot region (resolution: magnitude of (pos.right-pos.left) - 4)
  res = 10^(ceiling(log(pos.right-pos.left, 10))-2)
  
  if(res==0){
    res=1
    span.left = pos.left-5
    span.right = pos.right+5
  } else{
    span.left = ceiling(pos.left/res)*res-res*2
    span.right = ceiling(pos.right/res)*res+res*2
  }
  
  xlim=c(span.left,span.right)
  ylim = c(min(dat.gram.sub$score)-0.1,max(dat.gram.sub$score)+0.1)
  
  
  #par(mfrow=c(length(sample.top.var)+2,1), xpd=TRUE, mar=c(4,6,4,6))
  
  
  ### Plot average GRAM scores
  pdf(paste("figure/trackplot/trackplot_",region,"_mean", ".pdf", sep=""), width=plot.w, height=plot.h, useDingbats = FALSE )
  dat.gram.mean.score = unique(dat.gram.sub[dat.gram.sub$snpid %in% names(gram.mean), 
                                            c("chr","ed","snpid")])
  dat.gram.mean.score$score = gram.mean[dat.gram.mean.score$snpid]
  plot_score_track(dat.gram.mean.score, xlim, ylim, main="Max GRAM score", ylab="GRAM score")
  dev.off()
  
  ### Plot GRAM scores
  for(sampleid in dat.gram.sub$sampleid){
    pdf(paste("figure/trackplot/trackplot_",region,"_", sampleid, ".pdf", sep=""), width=plot.w, height=plot.h, useDingbats = FALSE )
    dat.gram.top.sample = dat.gram.sub[dat.gram.sub$sampleid == sampleid,]
    plot_score_track(dat.gram.top.sample, xlim, ylim, main=sampleid, ylab="GRAM score")
    dev.off()
  }
  
  ### Plot funseq scores
  pdf(paste("figure/trackplot/trackplot_",region,"_funseq", ".pdf", sep=""), width=plot.w, height=plot.h, useDingbats = FALSE )
  dat.funseq.sub.plot = dat.funseq.sub
  dat.funseq.sub.plot$snpid = dat.gram.sub$snpid
  dat.funseq.sub.plot = unique(dat.funseq.sub.plot)
  plot_score_track(dat.funseq.sub.plot, xlim, c(0,5), main="Funseq score", ylab="Funseq score")
  dev.off()
  
}

for(region in region.list){
  dat.gram.sub = dat.gram[dat.gram$region==region,]
  dat.funseq.sub = dat.funseq[dat.gram$region==region,]
  #pdf(paste("figures/trackplot/trackplot_",region,".pdf", sep=""), width=6, height=20)
  plot_tracks(dat.gram.sub, dat.funseq.sub)
  #dev.off()
}
