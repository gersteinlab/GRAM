

rm(list=ls())
#load("~/scratch/funseq3/mut-list/00final_ds_link/var.f14figure4.rdata")
options(stringsAsFactors=FALSE)

args<-commandArgs(TRUE)

db.ref.file=args[1]
db.mut.file=args[2]
snp.input.file=args[3]  ###
snp.uniq.file=args[4] ### snp.uniq  chr, st, ed, ref, mut, id
test.expr.file=args[5]
model.rd=args[6]
outfile=args[7]

if (length(args)!=7){

  print("Usage: Rscripts gram.predict.r ref-allele-deepbind mut-allele-deepbind snp-filter-cds-bed7 snp-unique-bed6 sample-expression model.rdata output\nbed7 format: chr, start, end, ref, alt, sample-id, rsid; bed6 format: chr, start, end, ref, alt, rsid\ngene expression: row name: gene symbol, colname: sample-id\n")
  
return(1);
}

print(paste("Rscript gram.predict.r", db.ref.file, db.mut.file, snp.input.file, snp.uniq.file, test.expr.file, model.rd, outfile,sep=" "))

require("randomForest")
require("glmnet")

load(model.rd)

step1.selex.featname=rownames(step1.rf.dpselex.model$importance)
step2.expr.featname=rownames(step2.expr.rf.model$importance) ##rownames(step2.var.refdb_selex.lasso.model$beta)
step2.selex.featname=rownames(step2.vodds.rf.model$importance) ##rownames(step2.var.refdb_selex.lasso.model$beta)

##################################################
## expression data and rank
#
test.snv=read.table(snp.input.file, sep="\t", header=F, stringsAsFactors=F)
colnames(test.snv)=c("chrome","start","end", "ref", "mut","sampleid","snpid")

test.expr.raw=read.table(test.expr.file,sep="\t", header=T,row.names=1,check.names = FALSE, stringsAsFactors=F)  ###

#selex=dp.cname

##gm.expr 20532 same dimension

com.gene=intersect(rownames(test.expr.raw), unique(genes))
expr.tf.uniqname=dp.comm

notfound=unique(setdiff(genes, rownames(test.expr.raw)), setdiff(expr.tf.uniqname, com.gene))

test.expr.na=as.data.frame(matrix(0, nrow=length(notfound), ncol=ncol(test.expr.raw)))
rownames(test.expr.na)=notfound
colnames(test.expr.na)=colnames(test.expr.raw)
##test.expr.na=0
##gene expre for rank
##gene name
test.expr.raw<-rbind(test.expr.raw[com.gene,], test.expr.na)

dim(test.expr.raw)

#[1] 20242     3

test.expr<-apply(test.expr.raw, 2, function(x){
  y=rank(x,ties.method="min")/length(x); y;}
 )

exp.profile=colnames(test.expr)


#cbind(test.expr, data.frame(rank=rank(test.expr.raw[,3])/nrow(test.expr)))
###rownames as hugo gene symbol
## same gene expression will be summed up
### row: gene, column is sample
#tfname=dp.cname[colnames(dp.ref.nodup),2]

test.tf.expr<-test.expr[expr.tf.uniqname,]

#test.tf.expr[duplicated(test.tf.expr$gname),]
#test.tf.expr<-test.tf.expr[which(!test.tf.expr$gid %in% c("ENSG00000258724.1", "ENSGR0000185960.8", "ENSGR0000214717.5")),]
#rownames(test.tf.expr)=test.tf.expr$gname

##################################################
###get cell line specific , get ref
test.db.ref0 = read.table(db.ref.file,sep="\t",header=T,check.names=F, stringsAsFactors=F)
test.db.mut0 = read.table(db.mut.file,sep="\t",header=T,check.names=F,stringsAsFactors=F)
snp.input=read.table(snp.input.file,sep="\t", header=F, stringsAsFactors=F)
snp.uniq = read.table(snp.uniq.file, sep="\t", header=F, stringsAsFactors=F)

colnames(snp.input)=c("chrome","start","end", "ref", "mut","sampleid","snpid")
####

if (nrow(test.db.ref0) != nrow(snp.uniq) || nrow(test.db.mut0) !=nrow(snp.uniq)){
  print("The uniq snp doesn't have the same of records from deepbind");
  return(1);
}
###snp.uniq[,6] == snp.input[,7]
rownames(test.db.ref0)=snp.uniq[,6]
rownames(test.db.mut0)=snp.uniq[,6]


##Filter based on the sample id of gene expression
print(paste("After filting, only ", length(intersect(unique(snp.input$sampleid), colnames(test.tf.expr))), "of ", length(unique(snp.input$sampleid))," sample left for the prediction",sep=""))

snp.input=snp.input[which(snp.input$sampleid %in% colnames(test.tf.expr)),, drop=F]

if (nrow(snp.input)==0){
  print("No snv left after filtering using matched gene expression data")
  return;
}


####get all deepbind results from unique deepbind results
##First filter the snp.input[,7] not  in snp.uniq[,6]
com.rsid=intersect(unique(snp.input[,7]), unique(snp.uniq[,6]))

if ( length(com.rsid) != length(unique(snp.input[,7])) ){
  ###it should be the same, since the unique snp is from snp.input
  print("The snp input file should have the same number of unique rsid as the uniuqe snp file");
  print("only the common rsid will be considered.");

}

snp.input=snp.input[which(snp.input[,7] %in% com.rsid),]

test.db.ref = test.db.ref0[snp.input[,7],]
test.db.mut = test.db.mut0[snp.input[,7],]


###require dp.ref.nodup, which is used for the rank estimation
test.db.ref1 = test.db.ref[, colnames(dp.ref.nodup)]
 test.num=nrow(test.db.ref1)
 test.db.ref1 <- rbind(test.db.ref1, dp.ref.nodup)   ###
  

 test.db.ref.rank = apply(test.db.ref1, 2, function(x, count){
   yy=unlist(lapply(x[1:count], function(xx, refcnt){
      sum(xx>=refcnt);
    },refcnt = x[-c(1:count)]))+1
}, count = test.num)

###gene expresion ref-rank based on a pre-defined expression file
test.tf.refrank = t(apply(cbind(data.frame(sample.id=snp.input[,6]), test.db.ref.rank), 1, function(x, tf, expr){
  sid=x[1];
  tt=unlist(lapply( split(x[-1], f=as.factor(tf)), function(xx){ max(xx);}))
  ttnn=names(tt)

  tt=as.numeric(tt)
  names(tt)=ttnn
    
  tt=sort(tt)
  #names(tt)
  y = expr[names(tt),sid]  ###here is the gene rank on sample_id

}, tf=expr.tf.uniqname, expr=test.tf.expr))

colnames(test.tf.refrank)<-paste("X",1:ncol(test.tf.refrank),sep="")



term.label=attr(step1.rf.dpselex.model$terms, "term.labels")


##################################################
###prediction 1. log FC 

test.db.ref.pred1=predict(step1.rf.dpselex.model, test.db.ref[,step1.selex.featname] , type="prob")
test.db.mut.pred1=predict(step1.rf.dpselex.model, test.db.mut[,step1.selex.featname], type="prob")


test.logfc_selex.pred=log2( test.db.mut.pred1[,2]/(1-test.db.mut.pred1[,2])) - log2(test.db.ref.pred1[,2]/(1-test.db.ref.pred1[,2]))

##################################################
###var prediction 2. var


test.var.refrank.pred= predict(step2.expr.rf.model, test.tf.refrank[, step2.expr.featname], type="prob")


test.var.refdb_selex.pred=predict(step2.vodds.rf.model,  test.db.ref[, step2.selex.featname], type="prob")
#####use selex: used for prediction , ref allele with selex scores

##################################################
##prediction: emvar vs non-emvar
test.step3.ds<-data.frame(logodds=abs(test.logfc_selex.pred), expr.rf = test.var.refrank.pred[,2], vodds.rf=test.var.refdb_selex.pred[,2])

#test.step3.ds<-data.frame(logodds=abs(test.logfc_selex.pred), vodds.rf = test.var.refrank.pred[,1], db.las=test.var.refdb_selex.pred[,1])

test.final_fc_db_rank_selex.pred=predict(step3.las.model, newx = as.matrix(test.step3.ds),type="response")


out=cbind(snp.input, data.frame(ref.enhAct=test.db.ref.pred1[,2],alt.enhAct=test.db.mut.pred1[,2], logodds=(test.logfc_selex.pred)), test.step3.ds[,-1], data.frame(gram.prob=test.final_fc_db_rank_selex.pred[,1]))


write.table(out, file=outfile, quote=F, row.names=F, col.names=T)





