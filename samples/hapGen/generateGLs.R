d <- read.delim("ex.haps", sep = " ", header=FALSE)
l <- read.delim("ex.leg", sep="", header=TRUE)

head(d)

g <- d[,seq(1,ncol(d), by = 2)]
g = g + d[,seq(2,ncol(d), by = 2)]
head(g)
?rbinom
?matrix
gCov <- matrix(rpois(nrow(g) * ncol(g),coverage), nrow=nrow(g))
head(gCov)

length(as.vector(gCov))
temp = rbinom(as.vector(gCov),as.vector(gCov),as.vector(g)/2)
head(temp)
coverage = 4

as.vector(gCov)
unlist(g)
gdf = cbind(gCov=as.vector(gCov), g=unlist(g))
head(gdf)

gdf = as.data.frame(gdf)
str(gdf)

### sample coverage of alt reads
gdf = cbind(gdf, fuzzedGT=sapply(gdf[,2], function(x) if(x==1) 1 else if(x==0) 0.01 else if (x==2) 1.99))

head(gdf)
temp = apply(gdf, 1, function(x) rbinom(1, x[1], x[3]/2))
head(temp)
gdf = cbind(gdf, altReads=temp)


temp = apply(gdf, 1, function(x) dbinom(x[4], x[1], 0.005))
temp2 = apply(gdf, 1, function(x) dbinom(x[4], x[1], 0.5))
temp3 = apply(gdf, 1, function(x) dbinom(x[4], x[1], 0.995))

gls = cbind(rr=temp, ra=temp2, aa=temp3)
head(gls[which(apply(gls, 1 , sum) >2),])

## gls look ok
out = cbind(gdf, gls)

### create bin from gls
?round
dim(d)
out = apply(gls, 1, function(x) paste(
    round(x[2]/sum(x), 3),
    round(x[3]/sum(x), 3)))
out = as.data.frame(matrix(out, nrow=1001), stringsAsFactors=FALSE)


names(out) = paste0("samp",1:60)
tmp = cbind(chr=rep("20",1000),pos=l[,2], allele=paste0(l[,3],l[,4]), out[1:1000,])

?write.table
write.table(tmp, gzfile("ex.bin"), sep="\t", quote=FALSE, row.names=FALSE)

### save in gen format
out = apply(gls, 1, function(x) paste(x[1], x[2], x[3]))
length(out)
out = as.data.frame(matrix(out, nrow=1001), stringsAsFactors=FALSE)
head(out)

tmp = cbind(chr=rep("20",1000),id=l[,2],pos=l[,2], a0=l[,3], a1=l[,4], out[1:1000,])
write.table(tmp, gzfile("ex.gen.gz"), sep=" ", quote=FALSE, row.names=FALSE, col.names=FALSE)
