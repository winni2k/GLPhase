library(pwr)
library(ggplot2)
# assuming effect sizes of 0.1, 0.2, 0.3 %

popsize <- c(0.5,1,2,4) * 10000

effectSize <- c(0.1,0.15,0,2,0.25,0.3,0.4,0.5,1, 2, 3)/1000
rsquare <- c(6,7,8,9,9.5,9.8)/10
rsquare
siglevel <- 5e-8

data <- expand.grid(effectSize,rsquare, siglevel,popsize)
names(data) <- c("effectSize", "rsquare", "sigLevel", "popSize")
data

res <- apply(data, 1, function(x) {
  pwr.f2.test(u=1, v= x[2]*x[4] -2, f2=x[1]/(1-x[1]), sig.level=x[2])
})      
length(res)

data$power <- sapply(res, function(x) x$power)

pdf(file="powerPlot.pdf")
ggplot(data=data, aes(x=effectSize, y=power,  color=factor(rsquare))) + geom_line() +facet_wrap(~popSize, nrow=2)
dev.off()

## now let's do cost analysis

genotyping <- 400
samplePrep <- 30
costPer1x <- 133
coverage <- (1:10)/10
budget <- 300000

names(data)
head(data)
popSize <- c(0.5,1,2,4)*10000

costData <- expand.grid(genotyping,samplePrep,costPer1x,budget,coverage)
names(costData) <- c("genotyping","samplePrep","costPer1x","budget","coverage")
costData$sampSize <- with(costData, (budget/( costPer1x*coverage + samplePrep)))

head(data)
names(data) <- c("effectSize", "rsquare", "sigLevel", "popSize")
pwr.f2.test(u=1, v= x[2]*x[4] -2, f2=x[1]/(1-x[1]), sig.level=x[2])

pdf(file="sampleSize.pdf")
ggplot(data=costData, aes(x=coverage, y=sampSize,  )) + geom_line() 
dev.off()
