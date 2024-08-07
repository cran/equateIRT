## -----------------------------------------------------------------------------
library("equateIRT")
data("data2pl", package = "equateIRT")

## ----message=FALSE, results='hide'--------------------------------------------
library("mirt")
m1 <- mirt(data2pl[[1]], SE = TRUE)
m2 <- mirt(data2pl[[2]], SE = TRUE)
m3 <- mirt(data2pl[[3]], SE = TRUE)
m4 <- mirt(data2pl[[4]], SE = TRUE)
m5 <- mirt(data2pl[[5]], SE = TRUE)

## -----------------------------------------------------------------------------
mlist<- list(m1, m2, m3, m4, m5)
test <- paste("test", 1:5, sep = "")
mod2pl <- modIRT(est.mods = mlist, names = test, display = FALSE)
coef(mod2pl$test1)[1:5]

## -----------------------------------------------------------------------------
lplan<-linkp(mods = mod2pl)
lplan

## ----message=FALSE, fig.cap="Linkage plan",fig.width=3.5,fig.height=3---------
library(sna)
par(mar=c(0, 0, 0, 0))
set.seed(6)
gplot(lplan, displaylabels = TRUE,  vertex.sides = 4, vertex.cex = 5, vertex.rot =45,  usearrows = FALSE, label.pos = 5, label.cex = 1, vertex.col = 0)

## -----------------------------------------------------------------------------
l12 <- direc(mods = mod2pl, which = c(1,2), method = "mean-mean")
l12
summary(l12)

## -----------------------------------------------------------------------------
direclist2pl <- alldirec(mods = mod2pl, method = "mean-mean")
direclist2pl

## -----------------------------------------------------------------------------


## -----------------------------------------------------------------------------
summary(direclist2pl, link="test1.test2")

## -----------------------------------------------------------------------------
cec4 <- chainec(r = 4, direclist = direclist2pl)
cec4

## -----------------------------------------------------------------------------
summary(cec4, path="test1.test2.test3.test4")

## -----------------------------------------------------------------------------
cec14 <- chainec(r = 4, direclist = direclist2pl, f1 = "test1", f2 = "test4")
cec14
summary(cec14)

## -----------------------------------------------------------------------------
pth <- paste("test", c(1,5,4), sep = "")
chainec154 <- chainec(direclist = direclist2pl, pths = pth)
summary(chainec154)

## -----------------------------------------------------------------------------
ecall <- c(cec14, chainec154)
fec <- bisectorec(ecall = ecall, weighted = TRUE, unweighted = TRUE)
fec
summary(fec)

## -----------------------------------------------------------------------------
 eqc(fec)

## -----------------------------------------------------------------------------
itm(fec, bistype = "weighted")

## -----------------------------------------------------------------------------
score(fec, bistype = "weighted")

## -----------------------------------------------------------------------------
score(fec, method = "OSE", bistype = "weighted")

## -----------------------------------------------------------------------------
score(chainec154, scores = 17)
score(cec4, path = "test1.test2.test3.test4", scores = 17)
score(fec, bistype = "unweighted", scores = 17)
score(fec, bistype = "weighted", scores = 17)

## -----------------------------------------------------------------------------
data(dataDIF)

## ----message=FALSE, results='hide'--------------------------------------------
library(mirt)
data1 <- dataDIF[dataDIF$group == 1, 1:20]
data2 <- dataDIF[dataDIF$group == 2, 1:20]
data3 <- dataDIF[dataDIF$group == 3, 1:20]
mod1 <- mirt(data1, SE = TRUE)
mod2 <- mirt(data2, SE = TRUE)
mod3 <- mirt(data3, SE = TRUE)

## -----------------------------------------------------------------------------
res_diftest2 <- dif.test(est.mods = list(mod1, mod2))
res_diftest2

## -----------------------------------------------------------------------------
res_diftest3 <- dif.test(est.mods = list(mod1, mod2, mod3))
res_diftest3

## -----------------------------------------------------------------------------
res_diftest3 <- dif.test(est.mods = list(mod1, mod2, mod3), reference = 2, 
                         method = "Haebara", purification = TRUE)
res_diftest3

## -----------------------------------------------------------------------------
data(est3pl)
test <- paste("test", 1:5, sep = "")
mod3pl <- modIRT(coef = est3pl$coef, var = est3pl$var, names = test, display = FALSE)
direclist3pl <- alldirec(mods = mod3pl, method = "Haebara")
pth3 <- paste("test", c(1:5,1), sep = "")
chainec_circle <- chainec(direclist = direclist3pl, pths = pth3)
summary(chainec_circle)
id.test(chainec_circle)

## -----------------------------------------------------------------------------
pth3 <- paste("test", 1:5, sep = "")
chainec3 <- chainec(direclist = direclist3pl, pths = pth3)
ecall <- c(chainec3, direclist3pl["test1.test5"])
summary(chainec3)
summary(direclist3pl$test1.test5)
sd.test(ecall)

