rm (list = ls ())

PATH = paste("...")
source(paste(PATH, 'Initioal function.r', sep=""))
source(paste(PATH, 'Mix-CMMNE.r', sep=""))
source(paste(PATH, 'Mix-CMMNW.r', sep=""))
source(paste(PATH, 'Mix-CrSN.r', sep=""))

TT = c("EII", "VII", "EEI", "VEI", 
       "EVI", "VVI", "EEE", "VEE",
       "EVE", "VVE", "EEV", "VEV",
       "EVV", "VVV")


panel.hist <- function(x, density = T, breaks = "Sturges", 
                       hist.col = "cyan", rug = TRUE, ...) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(usr[1], usr[2], 0, 1.5))
  tax <- table(x)
  if (length(tax) < 11) {
    breaks <- as.numeric(names(tax))
    y <- tax/max(tax)
    interbreak <- min(diff(breaks)) * (length(tax) - 
                                         1)/41
    rect(breaks - interbreak, 0, breaks + interbreak, 
         y, col = hist.col)
  }
  else {
    h <- hist(x, breaks = breaks, plot = FALSE)
    breaks <- h$breaks
    nB <- length(breaks)
    y <- h$counts
    y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col = hist.col)
  }
  if (density) {
    tryd <- try(d <- density(x, na.rm = TRUE, bw = "nrd", 
                             adjust = 1.2), silent = TRUE)
    if (!inherits(tryd, "try-error")) {
      d$y <- d$y/max(d$y)
      lines(d)
    }
  }
  if (rug) 
    rug(x)
}

library(gclus)
data("wine")
Data = wine

Golub = Data[, 2:14]
class = Data[, 1]

pairs(Golub, col = class, pch = 21,
      bg = c("deeppink3", "lemonchiffon3", "deepskyblue3")[class])
####----------- Set 1 --------------#
Y = Golub
para  = rgpar(Y, 3, clus = class)  # method = "trim-kmeans"


####------------- Fi####------------- Fi####------------- Fit FM-CN ---------------------##
###-----------------------------------------------## Atiken

CN.EM1 = ContaminatedMixt::CNmixt(
  Y,
  3,
  contamination = T,
  model = TT[3], initialization = "mixt", 
  parallel = FALSE,
  verbose = F
)
CN.EM1$models[[1]]$IC$BIC
km.clus = CN.EM1$models[[1]]$group

rand.index(class, km.clus)[2]
error.rate(km.clus, class)
table(CN.EM1$models[[1]]$detection)

####------------- Fit FM-rSN ---------------------##
###-----------------------------------------------## Atiken
rSN.EM1 = MIX.CrSN.EM(Y, para = para, tol = 1e-5, Family = "rSN",
                      type = TT[1], Stp.rule = "Log.like",
                      max.in.iter = 200, max.iter = 2000, per = 1,
                      print = F, class = class)

panel.lower <- function(x, y){
  points(x,y, col = class, pch = 21,
         bg = c("deeppink3", "lemonchiffon3", "deepskyblue3")[class])
}
# Customize upper panel
upper.panel <- function(x, y){
  points(x,y, col = class, pch = 21,
         bg = c("deeppink3", "lemonchiffon3", "deepskyblue3")[rSN.EM1$post.cate])
}
# Create the plots
pairs(Golub, diag.panel = panel.hist,
      lower.panel = panel.lower,
      upper.panel = upper.panel)
table(rSN.EM1$outlier)

rSN.EM2 = MIX.CrSN.EM(Y, para = para, tol = 1e-5, Family = "rSN",
                      type = TT[2], Stp.rule = "Log.like",
                      max.in.iter = 200, max.iter = 2000, per = 1,
                      print = F, class = class)

rSN.EM3 = MIX.CrSN.EM(Y, para = para, tol = 1e-5, Family = "rSN",
                      type = TT[3], Stp.rule = "Log.like",
                      max.in.iter = 200, max.iter = 2000, per = 1,
                      print = F, class = class)

rSN.EM4 = MIX.CrSN.EM(Y, para = para, tol = 1e-5, Family = "rSN",
                      type = TT[4], Stp.rule = "Log.like",
                      max.in.iter = 200, max.iter = 2000, per = 1,
                      print = F, class = class)

rSN.EM5 = MIX.CrSN.EM(Y, para = para, tol = 1e-5, Family = "rSN",
                      type = TT[5], Stp.rule = "Log.like",
                      max.in.iter = 200, max.iter = 2000, per = 1,
                      print = F, class = class)

rSN.EM6 = MIX.CrSN.EM(Y, para = para, tol = 1e-5, Family = "rSN",
                      type = TT[6], Stp.rule = "Log.like",
                      max.in.iter = 200, max.iter = 2000, per = 1,
                      print = F, class = class)

rSN.EM7 = MIX.CrSN.EM(Y, para = para, tol = 1e-5, Family = "rSN",
                      type = TT[7], Stp.rule = "Log.like",
                      max.in.iter = 200, max.iter = 2000, per = 1,
                      print = F, class = class)

rSN.EM8 = MIX.CrSN.EM(Y, para = para, tol = 1e-5, Family = "rSN",
                      type = TT[8], Stp.rule = "Log.like",
                      max.in.iter = 200, max.iter = 2000, per = 1,
                      print = F, class = class)

rSN.EM9 = MIX.CrSN.EM(Y, para = para, tol = 1e-5, Family = "rSN",
                      type = TT[9], Stp.rule = "Log.like",
                      max.in.iter = 200, max.iter = 2000, per = 1,
                      print = F, class = class)

rSN.EM10 = MIX.CrSN.EM(Y, para = para, tol = 1e-5, Family = "rSN",
                       type = TT[10], Stp.rule = "Log.like",
                       max.in.iter = 200, max.iter = 2000, per = 1,
                       print = F, class = class)

rSN.EM11 = MIX.CrSN.EM(Y, para = para, tol = 1e-5, Family = "rSN",
                       type = TT[11], Stp.rule = "Log.like",
                       max.in.iter = 200, max.iter = 2000, per = 1,
                       print = F, class = class)

rSN.EM12 = MIX.CrSN.EM(Y, para = para, tol = 1e-5, Family = "rSN",
                       type = TT[12], Stp.rule = "Log.like",
                       max.in.iter = 200, max.iter = 2000, per = 1,
                       print = F, class = class)

rSN.EM13 = MIX.CrSN.EM(Y, para = para, tol = 1e-5, Family = "rSN",
                       type = TT[13], Stp.rule = "Log.like",
                       max.in.iter = 200, max.iter = 86, per = 1,
                       print = F, class = class)

rSN.EM14 = MIX.CrSN.EM(Y, para = para, tol = 1e-5, Family = "rSN",
                       type = TT[14], Stp.rule = "Log.like",
                       max.in.iter = 200, max.iter = 2000, per = 1,
                       print = F, class = class)

rSN.EM.out = rbind(
  c(
    rSN.EM1$para.len,
    rSN.EM1$logli,
    rSN.EM1$AIC,
    rSN.EM1$BIC,
    rSN.EM1$MCR,
    rSN.EM1$Mallows.Index,
    rSN.EM1$variation.information,
    rSN.EM1$adjusted.Rand.Index,
    rSN.EM1$time
  ),
  c(
    rSN.EM2$para.len,
    rSN.EM2$logli,
    rSN.EM2$AIC,
    rSN.EM2$BIC,
    rSN.EM2$MCR,
    rSN.EM2$Mallows.Index,
    rSN.EM2$variation.information,
    rSN.EM2$adjusted.Rand.Index,
    rSN.EM2$time
  ),
  c(
    rSN.EM3$para.len,
    rSN.EM3$logli,
    rSN.EM3$AIC,
    rSN.EM3$BIC,
    rSN.EM3$MCR,
    rSN.EM3$Mallows.Index,
    rSN.EM3$variation.information,
    rSN.EM3$adjusted.Rand.Index,
    rSN.EM3$time
  ),
  c(
    rSN.EM4$para.len,
    rSN.EM4$logli,
    rSN.EM4$AIC,
    rSN.EM4$BIC,
    rSN.EM4$MCR,
    rSN.EM4$Mallows.Index,
    rSN.EM4$variation.information,
    rSN.EM4$adjusted.Rand.Index,
    rSN.EM4$time
  ),
  c(
    rSN.EM5$para.len,
    rSN.EM5$logli,
    rSN.EM5$AIC,
    rSN.EM5$BIC,
    rSN.EM5$MCR,
    rSN.EM5$Mallows.Index,
    rSN.EM5$variation.information,
    rSN.EM5$adjusted.Rand.Index,
    rSN.EM5$time
  ),
  c(
    rSN.EM6$para.len,
    rSN.EM6$logli,
    rSN.EM6$AIC,
    rSN.EM6$BIC,
    rSN.EM6$MCR,
    rSN.EM6$Mallows.Index,
    rSN.EM6$variation.information,
    rSN.EM6$adjusted.Rand.Index,
    rSN.EM6$time
  ),
  c(
    rSN.EM7$para.len,
    rSN.EM7$logli,
    rSN.EM7$AIC,
    rSN.EM7$BIC,
    rSN.EM7$MCR,
    rSN.EM7$Mallows.Index,
    rSN.EM7$variation.information,
    rSN.EM7$adjusted.Rand.Index,
    rSN.EM7$time
  ),
  c(
    rSN.EM8$para.len,
    rSN.EM8$logli,
    rSN.EM8$AIC,
    rSN.EM8$BIC,
    rSN.EM8$MCR,
    rSN.EM8$Mallows.Index,
    rSN.EM8$variation.information,
    rSN.EM8$adjusted.Rand.Index,
    rSN.EM8$time
  ),
  c(
    rSN.EM9$para.len,
    rSN.EM9$logli,
    rSN.EM9$AIC,
    rSN.EM9$BIC,
    rSN.EM9$MCR,
    rSN.EM9$Mallows.Index,
    rSN.EM9$variation.information,
    rSN.EM9$adjusted.Rand.Index,
    rSN.EM9$time
  ),
  c(
    rSN.EM10$para.len,
    rSN.EM10$logli,
    rSN.EM10$AIC,
    rSN.EM10$BIC,
    rSN.EM10$MCR,
    rSN.EM10$Mallows.Index,
    rSN.EM10$variation.information,
    rSN.EM10$adjusted.Rand.Index,
    rSN.EM10$time
  ),
  c(
    rSN.EM11$para.len,
    rSN.EM11$logli,
    rSN.EM11$AIC,
    rSN.EM11$BIC,
    rSN.EM11$MCR,
    rSN.EM11$Mallows.Index,
    rSN.EM11$variation.information,
    rSN.EM11$adjusted.Rand.Index,
    rSN.EM11$time
  ),
  c(
    rSN.EM12$para.len,
    rSN.EM12$logli,
    rSN.EM12$AIC,
    rSN.EM12$BIC,
    rSN.EM12$MCR,
    rSN.EM12$Mallows.Index,
    rSN.EM12$variation.information,
    rSN.EM12$adjusted.Rand.Index,
    rSN.EM12$time
  ),
  c(
    rSN.EM13$para.len,
    rSN.EM13$logli,
    rSN.EM13$AIC,
    rSN.EM13$BIC,
    rSN.EM13$MCR,
    rSN.EM13$Mallows.Index,
    rSN.EM13$variation.information,
    rSN.EM13$adjusted.Rand.Index,
    rSN.EM13$time
  ),
  c(
    rSN.EM14$para.len,
    rSN.EM14$logli,
    rSN.EM14$AIC,
    rSN.EM14$BIC,
    rSN.EM14$MCR,
    rSN.EM14$Mallows.Index,
    rSN.EM14$variation.information,
    rSN.EM14$adjusted.Rand.Index,
    rSN.EM14$time
  )
)
row.names(rSN.EM.out) = TT
colnames(rSN.EM.out) = c("m",
                         "lk",
                         "AIC",
                         "BIC",
                         "MCR",
                         "Mallows",
                         "VI",
                         "ARI",
                         "Time")
rSN.EM.out


####------------- Fit FM-MMNE ---------------------##
###-----------------------------------------------##

MMNE.EM1 = MIX.CMMNE.EM(Y, para = para, tol = 1e-5,
                        type = TT[1], Stp.rule = "Atiken",
                        max.in.iter = 200, max.iter = 2000, per = 1,
                        print = F, class = class)
MMNE.EM2 = MIX.CMMNE.EM(Y, para = para, tol = 1e-5,
                        type = TT[2], Stp.rule = "Atiken",
                        max.in.iter = 200, max.iter = 2000, per = 1,
                        print = F, class = class)
MMNE.EM3 = MIX.CMMNE.EM(Y, para = para, tol = 1e-5,
                        type = TT[3], Stp.rule = "Atiken",
                        max.in.iter = 200, max.iter = 2000, per = 1,
                        print = F, class = class)
MMNE.EM4 = MIX.CMMNE.EM(Y, para = para, tol = 1e-5,
                        type = TT[4], Stp.rule = "Atiken",
                        max.in.iter = 200, max.iter = 2000, per = 1,
                        print = F, class = class)
MMNE.EM5 = MIX.CMMNE.EM(Y, para = para, tol = 1e-5,
                        type = TT[5], Stp.rule = "Atiken",
                        max.in.iter = 200, max.iter = 2000, per = 1,
                        print = F, class = class)
MMNE.EM6 = MIX.CMMNE.EM(Y, para = para, tol = 1e-5,
                        type = TT[6], Stp.rule = "Atiken",
                        max.in.iter = 200, max.iter = 2000, per = 1,
                        print = F, class = class)
MMNE.EM7 = MIX.CMMNE.EM(Y, para = para, tol = 1e-5,
                        type = TT[7], Stp.rule = "Atiken",
                        max.in.iter = 200, max.iter = 2000, per = 1,
                        print = F, class = class)
MMNE.EM8 = MIX.CMMNE.EM(Y, para = para, tol = 1e-5,
                        type = TT[8], Stp.rule = "Atiken",
                        max.in.iter = 200, max.iter = 2000, per = 1,
                        print = F, class = class)
MMNE.EM9 = MIX.CMMNE.EM(Y, para = para, tol = 1e-5,
                        type = TT[9], Stp.rule = "Atiken",
                        max.in.iter = 200, max.iter = 2000, per = 1,
                        print = F, class = class)
MMNE.EM10 = MIX.CMMNE.EM(Y, para = para, tol = 1e-5,
                         type = TT[10], Stp.rule = "Atiken",
                         max.in.iter = 200, max.iter = 2000, per = 1,
                         print = F, class = class)
MMNE.EM11 = MIX.CMMNE.EM(Y, para = para, tol = 1e-5,
                         type = TT[11], Stp.rule = "Atiken",
                         max.in.iter = 200, max.iter = 2000, per = 1,
                         print = F, class = class)
MMNE.EM12 = MIX.CMMNE.EM(Y, para = para, tol = 1e-5,
                         type = TT[12], Stp.rule = "Atiken",
                         max.in.iter = 200, max.iter = 2000, per = 1,
                         print = F, class = class)
MMNE.EM13 = MIX.CMMNE.EM(Y, para = para, tol = 1e-5,
                         type = TT[13], Stp.rule = "Atiken",
                         max.in.iter = 200, max.iter = 2000, per = 1,
                         print = T, class = class)
MMNE.EM14 = MIX.CMMNE.EM(Y, para = para, tol = 1e-5,
                         type = TT[14], Stp.rule = "Atiken",
                         max.in.iter = 200, max.iter = 2000, per = 1,
                         print = T, class = class)

MMNE.EM.out = rbind(
  c(
    MMNE.EM1$para.len,
    MMNE.EM1$logli,
    MMNE.EM1$AIC,
    MMNE.EM1$BIC,
    MMNE.EM1$MCR,
    MMNE.EM1$Mallows.Index,
    MMNE.EM1$variation.information,
    MMNE.EM1$adjusted.Rand.Index,
    MMNE.EM1$time
  ),
  c(
    MMNE.EM2$para.len,
    MMNE.EM2$logli,
    MMNE.EM2$AIC,
    MMNE.EM2$BIC,
    MMNE.EM2$MCR,
    MMNE.EM2$Mallows.Index,
    MMNE.EM2$variation.information,
    MMNE.EM2$adjusted.Rand.Index,
    MMNE.EM2$time
  ),
  c(
    MMNE.EM3$para.len,
    MMNE.EM3$logli,
    MMNE.EM3$AIC,
    MMNE.EM3$BIC,
    MMNE.EM3$MCR,
    MMNE.EM3$Mallows.Index,
    MMNE.EM3$variation.information,
    MMNE.EM3$adjusted.Rand.Index,
    MMNE.EM3$time
  ),
  c(
    MMNE.EM4$para.len,
    MMNE.EM4$logli,
    MMNE.EM4$AIC,
    MMNE.EM4$BIC,
    MMNE.EM4$MCR,
    MMNE.EM4$Mallows.Index,
    MMNE.EM4$variation.information,
    MMNE.EM4$adjusted.Rand.Index,
    MMNE.EM4$time
  ),
  c(
    MMNE.EM5$para.len,
    MMNE.EM5$logli,
    MMNE.EM5$AIC,
    MMNE.EM5$BIC,
    MMNE.EM5$MCR,
    MMNE.EM5$Mallows.Index,
    MMNE.EM5$variation.information,
    MMNE.EM5$adjusted.Rand.Index,
    MMNE.EM5$time
  ),
  c(
    MMNE.EM6$para.len,
    MMNE.EM6$logli,
    MMNE.EM6$AIC,
    MMNE.EM6$BIC,
    MMNE.EM6$MCR,
    MMNE.EM6$Mallows.Index,
    MMNE.EM6$variation.information,
    MMNE.EM6$adjusted.Rand.Index,
    MMNE.EM6$time
  ),
  c(
    MMNE.EM7$para.len,
    MMNE.EM7$logli,
    MMNE.EM7$AIC,
    MMNE.EM7$BIC,
    MMNE.EM7$MCR,
    MMNE.EM7$Mallows.Index,
    MMNE.EM7$variation.information,
    MMNE.EM7$adjusted.Rand.Index,
    MMNE.EM7$time
  ),
  c(
    MMNE.EM8$para.len,
    MMNE.EM8$logli,
    MMNE.EM8$AIC,
    MMNE.EM8$BIC,
    MMNE.EM8$MCR,
    MMNE.EM8$Mallows.Index,
    MMNE.EM8$variation.information,
    MMNE.EM8$adjusted.Rand.Index,
    MMNE.EM8$time
  ),
  c(
    MMNE.EM9$para.len,
    MMNE.EM9$logli,
    MMNE.EM9$AIC,
    MMNE.EM9$BIC,
    MMNE.EM9$MCR,
    MMNE.EM9$Mallows.Index,
    MMNE.EM9$variation.information,
    MMNE.EM9$adjusted.Rand.Index,
    MMNE.EM9$time
  ),
  c(
    MMNE.EM10$para.len,
    MMNE.EM10$logli,
    MMNE.EM10$AIC,
    MMNE.EM10$BIC,
    MMNE.EM10$MCR,
    MMNE.EM10$Mallows.Index,
    MMNE.EM10$variation.information,
    MMNE.EM10$adjusted.Rand.Index,
    MMNE.EM10$time
  ),
  c(
    MMNE.EM11$para.len,
    MMNE.EM11$logli,
    MMNE.EM11$AIC,
    MMNE.EM11$BIC,
    MMNE.EM11$MCR,
    MMNE.EM11$Mallows.Index,
    MMNE.EM11$variation.information,
    MMNE.EM11$adjusted.Rand.Index,
    MMNE.EM11$time
  ),
  c(
    MMNE.EM12$para.len,
    MMNE.EM12$logli,
    MMNE.EM12$AIC,
    MMNE.EM12$BIC,
    MMNE.EM12$MCR,
    MMNE.EM12$Mallows.Index,
    MMNE.EM12$variation.information,
    MMNE.EM12$adjusted.Rand.Index,
    MMNE.EM12$time
  ),
  c(
    MMNE.EM13$para.len,
    MMNE.EM13$logli,
    MMNE.EM13$AIC,
    MMNE.EM13$BIC,
    MMNE.EM13$MCR,
    MMNE.EM13$Mallows.Index,
    MMNE.EM13$variation.information,
    MMNE.EM13$adjusted.Rand.Index,
    MMNE.EM13$time
  ),
  c(
    MMNE.EM14$para.len,
    MMNE.EM14$logli,
    MMNE.EM14$AIC,
    MMNE.EM14$BIC,
    MMNE.EM14$MCR,
    MMNE.EM14$Mallows.Index,
    MMNE.EM14$variation.information,
    MMNE.EM14$adjusted.Rand.Index,
    MMNE.EM14$time
  )
)
row.names(MMNE.EM.out) = TT
colnames(MMNE.EM.out) = c("m",
                          "lk",
                          "AIC",
                          "BIC",
                          "MCR",
                          "Mallows",
                          "VI",
                          "ARI",
                          "Time")
MMNE.EM.out

####------------- Fit FM-MMNW ---------------------##
###-----------------------------------------------##

MMNW.EM1 = MIX.CMMNW.EM(Y, para = para, tol = 1e-5,
                        type = TT[1], Stp.rule = "Log.like",
                        max.in.iter = 200, max.iter = 2000, per = 1,
                        print = F, class = class)
MMNW.EM2 = MIX.CMMNW.EM(Y, para = para, tol = 1e-5,
                        type = TT[2], Stp.rule = "Log.like",
                        max.in.iter = 200, max.iter = 2000, per = 1,
                        print = F, class = class)
MMNW.EM3 = MIX.CMMNW.EM(Y, para = para, tol = 1e-5,
                        type = TT[3], Stp.rule = "Log.like",
                        max.in.iter = 200, max.iter = 2000, per = 1,
                        print = F, class = class)
MMNW.EM4 = MIX.CMMNW.EM(Y, para = para, tol = 1e-5,
                        type = TT[4], Stp.rule = "Log.like",
                        max.in.iter = 200, max.iter = 2000, per = 1,
                        print = F, class = class)
MMNW.EM5 = MIX.CMMNW.EM(Y, para = para, tol = 1e-5,
                        type = TT[5], Stp.rule = "Log.like",
                        max.in.iter = 200, max.iter = 2000, per = 1,
                        print = F, class = class)
MMNW.EM6 = MIX.CMMNW.EM(Y, para = para, tol = 1e-5,
                        type = TT[6], Stp.rule = "Log.like",
                        max.in.iter = 200, max.iter = 2000, per = 1,
                        print = F, class = class)
MMNW.EM7 = MIX.CMMNW.EM(Y, para = para, tol = 1e-5,
                        type = TT[7], Stp.rule = "Log.like",
                        max.in.iter = 200, max.iter = 2000, per = 1,
                        print = F, class = class)
MMNW.EM8 = MIX.CMMNW.EM(Y, para = para, tol = 1e-5,
                        type = TT[8], Stp.rule = "Log.like",
                        max.in.iter = 200, max.iter = 2000, per = 1,
                        print = F, class = class)
MMNW.EM9 = MIX.CMMNW.EM(Y, para = para, tol = 1e-5,
                        type = TT[9], Stp.rule = "Log.like",
                        max.in.iter = 200, max.iter = 2000, per = 1,
                        print = F, class = class)
MMNW.EM10 = MIX.CMMNW.EM(Y, para = para, tol = 1e-5,
                         type = TT[10], Stp.rule = "Log.like",
                         max.in.iter = 200, max.iter = 2000, per = 1,
                         print = F, class = class)
MMNW.EM11 = MIX.CMMNW.EM(Y, para = para, tol = 1e-5,
                         type = TT[11], Stp.rule = "Log.like",
                         max.in.iter = 200, max.iter = 2000, per = 1,
                         print = F, class = class)
MMNW.EM12 = MIX.CMMNW.EM(Y, para = para, tol = 1e-5,
                         type = TT[12], Stp.rule = "Log.like",
                         max.in.iter = 200, max.iter = 2000, per = 1,
                         print = F, class = class)
MMNW.EM13 = MIX.CMMNW.EM(Y, para = para, tol = 1e-5,
                         type = TT[13], Stp.rule = "Log.like",
                         max.in.iter = 200, max.iter = 2000, per = 1,
                         print = F, class = class)
MMNW.EM14 = MIX.CMMNW.EM(Y, para = para, tol = 1e-5,
                         type = TT[14], Stp.rule = "Log.like",
                         max.in.iter = 200, max.iter = 2000, per = 1,
                         print = F, class = class)
MMNW.EM.out = rbind(
  c(
    MMNW.EM1$para.len,
    MMNW.EM1$logli,
    MMNW.EM1$AIC,
    MMNW.EM1$BIC,
    MMNW.EM1$MCR,
    MMNW.EM1$Mallows.Index,
    MMNW.EM1$variation.information,
    MMNW.EM1$adjusted.Rand.Index,
    MMNW.EM1$time
  ),
  c(
    MMNW.EM2$para.len,
    MMNW.EM2$logli,
    MMNW.EM2$AIC,
    MMNW.EM2$BIC,
    MMNW.EM2$MCR,
    MMNW.EM2$Mallows.Index,
    MMNW.EM2$variation.information,
    MMNW.EM2$adjusted.Rand.Index,
    MMNW.EM2$time
  ),
  c(
    MMNW.EM3$para.len,
    MMNW.EM3$logli,
    MMNW.EM3$AIC,
    MMNW.EM3$BIC,
    MMNW.EM3$MCR,
    MMNW.EM3$Mallows.Index,
    MMNW.EM3$variation.information,
    MMNW.EM3$adjusted.Rand.Index,
    MMNW.EM3$time
  ),
  c(
    MMNW.EM4$para.len,
    MMNW.EM4$logli,
    MMNW.EM4$AIC,
    MMNW.EM4$BIC,
    MMNW.EM4$MCR,
    MMNW.EM4$Mallows.Index,
    MMNW.EM4$variation.information,
    MMNW.EM4$adjusted.Rand.Index,
    MMNW.EM4$time
  ),
  c(
    MMNW.EM5$para.len,
    MMNW.EM5$logli,
    MMNW.EM5$AIC,
    MMNW.EM5$BIC,
    MMNW.EM5$MCR,
    MMNW.EM5$Mallows.Index,
    MMNW.EM5$variation.information,
    MMNW.EM5$adjusted.Rand.Index,
    MMNW.EM5$time
  ),
  c(
    MMNW.EM6$para.len,
    MMNW.EM6$logli,
    MMNW.EM6$AIC,
    MMNW.EM6$BIC,
    MMNW.EM6$MCR,
    MMNW.EM6$Mallows.Index,
    MMNW.EM6$variation.information,
    MMNW.EM6$adjusted.Rand.Index,
    MMNW.EM6$time
  ),
  c(
    MMNW.EM7$para.len,
    MMNW.EM7$logli,
    MMNW.EM7$AIC,
    MMNW.EM7$BIC,
    MMNW.EM7$MCR,
    MMNW.EM7$Mallows.Index,
    MMNW.EM7$variation.information,
    MMNW.EM7$adjusted.Rand.Index,
    MMNW.EM7$time
  ),
  c(
    MMNW.EM8$para.len,
    MMNW.EM8$logli,
    MMNW.EM8$AIC,
    MMNW.EM8$BIC,
    MMNW.EM8$MCR,
    MMNW.EM8$Mallows.Index,
    MMNW.EM8$variation.information,
    MMNW.EM8$adjusted.Rand.Index,
    MMNW.EM8$time
  ),
  c(
    MMNW.EM9$para.len,
    MMNW.EM9$logli,
    MMNW.EM9$AIC,
    MMNW.EM9$BIC,
    MMNW.EM9$MCR,
    MMNW.EM9$Mallows.Index,
    MMNW.EM9$variation.information,
    MMNW.EM9$adjusted.Rand.Index,
    MMNW.EM9$time
  ),
  c(
    MMNW.EM10$para.len,
    MMNW.EM10$logli,
    MMNW.EM10$AIC,
    MMNW.EM10$BIC,
    MMNW.EM10$MCR,
    MMNW.EM10$Mallows.Index,
    MMNW.EM10$variation.information,
    MMNW.EM10$adjusted.Rand.Index,
    MMNW.EM10$time
  ),
  c(
    MMNW.EM11$para.len,
    MMNW.EM11$logli,
    MMNW.EM11$AIC,
    MMNW.EM11$BIC,
    MMNW.EM11$MCR,
    MMNW.EM11$Mallows.Index,
    MMNW.EM11$variation.information,
    MMNW.EM11$adjusted.Rand.Index,
    MMNW.EM11$time
  ),
  c(
    MMNW.EM12$para.len,
    MMNW.EM12$logli,
    MMNW.EM12$AIC,
    MMNW.EM12$BIC,
    MMNW.EM12$MCR,
    MMNW.EM12$Mallows.Index,
    MMNW.EM12$variation.information,
    MMNW.EM12$adjusted.Rand.Index,
    MMNW.EM12$time
  ),
  c(
    MMNW.EM13$para.len,
    MMNW.EM13$logli,
    MMNW.EM13$AIC,
    MMNW.EM13$BIC,
    MMNW.EM13$MCR,
    MMNW.EM13$Mallows.Index,
    MMNW.EM13$variation.information,
    MMNW.EM13$adjusted.Rand.Index,
    MMNW.EM13$time
  ),
  c(
    MMNW.EM14$para.len,
    MMNW.EM14$logli,
    MMNW.EM14$AIC,
    MMNW.EM14$BIC,
    MMNW.EM14$MCR,
    MMNW.EM14$Mallows.Index,
    MMNW.EM14$variation.information,
    MMNW.EM14$adjusted.Rand.Index,
    MMNW.EM14$time
  )
)
row.names(MMNW.EM.out) = TT
colnames(MMNW.EM.out) = c("m",
                          "lk",
                          "AIC",
                          "BIC",
                          "MCR",
                          "Mallows",
                          "VI",
                          "ARI",
                          "Time")
MMNW.EM.out



