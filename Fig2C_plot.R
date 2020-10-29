##### Jordan's YSI, 09/19/19, Surface-90 m, Cast 2

##read in data set
data <- read.csv("Fig1CD_chemistry.csv")
df <- data.frame(Fig1CD_chemistry)

## DO, Sal, Temp, Chl, Density overlayed
par(mar = c(10,8,2,8))

with(df, plot(df$Depth, df$NOx, type="o", col="red", ylab=NA, pch=16, cex=0.7, 
              xaxs="i", yaxs="i",xlab=NA, ylim=c(0,18), lwd=1, xlim=c(0,100),
              xaxt="n"))
x=c(0,10,20,30,40,50,60,70,80,90,100)
axis(side = 1, at=x, labels=x, las=2)
mtext(side = 2, line = 2.1, expression(paste("[NO"[x]^{},"] (",mu,"M)")))
#mtext(side = 1, line = 2.3, 'Depth (m)')
par(new = T)

with(df, plot(df$Depth, df$PO4, type="o", axes=F, xlab=NA, ylab=NA, pch=16, 
              cex=0.7, lwd=1, xaxs="i", yaxs="i", ylim=c(0,6), xlim=c(0,100)))
axis(side = 4)
mtext(side = 4, line = 2.4, expression(paste("[PO"[4]^3^-{}, "] (",mu,"M)")))
par(new = T)

with(df, plot(df$Depth, df$NH4, type="o", axes=F, xlab=NA, ylab=NA, pch=16, 
              cex=0.7, ylim=c(0,50), xaxs="i", yaxs="i", xlim=c(0,100), 
              col="blue", lwd=1))
axis(side = 4, line = 3.6)
mtext(side = 4, line = 6.0, expression(paste("[NH"[4]^+{},  "] (",mu,"M)")))
par(new = T)

## DIC
with(df, plot(df$Depth, df$DIC, type="o", axes=F, xlab=NA, ylab=NA, pch=16, 
              cex=0.7, ylim=c(2.0,2.8), xlim=c(0,100), xaxs="i", yaxs="i",
              col="forestgreen", lwd=1))
axis(side = 2, line = 4.0)
mtext(side = 2, line = 6.0, expression(paste("DIC (mM)")))


