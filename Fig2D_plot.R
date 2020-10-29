##### Jordan's YSI, 09/19/19, Surface-90 m, Cast 2

##read in data set
data <- read.csv("Fig1CD_chemistry.csv")
df <- data.frame(Fig1CD_chemistry)

## DO, Sal, Temp, Chl, Density overlayed
par(mar = c(6,6,2,6))

# Total Fe
with(df, plot(df$Depth, df$FeT, type="o", col="red", ylab=NA, pch=16, xaxs="i", 
              yaxs="i", cex=0.7, xlab=NA, lwd=1, ylim=c(0,0.3), xlim=c(0,100),
              xaxt="n", yaxt="n"))
x=c(0,10,20,30,40,50,60,70,80,90, 100)
y=c(0, 0.1, 0.2, 0.3)
axis(side = 1, at=x, labels=x, las=2)
axis(side = 2, at=y, labels=y)
mtext(side = 2, line = 2.1, adj=0.85, expression(paste("[Fe"[d],"]"))) #(nM)
mtext(side = 1, line = 2.3)
par(new = T)

# Fe(II)
with(df, plot(df$Depth, df$Fe2, type="o", axes=F, xlab=NA, ylab=NA, pch=16, 
              xaxs="i", yaxs="i", cex=0.7, col="purple", lwd=1, xlim=c(0,100), ylim=c(0,0.3)))
#axis(side = 2, line = 4.0)
mtext(side = 2, line = 2.1, adj=0.15, expression(paste("[Fe(II)"[d],"]"))) #(nM)
par(new = T)

# total polyS (elemental sulfur)
with(df, plot(df$Depth, df$polySBOTH, type="o", axes=F, pch=16, cex=0.7, 
              xaxs="i", yaxs="i", xlab=NA, ylab=NA, xlim=c(0,100), 
              ylim=c(0,775), col="forestgreen", lwd=1))
axis(side=4)
mtext(side=4, line=2.4, adj=0.85, expression(paste("[S(0)]"))) # (",mu,"M)"
par(new = T)

# thiosulfate
with(df, plot(df$Depth, df$S2O3, type="o", axes=F, xlab=NA, pch=16, cex=0.7, 
              xaxs="i", yaxs="i", ylab=NA, xlim=c(0,100), ylim=c(0,775), col="blue1", lwd=1))
#axis(side = 4, line = 3.5)
mtext(side=4, line=2.4, adj=0.45, expression(paste("[S"[2],"O"[3]^2^-{},"]"))) # (",mu,"M)" )
par(new = T)

## S(-II)
with(df, plot(df$Depth, df$S, type="o", axes=F, pch=16, cex=0.7, xlab=NA, 
              xaxs="i", yaxs="i", ylab=NA, xlim=c(0,100), ylim=c(0,775), col="cornflowerblue", lwd=1))
#axis(side = 4, line = 7.3)
mtext(side=4, line=2.4, adj=0.05, expression(paste("[",sum(S(-II)),"]"))) # (",mu,"M)"
par(new = T)


##################################3
# polyS > 0.7 um
with(df, plot(df$Depth, df$polySLARGE, type="o", axes=F, pch=16, cex=0.8, xlab=NA, ylab=NA, ylim=c(0,600), col="darkolivegreen3", lwd=2))
axis(side = 2, line = 4.0)
mtext(side = 2, line = 6.0, expression(paste("polyS > 0.7",mu,"m (",mu,"M)")))
par(new = T)

# polyS < 0.7 um
with(df, plot(df$Depth, df$polySSMALL, type="o", axes=F, pch=16, cex=0.8, xlab=NA, ylab=NA, ylim=c(0,130), col="chartreuse4", lwd=2))
axis(side = 2, line = 7.5)
mtext(side = 2, line = 9.5, expression(paste("polyS < 0.7",mu,"m (",mu,"M)")))
