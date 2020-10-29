##### Emily's SeaBird CTD, 09/19/19, Surface-90 m

##read in data set
data <- read.csv("SeaBird_091919_1.csv")
df <- data.frame(SeaBird_091919_1)

## DO, Sal, Temp, Chl, Density overlayed
par(mar = c(8,7,2,7))

with(df, plot(df$Depth, df$OxygenPS, type="l", col="red", ylab=NA, xaxs="i", 
              yaxs="i",,xlab=NA, ylim=c(0,100), lwd=1, xlim=c(0,100),
              xaxt="n"))
x=c(0,10,20,30,40,50,60,70,80,90,100)
axis(side = 1, at=x, labels=x, las=2)
mtext(side = 2, line = 2.3, 'Dissolved Oxygen (% Sat.)')
#mtext(side = 1, line = 2.3, 'Depth (m)')
par(new = T)

with(df, plot(df$Depth, df$Turbidity, type="l", axes=F, xlab=NA, ylab=NA, ylim=c(0,2.5), 
              xlim=c(0,100), col="forestgreen", xaxs="i", yaxs="i", lwd=1))
axis(side = 4)
mtext(side = 4, line = 2, expression(paste("Turbidity (NTU)")))
par(new = T)

with(df, plot(df$Depth, df$Temp, type="l", axes=F,xlab=NA, ylab=NA, ylim=c(16,30), 
              col="blue", lwd=1, xaxs="i", yaxs="i"))
axis(side = 4, line = 3.5)
mtext(side = 4, line = 5.7, expression(paste("Temperature (",degree,"C)")))
par(new = T)

with(df, plot(df$Depth, df$Density, type="l", axes=F, xlab=NA, xaxs="i", 
              yaxs="i", ylab=NA, ylim=c(21,27), col="purple", lwd=1))
axis(side = 2, line = 4.0)
mtext(side = 2, line = 6.0, expression(paste(sigma ["T"])))

