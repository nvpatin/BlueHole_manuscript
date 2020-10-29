##### Emily's SeaBird CTD, 09/19/19, Surface-90 m

##read in data set
data <- read.csv("SeaBird_091919_1.csv")
df <- data.frame(SeaBird_091919_1)

## DO, Sal, Temp, Chl, Density overlayed
par(mar = c(4,6,4,6))
#par(las=1)

with(df, plot(df$Depth, df$Ph, type="l", col="darkorange3", ylab=NA, xaxs="i", 
              yaxs="i", xlab=NA, ylim=c(7.55,7.8), lwd=1, xlim=c(0,100),
              xaxt="n"))
x=c(0,10,20,30,40,50,60,70,80,90,100)
axis(side = 1, at=x, labels=x, las=2)
mtext(side = 2, line = 2.3, 'pH')
#mtext(side = 1, line = 2.3, 'Depth (m)', las=1)
par(new = T)

with(df, plot(df$Depth, df$Salinity, type="l", axes=F, xlab=NA, xlim=c(0,100), 
              ylab=NA, ylim=c(35,36.2), xaxs="i", yaxs="i", lwd=1))
axis(side = 4)
mtext(side = 4, line = 2, 'Salinity (PSU)')
par(new = T)

