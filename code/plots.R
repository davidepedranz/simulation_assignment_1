PlotAutocorrelation <- function(autocorrelation) {
  
  # analyze data
  data.len <- length(autocorrelation)
  data.max <- max(autocorrelation)
  data.min <- min(autocorrelation)
  data.h <- max(data.max, abs(data.min))
  
  # output file
  pdf('./autocorrelation.pdf', width=16, height=9)
  
  # margins
  par(mar=c(7, 7, 1, 1))
  par(mgp=c(4, 1, 0))
  
  # plot
  plot.new()
  plot.window(xlim=c(1, data.len), ylim=c(-data.h, data.h))

  # plot autocorrelation
  x <- 1:data.len
  y <- autocorrelation
  lines(x, y, lwd=2)
  
  # axis titles
  cex.title <- 2.4
  title(xlab="samples", cex.lab=cex.title)
  title(ylab="autocorrelation", cex.lab=cex.title)
  
  # axis numbers
  cex.asis <- 2.3
  axis(1, cex.axis=cex.asis)
  axis(2, cex.axis=cex.asis)
  
  # box around the plot
  box()
  
  # end
  dev.off()
}
