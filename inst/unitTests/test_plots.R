library(RUnit)

test_genoPlot <- function() {
    data(genoset)
    checkTrue( genoPlot( x=genoset.ds,y=genoset.ds[,1,"lrr"] ) )
    checkTrue( genoPlot( genoPos(genoset.ds), genoset.ds[,1,"lrr"], locs=locData(genoset.ds) ) )
    checkTrue( genoPlot( 1:10, Rle(c(rep(0,5),rep(3,4),rep(1,1))) ) )
}
