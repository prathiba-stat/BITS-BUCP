#set working directory here
#setwd("C:/foo/")
plots = function(x, compVal, ropeRad, maintitle, HDImass = .95, plotname){
source("./openGraphSaveGraph.R")
source("./plotPost.R")
source("./HDIofMCMC.R") 

paramMCMCsample = x
# Plot the posterior with HDI and ROPE:
postInfo = plotPost( paramMCMCsample , compVal=compVal ,
                     ROPE=compVal+c(-ropeRad,ropeRad) , showMode=TRUE ,
                     credMass=HDImass , xlab="effect size", main = 
                       maintitle, plotname = plotname)
}


# plotAreaInROPE = function( mcmcChain , maxROPEradius , compVal=0.0 ,
#                            HDImass=0.95 , ... ) {
#   ropeRadVec = seq( 0 , maxROPEradius , length=201 ) # arbitrary comb
#   areaInRope = rep( NA , length(ropeRadVec) )
#   for ( rIdx in 1:length(ropeRadVec) ) {
#     areaInRope[rIdx] = ( sum( mcmcChain > (compVal-ropeRadVec[rIdx])
#                               & mcmcChain < (compVal+ropeRadVec[rIdx]) )
#                          / length(mcmcChain) )
#   }
#   plot( ropeRadVec , areaInRope ,
#         xlab=bquote("Radius of ROPE around "*.(compVal)) ,
#         ylab="Posterior in ROPE" ,
#         type="l" , lwd=4 , col="darkred" , cex.lab=1.5 , ... )
#   # From http://www.indiana.edu/~kruschke/DoingBayesianDataAnalysis/Programs/
#   # get HDIofMCMC.R. Put it in the working directory of this script, or specify
#   # the path in the next line's source() command:
#   # source("./HDIofMCMC.R")
#   HDIlim = HDIofMCMC( mcmcChain , credMass=HDImass )
#   farHDIlim = HDIlim[which.max(abs(HDIlim-compVal))]
#   ropeRadHDI = abs(compVal-farHDIlim)
#   areaInFarHDIlim = ( sum( mcmcChain > (compVal-ropeRadHDI)
#                            & mcmcChain < (compVal+ropeRadHDI) )
#                       / length(mcmcChain) )
#   lines( c(ropeRadHDI,ropeRadHDI) , c(-0.5,areaInFarHDIlim) ,
#          lty="dashed" , col="darkred" )
#   text( ropeRadHDI , 0 ,
#         bquote( atop( .(100*HDImass)*"% HDI limit" ,
#                       "farthest from "*.(compVal) ) ) , adj=c(0.5,0) )
#   lines( c(-0.5,ropeRadHDI) ,c(areaInFarHDIlim,areaInFarHDIlim) ,
#          lty="dashed" , col="darkred" )
#   text( 0 , areaInFarHDIlim , bquote(.(signif(areaInFarHDIlim,3))) ,
#         adj=c(0,1.1) )
#   return( list( ROPEradius=ropeRadVec , areaInROPE=areaInRope ) )
# }

# Generate a fictitious MCMC posterior:
# m = 0.03
# s = (0.1-m)/3
# 
# # Plot the area in ROPE:
# ropeInfo = plotAreaInROPE( mcmcChain=paramMCMCsample , maxROPEradius=0.15 ,
#                            compVal=compVal , HDImass=HDImass )
