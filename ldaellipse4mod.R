# ldaellipse4mod.R Copyright (C) 2013 Jose Sergio Hleap
#
# R script for Linear discriminat analysis and 95% ellipse collision test for GM data. 
# This script will classify the variables in a GM file given a cluster file with an LDA,
# compute the 95% confidence ellipses around the groups, and compute the collisions among
# all compured ellipses
# Usage: Rscript featureSel.R prefix dimensions
#        Arguments:
#        1) Prefix : The name scheme use in the analysis (normally modularity using Moduler).
#                    Two files must be available an [prefix].lms with the agglomerated correlation
#                    vector magnitudes of a shape, and [prefix].community, with a membership vector
#                    as a single line text file with the membership vector as string splitted by
#                    space.
#        2) Number of dimensions of the structures (must be multiple of the 
#           total of variables)
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#E-mail: jshleap@squalus.org
###############################################################################

# TODO: Include MANOVA?? Include a boolean when an LD is repeated and create a 
#function for overlapping lines. 
# Author: Jose Sergio Hleap; jshleap@dal.ca


# Libraries ---------------------------------------------------------------
require(MASS)
require(ellipse)
require(fields)
require(ggplot2)

# Functions ---------------------------------------------------------------
# if one is an ellipse and the other a line, find if at least one point is inside the ellipse
# migth be used for two ellipses?
is.overlap<- function(m1,m2){
  # Check is any are points and therefore na in ellipse. if so, put the center in place of x, y
  if (any(is.na(m1$x) | is.na(m1$y))){m1$x <- m1$centerx; m1$y <- m1$centery}
  if (any(is.na(m2$x) | is.na(m2$y))){m2$x <- m2$centerx; m2$y <- m2$centery}
  A1<-findaxis(cbind(m1$x,m1$y),c(m1$centerx[1],m1$centery[1]))
  A2<-findaxis(cbind(m2$x,m2$y),c(m2$centerx[1],m2$centery[1]))
  # check if both are lines
  if(all(isline(A1),isline(A2))){
    # check orientation of line if could be done by computing the slope. but this way is a good proxy and faster
    if(var(m1$x)>=var(m1$y)){f1<-m1$x}else{f1<-m1$y}
    if(var(m2$x)>=var(m2$y)){f2<-m2$x}else{f2<-m2$y}
    # find the aproximation of maxium and minimum points
    minindex1 <- which(f1 == min(f1)) ; maxindex1 <- which(f1 == max(f1))
    minindex2 <- which(f2 == min(f2)) ; maxindex2 <- which(f2 == max(f2))
    p1<- c(m1$x[minindex1[1]],m1$y[minindex1[1]]) ; p2 <- c(m1$x[maxindex1[1]],m1$y[maxindex1[1]])
    p3<- c(m2$x[minindex2[1]],m2$y[minindex2[1]]) ; p4 <- c(m2$x[maxindex2[1]],m2$y[maxindex2[1]])
    #check if lines overlap parallely (one on the other)
    if (any(any(m2$x > p1[1] & m2$x < p2[1] ) & any(m2$y > p1[2] & m2$y < p2[2]) | 
        any(m1$x > p3[1] & m1$x < p4[1] ) & any(m1$y > p3[2] & m1$y < p4[2]))){ovr<-TRUE}
    # else check if the lines cross
    else{ovr<-doLinesIntersect(p1,p2,p3,p4)}
  }
  # if they are not both lines, at least one has to be an ellipse. find if any points in ellipse
  else{
    if(isline(A1)){ovr<- storeIneq(cbind(m2$x,m2$y),A2,cbind(m1$x,m1$y))}
    else{ovr<- storeIneq(cbind(m1$x,m1$y),A1,cbind(m2$x,m2$y))}}
  return(ovr)}

# Find if any point in an ellipse, using the ellipse inequality
storeIneq <- function(ellipse,axis,points){
  v<-numeric()
  #   v<-0
  for (p in seq(dim(points)[1])){
    x <- points[p,1] ; y<- points[p,2] ; h<-axis$center[1] ; k<- axis$center[2]; rx <- axis$minor ; ry<- axis$major
    i<-((((x-h)^2)/rx^2) + (((y-k)^2)/ry^2)); v<-c(v,i)}
  ovr <- any(v <= 1)
  return(ovr)
}

# Find the major and minor axis of an ellipse.
findaxis<-function(coor,center){
  A <- rbind(center,coor)
  D<-rdist(A)
  minor <- D[1,][which(D[1,] == min(D[1,2:dim(D)[2]]))] 
  major <- D[1,][which(D[1,] == max(D[1,2:dim(D)[2]]))]
  return(list(major=major,minor=minor,center=center))
}

# code translated from the http://stackoverflow.com/users/1998582/kris python version
doLinesIntersect<- function(p1,p2,p3,p4) {
  # all ps are vectors with x and y components
  s10_x = p2[1] - p1[1] ;  s10_y = p2[2] - p1[2] 
  s32_x = p4[1] - p3[1] ;  s32_y = p4[2] - p3[2]
  denom = s10_x * s32_y - s32_x * s10_y
  # check if collinear
  if(denom == 0){return(FALSE)}else{
    denom_is_positive <- denom > 0
    s02_x = p1[1] - p3[1] ; s02_y = p1[2] - p3[2]
    s_numer = s10_x * s02_y - s10_y * s02_x
    # no collision
    if (s_numer < 0 & denom_is_positive){return(FALSE)}else{
      t_numer = s32_x * s02_y - s32_y * s02_x
      # no collision
      if(t_numer < 0 & denom_is_positive){return(FALSE)}else{
        # no collision
        if(any((s_numer > denom & denom_is_positive), (t_numer > denom & denom_is_positive))){return(FALSE)}else{
          # collision detected
          #t = t_numer / denom
          #intersection_point = [ p0[0] + (t * s10_x), p0[1] + (t * s10_y) ]
          return(TRUE)} 
      }    
    }
  }
}

# Given an axis objects, find if is a line or not
isline <- function(A){
  if(any(A$minor < 0.01 , A$major < 0.01)){return(TRUE)}else{return(FALSE)}
}


# Application of the code -------------------------------------------------

#command line argumants
args <- commandArgs(TRUE)
fn<-args[1]
dimensions<-as.numeric(args[2])

D <- read.table(paste(substr(fn,1,nchar(fn)-3),'.lms',sep=''),sep=';')
#only one dimension
X <- as.matrix(D)
#get classifier
gc <- paste(substr(fn,1,nchar(fn)-3),'.community',sep='')#,'.graphcluster',sep='')
cls <-scan(gc,what=character(0))
cls<-as.vector(unlist(cls))
# clean the data so it does not contain singletons
X <- X[which(cls != '?'),which(cls != '?')]
fitX <-lda(as.factor(cls[which(cls!='?')]) ~ X, data = as.data.frame(X), tol = 1.0e-10) # Fit the data to LD
predX <- data.frame(cls = predict(fitX)$class,predict(fitX)$x)# get the prediction power
#fitX <-lda(as.factor(cls[which(cls!='?')]) ~ t(X), data = as.data.frame(t(X))) # Fit the data to LD
if(dim(predX)[2] <=2){ predX<- cbind(predX,LD2=predX$LD1)}
#prepare vars to store
LX<-list()
ellX <- data.frame()
#loop over each level (module) and plot them with a 95% confidence ellipse
for(g in levels(predX$cls)){
  if(all(dim(predX[predX$cls==g,]) != 0)){
    if (length(levels(predX$cls))==2){predX<-cbind(predX,LD2=predX$LD1)}
    v <- cbind(as.data.frame(with(predX[predX$cls==g,], ellipse(cor(LD1, LD2), scale=c(sd(LD1),sd(LD2)),centre=c(mean(LD1),mean(LD2)),#level = 0.99,
                                                                npoints=1000))),Module=g, centerx=mean(predX[predX$cls==g,]$LD1),centery = mean(predX[predX$cls==g,]$LD2))
    LX[[g]]<-v
    ellX <- rbind(ellX,  v)
  }}
Xplot <- ggplot(predX, aes(x=LD1, y=LD2, col=cls) ) + geom_point( size = 4, aes(color = cls))+theme_bw()+
  geom_path(data=ellX, aes(x=x,y=y,color=Module),size=1,linetype=2) + ggtitle(paste('LDA of the correlation vector maginitude matrix of dataset',substr(fn,1,nchar(fn)-3)))

#get overlaps
merges<-list()
for(m1 in LX){
  for(m2 in LX){
    if(levels(m1$Module) != levels(m2$Module)){
      if(is.overlap(m1,m2)){
        merges<-c(merges,sprintf('%s,%s',as.character(levels(m1$Module)),as.character(levels(m2$Module))))}
      else{merges<-c(merges,NA)}
    }
  }
}

lapply(merges,write,'merges.txt',append=T)
png(file=paste('LDA',substr(fn,1,nchar(fn)-3),'.png',sep=''))
Xplot
dev.off()