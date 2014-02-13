# featuresel Copyright (C) 2013 Jose Sergio Hleap
#
# R script for feature selection for GM data. Will pick variables that have 
# a significant kuskal-wallis difference given a classification scheme.
# Usage: Rscript featureSel.R <gmfile> <membershipfile> <dimensions> <FDM>
#        Arguments:
#        1) File name of the GM data
#        2) File name of the membership (must include the membership for all 
#           structures in the gm file)
#        3) Number of dimensions of the structures (must be multiple of the 
#           total of variables)
#        4) a boolean (TRUE/FALSE) if the feature selection should be made on 
#           the raw coordinates (FALSE) or in the form difference matrix
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
read.gm<- function(filename){D <- read.table(filename,sep=';', row.names=1)
                       D <- D[!sapply(D,function(k) all(is.na(k)))]
                       X <- as.matrix(D)
                       return(X)}

sig.var <- function(M,cls){
  V<-vector()
  for (i in seq(dim(M)[2])){
    uav <- kruskal.test(M[,i] ~ as.factor(cls))
    if(!is.nan(uav$p.value)){
      if(uav$p.value < 0.05){
        V<-c(V,i)
      }
    }
  }
return(V)}

## Convert GMfile into 'shapes'-readable array
## A is the read.table variable of a GMfile
## s is the sample size
## k is the number of landmarks
## d is the number of dimensions (2D,2;3D,3)

arr <- function(A,d){s<-dim(A)[1]; k<-(dim(A)[2]/d); Data <- as.numeric(A);
                     arr<-array(Data,c(s, d, k)); arr<-aperm(arr,c(3,2,1)); return(arr)}

## Mean shape coordinates##
mshape <- function(A){apply(A, c(1,2), mean)}

## Interlandmark distances between configurations #
## M1: First configuration matrix of k dimensions and p landmarks.
## M2: Second configuration matrix of k dimensions and p landmarks.
ild2 <- function(M1, M2){sqrt(apply((M1-M2)^2, 1, sum))}


## Euclidean distance between mean shape and a given configuration##
##returns the diagonal vector of euclidean distances between landmarks
euc.dist <- function(A){ n <- dim(A)[3]; k <- dim(A)[1];	M<-matrix(NA,n,k); msh <- mshape(A);
                         for (i in 1:n){ for (j in 1:k){
                           M[i,j] <- sqrt(sum((msh[j,1]-A[j,1,i])^2,(msh[j,2]-A[j,2,i])^2,
                                              (msh[j,3]-A[j,3,i])^2))}}
                         return(M)}

## full matrix of distances FM, and the code of the fm function that returns the vectorized form
## fm of the upper ## diagonal entry of p(p-1)/2 elements.
FM<-function(M){as.matrix(dist(M))}
fm<-function(M){mat<-FM(M); mat[col(mat)<row(mat)]}

## Form difference matrix (FDM) between two configurations; FDM_M1/M2 = FM1/FM2
FDM<- function(M1,M2){ FDM<-FM(M1)/FM(M2)}


## Convert 'shapes'-readable array into GMfile 
## A is an array Shapes-type
## d is the dimensions
A2GM<- function(A,d){m<-matrix(NA,dim(A)[3],dim(A)[1]*d)
                     for (i in 1:dim(A)[3]){ 
                       for (j in 1:d){
                         m[i,seq(j,dim(m)[2],d)]<-A[,j,i]}}
                     return(as.data.frame(m))}

normalize <- function(x) {(x - min(x, na.rm=TRUE)+1)/((max(x,na.rm=TRUE) -
                                                 min(x, na.rm=TRUE))+1)} 
# Application of the code -------------------------------------------------

#command line arguments
args <- commandArgs(TRUE)
fn<- args[1]
dims <- as.numeric(args[3])
fdm <- args[4]
cls <- scan(args[2],what=character(0),sep=';')
cls<-as.vector(unlist(cls))
print(length(cls))
if(length(args)>4){testset <- read.gm(args[5])
                   ntest <- dim(testset)[1]
                   cls <- c(cls,rep('test',ntest))
                   }
print(length(cls))
if (fdm == TRUE){
  temp <- read.gm(fn)
  A<-arr(temp,dims) ; ms <- mshape(A) ; p <- dim(A)[1]
  M <- matrix(NA,dim(A)[3], (p*(p-1))/2); rownames(M)<-rownames(temp)
  for(i in seq(dim(A)[3])){
    v <- FDM(A[,,i],ms) ; M[i,] <- v[upper.tri(v)]}
}else{  
M <- read.gm(fn)}
print(dim(testset))
print(dim(M))
DF <- data.frame(cls=cls,data=rbind(M,testset))
#DF <- data.frame(cls=cls,data=apply(rbind(M,testset),2,normalize))
TD <- as.matrix(DF[DF$cls!='test',][,2:dim(DF)[2]])
V <- sig.var(TD,cls[which(cls != 'test')])
print('Total variables before selection:')
print(dim(M)[2])
print("Number of variables after selection:")
print(length(V))
write(V, file=paste(substr(fn,1,nchar(fn)-2),'selectedvars',sep=''),sep=';',ncolumns=length(V))
for(g in levels(DF$cls)){G <- DF[DF$cls==g,V]
                         outf<-paste(g,'.sv.gm',sep='')
                         write.table(G,file=outf,sep=';',row.names=T,col.names=F)}
# 
# if(length(args)>4){
#   testset <- read.gm(args[5])
#   if(any(is.na(testset))){
#     print("WARNING: NAs FOUND IN TEST FILE.")
#   }
#   write.table(testset[,V],file=args[5],sep=';',row.names=T,col.names=F)
# }
# print(dim(M[,V]))
# print(dim(testset[,V]))
                       