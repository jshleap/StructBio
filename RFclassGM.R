# TODO: 
# classifyGM in R
# Author: jshleap
#################################################################################
require(randomForest)
require(MASS)
require(outliers)
require(ellipse)
require(fields)
require(ggplot2)
require(varSelRF)
require(smacof)

###################################################################################################
## Convert GMfile into 'shapes'-readable array
## A is the read.table variable of a GMfile
## s is the sample size
## k is the number of landmarks
## d is the number of dimensions (2D,2;3D,3)
arr <- function(A,d){s<-dim(A)[1]; k<-(dim(A)[2]/d); Data <- as.numeric(A);
                     arr<-array(Data,c(s, d, k)); arr<-aperm(arr,c(3,2,1)); 
                     return(arr)} 
###################################################################################################
## Mean shape coordinates##
mshape <- function(A){apply(A, c(1,2), mean)}
###################################################################################################
## Interlandmark distances between configurations #
## M1: First configuration matrix of k dimensions and p landmarks.
## M2: Second configuration matrix of k dimensions and p landmarks.
ild2 <- function(M1, M2){sqrt(apply((M1-M2)^2, 1, sum))}
###################################################################################################
## Euclidean distance between mean shape and a given configuration##
##returns the diagonal vector of euclidean distances between landmarks
euc.dist <- function(A){ n <- dim(A)[3]; k <- dim(A)[1];  M<-matrix(NA,n,k); msh <- mshape(A);
                         for (i in 1:n){ for (j in 1:k){ M[i,j] <- sqrt(sum((msh[j,1]-A[j,1,i])^2,(msh[j,2]-A[j,2,i])^2,
                                                                            (msh[j,3]-A[j,3,i])^2))}}; 
                         return(M)}
##################################################################################################
## full matrix of distances FM, and the code of the fm function that returns the vectorized form
## fm of the upper ## diagonal entry of p(p-1)/2 elements.
FM<-function(M){as.matrix(dist(M))}
fm<-function(M){mat<-FM(M); mat[col(mat)<row(mat)]}
##################################################################################################
## Form difference matrix (FDM) between two configurations; FDM_M1/M2 = FM1/FM2
FDM<- function(M1,M2){ FDM<-FM(M1)/FM(M2)}
##################################################################################################
## Convert 'shapes'-readable array into GMfile 
## A is an array Shapes-type
## d is the dimensions
A2GM<- function(A,d){m<-matrix(NA,dim(A)[3],dim(A)[1]*d) ;
                     for (i in 1:dim(A)[3]){ for (j in 1:d){ m[i,seq(j,dim(m)[2],d)]<-A[,j,i]}};
                     as.data.frame(m)}
##################################################################################################
#Compute the FDM and deviations from the mean for each elemenmt in the gm file
FDmatrix <- function(A){v <- vector();n=dim(A)[3]; msh<-mshape(A)
                              for(a in seq(n)){
                                fdm  <- FDM(A[,,a],msh)
                                mfdm <- abs(fdm-median(fdm, na.rm=T))
                                rownames(mfdm)<-1:dim(msh)[1]
                                rfdm <- round(apply(mfdm,2,sum,na.rm=T),2)
                                v <- rbind(v,rfdm)
                              }
                              LIS<- round(((apply(v,2,sum,na.rm=T))/n),3)
                              return(LIS)}

FDmatrix2 <- function(A){
  v <- vector();v <- vector();n=dim(A)[3]; msh<-mshape(A)
  for(a in seq(n)){
    fdm  <- FDM(A[,,a],msh)
    mfdm <- abs(fdm-median(fdm, na.rm=T))
    rownames(mfdm)<-1:dim(msh)[1]
    rfdm <- round(apply(mfdm,2,sum,na.rm=T),2)
    v <- rbind(v,rfdm)}
  return(v)}
##################################################################################################
# transform data in a range from -1 to 1
# x is a list or vector
scal <- function(x){2*(x - min(x))/(max(x) - min(x)) - 1}
##################################################################################################


read.gm<- function(filename){D <- read.table(filename,sep=';', row.names=1)
  D <- D[!sapply(D,function(k) all(is.na(k)))]; X <- as.matrix(D)
  return(X)}

sig.var <- function(M,cls,fs){
  if(fs == 'KW'){
  #feature selection using kurskall-wallis test
  V<-vector()
  for (i in seq(dim(M)[2])){
    uav <- kruskal.test(M[,i] ~ as.factor(cls))
    
    if(!is.nan(uav$p.value)){
      if(uav$p.value < 0.05){
        V<-c(V,i)}}}
  }else{
    V <- varSelRF(M,as.factor(cls),returnFirstForest = F)
    V <- V$selected.vars
  }
  print("Selected variables:"); print(V)
  return(V)}

cM <- function(pv,tv,n=levels(as.factor(tv))){
  # giving a true vector and a predicted vector return the confusion matrix per label
  cMs<-list()
  lab <- n
  for(l in seq(lab)){tp = 0 ; fp = 0 ; fn = 0 ; tn = 0; mat <- matrix(NA,2,2); colnames(mat)<- c('R+','R-') ; row.names(mat)<- c('T+','T-')
                     TV <- tr.memvec(tv,lab[l]); PV<- tr.memvec(pv,lab[l])
                     for(t in seq(TV)){
                       if(TV[t] == '+'){
                         if(TV[t] == PV[t]){tp = tp+1}else{fn = fn+1}
                       }
                       else if(TV[t] == '-'){
                         if(TV[t] == PV[t]){tn = tn+1}else{fp = fp+1}}}
                     mat['T+','R+'] <- tp ; mat['T+','R-'] <- fp ; mat['T-','R+'] <- fn ; mat['T-','R-'] <- tn
                     cMs[[l]]<-mat}
  return(cMs)}

tr.memvec<-function(vec,label){
  # binarize a vector given a label
  nv <- vector()
  for(i in vec){
    if(i == label){nv <- c(nv,'+')}else{nv <- c(nv,'-')}
  }
  return(nv)}

fscore<- function(cMat){
  # compute the sensitivity, specificity and F-score from a confusion matrix
  sn <- cMat['T+','R+']/(cMat['T+','R+']+cMat['T-','R+'])
  sp <- cMat['T+','R+']/(cMat['T+','R+']+cMat['T+','R-'])
  Fs <- 2*((sn*sp)/(sn+sp))
  print(paste('Sn',sn,"Sp",sp,'FS',Fs,sep=' '))
  if(is.nan(sn)){sn<-0}else if(is.nan(sp)){sp<-0}else if(is.nan(Fs)){Fs<-0}
  return(Fs)}

meanFS<-function(cMs){fs<-vector()
                      for(i in seq(length(cMs))){
                        f <- fscore(cMs[[i]])
                        print(f)
                        fs <- c(fs,f)}
                      mf <- mean(fs, na.rm=TRUE)
                      #if(mf == NA){return(0)}else{return(mf)}}
                      return(mf)}

split.f<- function(mat,nspl,cls){nspl<- dim(mat)[1]*(1/nspl)
                                 # randomize rows and split the dataset in nspl
                                 spl <- list()
                                 ve <- sample(nrow(mat))
                                 shuf <- mat[ve,]
                                 shcls <-cls[ve]
                                 j<-0
                                 for(i in seq(1,dim(mat)[1],nspl)){j<-j+1
                                                                   temp <- cbind(shcls[i:(i+(nspl-1))],shuf[i:(i+(nspl-1)),])
                                                                   spl[[j]]<-temp}
                                 return(spl)}

sumCM<- function(CM1,CM2){
  L<-list()
  for(i in seq(length(CM1))){
    L[[i]]<-CM1[[i]]+CM2[[i]]
  }                       
  return(L)}

CV<- function(mat,cls,nspl=3,ntrees=500,model=''){
  # perform the crossvalidation
  crossval<-list()
  spl <- split.f(mat,nspl,cls) 
  sp <- seq(nspl)
  NFS<-vector()
  #first<-T
  instancelist<-list()
  for(i in sp){
    test<- spl[[i]] ;tr <- setdiff(sp,i);train<-vector()
    for(j in tr){train<-rbind(train,spl[[j]])}; labTrain <- train[,1]; train<- as.matrix(train[,(2:dim(train)[2])])
    class(train)<-'numeric'; labTest <- as.vector(test[,1]); test <- as.matrix(test[,(2:dim(test)[2])]); class(test)<-"numeric"
    if(model != ''){
      load(model)
      }else{
        rfr<- randomForest(train,as.factor(labTrain),ntree=ntrees,xtest=test,ytest=as.factor(labTest))
      }
    Te<-rfr$test
    instancelist[[i]]<-Te
    #FS<-rfFScore(Te)
    #if(first){SC <- FS ; first<-F}else{SC<- SC + FS}
  }
  accuracy<-rfFScore(instancelist,cls)
  #accuracy <- SC/nspl
  print(accuracy)
  accuracy<-rbind(colnames(accuracy),accuracy)
  accuracy<-cbind(rownames(accuracy),accuracy)
  write.table(accuracy,file='CrossValidation',sep='\t', row.names = F, col.names = F,quote=F)
  crossval[[1]]<-Te$confusion ; crossval[[2]]<-test ; crossval[[3]]<-labTest
  return(crossval)}


rfFScore<- function(testinstancelist,cls){
  pred = c() ; first=T
  for(i in seq(length(testinstancelist))){pred = c(pred,levels(testinstancelist[[i]]$predicted))}
  lev=unique(cls)
  for(i in seq(length(testinstancelist))){
    CM <- testinstancelist[[i]]$confusion
    if(dim(CM)[1] < length(lev)){missing = setdiff(lev,row.names(CM)) ; v <- matrix(0,1,length(lev)+1); row.names(v)<- missing; CM<-rbind(CM,v)}
    CM <- CM[order(rownames(CM)),]
    if(first){Mat <- CM; first=F}else{Mat<- Mat + CM}
  }
  L<-matrix(NA,length(unique(pred)),3)
  rownames(L)<-as.vector(unique(pred))
  colnames(L)<-c("Sn","Sp","F-score")
  #Mat<- testinstance$confusion
  for(cl in unique(pred)){
    X <- Mat[cl,cl]/sum(Mat[,cl])#sn
    Y <- Mat[cl,cl]/sum(Mat[cl,])#sp
    FS<- 2*((X*Y)/(X+Y))
    L[cl,]<-c(X,Y,FS)}
  return(L)}

storeIneq <- function(ellipse,axis,points){
  v<-numeric()
  #   v<-0
  for (p in seq(dim(points)[1])){
    x <- points[p,1] ; y<- points[p,2] ; h<-axis$center[1] ; k<- axis$center[2]; rx <- axis$minor ; ry<- axis$major
    i<-((((x-h)^2)/rx^2) + (((y-k)^2)/ry^2)); v<-c(v,i)}
  ovr <- which(v <= 1)
  return(ovr)
}

findaxis<-function(coor,center){
  A <- rbind(center,coor)
  D<-rdist(A)
  minor <- D[1,][which(D[1,] == min(D[1,2:dim(D)[2]]))] 
  major <- D[1,][which(D[1,] == max(D[1,2:dim(D)[2]]))]
  return(list(major=major,minor=minor,center=center))
}

isline <- function(A){
  if(any(A$minor < 0.01 , A$major < 0.01)){return(TRUE)}else{return(FALSE)}
}

PsinLine<- function(X,Y,Px,Py){
  line <- coef(lm(Y~X,as.data.frame(cbind(X,Y))))
  sel <- which(round(Py,3) == round(((line[2]*Px)+line[1]),3))
return(sel)}

ifpoint<- function(ell){
  #create a fake little ellipse (i.e add noise)
  x <- ell$x[1] ; y<- ell$y[1] ; dimen<- dim(ell)[1]
  w <- 0.01 ; h <-0.02
  newx = x + rnorm(dimen,sd=0.0001)
  newy = y + rnorm(dimen,sd=0.0001)
  c=cor(cbind(newx,newy))
  e<- ellipse(c,centre=c(x,y),npoints = dimen,level=1e-7)
  ell$x<-e[,1]; ell$y<-e[,2]
return(ell)}

LDAfs <- function(matrix,memvec){
  #This function will use LDA to select features and structures
  fitX <-lda(as.factor(memvec) ~ matrix, data = as.data.frame(matrix), tol = 1.0e-10) # Fit the data to LD
  tt<-t.test(fitX$scaling[,1])
  selvar <- which(!fitX$scaling[,1]<tt$conf.int[1] & !fitX$scaling[,1] > tt$conf.int[2]) # feature selection based on loadings
  predX <- data.frame(cls = predict(fitX)$class,predict(fitX)$x)# get the prediction power
  #prepare vars to store
  LX<-list()
  ellX <- data.frame()
  #loop over each level (module) and plot them with a 95% confidence ellipse
  for(g in levels(predX$cls)){
    print(paste("evaluating",g,sep=' '))
    if(all(dim(predX[predX$cls==g,]) != 0)){
      if (length(levels(predX$cls))==2){predX<-cbind(predX,LD2=predX$LD1)}
      v <- cbind(as.data.frame(with(predX[predX$cls==g,], ellipse(cor(LD1, LD2), scale=c(sd(LD1),sd(LD2)),centre=c(mean(LD1),mean(LD2)),#level = 0.99,
                                                                  npoints=1000))),Module=g, centerx=mean(predX[predX$cls==g,]$LD1),centery = mean(predX[predX$cls==g,]$LD2))
      LX[[g]]<-v
      ellX <- rbind(ellX,  v)
    }}
  Xplot <- ggplot(predX, aes(x=LD1, y=LD2, col=cls) ) + geom_point( size = 4, aes(color = cls))+theme_bw()+
    geom_path(data=ellX, aes(x=x,y=y,color=Module),size=1,linetype=2) + ggtitle('Training Set LDA')
  #png('TrainingLDA.png')
  ggsave(filename="TrainingLDA.jpg", plot=Xplot)
  #dev.off()
  
  selected<-list()
  for(h in levels(as.factor(cls))){
    #print('Selecting structures inside 95% confidence ellipses')
    gr<-LX[[h]]
    tv <- which(cls==h) ; pv <- which(predX$cls == h) ; wrong <- setdiff(pv,tv)
    dp<-with(predX[predX$cls==levels(gr$Module),],cbind(LD1,LD2))
    if((var(gr$x) < 1e-10) & (var(gr$y) < 1e-10)){gr<-ifpoint(gr)}
    A <- findaxis(cbind(gr$x,gr$y),c(gr$centerx[1],gr$centery[1]))
    if(isline(A)){
      SS<-PsinLine(gr$x,gr$y,dp[,1],dp[,2])}else{SS <- storeIneq(cbind(gr$x,gr$y),A,dp)}
    retain<- setdiff(pv[SS],wrong)
    selected[[h]]<-retain}
  return(list(selected,selvar))}

# Application of the code -------------------------------------------------
print('USAGE:')
print('Rscript RFclassGM.R gmfile clsfile [options]')
print('Options:')
print("-FS  TYPE: Perform feature selection, either kuskall-wallis-based (KW) or random forest feature selection (RF)")
print('-ntrees NTREES  : provide the number of trees (default:500)')
print('-t TESTFILE: You need to provide a GM file with the testsets aligned to the same frame than the training set')
print('-LDfs : This will perform FS as before (KW or RF should be passed), and also will curate the training set getting rid of outlier structures.')
print('-varsel: Have to be passed in conjuction with LDfs. Will do feature selection by LDA loadings.')
print('-folds FOLDS  : Provide the number of folds to do the cross validation with (Default: 3)')
print('-model MODELFILE : If a model has already been done, load it crossvalidate and predict if testfile provided')
print('-FD : Use form difference matrix instead of raw coordinates')
print('-MDS : use multidimensional scaling transformation')

args <- commandArgs(TRUE)
if('-FD' %in% args){fd<-T}else{fd<-F}
if('-FS' %in% args){fs<-args[match('-FS',args)+1]}else{fs<-F}
if('-t' %in% args){TeSt<-T}else{TeSt<-F; ntest<-0; extraT<-T }
fn<- args[1]
if(grepl('comb',fn,fixed=T) & fd){
  st<-read.gm(sub('.comb.gm','.st.gm',fn,fixed=T)); seq<-read.gm(sub('.comb.gm','.seq.gm',fn,fixed=T))
  A<-arr(st,3); msh<-mshape(A);FD<-FDmatrix2(A,dim(A)[3],msh);nam<-row.names(st); st<-as.matrix(FD); row.names(st)<-nam;st<-scal(st)
  M<-cbind(st,seq); fd<-F
  if(TeSt){Tst<-read.gm(sub('.comb.gm','.st.gm',args[match('-t',args)+1],fixed=T))
           Tseq<-read.gm(sub('.comb.gm','.seq.gm',args[match('-t',args)+1],fixed=T)) ; Te<-arr(Tst,3);FDt<-FDmatrix2(Te,dim(Te)[3],msh); 
           tnam<-row.names(Tst); Tst<-as.matrix(scal(FDt)); row.names(Tst)<-tnam; testset<-cbind(Tst,Tseq); extraT<-F}}else{M <- read.gm(fn)}
cls <- scan(args[2],what=character(0),sep=';')
cls<-as.vector(unlist(cls))
#print("Original variables:"); print(colnames(M))
if('-model' %in% args){model<-args[match('-model',args)+1]}else{model<-''}
if('-LDfs' %in% args){fs<-F;LDA<-args[match('-LDfs',args)+1]}else{LDA<-F}
if('-varsel' %in% args){vs<-T}else{vs<-F}
if('-ntrees' %in% args){ntrees <- as.numeric(args[match('-ntrees',args)+1])}else{ntrees<-500}
if('-folds' %in% args){folds <- as.numeric(args[match('-folds',args)+1])}else{folds<-3}
if(TeSt & extraT){testfile <- args[match('-t',args)+1] ;testset <- read.gm(testfile); ntest <- dim(testset)[1]}
if(fd){A<-arr(M,3); msh<-mshape(A);FD<-FDmatrix2(A,dim(A)[3],msh);nam<-row.names(M); M<-as.matrix(FD); row.names(M)<-nam;M<-scal(M)
                    if(TeSt){Te<-arr(testset,3);FDt<-FDmatrix2(Te,dim(Te)[3],msh); tnam<-row.names(testset); testset<-as.matrix(scal(FDt)); row.names(testset)<-tnam}}
if('-MDS' %in% args){
  mds<-T;if(TeSt){M<-rbind(M,testset)}; if(file.exists('fit.Rsave')){load('fit.Rsave')}else{
  #shape <- cmdscale(dist(M), k=(dim(M)[1]-ntest-1)); fit <- isoMDS(dist(M),y=shape,k=dim(shape)[2]); fit<-fit$points
  fit <- smacofSym(dist(M),ndim=(dim(M)[1]-ntest-1),type="ordinal"); save(fit,file='fit.Rsave'); fit<-fit$conf
  #fit <- cmdscale(dist(M), k=(dim(M)[1]-ntest-1)) ; 
   M<-fit[1:(dim(M)[1]-ntest),]; testset<- fit[(dim(M)[1]-ntest):dim(M)[1],]}
  }                    

if(fs != F){V<- sig.var(M,cls,fs); M<- M[,V] ; if(TeSt == T){testset<-testset[,V]}}
if(LDA != F){ #if(LDA == 'KW'){V<- sig.var(M,cls,'KW'); M<- M[,V]} else if(LDA == 'RF'){
  V<- sig.var(M,cls,LDA); M<- M[,V]#}
     if(TeSt == T){testset<-testset[,V]} ; S <- LDAfs(M,cls); sel <- S[[1]]; selvar<- S[[2]]; SV <- unlist(sel)
     if(vs){M<- M[SV,selvar]; cls <- cls[SV]; if(TeSt){testset<-testset[,selvar]}}else{M<- M[SV,]; cls <- cls[SV]}}
if(model != ''){load(model)
                }else{
                  write(colnames(M),file="usedVars.txt")
                  crossVal <- CV(M,cls,nspl=folds,ntrees)
                  rfr<- randomForest(M,as.factor(cls),ntree=ntrees)
                  save(rfr,file=paste(substr(fn,1,nchar(fn)-2),'model',sep=''))
}
if(TeSt){p <- predict(rfr,testset); write.table(p,file='prediction',sep='\t',col.names = F,quote=F)}
