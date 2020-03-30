library(factoextra)
source('periodo.r')
source('NCSDE2.r')
library(plot3D)
library(fields)
library(animation)
library(mclust)
library(fda)
library(clusteval)
###################################################################################################
#############  SIMULATION 1 #######################################################################

##m=30 p=p1  change p=60, 480, 960 and p=p2
m = 30
p = 1

time_start = Sys.time()
####################################
mm = m/3
grid = list( x= seq( 0,40,,40), y= seq(0,40,,40)) 
tru = c(rep(1,mm),rep(2,mm),rep(3,mm))

set.seed(1)
sim_ncsdea = 0   ## weighted a
sim_ncsdes = 0    
sim_ncsdea_no = 0 ## non-weighted a 
sim_sp = 0         ## sdf
sim_kernel = 0  
sim_sepa = 0    ## sep a

jac_ncsdea = 0   ## weighted a
jac_ncsdes = 0    
jac_ncsdea_no = 0 ## non-weighted a 
jac_sp = 0        ## sdf
jac_kernel = 0  
jac_sepa = 0   ## sep a
N_basis = 6

for(k in 1:100){
  XXs = list()
  for(i in 1:3){    #generate random fields
    for (j in 1:mm){
      if (p == 1) {obj = matern.image.cov(grid=grid, theta=0.4*i, smoothness =0.4*i, setup=TRUE)}
      if (p == 2) {obj = matern.image.cov(grid=grid, theta=0.4*i, smoothness =0.4*(4-i), setup=TRUE)}
      XXs[[(i-1)*mm+j]] = sim.rf(obj)
    }}
  init = initial2(XXs, d=4, K=N_basis, m=2, pen='diff', dim=length(Xs), KK=3) 
  a = init$A
  weight = init$sv[1:3]
  weight = weight/sum(weight)
  iter = iterate(lambda=0, eps=0.01, n.iter=30, update.lam=TRUE, combine.lam=FALSE, which.spec=NULL, trace.graph=FALSE)
  n = length(iter$A)
  B = init$B
  A = iter$A[[n]]
  Theta = iter$Theta[[n]]
  sdf = exp(B%*%Theta%*%t(A))
  
  dist0=dist(A)
  cluster0=hclust(dist0, method = 'ward.D')
  clust0=cutree(cluster0,3)
  
  for(i in 1:3){
    A[,i]=A[,i]*weight[i]
  }
  dist1 = dist(A)
  cluster1 = hclust(dist1, method = 'ward.D')
  clust1 = cutree(cluster1,3)
  
  sdf = exp(B%*%Theta%*%t(A))
  dist5 = dist(t(sdf))
  cluster5 = hclust(dist5, method = 'ward.D')
  clust5 = cutree(cluster5,3)
  
  I = matrix(0,1600,m)
  for (i in 1:m) {I[,i] = as.vector(fftshift(Re(periodogram2(XXs[[i]]))))}
  smooth1 = B%*%solve(t(B)%*%B)%*%t(B)%*%log(I)
  dist2 = dist(t(exp(smooth1)))
  cluster2 = hclust(dist2,method = 'ward.D')
  clust2 = cutree(cluster2,3)
  
  smooth2 = smooth1
  for(i in 1:m){
    smooth2[,i] = as.vector(image.smooth(matrix(I[,i], 40, 40))$z)
  }
  dist22 = dist(t(smooth2))
  cluster22 = hclust(dist22,method = 'ward.D')
  clust22 = cutree(cluster22,3)
  
  init = initial2(XXs, d=4, K=N_basis, m=2, pen='diff', dim=length(XXs), KK=3)
  ansi = list()
  ansi.Theta.tA = NULL
  for (i in 1:m) { 
    ansi[[i]] = iterate(lambda=0, 0.01, n.iter=30, update.lam=T, combine.lam=F, which.spec=i, trace.graph=F)
    ansi.Theta.tA = cbind(ansi.Theta.tA, ansi[[i]]$Theta[[length(ansi[[i]]$Theta)]]%*%t(ansi[[i]]$A[[length(ansi[[i]]$A)]])   )
  } 
  ans.A = svd(ansi.Theta.tA)$v[,1:3]%*%diag(svd(ansi.Theta.tA)$d[1:3])	
  ans.Theta = svd(ansi.Theta.tA)$u[,1:3]	
  ans.Theta.tA = ans.Theta%*%t(ans.A)
  
  dist4 = dist(ans.A)
  cluster4 = hclust(dist4,method = 'ward.D')
  clust4 = cutree(cluster4,3)
  
  sim_ncsdea = sim_ncsdea + adjustedRandIndex(clust1,tru)/100
  sim_ncsdea_no = sim_ncsdea_no + adjustedRandIndex(clust0,tru)/100
  sim_ncsdes = sim_ncsdes + adjustedRandIndex(clust5,tru)/100
  sim_sp = sim_sp + adjustedRandIndex(clust2,tru)/100
  sim_kernel = sim_kernel + adjustedRandIndex(clust22,tru)/100
  sim_sepa = sim_sepa + adjustedRandIndex(clust4,tru)/100
  jac_ncsdea = jac_ncsdea+cluster_similarity(clust1,tru)/100
  jac_ncsdea_no = jac_ncsdea_no + cluster_similarity(clust0,tru)/100
  jac_ncsdes = jac_ncsdes + cluster_similarity(clust5,tru)/100
  jac_sp = jac_sp+ cluster_similarity(clust2,tru)/100
  jac_kernel = jac_kernel + cluster_similarity(clust22,tru)/100
  jac_sepa = jac_sepa + cluster_similarity(clust4,tru)/100
  print(k)
}
time_end = Sys.time()
print(sim_ncsdea)
print(sim_ncsdea_no)
print(sim_ncsdes)
print(sim_sp)
print(sim_kernel)
print(sim_sepa)
print('----------')
print(jac_ncsdea)
print(jac_ncsdea_no)
print(jac_ncsdes)
print(jac_sp)
print(jac_kernel)
print(jac_sepa)
print('Time cost:')
print(time_end - time_start)
####################################


###################################################################################################
############# Simulation 2#########################################################################
set.seed(2)
grid = list(x=seq(0, 40, , 40), y=seq(0, 40, ,40)) 
XXs = list() 

for(i in 1:50){    #generate random fields
for (j in 1:20){
  obj=matern.image.cov( grid=grid, theta=0.5+0.05*i,smoothness =0.5+0.05*i, setup=TRUE)
  XXs[[(i-1)*20+j]]=sim.rf(obj)
}}

####indentify K
freq = get_Fourier_freqs(nrow(XXs[[1]]), ncol(XXs[[1]]))/2/pi; 
indx = 1:nrow(freq)
I = matrix(0,nr=nrow(freq), nc=1000)
for (i in 1:1000) I[,i] = as.vector(fftshift(Re(periodogram2(XXs[[i]]))))[indx]
B1 = fda::bsplineS(freq[,1], breaks=seq(min(freq[,1]), max(freq[,1]), length.out=4), norder=4);
B2 = fda::bsplineS(freq[,2], breaks=seq(min(freq[,2]), max(freq[,2]), length.out=4), norder=4);
B = mgcv::tensor.prod.model.matrix(X=list(B1,B2))
smooth = solve(t(B)%*%B)%*%t(B)
Theta.tA = (smooth%*%log(I))
den = B%*%Theta.tA
fviz_nbclust(t(den),  FUNcluster = hcut, method = "wss", k.max = 12)+
  geom_vline(xintercept = 4, linetype = 2) 


init = initial2(XXs, d=4, K=6, pen='diff', dim=length(XXs), KK=4)
itersim = iterate(lambda=0, eps=0.01, n.iter=50, update.lam=TRUE, combine.lam=FALSE, which.spec=NULL, trace.graph=FALSE)
itersimsmo = iterate13(lambda=0,lambda2=0, eps=0.01, n.iter=50, update.lam=TRUE, combine.lam=FALSE, which.spec=NULL, trace.graph=FALSE, dima=c(20,50))
n = length(itersimsmo$A)

##generate gif## need ImageMagick
for (i in 1:n)
{
  B = init$B
  Theta = itersimsmo$Theta[[i]]
  A = itersimsmo$A[[i]]
  sdf = B%*%Theta%*%t(A)
  sdf1[,i]=sdf[,i]
}
B = init$B
Theta =init$Theta
A=init $Az
sdf = B%*%Theta%*%t(A)
u = rep(0,n)
l = rep(0,n)
for(i in 1:n){
  u[i] = max(itersimsmo$A[[i]][,1])
  l[i] = min(itersimsmo$A[[i]][,1])
}
u = max(u)
l = min(l)
saveGIF(
  interval=.25,convert='convert',
  movie.name = "asim.gif",
  expr={
    for (i in 1:(n-1)) {
      for(i in 1:4){
        asmo[,i]=asmo[,i]*weight[i]
      }
      dist3 = dist(asmo)
      cluster3 = hclust(dist3,method = 'ward.D')
      clust3 = cutree(cluster3,4)
      image(t(matrix(clust3, 20, 50)), xaxt='n', yaxt='n',col=tim.colors())
    }}
)

init = initial2(XXs,d=4, K=6, pen='diff', dim=length(XXs), KK=4)
a = itersim$A[[length(itersim$A)]]
asmo = itersimsmo$A[[length(itersimsmo$A)]]
weight = init$sv[1:4]
weight = weight/sum(weight)
for(i in 1:4){
  a[,i] = a[,i]*weight[i]
  asmo[,i] = asmo[,i]*weight[i]
}
set.panel(2, 1)
dist3 = dist(a)
cluster3 = hclust(dist3,method = 'ward.D')
clust3 = cutree(cluster3, 4)
image(t(matrix(clust3, 20, 50)), xaxt='n', yaxt='n', col=tim.colors())

dist3 = dist(asmo)
cluster3 = hclust(dist3, method = 'ward.D')
clust3 = cutree(cluster3, 4)
image(t(matrix(clust3, 20, 50)), xaxt='n', yaxt='n', col=tim.colors())

