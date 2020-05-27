library("Rcpp")

set.seed(1)
nthreads=1
source("simulate_data_fxns.R")
source("source.R")


############An illustrated example

p=5e3 # ovserved dimension
d=35  #unobservec latent dimension

nsamp=2e3 #number of samples

k=7 ####number of true cluster


probs= sample.int(10,k,replace=T); lab= sample.int(n=k,size=nsamp, prob=probs  ,replace = T)

eta=matrix(rnorm(nsamp*d), nrow=nsamp,ncol=d)

##########Generate from case (iii) Lamb model ########################
lambda=simulate_lambda(d,p,1) # simulate a sparse loading matrix lambda
for(i in 1:k){
  inds=which(lab==i)
  vars=sqrt(rchisq(d,1))
  m=  diag(vars)
  eta[inds,]= eta[inds,]%*% t(m) +i*2
}
Y=tcrossprod(eta,lambda) + matrix(rnorm(p*nsamp,sd=2),nrow=nsamp,ncol=p)
######################################################################


###############Initial values for MCMC input###########
Y=scale(Y) 

##### PCA on the data #####
pca.results=irlba::irlba(Y,nv=50)
# pca.results=irlba::irlba(Y,nv=50,center = Matrix::colMeans(Y))
cum_eigs= cumsum(pca.results$d)/sum(pca.results$d)
d=min(which(cum_eigs>.9)) ## settinf the latent dimension
###########################

eta= pca.results$u [,1:d]%*% diag(pca.results$d [1:d]) ##initial value for latent variables
lambda=pca.results$v[,1:d] #initial value for factor loading


nmix=80
kmeans_clust=kmeans(Y,nmix)
labels2=kmeans_clust$cluster-1 ##initial value for cluster indicators
rm(kmeans_clust,pca.results)
########################################################

as = 1                          # gamma hyperparameters for residual precision
bs = 0.3
a=0.5



result=DL_mixture( dir_prec=1, #Dirichlet-process precision parameter
                   diag_psi_iw=20, # the inverse scale matrix of NIW prior is diag_psi_iw \times Identity matrix
                   niw_kap=.05, # precision parameter in NIW prior
                   as=1, #shape paramter of the IG prior on residual error variance
                   bs=0.3, #scale paramter of the IG prior on residual error variance
                   a=0.5, #Dirichlet-Laplace parameter 
                   nrun=2e3, #number of efective MCMC samples
                   burn=1e2, #number of burn-in samples
                   thin=5, #thinning value
                   nstep=5, #number of Gibbs scans to perform in the split-merge sampler
                   prob=0.5, # probability of switiching between split-merge sampler and Gibbs sampler while updating the cluster indicators. 0 implies Gibbs sampler will always be selected
                   lambda, #initil value for the factor loading matrix
                   eta, #initial value of the latent variables
                   Y, #the (scaled) observed data matrix
                   del=labels2, #initial values for the cluster indicator variables
                   1,1 #two stabilizing constants rregarding the code
                   )

