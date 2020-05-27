library("Rcpp")

nthreads=1
source("simulate_data_fxns.R")

#### Optimizing RcpArmadillo ####
library(OpenMPController )
omp_set_num_threads(1)

Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
Sys.setenv("PKG_LIBS"="-fopenmp")


#########################set OPENBLAS nthreads
require(inline)
openblas.set.num.threads <- cfunction( signature(ipt="integer"),
                                       body = 'openblas_set_num_threads(*ipt);',
                                       otherdefs = c ('extern void openblas_set_num_threads(int);'),
                                       libargs = c ('-L/opt/openblas/lib -lopenblas'),
                                       language = "C",
                                       convention = ".C"
)
openblas.set.num.threads(4)
#####################################


sourceCpp("DL_linear_split_merge_package.cpp")