# Source from the project root; no dependency on a thesis checkout.
stopifnot(file.exists('R/kennedy_ratio.R'))
.libPaths(unique(c(normalizePath('library'),.libPaths())))
Sys.setenv(OMP_NUM_THREADS=1,OPENBLAS_NUM_THREADS=1,MKL_NUM_THREADS=1)
for(p in c('drcmd','SuperLearner','gam'))stopifnot(requireNamespace(p,quietly=TRUE))
stopifnot(packageDescription('drcmd')$RemoteSha=='76dc26abdb921025a2710b2826e96efe9b9df800',
          packageVersion('SuperLearner')==package_version('2.0-40'),
          packageVersion('gam')==package_version('1.22-4'))
for(f in c('utils','dgp','kennedy_ratio','kennedy_direct','drcmd','audit','runner')) source(paste0('R/',f,'.R'))
RNGkind("L'Ecuyer-CMRG")
