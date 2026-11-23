source_fingerprints<-function() {
  f<-c(list.files('R',pattern='\\.R$',full.names=TRUE),file.path(scenario_paths,'dgp.R'))
  tools::md5sum(f)
}
initialize_study<-function() {
  path<-'results/manifest.rds';if(file.exists(path))return(readRDS(path))
  s<-population_settings(80);hi<-population_settings(120)
  close_numeric(s$c,hi$c,1e-11);close_numeric(s$truth,hi$truth,1e-12)
  stopifnot(all(abs(s$observation_rate-.2)<1e-12))
  checks<-data.frame(scenario=s$scenario,c_error=s$c-hi$c,truth_error=s$truth-hi$truth)
  for(g in 1:3){ref<-read.csv(sprintf('reference/scenario_%d_settings.csv',g));z<-s[s$group==g,]
    close_numeric(z$c,ref$c);close_numeric(z$truth,ref$truth)
    write.csv(z,file.path(scenario_paths[g],'settings.csv'),row.names=FALSE)
  }
  write.csv(checks,'validation/integration_checks.csv',row.names=FALSE)
  jobs<-do.call(rbind,lapply(seq_len(nrow(s)),function(i)data.frame(setting=i,rep=1:200)))
  jobs$key<-sprintf('%s_r%04d',s$scenario[jobs$setting],jobs$rep)
  jobs$group<-s$group[jobs$setting];jobs$task<-ceiling(seq_len(nrow(jobs))/10)
  set.seed(20261123L);stream<-.Random.seed;seeds<-vector('list',nrow(jobs))
  for(i in seq_len(nrow(jobs))){
    ds<-parallel::nextRNGSubStream(stream);ks<-parallel::nextRNGSubStream(ds)
    seeds[[i]]<-list(data=stream,drcmd=ds,kennedy=ks);stream<-parallel::nextRNGStream(stream)
  }
  software<-software_info()
  m<-list(settings=s,jobs=jobs,seeds=seeds,seed=20261123L,sources=source_fingerprints(),
          software=software,R=R.version.string,created=as.character(Sys.time()))
  write.csv(jobs,'results/attempt_plan.csv',row.names=FALSE)
  write.csv(software,'software/versions.csv',row.names=FALSE)
  writeLines(capture.output(sessionInfo()),'software/session_info.txt')
  writeLines(capture.output(SuperLearner::SL.gam),'software/SL.gam_source.txt')
  writeLines(capture.output(args(drcmd::drcmd)),'software/drcmd_interface.txt')
  atomic_save(m,path);m
}
check_study<-function() {
  m<-readRDS('results/manifest.rds');stopifnot(identical(m$sources,source_fingerprints()))
  p<-software_info();stopifnot(identical(p[,1:3],m$software[,1:3]),identical(R.version.string,m$R));m
}
start_observers<-function() {study_start_audit();dmg_start_audit()}
checkpoint_path<-function(job,part)file.path(scenario_paths[job$group],'results/checkpoints',job$key,paste0(part,'.rds'))
kennedy_direct_fit<-function(dat,folds,upstream) {
  if(is.null(upstream$value))return(list(result=NULL,error=upstream$error,
    failure_type='upstream_kennedy_failure',warnings=upstream$warnings,seconds=0))
  t<-proc.time()[3];f<-upstream$value;pieces<-list();scores<-matrix(NA_real_,nrow(dat),2)
  for(k in seq_along(folds)) {
    ff<-folds[[k]];p<-f$fold_details[[k]];stopifnot(identical(p$train_ids,ff$train),identical(p$test_ids,ff$test))
    study_context$method<-'direct_m_gam';study_context$fold<-k
    z<-dmg_fit_fold(dat[ff$train,],dat[ff$test,],p$generated_training$lambda1)
    z$test_ids<-ff$test;pieces[[k]]<-z
    if(z$success){pr<-p$predictions
      scores[ff$test,]<-kennedy_score(dat[ff$test,],pr$kappa,pr$lambda1,pr$gamma1,(1-pr$gamma1)*z$m0,pr$gamma1*z$m1)
    }
  }
  fits<-unlist(lapply(pieces,`[[`,'fits'),recursive=FALSE)
  ok<-all(vapply(pieces,`[[`,TRUE,'success'))
  list(result=if(ok)kennedy_finish(scores)else NULL,scores=scores,fold_details=pieces,
       failure_type=if(ok)''else'weighted_gam_failure',
       error=paste(unique(Filter(nzchar,vapply(fits,`[[`,'','error'))),collapse=' | '),
       warnings=unique(unlist(lapply(fits,`[[`,'warnings'))),seconds=unname(proc.time()[3]-t))
}
run_attempt<-function(index,m=check_study()) {
  job<-m$jobs[index,];s<-m$settings[job$setting,];seeds<-m$seeds[[index]]
  path<-function(part)checkpoint_path(job,part)
  if(file.exists(path('complete')))return(invisible(NULL))
  dir.create(dirname(path('data')),recursive=TRUE,showWarnings=FALSE)
  cat(format(Sys.time()),'start',job$key,'\n');flush.console()
  # Resumption uses this exact stored dataset and method streams, never a replacement.
  if(!file.exists(path('data'))){
    set_stream(seeds$data);gen<-generate_dataset(s)
    atomic_save(gen$evaluation,path('evaluation'))
    atomic_save(list(observed=gen$observed,seed=seeds$data,hash=object_md5(gen$observed),setting=s),path('data'))
  }
  d<-readRDS(path('data'));dat<-d$observed
  stopifnot(identical(d$seed,seeds$data),identical(d$hash,object_md5(dat)),all(is.na(dat$A_obs)==(dat$R==0)))
  audit_event<-function(method)cat(format(Sys.time()),Sys.getenv('SLURM_JOB_ID'),method,'start\n',file=file.path(dirname(path('data')),'execution.log'),append=TRUE)
  if(!file.exists(path('drcmd'))){
    audit_event('drcmd');study_context$method<-'drcmd';study_audit_reset()
    e<-get('.gam_audit',.GlobalEnv);e$active<-TRUE;set_stream(seeds$drcmd)
    dr<-drcmd_point(dat);e$active<-FALSE;atomic_save(dr,path('drcmd'))
  }
  dr<-readRDS(path('drcmd'));folds<-dr$audit$folds
  # A failure before native fold creation can still leave Kennedy eligible.
  if(is.null(folds)){set_stream(seeds$drcmd);folds<-get('create_folds',asNamespace('drcmd'))(nrow(dat),5)}
  verify_folds(folds,nrow(dat));atomic_save(folds,path('folds'))
  stopifnot(dr$audit$counts['drcmd']==1,dr$audit$counts['est_varphi_eem']==0,dr$audit$counts['est_psi_tml']==0)
  if(!file.exists(path('original'))){
    audit_event('original');study_context$method<-'kennedy';study_context$fold<-0L;set_stream(seeds$kennedy)
    k<-capture_fit(kennedy_fit(dat,folds));atomic_save(k,path('original'))
  }
  k<-readRDS(path('original'))
  if(!file.exists(path('direct_m_gam'))){
    audit_event('direct_m_gam');set_stream(seeds$kennedy);before<-.Random.seed
    direct<-kennedy_direct_fit(dat,folds,k);stopifnot(identical(before,.Random.seed));atomic_save(direct,path('direct_m_gam'))
  }
  atomic_save(list(key=job$key,data_hash=d$hash,fold_hash=object_md5(folds),seeds=seeds,
    host=Sys.info()[['nodename']],job=Sys.getenv('SLURM_JOB_ID'),time=as.character(Sys.time())),path('complete'))
  cat(format(Sys.time()),'complete',job$key,'\n');flush.console();invisible(NULL)
}
