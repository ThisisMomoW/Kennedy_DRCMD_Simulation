# Loaded by run.R after the original, unchanged runtime.
for(f in c('estimator','runtime','summary','verify'))source(paste0('extensions/drcmd_eem/',f,'.R'))
eem_cli<-function(args){
 if(!length(args))stop('Usage: Rscript run.R eem init|verify|all|1|2|3|task ID|summarize [--partial]|render')
 mode<-args[1]
 if(mode=='init'){eem_initialize();return(invisible(NULL))}
 if(mode=='verify'){eem_verify();return(invisible(NULL))}
 if(mode=='summarize'){eem_summarize(partial='--partial'%in%args);return(invisible(NULL))}
 if(mode=='render'){
  status<-readLines(eem_path('comparison/STATUS.txt'));stopifnot(status[1]=='COMPLETE')
  rmarkdown::render('extensions/drcmd_eem/report.Rmd',output_format='pdf_document',
     output_file='drcmd_eem_comparison.pdf',output_dir=normalizePath(eem_path('report')),
     envir=new.env(parent=globalenv()),quiet=TRUE)
  eem_guard();return(invisible(NULL))
 }
 x<-eem_guard();b<-check_study();eem_start_observers()
 ids<-if(mode=='all')seq_len(nrow(x$jobs))else if(mode=='task'){
  if(length(args)!=2||!grepl('^[0-9]+$',args[2]))stop('Specify task ID 1..320')
  which(x$jobs$task==as.integer(args[2]))
 }else if(mode%in%c('1','2','3'))which(x$jobs$group==as.integer(mode))else integer()
 if(!length(ids))stop('Unknown EEM selection: ',paste(args,collapse=' '))
 for(i in ids)eem_one(i,b,x)
 eem_guard()
}
