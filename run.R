# Usage: Rscript run.R init|smoke|all|1|2|3|task TASK_ID
# EEM extension: Rscript run.R eem init|verify|all|1|2|3|task ID|summarize [--partial]|render
source('R/load.R')
a<-commandArgs(TRUE);mode<-if(length(a))a[1]else'all'
if(mode=='eem'){
 source('extensions/drcmd_eem/cli.R')
 eem_cli(a[-1]);quit(status=0)
}
if(mode=='init'){initialize_study();quit(status=0)}
m<-check_study();start_observers()
ids<-if(mode=='all')seq_len(nrow(m$jobs))else if(mode=='smoke')which(m$jobs$rep==1)else if(mode=='task')which(m$jobs$task==as.integer(a[2]))else which(m$jobs$group==as.integer(mode))
stopifnot(length(ids)>0)
for(i in ids)run_attempt(i,m)
if(mode=='smoke')writeLines('16 settings completed end to end.','validation/SMOKE_COMPLETE')
