pdf<-'reports/Kennedy_DR_CMD_Simulation_Report.pdf';txt<-tempfile(fileext='.txt')
stopifnot(system2('pdftotext',c('-layout',shQuote(pdf),shQuote(txt)))==0)
text<-readChar(txt,file.info(txt)$size);unlink(txt)
stopifnot(all(vapply(c('Observation-probability heterogeneity','Outcome slope','Quadratic outcome logit','Appendix. Numerical results'),function(x)grepl(x,text,fixed=TRUE),TRUE)))
lines<-strsplit(text,'\n',fixed=TRUE)[[1]]
num<-'(-?[0-9]+[.][0-9]+)';pair<-paste0(num,'\\s+\\(',num,'\\)')
pattern<-paste0('^\\s*(500|1,000)\\s+([0-9.]+)\\s+(Ratio GAM|Direct GAM|DR-CMD)\\s+',pair,'\\s+',pair,'\\s+',pair,'\\s+([0-9]+)\\s*$')
a<-regmatches(lines,regexec(pattern,lines,perl=TRUE));a<-a[lengths(a)>0];stopifnot(length(a)==48)
a<-do.call(rbind,lapply(a,function(x)x[-1]));s<-read.csv('results/summary.csv')
s<-s[order(s$group,s$n,s$level,match(s$method,c('original','direct_m_gam','drcmd'))),]
stopifnot(all(as.numeric(gsub(',','',a[,1]))==s$n),all(as.numeric(a[,2])==s$level),all(as.numeric(a[,10])==s$ci_n))
for(i in seq_along(fields<-c('bias','bias_mcse','rmse','rmse_mcse','coverage','coverage_mcse')))stopifnot(all(a[,i+3]==sprintf('%.3f',s[[fields[i]]])))
writeLines('PASS: all three scenarios and all 48 printed numeric rows, MCSEs and CI denominators agree with saved CSVs.','validation/PDF_VERIFIED')
cat(readLines('validation/PDF_VERIFIED'),'\n')
