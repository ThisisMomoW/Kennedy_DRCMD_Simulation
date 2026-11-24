# Render the four-method report from saved results. Never launches a simulation.
source('R/load.R');source('extensions/drcmd_eem/cli.R')
invisible(eem_guard())
stopifnot(readLines(eem_path('comparison/STATUS.txt'))[1]=='COMPLETE')
rmarkdown::render('reports/Kennedy_DRCMD_Simulation_Report.Rmd',output_format='pdf_document',
 envir=new.env(parent=globalenv()),quiet=TRUE)
pdf<-'reports/Kennedy_DRCMD_Simulation_Report.pdf'
stopifnot(file.info(pdf)$size>10000)
invisible(file.copy(pdf,eem_path('report/Kennedy_DRCMD_Simulation_Report.pdf'),overwrite=TRUE))
invisible(file.copy(pdf,eem_path('report/drcmd_eem_comparison.pdf'),overwrite=TRUE))
invisible(eem_guard())
cat('Rendered',pdf,'\n')
