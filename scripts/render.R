stopifnot(file.exists('reports/Kennedy_DR_CMD_Simulation_Report.Rmd'),file.exists('validation/AGGREGATION_VERIFIED'))
rmarkdown::render('reports/Kennedy_DR_CMD_Simulation_Report.Rmd',output_format='pdf_document',envir=new.env(parent=globalenv()),quiet=TRUE)
stopifnot(file.info('reports/Kennedy_DR_CMD_Simulation_Report.pdf')$size>10000)
cat('Rendered reports/Kennedy_DR_CMD_Simulation_Report.pdf\n')
