readme.html: EFTrajectories_Demo.Rmd Simdataset1_EF.Rdata Simdataset2_EF.Rdata
	Rscript -e 'if(!require(rmarkdown)){install.packages("rmarkdown");library(rmarkdown)}; render("EFTrajectories_Demo.Rmd",output_format=github_document(), output_file="readme.md")'
