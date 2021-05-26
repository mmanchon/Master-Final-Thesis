Rscript GEO_extraction.R GSE12059_series_matrix GPL887-20438 GENE_SYMBOL
Rscript GEO_extraction.R GSE40349_series_matrix GPL570-55999
Rscript GEO_extraction.R GSE50911_series_matrix GPL4133-12599 GENE_SYMBOL
Rscript GEO_extraction.R GSE55942_series_matrix GPL571-17391
Rscript GEO_extraction.R GSE115406_series_matrix GPL570-55999 GENE_SYMBOL
Rscript GEO_extraction.R GSE12764_series_matrix GPL570-55999 GENE_SYMBOL

# Ejecutar Overexpression
Rscript GSE12059_fold_change_over.R 2000
Rscript GSE12764_fold_change_over.R 2000
Rscript GSE50911_fold_change_over.R 2000
Rscript GSE115406_fold_change_over.R 2000


Rscript GEO_extraction.R GSE40349_series_matrix GPL570-55999
Rscript GEO_extraction.R GSE55942_series_matrix GPL571-17391
Rscript GEO_extraction.R GSE115406_series_matrix GPL570-55999 GENE_SYMBOL
Rscript GEO_extraction.R GSE46806_series_matrix GPL96-57554
#Rscript GEO_extraction.R GSE85170_series_matrix GPL18573_family GENE_SYMBOL

Rscript GSE46806_fold_change_down.R 2000
Rscript GSE40349_fold_change_down.R 2000
Rscript GSE115406_fold_change_down.R 2000
Rscript GSE55942_fold_change_down.R 2000

