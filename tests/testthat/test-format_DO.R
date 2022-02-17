#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
library(DOSE)
data(DO2EG)
dict_DO <- enrichDO(unlist(DO2EG), ont = "DO", pvalueCutoff = 1,
                    pAdjustMethod = "BH", minGSSize = 0, maxGSSize = 1e+10,
                    qvalueCutoff = 1, readable = FALSE)
human_DO <- format_DO(dict = dict_DO@result, all_geneIDs = dict_DO@gene,
                      orgdb = org.Hs.eg.db::org.Hs.eg.db)

test_that("format_DO() runs properly.", {
  expect_equal(length(human_DO), 2)
  expect_equal(ncol(human_DO$disease), 6)
})
