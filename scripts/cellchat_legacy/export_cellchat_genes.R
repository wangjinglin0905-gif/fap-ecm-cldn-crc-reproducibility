library(CellChat)

data(CellChatDB.human)
genes <- sort(unique(extractGeneSubsetFromPair(
  CellChatDB.human$interaction,
  complex_input = CellChatDB.human$complex,
  geneInfo = CellChatDB.human$geneInfo
)))
writeLines(genes, "work/cellchat_db_genes.txt")
cat("CellChat genes:", length(genes), "\n")
