for (pkg in c("crumblr", "dreamlet", "variancePartition", "limma", "edgeR")) {
  installed <- require(pkg, character.only=TRUE, quietly=TRUE)
  cat(pkg, ":", ifelse(installed, "YES", "NO"), "\n")
}
