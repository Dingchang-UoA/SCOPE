annotate_PSEA <- function(mSet, mummichog.lib.new) {
  lapply(mSet[["path.hits"]], function(x) {
    data.frame(
      cpd_id = paste0(x, collapse = ";"),
      cpd_names = paste0(mummichog.lib.new$cpd.lib$name[match(x, mummichog.lib.new$cpd.lib$id)], collapse = ";")
    )
  }) -> res
  res <- do.call(rbind, res)
  cbind(mSet[["mummi.resmat"]], res) -> resdf
  return(resdf)
}