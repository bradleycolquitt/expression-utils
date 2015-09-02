prep_sql = function(fname) {
  sql = scan(fname, what=character(), sep="\n", quote = "")
  sql = paste(sql, collapse=" ")
  return(sql)
}