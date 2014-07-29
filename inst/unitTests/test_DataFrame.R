#### Tests for RleDataFrame and  additional methods for DataFrame

test_RleDataFrame <- function() {
  foo = new("RleDataFrame", listData=list(A=Rle(1:5, rep(2, 5)), B=Rle(6:10,rep(2, 5))), nrows=10L)
  foo2 = new("RleDataFrame", listData=list(A=Rle(1, 2), B=Rle(6,2)), nrows=2L)
  foo3 = RleDataFrame(A=Rle(1:5, rep(2, 5)), B=Rle(6:10,rep(2, 5)))
  foo4 = DataFrame(A=Rle(1:5, rep(2, 5)), B=Rle(6:10,rep(2, 5)))
  checkEquals( foo3, as(foo3, "RleDataFrame"), "Coercion from DataFrame to RleDataFrame")
  checkEquals( RleList(foo3@listData, compress=FALSE), as(foo3, "RleList"), "Coercion to a List")
  checkEquals( as.list(foo3), foo3@listData, "Coercion from RleDataFrame to list")
  checkEquals( as.data.frame(foo3), data.frame(lapply(foo3, as.vector)), "data.frame coercion")
  checkEquals( as.matrix(as.data.frame(foo3)), as(foo3, "matrix"), "Matrix coercion")
  checkEquals( as.matrix(as.data.frame(foo3)), as.matrix(foo3), "Matrix coercion again")
  checkEquals( data.frame(A=foo3[[1]], B=foo3[[2]]), as.data.frame(foo3), "Coercion to data.frame")
  checkEquals( foo[1:2, ], foo2, "Row subset")
  checkEquals( foo, foo3, "Create RleDataFrame with new or with RleDataFrame" )
  checkException(new("RleDataFrame", listData=list(1:10, 1:10), nrows=10L), silent=TRUE)
}

test_rowMeans_and_rowSums <- function() {
  foo = new("RleDataFrame", listData=list(A=Rle(c(NA, 2:3, NA, 5), rep(2, 5)), B=Rle(c(6:7, NA, 8:10),c(3,2,1,2,1,1))), nrows=10L)
  mat = do.call(cbind, lapply(foo, as.numeric))
  checkEquals(rowMeans(mat), as.numeric(rowMeans(foo)))
  checkEquals(rowMeans(mat, na.rm=TRUE), as.numeric(rowMeans(foo, na.rm=TRUE)))
  checkEquals(rowSums(mat), as.numeric(rowSums(foo)))
  checkEquals(rowSums(mat, na.rm=TRUE), as.numeric(rowSums(foo, na.rm=TRUE)))
}

test_colMeans_and_colSums <- function() {
  df.ds = DataFrame( a = Rle(c(5,4,3),c(2,2,2)), b = Rle(c(3,6,9),c(1,1,4)) )
  rle.df = new("RleDataFrame", listData=df.ds@listData, nrows=nrow(df.ds))
  mat.ds = matrix( c(5,5,4,4,3,3,3,6,9,9,9,9), ncol=2, dimnames=list(NULL,c("a","b")))
  checkEquals( suppressWarnings(colMeans(df.ds)), colMeans(mat.ds) )
  checkEquals( colMeans(rle.df), colMeans(mat.ds) )
}
