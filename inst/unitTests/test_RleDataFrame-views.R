### Tests for summaries of views on RleDataFrame
library(genoset)
library(RUnit)

test_RleDataFrame_views <- function() {
  # input
  rle_list = list(a = Rle(1:5, rep(2, 5)), b=Rle(6:10, rep(2, 5)))
  rle_df = RleDataFrame(rle_list, row.names=LETTERS[1:10])
  rle_list = RleList(rle_list, compress=FALSE)
  ir = IRanges(start=c(1, 6), end=c(3, 9), names=c("GENE1", "GENE2"))
  mat = as.matrix(ir)
  # output
  matricize <- function(x, rownames=c("GENE1", "GENE2")) {
    y = do.call(cbind, x)
    rownames(y) = rownames
    y
  }
 hackRleViewsList <- function(rle_list, ir) {
      RleViewsList(
          lapply( rle_list, Views, start=ir )
          )
  }
  rle_view_list = hackRleViewsList(rle_list, ir)
  sum_list       = as.list(viewSums(rle_view_list))
  mean_list      = as.list(viewMeans(rle_view_list))
  min_list       = as.list(viewMins(rle_view_list))
  max_list       = as.list(viewMaxs(rle_view_list))
  which_min_list = as.list(viewWhichMins(rle_view_list))
  which_max_list = as.list(viewWhichMaxs(rle_view_list))
  # tests
  checkEquals(rangeSums(rle_df, ir, simplify=FALSE), sum_list)
  checkEquals(rangeMeans(rle_df, ir, simplify=FALSE, na.rm=TRUE), mean_list)
  checkEquals(rangeMins(rle_df, ir, simplify=FALSE), min_list)
  checkEquals(rangeMaxs(rle_df, ir, simplify=FALSE), max_list)
  checkEquals(rangeWhichMins(rle_df, ir, simplify=FALSE), which_min_list)
  checkEquals(rangeWhichMaxs(rle_df, ir, simplify=FALSE), which_max_list)
  checkEquals(rangeSums(rle_df, ir, simplify=TRUE),      matricize(sum_list))
  checkEquals(rangeMeans(rle_df, ir, simplify=TRUE, na.rm=TRUE),     matricize(mean_list))
  checkEquals(rangeMins(rle_df, ir, simplify=TRUE),      matricize(min_list))
  checkEquals(rangeMaxs(rle_df, ir, simplify=TRUE),      matricize(max_list))
  checkEquals(rangeWhichMins(rle_df, ir, simplify=TRUE), matricize(which_min_list))
  checkEquals(rangeWhichMaxs(rle_df, ir, simplify=TRUE), matricize(which_max_list))

}
