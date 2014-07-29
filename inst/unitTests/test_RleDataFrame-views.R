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
  sum_list       = as.list(viewSums(Views(rle_list, ir)))
  mean_list      = as.list(viewMeans(Views(rle_list, ir)))
  min_list       = as.list(viewMins(Views(rle_list, ir)))
  max_list       = as.list(viewMaxs(Views(rle_list, ir)))
  which_min_list = as.list(viewWhichMins(Views(rle_list, ir)))
  which_max_list = as.list(viewWhichMaxs(Views(rle_list, ir)))
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
