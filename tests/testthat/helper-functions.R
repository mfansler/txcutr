expect_equal_applied <- function (x, y, fns) {
  for (fn in fns) {
    expect_equal(fn(x), fn(y))
  }
}

get_distal_tx_name <- function (gr) {
  if (all(strand(gr) == '+')) {
    gr$tx_name[which.max(end(gr))]
  } else {
    gr$tx_name[which.min(start(gr))]
  }
}
