expect_equal_applied <- function (x, y, fns) {
  for (fn in fns) {
    expect_equal(fn(x), fn(y))
  }
}
