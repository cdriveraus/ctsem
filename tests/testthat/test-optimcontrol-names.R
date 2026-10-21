# `$` on a list falls back to the one element whose name begins with what was
# asked for, so a read of an optimcontrol field that is a prefix of another
# returns the other whenever the first is absent: `optimcontrol$uncertainty`
# read a user's `uncertaintyControl` and stopped the fit in `match.arg()`,
# `restarts` would have read `restartsd` as a count, `stalltighten` the
# count `stalltightenings`. Every such read takes `[[`; this keeps it so.

test_that("no optimcontrol field is read by a prefix of another", {
  files <- list.files(test_path("..", "..", "R"), pattern = "[.]R$", full.names = TRUE)
  skip_if(!length(files), "package sources not available")
  # Code only: string literals (messages name fields as text) and comments
  # out.
  src <- unlist(lapply(files, function(f) {
    lines <- readLines(f, warn = FALSE)
    lines <- gsub('"([^"\\\\]|\\\\.)*"', '""', lines)
    lines <- gsub("'([^'\\\\]|\\\\.)*'", "''", lines)
    sub("#.*$", "", lines[!grepl("^\\s*#", lines)])
  }))
  names <- unique(sub("^optimcontrol[$]", "", unlist(regmatches(src,
    gregexpr("optimcontrol[$][A-Za-z_.][A-Za-z0-9_.]*", src)))))
  # A read: the field taken without being assigned on that line. The name is
  # matched possessively, so it cannot back off to a shorter prefix that the
  # assignment test would then pass.
  reads <- unlist(regmatches(src, gregexpr(
    "optimcontrol[$][A-Za-z_.][A-Za-z0-9_.]*+(?!\\s*(<-|=[^=]))", src, perl = TRUE)))
  reads <- unique(sub("^optimcontrol[$]", "", reads))
  shadowed <- reads[vapply(reads, function(r) any(startsWith(setdiff(names, r), r)), TRUE)]
  expect_identical(shadowed, character())
})
