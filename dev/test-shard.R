# Split the test suite into balanced shards, for CI.
#
#   Rscript dev/test-shard.R 2 4     # the files in shard 2 of 4
#
# Why this exists. Tier 2 (julia-tests.yaml: NOT_CRAN=true, Julia installed)
# does not finish inside a GitHub runner's budget as one job. Seven
# consecutive runs were killed by the 120-minute cap having spent 105 minutes
# in the suite without reaching the end, and the workflow had never been green
# once. The same files measure about 70 minutes on the development laptop, and
# a runner core is roughly half its speed, so this is slowness rather than a
# hang. Four shards put each job well under the cap and the tier finishes in
# the wall clock of its slowest shard rather than the sum.
#
# How. Longest-processing-time-first: take the files in descending order of
# measured cost and give each to whichever shard is lightest so far. That is
# the standard greedy partition and it lands within a few percent of balanced
# here.
#
# dev/test-timings.csv is the measurement and it goes stale -- it is a
# snapshot, files get added, and eleven on disk are already missing from it. A
# file it does not know about is weighted at the median of the ones it does.
# That is only ever a balance question: every test file on disk lands in
# exactly one shard whatever the csv says, and the workflow asserts that the
# files it was handed are the files that ran.

ctTestShard <- function(shard, shards, dir = 'tests/testthat',
  timings = 'dev/test-timings.csv') {

  shard <- as.integer(shard); shards <- as.integer(shards)
  stopifnot(length(shard) == 1L, length(shards) == 1L, !is.na(shard),
    !is.na(shards), shards >= 1L, shard >= 1L, shard <= shards)

  files <- sort(basename(Sys.glob(file.path(dir, 'test*.[rR]'))))
  if(!length(files)) stop('no test files found under ', dir)

  secs <- rep(NA_real_, length(files))
  if(file.exists(timings)) {
    tm <- utils::read.csv(timings, stringsAsFactors = FALSE)
    secs <- tm$secs[match(files, tm$file)]
  }
  secs[is.na(secs)] <- if(all(is.na(secs))) 1 else
    stats::median(secs, na.rm = TRUE)

  load <- numeric(shards)
  bucket <- integer(length(files))
  for(i in order(secs, decreasing = TRUE)) {
    b <- which.min(load)
    bucket[i] <- b
    load[b] <- load[b] + secs[i]
  }

  mine <- files[bucket == shard]

  # testthat's `filter` is an unanchored grepl against `context_name(files)`,
  # and it hands that the full path rather than the base name -- so the
  # `^test[-_]` strip inside it does not fire and the subject is
  # "tests/testthat/test-julia-sample", not "julia-sample". Match either
  # shape: end-anchored, preceded by a path separator or the start of the
  # string, with the prefix optional. Metacharacters are escaped and the
  # result is checked against testthat's own selection in
  # `Rscript dev/test-shard.R` with no arguments, because a filter that
  # quietly matches fewer files than it was given still looks exactly like a
  # passing shard. The workflow makes the same check on the live run, by
  # comparing the files that reported results against the ones it was handed.
  nm <- sub('[.][rR]$', '', sub('^test[-_]', '', mine))
  nm <- gsub('([][{}()*+?.^$|\\\\])', '\\\\\\1', nm)
  filter <- paste0('(^|/)(test[-_])?(', paste(nm, collapse = '|'), ')$')

  sel <- sort(basename(testthat:::find_test_scripts(dir, filter = filter)))
  if(!identical(sel, sort(mine))) stop('shard ', shard, ' of ', shards,
    ': the filter selects ', length(sel), ' of ', length(mine), ' files. ',
    'Missing: ', paste(setdiff(mine, sel), collapse = ' '), '. ',
    'Extra: ', paste(setdiff(sel, mine), collapse = ' '), '.')

  list(files = mine, secs = load[shard], filter = filter)
}

# With no arguments, check the partition instead of printing one: that every
# file on disk lands in exactly one shard, across a range of shard counts.
# (Each call already checks its own filter against testthat's selection.)
ctTestShardCheck <- function(upto = 8L, dir = 'tests/testthat', ...) {
  disk <- sort(basename(Sys.glob(file.path(dir, 'test*.[rR]'))))
  for(n in seq_len(upto)) {
    got <- character()
    for(i in seq_len(n)) {
      s <- ctTestShard(i, n, dir = dir, ...)   # verifies its own selection
      if(!length(s$files)) stop('shards = ', n, ', shard ', i, ' is empty')
      got <- c(got, s$files)
    }
    if(!identical(sort(got), disk)) stop('shards = ', n,
      ': partition is not the file list; ',
      'missing ', paste(setdiff(disk, got), collapse = ' '), '; ',
      'duplicated ', paste(unique(got[duplicated(got)]), collapse = ' '))
  }
  cat(length(disk), 'files partition cleanly for 1 to', upto, 'shards\n')
  invisible(TRUE)
}

if(sys.nframe() == 0L && !interactive()) {
  a <- commandArgs(trailingOnly = TRUE)
  if(!length(a)) ctTestShardCheck() else {
    if(length(a) != 2L) stop('usage: Rscript dev/test-shard.R [<shard> <shards>]')
    s <- ctTestShard(a[1], a[2])
    cat(sprintf('shard %s/%s: %d files, %.0f s by the last measurement\n',
      a[1], a[2], length(s$files), s$secs))
    cat(paste0('  ', s$files, collapse = '\n'), '\n', sep = '')
  }
}
