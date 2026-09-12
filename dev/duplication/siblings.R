# Two detectors for CONCEPTUAL duplication -- functions that answer the same
# question. Neither looks at code text, so both find independently-written
# pairs a clone detector cannot see.
#
#   (1) INTERFACE SIBLINGS. Functions whose formal-argument sets overlap
#       heavily, where neither calls the other. In a codebase with consistent
#       argument naming, "takes (fit, samples, cells, layout)" is a strong
#       statement about what a function is FOR -- two functions with the same
#       inputs and no wrapper relation are two answers to one question.
#
#   (2) ROUTE SIBLINGS. ctsem's fork idiom is an early return:
#           f <- function(fit, ...) {
#             if (inherits(fit,'ctJuliaFit')) return(g(fit, ...))
#             ...the other implementation, inline...
#           }
#       By construction f-minus-the-guard and g answer the same question: the
#       caller asked f, and either g or f's own body supplies the answer. This
#       finds every such pair and names the predicate that selects between them.

files <- list.files("R", pattern = "\\.[Rr]$", full.names = TRUE)

defs <- list(); asts <- list()
for (f in files) {
  ex <- try(parse(f, keep.source = TRUE), silent = TRUE)
  if (inherits(ex, "try-error")) next
  asts[[f]] <- ex
  refs <- attr(ex, "srcref")
  for (i in seq_along(ex)) {
    e <- ex[[i]]
    if (is.call(e) && length(e) >= 3 &&
        as.character(e[[1]])[1] %in% c("<-", "=") &&
        is.call(e[[3]]) && as.character(e[[3]][[1]])[1] == "function") {
      # e[[3]] is `function`(formals, body); formals are element 2.
      fm <- names(as.list(e[[3]][[2]]))
      defs[[length(defs) + 1L]] <- list(name = paste(deparse(e[[2]]), collapse=""),
        file = basename(f), line = refs[[i]][1],
        formals = setdiff(fm, c("", "...")),
        fn = e[[3]],
        body = paste(deparse(e[[3]]), collapse = "\n"))
    }
  }
}
cat("function definitions found:", length(defs), "\n")
cat("with >= 3 named formals:", sum(vapply(defs, function(d) length(d$formals) >= 3, TRUE)), "\n")

# ---- (1) interface siblings ------------------------------------------------
# fixed=TRUE rather than a regex: ctsem names are full of dots and brackets,
# and escaping them is how this went wrong the first time.
calls_inside <- function(d, other)
  grepl(paste0(other, "("), d$body, fixed = TRUE) ||
  grepl(paste0(other, " ("), d$body, fixed = TRUE)

cat("\n== INTERFACE SIBLINGS (same inputs, neither wraps the other) ==\n")
cand <- Filter(function(d) length(d$formals) >= 3, defs)
rows <- list()
for (i in seq_along(cand)) for (j in seq_len(i - 1L)) {
  a <- cand[[i]]; b <- cand[[j]]
  inter <- intersect(a$formals, b$formals); uni <- union(a$formals, b$formals)
  jac <- length(inter) / length(uni)
  if (jac < 0.6 || length(inter) < 3) next
  if (calls_inside(a, b$name) || calls_inside(b, a$name)) next
  rows[[length(rows)+1L]] <- data.frame(
    A = sprintf("%s (%s:%d)", a$name, a$file, a$line),
    B = sprintf("%s (%s:%d)", b$name, b$file, b$line),
    J = jac, shared = paste(inter, collapse=","),
    cross = a$file != b$file, stringsAsFactors = FALSE)
}
if (length(rows)) {
  r <- do.call(rbind, rows); r <- r[order(-r$cross, -r$J), ]
  cat(sprintf("%-40s %-40s %5s  %s\n","function A","function B","J","shared formals"))
  cat(strrep("-",120),"\n")
  for (k in seq_len(min(nrow(r), 28L)))
    cat(sprintf("%-40s %-40s %5.2f  %s\n", r$A[k], r$B[k], r$J[k], substr(r$shared[k],1,52)))
  cat("\ntotal pairs:", nrow(r), "  cross-file:", sum(r$cross), "\n")
} else cat("none\n")

# ---- (2) route siblings ----------------------------------------------------
# Find `if (<pred>) return(<call>)` anywhere in a function body.
guards <- list()
scan_guard <- function(e, owner) {
  if (!is.call(e)) return(invisible(NULL))
  if (identical(as.character(e[[1]])[1], "if") && length(e) >= 3) {
    pred <- paste(deparse(e[[2]]), collapse = " ")
    arm <- e[[3]]
    # unwrap a one-statement brace
    if (is.call(arm) && identical(as.character(arm[[1]])[1], "{") && length(arm) == 2L)
      arm <- arm[[2]]
    if (is.call(arm) && identical(as.character(arm[[1]])[1], "return") && length(arm) == 2L) {
      inner <- arm[[2]]
      if (is.call(inner) && is.name(inner[[1]])) {
        guards[[length(guards)+1L]] <<- list(owner = owner, pred = pred,
          target = as.character(inner[[1]]))
      }
    }
  }
  for (k in seq_along(e)[-1]) if (is.call(e[[k]])) scan_guard(e[[k]], owner)
  invisible(NULL)
}
for (d in defs) scan_guard(d$fn, sprintf("%s (%s:%d)", d$name, d$file, d$line))

known <- vapply(defs, function(d) d$name, character(1))
cat("\n== ROUTE SIBLINGS: `if (pred) return(g(...))` -- f's own body is the other answer ==\n")
cat(sprintf("%-44s %-34s %s\n", "entry point f", "delegates to g", "selecting predicate"))
cat(strrep("-", 132), "\n")
n <- 0L
for (g in guards) {
  if (!(g$target %in% known)) next            # g must be a ctsem function
  if (nchar(g$pred) > 62) g$pred <- paste0(substr(g$pred, 1, 59), "...")
  cat(sprintf("%-44s %-34s %s\n", g$owner, g$target, g$pred))
  n <- n + 1L
}
cat("\nroute-sibling pairs:", n, "\n\nDONE\n")
