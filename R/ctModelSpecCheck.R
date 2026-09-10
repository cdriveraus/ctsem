# What a model matrix cell is allowed to say, checked once, where the names are
# still the user's own.
#
# A cell can be four different things, and ctsem tells them apart by what the
# text looks like rather than by anything the user declares:
#
#   a number            fixed, not estimated
#   a bare name         a free parameter; the same name twice is one parameter
#   `name|transform|…`  a free parameter with fields (see R/ctParSpec.R)
#   an expression       computed from latent states, tdpreds and PARS cells
#
# The classification is silent, so a cell that means one of these and reads as
# another produces a different model without complaint. The checks below are the
# cases where that has happened, or can. They run in `ctModelConvertOMX()` after
# the field parsing and before the cycle check, while `ctspec$param` still holds
# the text the user wrote -- `ctModelStatesAndPARS()` rewrites latent and tdpred
# names into `state[i]` and `tdpreds[rowi, j]` later, in ctFit(), by which point
# the collisions are no longer visible.

# The labels PARS cells declare. A PARS cell is how a user names a parameter
# that an expression elsewhere can reference, so these are legitimate names in
# any other cell.
.ctSpecParsLabels <- function(ctspec) {
  lab <- ctspec$param[ctspec$matrix %in% 'PARS' & !is.na(ctspec$param)]
  if (!length(lab)) return(character(0))
  lab <- trimws(vapply(lab, function(x) strsplit(x, '|', fixed = TRUE)[[1]][1],
    character(1), USE.NAMES = FALSE))
  unique(lab[nzchar(lab) & !grepl('\\W', lab)])
}

.ctSpecCellName <- function(ctspec, i) {
  paste0(ctspec$matrix[i], '[', ctspec$row[i], ',', ctspec$col[i], ']')
}

# ---------------------------------------------------------------------------

.ctCheckModelSpec <- function(ctspec, latentNames, manifestNames, TDpredNames,
  TIpredNames) {

  for (i in seq_len(nrow(ctspec))) {
    p <- ctspec$param[i]
    if (is.na(p) || !nzchar(trimws(p))) next
    p <- trimws(p)
    cell <- .ctSpecCellName(ctspec, i)

    # --- a bare latent or tdpred name -------------------------------------
    #
    # The whole cell is one name and that name belongs to a state or a tdpred.
    # Read as a parameter name it is a free parameter; read as a reference it
    # is that state's current value, and ctsem takes the second reading. The
    # two are entirely different models -- on a one-latent example, 27
    # parameters against 11 -- and nothing said which was built. An expression
    # (`eta1 * dr`) is unambiguous and says nothing here; only a cell that is
    # exactly a name is caught.
    if (p %in% latentNames) {
      warning(cell, ' is the latent state name "', p, '", so the cell takes ',
        'that state\'s value rather than being a free parameter. Write ',
        'state[', match(p, latentNames), '] if that is what you meant -- it ',
        'builds the same model without this warning -- or give the parameter ',
        'a name of its own.', call. = FALSE)
      next
    }
    if (p %in% TDpredNames) {
      warning(cell, ' is the time dependent predictor name "', p, '", so the ',
        'cell takes that predictor\'s value rather than being a free ',
        'parameter. Write tdpreds[rowi, ', match(p, TDpredNames), '] if that ',
        'is what you meant, or give the parameter a name of its own.',
        call. = FALSE)
      next
    }

    # --- a bare name for something a cell cannot reference ----------------
    #
    # A latent or tdpred name at least becomes the reference the user probably
    # wanted. These do not: the cell stays an ordinary free parameter that
    # merely happens to share a name, so the model builds, fits, and reports a
    # parameter called `age` next to a time independent predictor called `age`.
    # Nothing is numerically wrong; nothing asked for is delivered either.
    if (p %in% TIpredNames) {
      warning(cell, ' is "', p, '", which is also a time independent ',
        'predictor of this model. The cell is a free parameter of that name, ',
        'not a reference to the predictor -- a cell cannot reference one. To ',
        'let "', p, '" predict this parameter, name it in the cell\'s tipreds ',
        'field: "', p, '_par, tipreds=', p, '".', call. = FALSE)
      next
    }
    if (p %in% manifestNames) {
      warning(cell, ' is "', p, '", which is also a manifest variable of this ',
        'model. The cell is a free parameter of that name, not a reference to ',
        'the observation -- a cell cannot reference one. Reference the latent ',
        'state it loads on instead.', call. = FALSE)
      next
    }
    if (p %in% c('dt', 'time')) {
      warning(cell, ' is "', p, '". The cell is a free parameter of that name; ',
        'the time interval is not available to a matrix cell, so nothing here ',
        'varies with it. Discrete time effects over an interval come from the ',
        'continuous time parameters, not from writing the interval into a ',
        'cell.', call. = FALSE)
      next
    }

    # --- an expression that references nothing computable -----------------
    #
    # A compound cell becomes a *calculation* in the generated program, and
    # what makes it one is that it references something the program computes:
    # a latent state, a time dependent predictor, or a PARS cell. A fresh
    # parameter name inside such an expression is fine -- `lbystate * eta2 + 1`
    # declares `lbystate` and works -- because the cell is a calculation and
    # the name is registered as it is rewritten.
    #
    # An expression that references none of those is not a calculation, and
    # nothing registers anything: the whole text becomes a free parameter's
    # *label*, and its transform is built from NA, reaching the backend as the
    # unfinished string "99999 + 99999*" and dying there. So
    # `DRIFT = matrix('-exp(dr)')` -- the way most software would have you
    # write it -- never worked, and said so only at the end, in the generated
    # program's own words.
    #
    # The test is therefore about the reference, not about the names: an
    # unknown name is only a problem when it is *all* the cell has.
    if (grepl('\\W', gsub('.', '', p, fixed = TRUE))) {
      computable <- grepl('[', p, fixed = TRUE) ||
        any(vapply(c(latentNames, TDpredNames, .ctSpecParsLabels(ctspec)),
          function(nm) nzchar(nm) &&
            grepl(paste0('\\b', nm, '\\b'), p), logical(1)))
      if (!computable) {
        # Which unreachable thing was reached for, when it is one of the ones
        # people reach for. `dt` in particular reads as "scaled by the
        # interval", and a cell has no interval.
        reached <- character(0)
        for (nm in c('dt', 'time')) {
          if (grepl(paste0('\\b', nm, '\\b'), p)) {
            reached <- c(reached, paste0('"', nm,
              '" is the time interval, which a cell cannot see'))
            break
          }
        }
        for (nm in TIpredNames) {
          if (nzchar(nm) && grepl(paste0('\\b', nm, '\\b'), p)) {
            reached <- c(reached, paste0('"', nm, '" is a time independent ',
              'predictor, which acts through the cell\'s tipreds field ',
              'rather than by being referenced'))
            break
          }
        }
        for (nm in manifestNames) {
          if (nzchar(nm) && grepl(paste0('\\b', nm, '\\b'), p)) {
            reached <- c(reached, paste0('"', nm, '" is a manifest variable, ',
              'which a cell cannot see -- reference the latent state instead'))
            break
          }
        }
        stop(cell, ' is the expression "', p, '", which references no latent ',
          'state, time dependent predictor or PARS cell, so there is nothing ',
          'for it to be computed from. ',
          if (length(reached)) paste0('Note that ',
            paste(reached, collapse = ', and '), '. ') else '',
          'To transform a single parameter, put the transform in the second ',
          'field and call the parameter param there: "dr|-exp(param)" rather ',
          'than "-exp(dr)". To build a cell from several parameters, name each ',
          'one in a PARS cell and reference it from here. For a fixed cell, ',
          'write the number itself.', call. = FALSE)
      }
    }
  }

  # --- an index that is outside the thing it indexes ---------------------
  #
  # `state[9]` in a one-latent model and `PARS[2,2]` against a 1x1 PARS both
  # reported "subscript out of bounds" from inside the matrix unfolding, which
  # names neither the cell nor the bound. The reference is the user's text, so
  # this is the place that can still say where it came from.
  nlatent <- length(latentNames)
  parsrows <- which(ctspec$matrix %in% 'PARS')
  parsdim <- if (length(parsrows)) c(max(ctspec$row[parsrows]),
    max(ctspec$col[parsrows])) else c(0L, 0L)
  for (i in seq_len(nrow(ctspec))) {
    p <- ctspec$param[i]
    if (is.na(p) || !grepl('[', p, fixed = TRUE)) next
    cell <- .ctSpecCellName(ctspec, i)
    for (ref in regmatches(p,
      gregexpr('\\bstate[[:space:]]*\\[[[:space:]]*[0-9]+[[:space:]]*\\]',
        p))[[1]]) {
      k <- as.integer(gsub('[^0-9]', '', ref))
      if (k > nlatent || k < 1L) {
        stop(cell, ' references ', ref, ', but the model has ', nlatent,
          ' latent state', if (nlatent == 1) '' else 's', ' (',
          paste(latentNames, collapse = ', '), ').', call. = FALSE)
      }
    }
    for (ref in regmatches(p,
      gregexpr(paste0('\\bPARS[[:space:]]*\\[[[:space:]]*[0-9]+[[:space:]]*,',
        '[[:space:]]*[0-9]+[[:space:]]*\\]'), p))[[1]]) {
      rc <- as.integer(strsplit(gsub('[^0-9,]', '', ref), ',', fixed = TRUE)[[1]])
      if (!length(parsrows)) {
        stop(cell, ' references ', ref, ', but the model has no PARS matrix. ',
          'Give ctModel() a PARS argument declaring the parameters the ',
          'expression needs.', call. = FALSE)
      }
      if (rc[1] > parsdim[1] || rc[2] > parsdim[2]) {
        stop(cell, ' references ', ref, ', but PARS is ', parsdim[1], 'x',
          parsdim[2], '.', call. = FALSE)
      }
    }
  }

  # --- one parameter, two different sdscales -----------------------------
  #
  # Repeating a name is ctsem's equality constraint, so the cells are one
  # parameter and it has one population sd prior. Two cells asking for
  # different scales is resolved by position -- whichever row is reached first
  # -- so the same model with its cells written in the other order gets a
  # different prior, and nothing reports which was used.
  free <- !is.na(ctspec$param) & is.na(ctspec$value)
  simple <- free & !grepl('\\W', gsub('.', '', ctspec$param, fixed = TRUE))
  for (nm in unique(ctspec$param[simple])) {
    rows <- which(simple & ctspec$param %in% nm)
    if (length(rows) < 2) next
    scales <- unique(ctspec$sdscale[rows][!is.na(ctspec$sdscale[rows])])
    if (length(scales) > 1) {
      warning('Parameter "', nm, '" is given sdscale ',
        paste(scales, collapse = ' and '), ' by different cells (',
        paste(vapply(rows, function(r) .ctSpecCellName(ctspec, r),
          character(1)), collapse = ', '),
        '). Those cells are one parameter, so it has one prior scale: ',
        ctspec$sdscale[rows[1]], ' is used, from the first cell. Set the same ',
        'sdscale in every cell that names it.', call. = FALSE)
    }
  }

  # --- sdscale where nothing varies --------------------------------------
  #
  # sdscale multiplies the prior on a population standard deviation, so it
  # means nothing for a parameter that has no population distribution. It was
  # accepted and dropped.
  indcols <- grep('^indvarying', names(ctspec), value = TRUE)
  for (i in which(free)) {
    sd <- ctspec$sdscale[i]
    if (is.na(sd) || sd == 1) next
    anyvary <- any(vapply(indcols, function(cl)
      isTRUE(as.logical(ctspec[[cl]][i])), logical(1)))
    if (!anyvary) {
      warning(.ctSpecCellName(ctspec, i), ' sets sdscale ', sd,
        ' but the parameter does not vary over subjects, so there is no ',
        'population standard deviation for it to scale. Set indvarying=TRUE ',
        'as well, or drop the sdscale.', call. = FALSE)
    }
  }

  invisible(NULL)
}
