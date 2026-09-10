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
# cases where that has happened, or can.
#
# They run from `ctModelConvertOMX()` -- which, despite the name, is where every
# ct/dt model is built: `ctModel()` calls it and `ctStanModel` is an alias for
# it, so it is the one chokepoint every fittable model passes through
# (`ctModel(type='omx')` alone skips it, and that object is not fittable until
# converted). Placed after the field parsing and before the cycle check, while
# `ctspec$param` still holds the text the user wrote: `ctModelStatesAndPARS()`
# rewrites latent and tdpred names into `state[i]` and `tdpreds[rowi, j]` later,
# in ctFit(), by which point the collisions are no longer visible.

# Bare identifiers in a cell's text: names that are not a function call and not
# an indexed reference. Those two are the only other things a name can be in a
# cell, so what is left has to resolve to something the model declares.
#
# Two details are load-bearing. The lookbehind stops a match starting inside a
# number or another identifier, so the `e` of `1e-5` and the `ARS` of `xPARS`
# are not names. And an identifier may not begin with `.`, so the `.5` of
# `0.5*PARS[1,1]` is not one either -- ctsem's own names are validated
# alphanumeric (R/ctModel.R), so nothing legitimate is lost.
.ctSpecNameRE <- '(?<![A-Za-z0-9._])[A-Za-z][A-Za-z0-9._]*[[:space:]]*[([]?'

.ctSpecBareNames <- function(txt) {
  if (is.null(txt) || length(txt) != 1L || is.na(txt) || !nzchar(txt)) {
    return(character(0))
  }
  got <- regmatches(txt, gregexpr(.ctSpecNameRE, txt, perl = TRUE))[[1]]
  if (!length(got)) return(character(0))
  got <- got[!grepl('[([]$', got)]
  unique(trimws(got))
}

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
    # In PARS the same text is refused outright. A PARS cell is there to
    # declare a parameter, so a bare name reads as the declaration -- taking
    # it as a state reference instead is not a model anyone wants, and the
    # cell's label then collides with the state name everywhere else too.
    inpars <- identical(as.character(ctspec$matrix[i]), 'PARS')
    if (p %in% latentNames) {
      hint <- paste0(cell, ' is the latent state name "', p, '", so the cell ',
        'takes that state\'s value rather than being a free parameter. Write ',
        'state[', match(p, latentNames), '] if that is what you meant',
        if (inpars) ', or as part of an expression' else
          ' -- it builds the same model without this warning',
        ', or give the parameter a name of its own.')
      if (inpars) stop(hint, call. = FALSE) else warning(hint, call. = FALSE)
      next
    }
    if (p %in% TDpredNames) {
      hint <- paste0(cell, ' is the time dependent predictor name "', p,
        '", so the cell takes that predictor\'s value rather than being a ',
        'free parameter. Write tdpreds[rowi, ', match(p, TDpredNames),
        '] if that is what you meant, or give the parameter a name of its own.')
      if (inpars) stop(hint, call. = FALSE) else warning(hint, call. = FALSE)
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

    # --- an expression naming something undeclared ------------------------
    #
    # A compound cell becomes a *calculation* in the generated program, and it
    # is resolved by name against the latent states, the time dependent
    # predictors and the PARS cells. Nothing else introduces a name: writing a
    # fresh one inside an expression does not declare it. Measured across the
    # matrices, `DRIFT = matrix('-exp(zz)*eta1')` with no PARS cell for `zz`
    # dies at ctFit() with `object 'zz' not found`, and so do the same shapes
    # in LAMBDA, DIFFUSION, MANIFESTMEANS, CINT and TDPREDEFFECT. The one that
    # does not is worse: `DRIFT = matrix('-log1p(exp(zz))*(1+TD1)')` builds and
    # `zz` is simply not a parameter, so the cell is not what it looks like.
    #
    # Declared in a PARS cell, every one of those works. So the rule is that
    # each name has to be declared, and PARS is where a parameter for an
    # expression gets declared.
    #
    # An expression with no names at all is caught separately, since there is
    # nothing to report as undeclared: `'2/3'` becomes a free parameter whose
    # label is the text and whose transform is built from NA, reaching the
    # backend as the unfinished string "99999 + 99999*".
    if (grepl('\\W', gsub('.', '', p, fixed = TRUE))) {
      declared <- c(latentNames, TDpredNames, .ctSpecParsLabels(ctspec),
        'state', 'tdpreds', 'PARS', 'param', 'rowi')
      bare <- .ctSpecBareNames(p)
      undeclared <- setdiff(bare, declared)
      if (length(undeclared)) {
        # Some of these are not names a PARS cell could supply. Saying so is
        # the point -- `dt` reads as "scaled by the interval", and a cell has
        # no interval, so telling the user to declare a parameter called dt
        # would send them the wrong way entirely.
        why <- character(0)
        for (nm in undeclared) {
          if (nm %in% c('dt', 'time')) {
            why <- c(why, paste0('"', nm, '" is the time interval, which a ',
              'cell cannot see -- an interval effect comes from the ',
              'continuous time parameters, not from a cell'))
          } else if (nm %in% TIpredNames) {
            why <- c(why, paste0('"', nm, '" is a time independent predictor, ',
              'which acts through the cell\'s tipreds field rather than by ',
              'being referenced'))
          } else if (nm %in% manifestNames) {
            why <- c(why, paste0('"', nm, '" is a manifest variable, which a ',
              'cell cannot see -- reference the latent state it loads on'))
          }
        }
        fixable <- setdiff(undeclared,
          c('dt', 'time', TIpredNames, manifestNames))
        stop(cell, ' is the expression "', p, '", which uses ',
          paste0('"', undeclared, '"', collapse = ', '),
          if (length(undeclared) > 1) ' -- none of which are declared'
            else ' -- which is not declared',
          ' by this model. Writing a name inside an expression does not ',
          'declare it: a name has to be a latent state, a time dependent ',
          'predictor, or a parameter given a PARS cell. ',
          if (length(why)) paste0('Note that ',
            paste(why, collapse = ', and '), '. ') else '',
          if (length(fixable)) paste0('Add PARS = c(\'', fixable[1],
            '\') and the expression works as written; to transform a single ',
            'parameter and nothing else, the second field is shorter: "',
            fixable[1], '|-exp(param)".') else '',
          call. = FALSE)
      }
      computable <- grepl('[', p, fixed = TRUE) || length(bare) > 0
      if (!computable) {
        stop(cell, ' is the expression "', p, '", which names nothing at all, ',
          'so there is nothing for it to be computed from. For a fixed cell, ',
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
          'expression needs, and refer to them by the names you give them.',
          call. = FALSE)
      }
      if (rc[1] > parsdim[1] || rc[2] > parsdim[2]) {
        stop(cell, ' references ', ref, ', but PARS is ', parsdim[1], 'x',
          parsdim[2], '.', call. = FALSE)
      }
    }
  }

  # --- which cells actually hold a parameter -----------------------------
  #
  # The three checks below aggregate over cells that name the same parameter,
  # so they need to know which cells *hold* one. A simple name in a cell is
  # not always a parameter name:
  #
  #   a latent state or tdpred name  is a reference to that state or predictor,
  #                                  and several matrices may reference the
  #                                  same one -- that is not one parameter
  #                                  appearing twice
  #   a PARS label, outside PARS     is a reference to the PARS cell that
  #                                  declares it; the referencing cell's own
  #                                  transform column is discarded when the
  #                                  reference is rewritten, so it says nothing
  #                                  about the parameter
  #
  # Only the declaring cell holds the parameter. Note that this still leaves
  # two PARS cells sharing one label as co-holders, which they are -- that is
  # the equality constraint inside PARS -- so a disagreement between them is a
  # real one. Exempting every PARS label instead, as a first attempt did, hid
  # exactly that case.
  parslabels <- .ctSpecParsLabels(ctspec)
  simplename <- !is.na(ctspec$param) & is.na(ctspec$value) &
    !grepl('\\W', gsub('.', '', ctspec$param, fixed = TRUE))
  isref <- ctspec$param %in% c(latentNames, TDpredNames) |
    (ctspec$param %in% parslabels & !(ctspec$matrix %in% 'PARS'))
  holds <- simplename & !is.na(isref) & !isref

  # --- one parameter, two different transforms ---------------------------
  #
  # Repeating a name is the equality constraint: the cells hold one parameter,
  # one raw value, one random effect. The transform, though, stays per cell --
  # so two cells naming the same parameter can turn that one raw value into
  # two different numbers, and there is no answer to "what is this parameter's
  # value" or "what does its population sd mean".
  #
  # The reachable ways in are wider than writing two transforms out. Only the
  # first of these states a transform at all:
  #
  #   MANIFESTMEANS = c('mm|exp(param)', 'mm|log(param)')   exp vs log
  #   MANIFESTMEANS = c('mm|exp(param)', 'mm')              exp vs the default
  #   DRIFT = matrix('a', 2, 2)                             the diagonal default
  #                                                         is negative-bounded
  #                                                         and the off-diagonal
  #                                                         is identity
  #   DRIFT[1,1] = 'p' with DIFFUSION[1,1] = 'p'            negative vs positive
  #
  # The third is a plausible thing to write -- "constrain every drift cell to
  # one parameter" -- and it silently produced a negative-bounded number in
  # the diagonal cells and the raw value in the off-diagonal ones.
  #
  # So the rule is about the resolved transform rather than about what was
  # written: cells sharing a name must resolve to the same one. Two drift
  # diagonals sharing a name still work, because their defaults agree, and so
  # does a T0MEANS cell sharing with a MANIFESTMEANS cell.
  for (nm in unique(ctspec$param[holds])) {
    rows <- which(holds & ctspec$param %in% nm)
    if (length(rows) < 2) next
    tf <- ctspec$transform[rows]
    if (length(unique(tf)) < 2) next

    # One cell per distinct transform, at most two. Showing the first two
    # *cells* instead would often show two that agree -- `DRIFT = matrix('a',
    # 2, 2)` has four cells and two transforms -- which reads as a
    # contradiction of the sentence above it. Two is enough to see the
    # disagreement, and a whole constrained matrix would otherwise fill the
    # message.
    shown <- rows[!duplicated(tf)][1:2]
    where <- vapply(shown, function(r) paste0(.ctSpecCellName(ctspec, r),
      ' as ', ctspec$transform[r]), character(1))
    ntf <- length(unique(tf))
    more <- if (length(rows) > 2) paste0('\n  (', length(rows),
      ' cells name it, with ', ntf, ' different transforms)') else ''

    # The worked example uses each shown cell's own transform, with the
    # parameter's name where `param` stood, so it is the transform that cell
    # actually had rather than a plausible-looking stand-in.
    asexpr <- function(r) paste0(.ctSpecCellName(ctspec, r), " = '",
      gsub('\\bparam\\b', nm, ctspec$transform[r]), "'")

    allpars <- all(ctspec$matrix[rows] %in% 'PARS')
    stop('Parameter "', nm, '" is declared by cells that transform it ',
      'differently:\n  ', paste(where, collapse = '\n  '), more,
      '\nThose cells are one parameter -- one raw value and one random ',
      'effect -- so a single transform has to apply to all of them; there is ',
      'otherwise no answer to what the parameter is, or what its population ',
      'standard deviation means. ',
      if (allpars) paste0('Declare "', nm, '" in one PARS cell rather than ',
        'two, and use its name in the cells that need it.') else
      paste0('Either give every cell the same transform explicitly, or ',
        'declare "', nm, '" in one additional PARS cell -- PARS = c(\'', nm,
        '\') -- and then use its name in an expression wherever it is needed, ',
        'e.g. ', asexpr(shown[1]), ' and ', asexpr(shown[2]),
        '. One parameter then serves cells that need different transforms, ',
        'and each transform is stated where it applies.'), call. = FALSE)
  }

  # --- one parameter, two different sdscales -----------------------------
  #
  # Repeating a name is ctsem's equality constraint, so the cells are one
  # parameter and it has one population sd prior. Two cells asking for
  # different scales is resolved by position -- whichever row is reached first
  # -- so the same model with its cells written in the other order gets a
  # different prior, and nothing reports which was used.
  for (nm in unique(ctspec$param[holds])) {
    rows <- which(holds & ctspec$param %in% nm)
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

  # --- one parameter, two different sets of predictor effects ------------
  #
  # A predictor effect displaces the *raw* parameter, and `TIPREDEFFECTsetup`
  # is indexed by parameter rather than by cell, so "age affects this cell but
  # not that one" cannot be represented for cells that are one parameter. It
  # was accepted anyway and resolved by position: measured,
  # `MANIFESTMEANS = c('mm||||age', 'mm||||sex')` and the same two cells in the
  # other order produce different sets of estimated effects.
  effectcols <- if (length(TIpredNames)) paste0(TIpredNames, '_effect') else
    character(0)
  effectcols <- intersect(effectcols, names(ctspec))
  if (length(effectcols)) {
    for (nm in unique(ctspec$param[holds])) {
      rows <- which(holds & ctspec$param %in% nm)
      if (length(rows) < 2) next
      stated <- vapply(rows, function(r) paste(
        .ctTipredEffectSpec(ctspec[r, effectcols]), collapse = ', '),
        character(1))
      if (length(unique(stated)) < 2) next
      shown <- rows[!duplicated(stated)][1:2]
      stop('Parameter "', nm, '" is declared by cells that give it different ',
        'time independent predictor effects:\n  ',
        paste(vapply(shown, function(r) paste0(.ctSpecCellName(ctspec, r),
          ' as ', paste0(sub('_effect$', '', effectcols), '=',
            .ctTipredEffectSpec(ctspec[r, effectcols]), collapse = ', ')),
          character(1)), collapse = '\n  '),
        if (length(rows) > 2) paste0('\n  (', length(rows),
          ' cells name it)') else '',
        '\nAn effect displaces the parameter itself, not one cell of it, so ',
        'the cells that are one parameter cannot differ here -- which of them ',
        'was used depended on the order they were written in. State the same ',
        'effects in every cell that names "', nm, '", or give the cells ',
        'different parameter names.', call. = FALSE)
    }
  }

  # --- one effect name, two different transforms -------------------------
  #
  # Naming an effect constrains the parameters that carry it to one
  # coefficient (`.ctTipredEffectKey()`). That coefficient displaces each
  # parameter's raw value equally, so what it does on the natural scale is the
  # same for all of them only when their transforms are the same. Where they
  # differ there is one number and two meanings, and nothing a summary or a
  # plot could report as "the effect".
  if (length(effectcols)) {
    for (ci in seq_along(effectcols)) {
      spec <- .ctTipredEffectSpec(ctspec[[effectcols[ci]]])
      label <- .ctTipredEffectLabel(spec)
      for (nm in unique(label[!is.na(label) & holds])) {
        rows <- which(holds & !is.na(label) & label %in% nm)
        pars <- unique(ctspec$param[rows])
        if (length(pars) < 2) next
        tf <- ctspec$transform[rows][!duplicated(ctspec$param[rows])]
        if (length(unique(tf)) < 2) next
        shown <- rows[!duplicated(ctspec$transform[rows])][1:2]
        stop('The ', sub('_effect$', '', effectcols[ci]), ' effect named "',
          nm, '" is shared by parameters that are transformed differently:\n  ',
          paste(vapply(shown, function(r) paste0('"', ctspec$param[r],
            '" in ', .ctSpecCellName(ctspec, r), ', as ',
            ctspec$transform[r]), character(1)), collapse = '\n  '),
          '\nOne coefficient displaces every parameter that carries it by the ',
          'same amount on the raw scale, so it only means one thing where the ',
          'transforms agree -- here it would mean something different for each, ',
          'and there is no single effect for a summary or a plot to report. ',
          'Share the name only between parameters with the same transform, or ',
          'let each have its own effect by writing a bare predictor name.',
          call. = FALSE)
      }
    }
  }

  # --- sdscale where nothing varies --------------------------------------
  #
  # sdscale multiplies the prior on a population standard deviation, so it
  # means nothing for a parameter that has no population distribution. It was
  # accepted and dropped.
  indcols <- grep('^indvarying', names(ctspec), value = TRUE)
  for (i in which(holds)) {
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
