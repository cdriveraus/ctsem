# Generate and optionally compile latex equation of subject level ctsem model.

Generate and optionally compile latex equation of subject level ctsem
model.

## Usage

``` r
ctModelLatex(
  x,
  matrixnames = TRUE,
  digits = 3,
  linearise = inherits(x, "ctStanFit"),
  textsize = "normalsize",
  folder = tempdir(),
  filename = paste0("ctsemTex", as.numeric(Sys.time())),
  tex = TRUE,
  equationonly = FALSE,
  compile = TRUE,
  open = TRUE,
  includeNote = TRUE,
  minimal = FALSE,
  splitDynamics = TRUE,
  splitMeasurement = TRUE,
  savepng = FALSE
)
```

## Arguments

- x:

  A model from
  [`ctModel`](https://cdriveraus.github.io/ctsem/reference/ctModel.md),
  or a fit from
  [`ctFit`](https://cdriveraus.github.io/ctsem/reference/ctFit.md) on
  either backend (`ctStanFit` or `ctJuliaFit`). Given a fit, the
  estimates are substituted for the parameter labels.

- matrixnames:

  Logical. If TRUE, includes ctsem matrix names such as DRIFT and
  DIFFUSION under the matrices.

- digits:

  Precision of decimals for numeric values.

- linearise:

  Logical. Show the linearised normal approximation for subject
  parameters and covariate effects, or the raw parameters? Only relevant
  where a subject parameter distribution is shown, and not available for
  a `ctJuliaFit`, which holds that distribution on the raw scale.

- textsize:

  Standard latex text sizes – tiny scriptsize footnotesize small
  normalsize large Large LARGE huge Huge. Useful if output overflows
  page.

- folder:

  Character string specifying folder to save to, defaults to temporary
  directory, use "./" for working directory.

- filename:

  filename, without suffix, to output .tex and .pdf files too.

- tex:

  Save .tex file? Otherwise latex is simply returned within R as a
  string.

- equationonly:

  Logical. If TRUE, output is only the latex relevant to the equation,
  not a compileable document.

- compile:

  Compile to .pdf? (Depends on `tex = TRUE`)

- open:

  Open after compiling? (Depends on `compile = TRUE`)

- includeNote:

  Include text describing matrix transformations and subject notation?
  triangular matrices (which results in a covariance or Cholesky matrix)
  is shown – the latter is a more direct representation of the model,
  while the former is often simpler to convey.

- minimal:

  if TRUE, outputs reduced form version displaying matrix dimensions and
  equation structure only.

- splitDynamics:

  Logical. If TRUE, split latent process dynamics across deterministic
  and random change lines. If FALSE, show the full dynamics equation on
  one line.

- splitMeasurement:

  Logical. If TRUE, split measurement equations across observation and
  observation noise lines. If FALSE, show the full measurement equation
  on one line.

- savepng:

  Logical. If TRUE, renders the equation as a png file instead of a pdf,
  viewing the png in RStudio viewer when available.

## Value

character string of latex code. Side effects include saving a .tex,
.pdf, and displaying the pdf.

## Examples

``` r
ctmodel <- ctModel(type='ct', 
n.latent=2, n.manifest=1, 
manifestNames='sunspots', 
latentNames=c('ss_level', 'ss_velocity'),
LAMBDA=matrix(c( 1, 'ma1' ), nrow=1, ncol=2),
DRIFT=matrix(c(0, 1,   'a21', 'a22'), nrow=2, ncol=2, byrow=TRUE),
MANIFESTMEANS=matrix(c('m1'), nrow=1, ncol=1),
CINT=matrix(c(0, 0), nrow=2, ncol=1),
DIFFUSION=matrix(c(
  0, 0,
  0, "diffusion"), ncol=2, nrow=2, byrow=TRUE))
  
l=ctModelLatex(ctmodel,compile=FALSE, open=FALSE)
cat(l)
```
