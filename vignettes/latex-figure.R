# Shows ctModelLatex()'s equations as an image in a rendered vignette or
# article. ctModelLatex() compiles them to a pdf with LaTeX, and the page is
# rendered, trimmed and written as a small png where knitr keeps figures.
latexFigure <- function(x, ...) {
  dir <- tempfile('latex')
  dir.create(dir)
  ctModelLatex(x, folder = dir, filename = 'model', open = FALSE, ...)
  page <- pdftools::pdf_render_page(file.path(dir, 'model.pdf'), dpi = 130)
  img <- magick::image_trim(magick::image_read(page))
  file <- knitr::fig_path('.png')
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  magick::image_write(magick::image_quantize(img, 16), file, format = 'png')
  knitr::include_graphics(file)
}
