## code to prepare `palclonepaletteette` dataset goes here
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package 'ggplot2' is required to build the data")
}

clonepalette <- vector(mode = "list", length = 3) |>
  setNames(c("hex", "colour", "fill"))

clonepalette$hex <- setNames(
  colours <- c(
    "#D3D3D3", "#E65C5C", "#E6765C", "#E6915C", "#E6AC5C", "#E6C65C",
    "#E6E15C", "#CFE65C", "#B5E65C", "#9AE65C", "#7FE65C", "#65E65C",
    "#5CE66E", "#5CE688", "#5CE6A3", "#5CE6BE", "#5CE6D8", "#5CD8E6",
    "#5CBEE6", "#5CA3E6", "#5C88E6", "#5C6EE6", "#655CE6", "#7F5CE6",
    "#9A5CE6", "#B55CE6", "#CF5CE6", "#E65CE1", "#E65CC6", "#E65CAC",
    "#E65C91", "#E65C76"
  ),
  c(
    0, 8, 28, 31, 21, 9, 4, 30, 29, 18, 6, 25, 2, 22, 24, 12, 7,
    19, 5, 1, 13, 14, 27, 26, 3, 23, 10, 11, 17, 16, 20, 15
  )
)

clonepalette$colour <- ggplot2::scale_colour_manual(values = clonepalette$hex)
clonepalette$fill <- scale_fill_manual(values = clonepalette$hex)

usethis::use_data(clonepalette, overwrite = TRUE)
