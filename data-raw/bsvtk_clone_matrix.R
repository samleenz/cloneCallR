## code to prepare `bsvtk_clone_matrix` dataset goes here
## This creates the canonical BSVTK LeGOteknik clone matrix
bsvtk_clone_matrix <- createCodeMatrix(c("eBFP2", "Sapphire", "Venus", "Tomato", "Katushka"))

usethis::use_data(bsvtk_clone_matrix, overwrite = TRUE)
