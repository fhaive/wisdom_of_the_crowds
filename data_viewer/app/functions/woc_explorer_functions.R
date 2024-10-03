#' Read All Sheets from an Excel File
#' @description
#' This function reads all sheets from an Excel file and returns a list of data frames or tibbles.
#' @param filename A character string. The path to the Excel file.
#' @param tibble A logical. If `TRUE`, the sheets are returned as tibbles. If `FALSE`, the sheets are returned as data frames. Default is `FALSE`.
#' @return A named list of data frames or tibbles, where each element corresponds to a sheet in the Excel file. The names of the list elements are the sheet names.
#' @details The function uses `readxl::excel_sheets` to get the names of all sheets in the Excel file. It then reads each sheet using `readxl::read_excel`. If `tibble` is set to `FALSE`, each sheet is converted to a data frame using `as.data.frame`.
#'
#' @import readxl
read_excel_allsheets <- function(filename, tibble = FALSE) { 
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X){ 
    y <- readxl::read_excel(filename, sheet = X)
  })
  
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  return(x)
}



#' Create a Hyperlink to Ensembl
#' @description 
#' This function generates an HTML hyperlink to the Ensembl database with a specified display text.
#' @param val1 A character string. The ID to be used in the Ensembl URL.
#' @param val2 A character string. The text to be displayed for the hyperlink.
#' @return A character string containing the HTML code for the hyperlink. The link points to the Ensembl database and opens in a new tab.
#' @details The function uses `sprintf` to format the input values into an HTML anchor tag. The `href` attribute is constructed using `val1`, and the link text is set to `val2`. The link will open in a new browser tab.
#' @examples
#' # Create a hyperlink to Ensembl with a specific ID and display text
#' link <- createLink("ENSG00000139618", "BRCA2")
#' cat(link)
#' @export
createLink <- function(val1, val2) { 
  sprintf('<a href="http://www.ensembl.org/id/%s" target="_blank" class="btn btn-primary">%s</a>', val1, val2)
}



#' Row Sum Ordering of a Matrix
#' @description 
#' This function orders the rows of a given matrix in decreasing order based on their row sums.
#' @param mat A numeric matrix. The matrix to be row sum ordered.
#' @return A numeric matrix with rows ordered by decreasing row sums.
#' @details The function computes the sum of each row in the input matrix, orders the rows based on these sums in decreasing order, and returns the reordered matrix. The column names and row names of the original matrix are preserved.
#' @examples
#' # Create a sample matrix
#' mat <- matrix(1:9, nrow = 3, byrow = TRUE)
#' colnames(mat) <- c("A", "B", "C")
#' rownames(mat) <- c("Row1", "Row2", "Row3")
#' print(mat)           # Original matrix
#' ordered_mat <- rowsum_ordering(mat)
#' print(ordered_mat)   # Ordered matrix
rowsum_ordering = function(mat) { 
  row_sums <- rowSums(mat)
  ordered_indices <- order(row_sums, decreasing = TRUE)
  ordered_mat <- mat[ordered_indices, , drop = FALSE]
  return(ordered_mat)
}
