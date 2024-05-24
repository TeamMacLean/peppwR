#' Find Fits for Nested Data
#'
#' @description This groups the data frame by specified columns, nests the data, and applies the fitting functions to each group.
#'
#' @param df A data frame containing the data to be fitted.
#' @param id_col A character string specifying the column name for the identifier. Default is "id".
#' @param group_col A character string specifying the column name for the group. Default is "group".
#' @param value_col A character string specifying the column name for the values to be fitted. Default is "value".
#'
#' @return A data frame with nested data and the fit results for each group.
#'
#' @export
find_fits <- function(df, id_col="id", group_col="group", value_col="value"){
    tidyr::nest(df, .by = tidyr::all_of(c(id_col, group_col )), data = {{value_col}} ) |>
    dplyr::mutate(
    fits = purrr::map(data, do_fits)
    )
}
