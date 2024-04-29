
find_fits <- function(df, id_col="id", group_col="group", value_col="value"){
    tidyr::nest(df, .by = tidyr::all_of(c(id_col, group_col )), data = {{value_col}} ) |>
    dplyr::mutate(
    fits = purrr::map(data, do_fits)
    )
}
