capture_warnings <- function(expr) {
  warnings <- character()
  value <- withCallingHandlers(
    expr,
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, warnings = warnings)
}

expect_captured_warning <- function(captured, pattern) {
  expect_true(
    any(grepl(pattern, captured$warnings)),
    info = paste0("Expected warning matching: ", pattern)
  )
}
