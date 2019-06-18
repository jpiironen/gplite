.onAttach <- function(...) {
  ver <- utils::packageVersion("gplite")
  packageStartupMessage("This is gplite version ", ver)
}