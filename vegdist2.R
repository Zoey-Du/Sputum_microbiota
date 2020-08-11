library(ape)
library(lattice)
library(matrixStats)
library(permutes)
library(lmPerm)
library(permute)
library(vegan)
library(GUniFrac)

vegdist2 <- function (x, method = "bray", tree = NULL, binary = FALSE, diag = FALSE, upper = FALSE, 
                      na.rm = FALSE, ...) 
{
  if (method != 17 & method != 18) x1 = x
  ZAP <- 1e-15
  if (!is.na(pmatch(method, "euclidian"))) 
    method <- "euclidean"
  METHODS <- c("manhattan", "euclidean", "canberra", "bray", 
               "kulczynski", "gower", "morisita", "horn", "mountford", 
               "jaccard", "raup", "binomial", "chao", "altGower", "cao", 
               "mahalanobis", "dw", "du")
  method <- pmatch(method, METHODS)
  inm <- METHODS[method]
  if (is.na(method)) 
    stop("invalid distance method")
  if (method == -1) 
    stop("ambiguous distance method")
  if (!method %in% c(1, 2, 6, 16) && any(rowSums(x, na.rm = TRUE) == 
                                         0)) 
    warning("you have empty rows: their dissimilarities may be meaningless in method ", 
            dQuote(inm))
  if (!method %in% c(1, 2, 3, 6, 16) && any(x < 0, na.rm = TRUE)) 
    warning("results may be meaningless because data have negative entries in method ", 
            dQuote(inm))
  if (method == 11 && any(colSums(x) == 0)) 
    warning("data have empty species which influence the results in method ", 
            dQuote(inm))
  if (method == 6) 
    x <- decostand(x, "range", 2, na.rm = TRUE, ...)
  if (method == 16) 
    x <- veganMahatrans(scale(x, scale = FALSE))
  if (binary) 
    x <- decostand(x, "pa")
  N <- nrow(x <- as.matrix(x))
  if (method %in% c(7, 13, 15) && !identical(all.equal(as.integer(x), 
                                                       as.vector(x)), TRUE)) 
    warning("results may be meaningless with non-integer data in method ", 
            dQuote(inm))
  if (method != 17 & method != 18) {
    d <- .C("veg_distance", x = as.double(x), nr = N, nc = ncol(x), 
            d = double(N * (N - 1)/2), diag = as.integer(FALSE), 
            method = as.integer(method), NAOK = na.rm, PACKAGE = "vegan")$d
    if (method == 10) 
      d <- 2 * d/(1 + d)
    d[d < ZAP] <- 0
    if (any(is.na(d))) 
      warning("missing values in results")
    attr(d, "Size") <- N
    attr(d, "Labels") <- dimnames(x)[[1]]
    attr(d, "Diag") <- diag
    attr(d, "Upper") <- upper
    attr(d, "method") <- paste(if (binary) 
      "binary ", METHODS[method], sep = "")
    attr(d, "call") <- match.call()
    class(d) <- "dist"
  }
  if (method == 17) {
    unifrac <- GUniFrac(x1, tree)
    unifrac <- unifrac$unifracs
    d <- as.dist(unifrac[, , 'd_1'])		# Weighted UniFrac
  }
  if (method == 18) {
    unifrac <- GUniFrac(x1, tree)
    unifrac <- unifrac$unifracs
    d <- as.dist(unifrac[, , 'd_UW'])		# Unweighted UniFrac
  }
  d
}
