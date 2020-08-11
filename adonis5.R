library(ape)
library(lattice)
library(matrixStats)
library(permutes)
library(lmPerm)
library(permute)
library(vegan)
library(GUniFrac)


adonis5 <-  function (formula, data = NULL, permutations = 999, method = "bray", tree = NULL,
                      strata = NULL, contr.unordered = "contr.sum", contr.ordered = "contr.poly", 
                      parallel = getOption("mc.cores"), ...) 
{
  EPS <- sqrt(.Machine$double.eps)
  TOL <- 1e-07
  Terms <- terms(formula, data = data)
  lhs <- formula[[2]]
  lhs <- eval(lhs, data, parent.frame())
  formula[[2]] <- NULL
  rhs.frame <- model.frame(formula, data, drop.unused.levels = TRUE)
  op.c <- options()$contrasts
  options(contrasts = c(contr.unordered, contr.ordered))
  rhs <- model.matrix(formula, rhs.frame)
  options(contrasts = op.c)
  grps <- attr(rhs, "assign")
  qrhs <- qr(rhs)
  rhs <- rhs[, qrhs$pivot, drop = FALSE]
  rhs <- rhs[, 1:qrhs$rank, drop = FALSE]
  grps <- grps[qrhs$pivot][1:qrhs$rank]
  u.grps <- unique(grps)
  nterms <- length(u.grps) - 1
  if (nterms < 1) 
    stop("right-hand-side of formula has no usable terms")
  H.s <- lapply(2:length(u.grps), function(j) {
    Xj <- rhs[, grps %in% u.grps[1:j]]
    qrX <- qr(Xj, tol = TOL)
    Q <- qr.Q(qrX)
    tcrossprod(Q[, 1:qrX$rank])
  })
  if (inherits(lhs, "dist")) {
    if (any(lhs < -TOL)) 
      stop("dissimilarities must be non-negative")
    dmat <- as.matrix(lhs^2)
  }
  else if ((is.matrix(lhs) || is.data.frame(lhs)) && isSymmetric(unname(as.matrix(lhs)))) {
    dmat <- as.matrix(lhs^2)
    lhs <- as.dist(lhs)
  }
  else {
    dist.lhs <- as.matrix(vegdist2(lhs, method = method,tree = tree, ...))
    dmat <- dist.lhs^2
  }
  n <- nrow(dmat)
  G <- -sweep(dmat, 1, rowMeans(dmat))/2
  SS.Exp.comb <- sapply(H.s, function(hat) sum(G * t(hat)))
  SS.Exp.each <- c(SS.Exp.comb - c(0, SS.Exp.comb[-nterms]))
  H.snterm <- H.s[[nterms]]
  tIH.snterm <- t(diag(n) - H.snterm)
  if (length(H.s) > 1) 
    for (i in length(H.s):2) H.s[[i]] <- H.s[[i]] - H.s[[i - 
        1]]
  SS.Res <- sum(G * tIH.snterm)
  df.Exp <- sapply(u.grps[-1], function(i) sum(grps == i))
  df.Res <- n - qrhs$rank
  if (inherits(lhs, "dist")) {
    beta.sites <- qr.coef(qrhs, as.matrix(lhs))
    beta.spp <- NULL
  }
  else {
    beta.sites <- qr.coef(qrhs, dist.lhs)
    beta.spp <- qr.coef(qrhs, as.matrix(lhs))
  }
  colnames(beta.spp) <- colnames(lhs)
  colnames(beta.sites) <- rownames(lhs)
  F.Mod <- (SS.Exp.each/df.Exp)/(SS.Res/df.Res)
  f.test <- function(tH, G, df.Exp, df.Res, tIH.snterm) {
    (sum(G * tH)/df.Exp)/(sum(G * tIH.snterm)/df.Res)
  }
  SS.perms <- function(H, G, I) {
    c(SS.Exp.p = sum(G * t(H)), S.Res.p = sum(G * t(I - H)))
  }
  p <- vegan:::getPermuteMatrix(permutations, n, strata = strata)
  permutations <- nrow(p)
  if (permutations) {
    tH.s <- lapply(H.s, t)
    if (is.null(parallel)) 
      parallel <- 1
    hasClus <- inherits(parallel, "cluster")
    isParal <- hasClus || parallel > 1
    isMulticore <- .Platform$OS.type == "unix" && !hasClus
    if (isParal && !isMulticore && !hasClus) {
      parallel <- makeCluster(parallel)
    }
    if (isParal) {
      if (isMulticore) {
        f.perms <- sapply(1:nterms, function(i) unlist(mclapply(1:permutations, 
                                                                function(j) f.test(tH.s[[i]], G[p[j, ], p[j, 
                                                                                                          ]], df.Exp[i], df.Res, tIH.snterm), mc.cores = parallel)))
      }
      else {
        f.perms <- sapply(1:nterms, function(i) parSapply(parallel, 
                                                          1:permutations, function(j) f.test(tH.s[[i]], 
                                                                                             G[p[j, ], p[j, ]], df.Exp[i], df.Res, tIH.snterm)))
      }
    }
    else {
      f.perms <- sapply(1:nterms, function(i) sapply(1:permutations, 
                                                     function(j) f.test(tH.s[[i]], G[p[j, ], p[j, 
                                                                                               ]], df.Exp[i], df.Res, tIH.snterm)))
    }
    if (isParal && !isMulticore && !hasClus) 
      stopCluster(parallel)
    P <- (rowSums(t(f.perms) >= F.Mod - EPS) + 1)/(permutations + 
                                                     1)
  }
  else {
    f.perms <- P <- rep(NA, nterms)
  }
  SumsOfSqs = c(SS.Exp.each, SS.Res, sum(SS.Exp.each) + SS.Res)
  tab <- data.frame(Df = c(df.Exp, df.Res, n - 1), SumsOfSqs = SumsOfSqs, 
                    MeanSqs = c(SS.Exp.each/df.Exp, SS.Res/df.Res, NA), F.Model = c(F.Mod, 
                                                                                    NA, NA), R2 = SumsOfSqs/SumsOfSqs[length(SumsOfSqs)], 
                    P = c(P, NA, NA))
  rownames(tab) <- c(attr(attr(rhs.frame, "terms"), "term.labels")[u.grps], 
                     "Residuals", "Total")
  colnames(tab)[ncol(tab)] <- "Pr(>F)"
  attr(tab, "heading") <- c(vegan:::howHead(attr(p, "control")), "Terms added sequentially (first to last)\n")
  class(tab) <- c("anova", class(tab))
  out <- list(aov.tab = tab, call = match.call(), coefficients = beta.spp, 
              coef.sites = beta.sites, f.perms = f.perms, model.matrix = rhs, 
              terms = Terms)
  class(out) <- "adonis"
  out
}

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

