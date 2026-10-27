# Exact weighted GLS for an intercept + binary condition and blocks of size <=2.
# This computes the same GLS normal equations as limma::gls.series, using the
# known 2x2 inverse for paired blocks instead of dense n-sample factorization.
# A deterministic independent comparison with limma is required before use.
fast_pairblock_gls <- function(E, weights, condition, block, rho, batch = 256L) {
  stopifnot(identical(dim(E), dim(weights)), all(condition %in% 0:1), length(condition) == ncol(E), abs(rho) < 1)
  groups <- split(seq_along(block), as.character(block))
  stopifnot(all(lengths(groups) <= 2L))
  pairs <- groups[lengths(groups) == 2L]
  ix1 <- vapply(pairs, function(z) z[which(condition[z] == 0L)], integer(1))
  ix2 <- vapply(pairs, function(z) z[which(condition[z] == 1L)], integer(1))
  stopifnot(length(ix1) == length(pairs), length(ix2) == length(pairs))
  tumor <- which(condition == 1L)
  coef <- su <- matrix(NA_real_, nrow(E), 2L, dimnames = list(rownames(E), c("(Intercept)", "conditionTumor")))
  sig <- numeric(nrow(E))
  for (st in seq.int(1L, nrow(E), by = batch)) {
    ids <- st:min(st + batch - 1L, nrow(E))
    y <- E[ids, , drop = FALSE]
    w <- weights[ids, , drop = FALSE]
    off <- -rho * sqrt(w[, ix1, drop = FALSE] * w[, ix2, drop = FALSE])/(1-rho^2)
    w[, c(ix1, ix2)] <- w[, c(ix1, ix2), drop = FALSE]/(1-rho^2)
    a <- rowSums(w) + 2 * rowSums(off)
    c <- rowSums(w[, tumor, drop = FALSE])
    b <- c + rowSums(off)
    d <- rowSums(w * y) + rowSums(off * (y[, ix1, drop = FALSE] + y[, ix2, drop = FALSE]))
    e <- rowSums(w[, tumor, drop = FALSE] * y[, tumor, drop = FALSE]) + rowSums(off * y[, ix1, drop = FALSE])
    den <- a*c-b*b
    stopifnot(all(den > 0))
    b0 <- (c*d-b*e)/den
    b1 <- (a*e-b*d)/den
    residual <- y - b0 - outer(b1, condition)
    rss <- rowSums(w*residual^2) + 2*rowSums(off*residual[,ix1,drop=FALSE]*residual[,ix2,drop=FALSE])
    stopifnot(all(rss >= -1e-8))
    coef[ids,] <- cbind(b0,b1)
    su[ids,] <- sqrt(cbind(c/den,a/den))
    sig[ids] <- sqrt(pmax(rss,0)/(ncol(E)-2L))
  }
  design <- cbind("(Intercept)" = 1, conditionTumor = condition)
  ans <- list(coefficients=coef, stdev.unscaled=su, sigma=sig,
              df.residual=rep(ncol(E)-2L,nrow(E)), Amean=rowMeans(E),
              design=design, rank=2L, cov.coefficients=solve(crossprod(design)),
              method="weighted_GLS_exact_2x2_patient_blocks", correlation=rho)
  methods::new("MArrayLM",ans)
}

validate_fast_pairblock_gls <- function() {
  set.seed(20261023)
  condition <- rep(c(0L,1L),c(8,16))
  block <- c(1:8,1:8,9:16)
  y <- matrix(rnorm(40*24),40,24)
  rownames(y) <- paste0("g",1:40)
  w <- matrix(runif(40*24,.1,4),40,24)
  design <- cbind("(Intercept)"=1,conditionTumor=condition)
  rbindlist(lapply(c(-.4,0,.25,.8),function(rho) {
    a<-fast_pairblock_gls(y,w,condition,block,rho)
    b<-limma::lmFit(y,design,weights=w,block=block,correlation=rho)
    err<-c(coefficients=max(abs(a$coefficients-b$coefficients)),
           stdev_unscaled=max(abs(a$stdev.unscaled-b$stdev.unscaled)),sigma=max(abs(a$sigma-b$sigma)))
    stopifnot(all(err<1e-9))
    data.table(rho=rho,quantity=names(err),max_absolute_error=unname(err))
  }))
}
