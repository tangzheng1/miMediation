# -------------------------------------------------------------------------
# Helper Functions for CAMRA
# -------------------------------------------------------------------------

fastCCLasso <- function(xx, isCnt = FALSE, pseudo = 0.5, k_cv = 3, 
                        lam_min_ratio = 1E-4, k_max = 20, n_boot=100, aa=NULL, bb=NULL) {
  n = nrow(xx)
  p = ncol(xx)
  if(isCnt) {
    xx = xx + pseudo
    xx = xx / rowSums(xx)
  }
  xx2 = log(xx) - rowMeans(log(xx))
  vxx2 = stats::var(xx2)
  
  if(is.null(aa)){
    aa = rep(1, p)
  }
  if(is.null(bb)){
    bb = 1 / diag(vxx2)
  }		
  
  #-Golden section method for the selection of lambda (log10 scale)
  xx = vxx2 * (aa * rep(bb, each = p) + bb * rep(aa, each = p))/2
  diag(xx) = 0
  lam_max = max(abs(xx))
  lam_int2 = log10(lam_max * c(lam_min_ratio, 1))
  a1 = lam_int2[1] 
  b1 = lam_int2[2]
  
  #-Store lambda and corresponding cross validation's loss
  lams = NULL 
  fvals = NULL
  #-Two trial points in first 
  a2 = a1 + 0.382 * (b1 - a1) 
  b2 = a1 + 0.618 * (b1 - a1)
  fb2 = cvfastCCLasso(lambda = 10^b2, k_cv = k_cv, xx2 = xx2, 
                      aa = aa, bb = bb)
  lams = c(lams, b2) 
  fvals = c(fvals, fb2)
  fa2 = cvfastCCLasso(lambda = 10^a2, k_cv = k_cv, xx2 = xx2, 
                      aa = aa, bb = bb)
  lams = c(lams, a2) 
  fvals = c(fvals, fa2)
  #Error tolerance for convergence
  err_lam2 = 1e-1 * max(1, lam_int2)
  err_fval = 1e-4
  err = b1 - a1
  k = 0
  
  while(err > err_lam2 && k < k_max) {
    fval_max = max(fa2, fb2)
    if(fa2 > fb2) {
      a1 = a2
      a2 = b2
      fa2 = fb2
      b2 = a1 + 0.618 * (b1 - a1)
      fb2 = cvfastCCLasso(lambda = 10^b2, k_cv = k_cv, xx2 = xx2, 
                          aa = aa, bb = bb)
      lams = c(lams, b2) 
      fvals = c(fvals, fb2)
    } else {
      b1 = b2
      b2 = a2
      fb2 = fa2
      a2 = a1 + 0.382 * (b1 - a1)
      fa2 = cvfastCCLasso(lambda = 10^a2, k_cv = k_cv, xx2 = xx2, 
                          aa = aa, bb = bb)
      lams = c(lams, a2)
      fvals = c(fvals, fa2)
    }
    fval_min = min(fa2, fb2)
    k = k + 1
    err = b1 - a1
    if(abs(fval_max - fval_min) / (1 + fval_min) <= err_fval) {
      break
    }
  }
  info_cv = list(lams = lams, fvals = fvals, k = k + 2, 
                 lam_int = 10^c(a1, b1))
  lambda = 10^((a2 + b2)/2)
  fit_res = fastcclasso_sub(lambda = lambda, SS2 = vxx2, aa = aa, bb = bb)
  
  sigma_mod <- boot_fastCCLasso(xx2=xx2,sigma_hat= fit_res$sigma,
                                lambda=lambda,aa=aa,bb=bb,
                                n_boot = n_boot, max_iter=200,
                                stop_eps=1e-6)
  
  return(list(rho=sigma_mod$cor_w,
              cov_diag=sigma_mod$var_w,
              lambda_best=lambda,
              info_cv=info_cv,
              p_vals=sigma_mod$p_vals))
}

cvfastCCLasso <- function(lambda, k_cv, xx2, aa, bb) {
  n = nrow(xx2)
  p = ncol(xx2)
  n_b = floor(n / k_cv)
  cv.loss = 0
  for(k in 1:k_cv) {
    ite = (n_b * (k - 1) + 1):(n_b * k)
    vxx2te = stats::var(xx2[ite, ])
    vxx2tr = stats::var(xx2[-ite, ])
    out = fastcclasso_sub(lambda = lambda, SS2 = vxx2tr, aa = aa, bb = bb)
    mm = out$sigma - out$ww - rep(out$ww, each = p) - vxx2te
    cv.loss = cv.loss + mean(mm^2 * aa * rep(bb, each = p))
  }
  return(cv.loss)
}

fastcclasso_sub <- function(lambda, SS2, aa, bb, k_max = 200, x_tol = 1E-4) {
  p = ncol(SS2)
  cc = 1 / (aa * sum(bb) + bb * sum(aa))
  aa2 = aa * cc
  bb2 = bb * cc
  cab1 = 1 + sum(aa * bb2)
  caa = sum(aa * aa2)
  cbb = sum(bb * bb2)
  aabb = aa * rep(bb, each = p) + bb * rep(aa, each = p)
  lambda2 = 2 * lambda / aabb
  ss2 = rowSums(SS2 * aabb)
  sigma = SS2
  ww = colMeans(sigma) - mean(sigma)/2
  k = 0
  err = 1
  while(err > x_tol && k < k_max) {
    xx = rowSums(sigma * aabb) - ss2
    ax1 = sum(aa2 * xx)
    bx1 = sum(bb2 * xx)
    ww2 = xx * cc + (aa2 * (cbb * ax1 - cab1 * bx1) + bb2 * (caa * bx1 -
                                                               cab1 * ax1)) / (cab1^2 - caa * cbb)
    sigma2 = SS2 + ww2 + rep(ww2, each = p)
    oo = diag(sigma2)
    sigma2 = (sigma2 > lambda2) * (sigma2 - lambda2) + 
      (sigma2 < - lambda2) * (sigma2 + lambda2)
    diag(sigma2) = oo
    err = max(abs(sigma2 - sigma)/(abs(sigma) + 1)) 
    k = k + 1
    sigma = sigma2
  }
  return(list(sigma = sigma, ww = ww2))
}

boot_fastCCLasso <- function(xx2,sigma_hat,lambda,aa,bb, 
                             n_boot = 100,
                             max_iter=200, 
                             stop_eps=1e-6) {
  n <- nrow(xx2)
  p <- ncol(xx2)
  
  cors_boot <- matrix(0, nrow = p * (p - 1)/2, ncol = n_boot + 1)
  vars_boot <- matrix(0, nrow = p, ncol = n_boot + 1)
  cors_mat <- matrix(0, p, p)
  ind_low <- lower.tri(cors_mat)
  
  sam_boot <- matrix(sample(1:n, size = n * n_boot, replace = T), 
                     ncol = n_boot)
  for(k in 1:n_boot) {
    ind_samp <- sam_boot[, k]
    S_samp <- stats::var(xx2[ind_samp,])
    cov_est<- fastcclasso_sub(lambda, SS2=S_samp, 
                              aa=aa, bb=bb,
                              k_max = 200, x_tol = stop_eps)
    vars_boot[, k] <- diag(cov_est$sigma)
    Is <- 1 / sqrt(vars_boot[, k])
    cor_est <- Is * cov_est$sigma * rep(Is, each = p)
    cors_boot[, k] <- cor_est[ind_low]
  }
  
  vars_boot[, n_boot + 1] <- diag(sigma_hat)
  Is <- 1 / sqrt(vars_boot[, n_boot + 1])
  cor_est <- Is * sigma_hat * rep(Is, each = p)  
  cors_boot[, n_boot + 1] <- cor_est[ind_low] 
  
  vars2 <- rowMeans(vars_boot)
  cors2mod <- rowMeans(cors_boot)
  cors2_mat <- diag(p)
  cors2_mat[ind_low] <- cors2mod
  cors2_mat <- t(cors2_mat)
  cors2_mat[ind_low] <- cors2mod
  p_vals <- pt(cors2mod * sqrt((n - 2) / (1 - cors2mod^2)), df = n - 2)
  p_vals <- ifelse(p_vals <= 0.5, p_vals, 1 - p_vals)
  pval_mat <- diag(p)
  pval_mat[ind_low] <- p_vals
  pval_mat <- t(pval_mat)
  pval_mat[ind_low] <- p_vals
  
  return(list(var_w = vars2, cor_w = cors2_mat, p_vals = pval_mat))   
}

Callallmethods <- function(method,xMat,cv_k, lambda_min_ratio=1e-4,
                           Edge_eps=1e-4){
  p <- dim(xMat)[2]
  S <- stats::var(log(xMat) - rowMeans(log(xMat)))
  lambda_max <- max(max(S - diag(p)), -min(S - diag(p)))
  lambda_min <- lambda_min_ratio * lambda_max
  lambda_int <- c(lambda_min,lambda_max)
  
  if(method=="fastCCLasso"){
    begin_time <- proc.time()
    result <- fastCCLasso(xx=xMat,lam_min_ratio = lambda_min_ratio, 
                          k_cv = cv_k, k_max = 20)
    end_time <- proc.time()
    result_cor <- result$rho
  }else if(method=="SparCC"){
    begin_time <- proc.time()
    result <- compute_corr_mod(fracs=xMat, iter=10, th=0.1)
    end_time <- proc.time()
    result_cor <- result$Cor.mat
  }else if(method=="CCLasso"){    
    begin_time <- proc.time()
    result <- cclasso(xMat, counts = FALSE, pseudo = 0.5, k_cv = cv_k, 
                      lam_int = lambda_int, k_max = 20)
    end_time <- proc.time()	
    result_cor <- result$cor_w
  }else if(method=="COAT"){
    begin_time <- proc.time()
    result <- coat(xMat, nFoler = cv_k, soft = 1)
    end_time <- proc.time()
    result_cor <- result$corr
  }else{
    return(message("This method is not exist, please check it! "))
  } 
  result_cor[abs(result_cor)< Edge_eps] <- 0 
  
  return(list( runtime=as.numeric((end_time-begin_time)[3]),
               est_lower=result_cor[lower.tri(result_cor)] ,
               cor_est=result_cor))
}

SparCC.count <- function(x, imax = 20, kmax = 10, alpha = 0.1, Vmin = 1e-4) {
  p <- ncol(x)
  n <- nrow(x)
  x <- x + 1
  y <- matrix(0, n, p)
  cov.w <- cor.w <- matrix(0, p, p)
  indLow <- lower.tri(cov.w, diag = T)
  covs <- cors <- matrix(0, p * (p + 1) / 2, imax)
  for(i in 1:imax) {
    y <- t(apply(x, 1, function(x) 
      gtools::rdirichlet(n = 1, alpha = x)))
    cov_cor <- SparCC.frac(x = y, kmax = kmax, alpha = alpha, Vmin = Vmin)
    covs[, i] <- cov_cor$cov.w[indLow]
    cors[, i] <- cov_cor$cor.w[indLow]
  }
  cov.w[indLow] <- apply(covs, 1, median) 
  cor.w[indLow] <- apply(cors, 1, median)
  cov.w <- cov.w + t(cov.w)
  diag(cov.w) <- diag(cov.w) / 2
  cor.w <- cor.w + t(cor.w)
  diag(cor.w) <- 1
  return(list(cov.w = cov.w, cor.w = cor.w))
}

SparCC.frac <- function(x, kmax = 10, alpha = 0.1, Vmin = 1e-4) {
  x <- log(x)
  p <- ncol(x)
  TT <- stats::var(x)
  T0 <- diag(TT) + rep(diag(TT), each = p) - 2 * TT
  rowT0 <- rowSums(T0)
  var.w <- (rowT0 - sum(rowT0) / (2 * p - 2))/(p - 2)
  var.w[var.w < Vmin] <- Vmin
  Is <- sqrt(1/var.w)
  cor.w <- (var.w + rep(var.w, each = p) - T0) * Is * rep(Is, each = p) * 0.5
  cor.w[cor.w <= - 1] <- - 1 
  cor.w[cor.w >= 1] <- 1
  Lmat <- diag(rep(p - 2, p)) + 1 
  rp <- NULL
  cp <- rep(TRUE, p)
  k <- 0  
  while(k < kmax && sum(cp) > 3) {
    T02 <- T0
    curr_cor.w <- cor.w
    diag(curr_cor.w) <- 0
    if(!is.null(rp)) {
      curr_cor.w[rp] <- 0
    }
    n_rp <- which.max(abs(curr_cor.w))
    if(abs(curr_cor.w[n_rp]) >= alpha) {
      t_id <- c(arrayInd(n_rp, .dim = c(p, p)))
      Lmat[t_id, t_id] <- Lmat[t_id, t_id] - 1
      n_rp <- c(n_rp, (p + 1) * sum(t_id) - 2 * p - n_rp)
      rp <- c(rp, n_rp)
      T02[rp] <- 0
      cp <- (diag(Lmat) > 0)
      var.w[cp] <- solve(Lmat[cp, cp], rowSums(T02[cp, cp]))
      var.w[var.w <= Vmin] <- Vmin
      Is <- sqrt(1/var.w)
      cor.w <- (var.w + rep(var.w, each = p) - T0) * Is * rep(Is, each = p) * 0.5
      cor.w[cor.w <= - 1] <- - 1
      cor.w[cor.w >= 1] <- 1
    }
    else {
      break
    }
    k <- k + 1
  }
  Is <- sqrt(var.w)
  cov.w <- cor.w * Is * rep(Is, each = p)
  return(list(cov.w = cov.w, cor.w = cor.w))
}

recover_l_PALM <- function(count_m, treat_cov, cov_ad = NULL,
                           prev.filter = 0, eps_p = 1e-10) {
  count_m <- as.matrix(count_m)
  storage.mode(count_m) <- "numeric"
  n <- nrow(count_m)
  p <- ncol(count_m)
  
  orig_taxa <- colnames(count_m)
  orig_taxa <- paste0("O", seq_len(p))
  colnames(count_m) <- orig_taxa
  
  stopifnot(length(treat_cov) == n)
  treat_cov <- matrix(as.numeric(treat_cov), ncol = 1)
  colnames(treat_cov) <- "treat"
  
  rn <- paste0("T", seq_len(n))
  rownames(count_m) <- rn
  rownames(treat_cov) <- rn
  
  if (!is.null(cov_ad)) {
    cov_ad <- data.frame(cov_ad)
    stopifnot(nrow(cov_ad) == n)
    rownames(cov_ad) <- rn
    colnames(cov_ad) <- paste0("Cov", seq_len(ncol(cov_ad)))
  }
  
  result1 <- PALM::palm(
    rel.abd            = count_m,
    covariate.interest = treat_cov,
    covariate.adjust   = cov_ad,
    prev.filter        = prev.filter
  )
  
  r1 <- result1$treat
  
  p_full    <- rep(1, p)  
  z_full    <- rep(0, p)  
  beta_full <- rep(0, p)  
  
  names(p_full) <- names(z_full) <- names(beta_full) <- orig_taxa
  
  feat <- as.character(r1$feature)
  if (length(feat) > 0) {
    idx <- match(feat, orig_taxa)
    ok  <- which(!is.na(idx))
    
    p_kept    <- as.numeric(r1$pval)
    beta_kept <- as.numeric(r1$coef)
    
    p_full[idx[ok]]    <- p_kept[ok]
    beta_full[idx[ok]] <- beta_kept[ok]
    
    p_adj <- pmax(p_kept[ok], eps_p)
    z_full[idx[ok]] <- stats::qnorm(1 - p_adj / 2) * sign(beta_kept[ok])
  }
  
  return(list(
    p      = p_full,
    z      = z_full,
    beta_l = beta_full,
    feature_kept = feat
  ))
}

recover_r <- function(count_m, treat_cov, y, sudo = 0.5, cov_ad = NULL,
                      CClasso = FALSE, cov_true = NULL) {
  logdata  <- log((count_m + sudo) / rowSums(count_m + sudo))
  por_data <- (count_m + sudo) / rowSums(count_m + sudo)
  
  CClasso_core <- function(count_m) {
    res_cov <- fastCCLasso(count_m, isCnt = TRUE)
    diag(sqrt(res_cov$cov_diag)) %*% res_cov$rho %*% diag(sqrt(res_cov$cov_diag))
  }
  
  if (CClasso) {
    est_cov <- CClasso_core(count_m)
  } else {
    res_l  <- SparCC.count(count_m)
    est_cov <- res_l$cov.w
  }
  
  if (!is.null(cov_true)) {
    est_cov <- cov_true
  }
  
  p <- ncol(count_m)
  n <- nrow(count_m)
  
  ilr_basis <- compositions::ilrBase(por_data)
  lasso_data_ilr <- as.matrix(compositions::ilr(por_data))
  R2 <- ilr_basis
  Z_ilr <- lasso_data_ilr %*% solve(t(R2) %*% est_cov %*% R2) %*% t(R2) %*% est_cov
  
  if (is.null(dim(treat_cov))) {
    treat_df <- data.frame(treat = as.numeric(treat_cov))
  } else {
    treat_df <- as.data.frame(treat_cov)
    if (nrow(treat_df) != n) stop("The number of rows in treat_cov is inconsistent with count_m.")
  }
  
  if (!is.null(cov_ad)) {
    cov_df <- as.data.frame(cov_ad)
    if (nrow(cov_df) != n) stop("The number of rows in cov_ad is inconsistent with count_m.")
    Zdf <- cbind(treat_df, cov_df)
  } else {
    Zdf <- treat_df
  }
  
  Z0 <- model.matrix(~ ., data = Zdf)
  qrZ <- qr(Z0)
  
  y_hat   <- qr.fitted(qrZ, y)
  y_tilde <- as.numeric(y - y_hat)
  
  X_hat   <- qr.fitted(qrZ, Z_ilr)
  X_tilde <- Z_ilr - X_hat
  
  outRidge <- hdi::ridge.proj(x = X_tilde, y = y_tilde)
  
  all_p  <- outRidge$pval
  beta_r <- outRidge$bhat
  z      <- stats::qnorm(1 - all_p / 2) * sign(beta_r)
  
  beta_r <- as.vector(beta_r)
  z <- as.vector(z)
  return(list(
    p      = all_p,     
    z      = z,
    beta_r = beta_r,    
    y_tilde = y_tilde,
    X_tilde = X_tilde,
    Z0      = Z0,
    X_doubel = Z_ilr
  ))
}

pre_filter_fun <- function(count_m, treat_cov, y,
                           const = 2,
                           seed = 42,
                           sudo = 0.5,
                           cov_ad = NULL,
                           adaptive_L = FALSE) {
  set.seed(seed)
  count_m <- as.matrix(count_m)
  storage.mode(count_m) <- "numeric"
  n <- nrow(count_m)
  p <- ncol(count_m)
  
  if (length(treat_cov) != n) stop("treat_cov length != nrow(count_m)")
  y <- as.numeric(y)
  if (length(y) != n) stop("y length != nrow(count_m)")
  
  if (!is.null(cov_ad)) {
    cov_ad <- as.matrix(cov_ad)
    storage.mode(cov_ad) <- "numeric"
    if (nrow(cov_ad) != n) stop("nrow(cov_ad) != nrow(count_m)")
  }
  
  logdata <- log((count_m + sudo) / rowSums(count_m + sudo))
  logdata[logdata < (-10)] <- (-10)
  
  res_l  <- SparCC.count(count_m)
  est_cov <- res_l$cov.w
  por_data <- (count_m + sudo) / rowSums(count_m + sudo)
  
  ilr_basis <- compositions::ilrBase(por_data)
  R2 <- ilr_basis
  
  Z_ilr <- (logdata %*% R2) %*%
    solve(t(R2) %*% est_cov %*% R2) %*%
    t(R2) %*% est_cov
  
  treat_vec <- as.numeric(treat_cov)
  Z0 <- cbind(treat = treat_vec)
  if (!is.null(cov_ad)) {
    Z0 <- cbind(Z0, cov_ad)
    colnames(Z0) <- make.names(colnames(Z0), unique = TRUE)
  }
  
  X <- cbind(Z_ilr, Z0)
  
  pZ <- ncol(Z_ilr)   
  p0 <- ncol(Z0)      
  pf <- c(rep(1, pZ), rep(0, p0))
  
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("Package 'glmnet' is required for pre_filter_fun().")
  }
  
  if (isTRUE(adaptive_L)) {
    cvfit <- glmnet::cv.glmnet(
      x = X, y = y,
      alpha = 1,
      penalty.factor = pf,
      nfolds = 5,
      type.measure = "mse",
      standardize = TRUE
    )
    b <- as.matrix(stats::coef(cvfit, s = "lambda.min"))  
    beta_all <- as.numeric(b)[-1]  
    beta_Z   <- beta_all[1:pZ]     
    selection_set <- which(beta_Z != 0)
    if (length(selection_set) == 0) {
      ord <- order(abs(beta_Z), decreasing = TRUE)
      selection_set <- ord[1]
    }
    return(sort(unique(as.integer(selection_set))))
  }
  
  K_raw <- floor(const * n / log(max(n, 3)))
  K     <- max(1L, min(pZ, K_raw))
  
  fit <- glmnet::glmnet(
    x = X, y = y,
    alpha = 1,
    penalty.factor = pf,
    dfmax = min(K + p0, pZ + p0),  
    nlambda = 500, lambda.min.ratio = 1e-6,
    standardize = TRUE
  )
  
  B  <- as.matrix(fit$beta)  
  dfZ <- colSums(B[1:pZ, , drop = FALSE] != 0)
  
  idx <- which(dfZ >= K)[1]
  if (is.na(idx)) idx <- ncol(B)
  
  beta_Z_full <- as.numeric(B[1:pZ, idx])
  ord  <- order(abs(beta_Z_full), decreasing = TRUE)
  keep <- ord[seq_len(min(K, length(ord)))]
  
  sort(unique(as.integer(keep)))
}

p_mediation_maxp <- function(p_alpha, p_beta,
                             pi_alpha0 = NULL, pi_beta0 = NULL,
                             pi_method = c("JC","cp4p"),
                             weight_method = c("maxp","product","indenp")) {
  stopifnot(length(p_alpha) == length(p_beta))
  pi_method     <- match.arg(pi_method)
  weight_method <- match.arg(weight_method)
  
  mix_weights_product <- function(pi_alpha0, pi_beta0) {
    eps <- 1e-8
    pa <- min(max(pi_alpha0, eps), 1 - 1e-6)
    pb <- min(max(pi_beta0,  eps), 1 - 1e-6)
    pi0 <- 1 - (1 - pa) * (1 - pb)        
    pi0 <- max(pi0, 1e-6)                 
    w00 <- (pa * pb) / pi0
    w10 <- ((1 - pa) * pb) / pi0
    w01 <- (pa * (1 - pb)) / pi0
    c(w00 = w00, w10 = w10, w01 = w01, pi0 = pi0)
  }
  
  from_miMed_estpi0 <- function(z) {
    xi = c(0:100)/100
    tmax = sqrt(log(length(z)))
    tt = seq(0, tmax, 0.05)
    epsest = NULL
    for (j in 1:length(tt)) {
      t = tt[j]
      f = t * xi
      f = exp(f^2/2)
      w = (1 - abs(xi))
      co = 0 * xi
      for (i in 1:101) {
        co[i] = mean(cos(t * xi[i] * z))
      }
      epshat = sum(w * f * co)/sum(w)
      epsest = c(epsest, epshat)
    }
    tmp = min(epsest)
    if (tmp > 1) 
      tmp = 1
    return(tmp)
  }
  
  mix_weights_maxp <- function(p_alpha, p_beta,
                               pi_alpha0 = NULL, pi_beta0 = NULL,
                               pi_method = c("cp4p","JC")) {
    pi_method <- match.arg(pi_method)
    if (is.null(pi_alpha0) || is.null(pi_beta0)) {
      if (pi_method == "cp4p") {
        pa <- cp4p::estim.pi0(p_alpha); pb <- cp4p::estim.pi0(p_beta)
        grab <- function(x) if (!is.null(x$pi0)) as.numeric(x$pi0) else mean(unlist(x), na.rm = TRUE)
        pi_alpha0 <- grab(pa); pi_beta0 <- grab(pb)
      } else {
        z1 <- stats::qnorm(1 - p_alpha); z2 <- stats::qnorm(1 - p_beta)
        pi_alpha0 <- from_miMed_estpi0(z1)
        pi_beta0  <- from_miMed_estpi0(z2)
      }
    }
    
    p_max <- pmax(p_alpha, p_beta)
    if (pi_method == "cp4p") {
      pi0_hat <- {
        obj <- cp4p::estim.pi0(p_max)
        if (!is.null(obj$pi0)) as.numeric(obj$pi0) else mean(unlist(obj))
      }
    } else {
      zmax   <- stats::qnorm(1 - p_max)
      pi0_hat <- from_miMed_estpi0(zmax)
    }
    
    clip01 <- function(x) min(max(x, 1e-6), 1 - 1e-6)
    pi_alpha0 <- clip01(pi_alpha0); pi_beta0 <- clip01(pi_beta0); pi0_hat <- clip01(pi0_hat)
    
    w00 <- (pi_alpha0 + pi_beta0 - pi0_hat) / pi0_hat
    w10 <- (pi0_hat - pi_alpha0) / pi0_hat
    w01 <- (pi0_hat - pi_beta0)  / pi0_hat
    
    w <- pmax(c(w00, w10, w01), 0); w <- w / sum(w)
    names(w) <- c("w00","w10","w01")
    w
  }
  
  p_alpha_pi0 <- p_alpha
  p_beta_pi0  <- p_beta
  
  bad_a <- !is.finite(p_alpha_pi0)
  bad_b <- !is.finite(p_beta_pi0)
  p_alpha_pi0[bad_a] <- 1
  p_beta_pi0 [bad_b] <- 1
  
  p_alpha_pi0 <- pmin(pmax(p_alpha_pi0, 0), 1)
  p_beta_pi0  <- pmin(pmax(p_beta_pi0,  0), 1)
  
  if (is.null(pi_alpha0) || is.null(pi_beta0)) {
    if (pi_method == "cp4p") {
      pa <- cp4p::estim.pi0(p_alpha_pi0)
      pb <- cp4p::estim.pi0(p_beta_pi0)
      grab <- function(x) {
        if (!is.null(x$pi0)) as.numeric(x$pi0) else mean(unlist(x), na.rm = TRUE)
      }
      pi_alpha0 <- grab(pa)
      pi_beta0  <- grab(pb)
    } else {  
      z1 <- stats::qnorm(1 - p_alpha_pi0)
      z2 <- stats::qnorm(1 - p_beta_pi0)
      pi_alpha0 <- from_miMed_estpi0(z1)
      pi_beta0  <- from_miMed_estpi0(z2)
    }
  }
  
  eps <- 1e-8
  pi_alpha0 <- min(max(pi_alpha0, eps), 1 - eps)
  pi_beta0  <- min(max(pi_beta0,  eps), 1 - eps)
  
  keep <- is.finite(p_alpha) & is.finite(p_beta)
  out  <- rep(NA_real_, length(p_alpha))
  if (!any(keep)) return(out)
  
  p_a <- p_alpha[keep]
  p_b <- p_beta[keep]
  p_a <- pmin(pmax(p_a, eps), 1 - eps)
  p_b <- pmin(pmax(p_b, eps), 1 - eps)
  
  if (weight_method == "maxp") {
    w <- mix_weights_maxp(p_a, p_b,
                          pi_alpha0 = pi_alpha0,
                          pi_beta0  = pi_beta0,
                          pi_method = pi_method)
  } else if (weight_method == "product") {
    w_raw <- mix_weights_product(pi_alpha0, pi_beta0) 
    w_vec <- pmax(w_raw[c("w00","w10","w01")], 0)
    w     <- w_vec / sum(w_vec)
    names(w) <- c("w00","w10","w01")
  } else if (weight_method == "indenp") {
    w_vec <- c(
      w00 = pi_alpha0 * pi_beta0,
      w10 = (1 - pi_alpha0) * pi_beta0,
      w01 = pi_alpha0 * (1 - pi_beta0)
    )
    w_vec <- pmax(w_vec, 0)
    w     <- w_vec / sum(w_vec)
  }
  
  w00 <- as.numeric(w["w00"])
  w10 <- as.numeric(w["w10"])
  w01 <- as.numeric(w["w01"])
  
  estimate_F1_grenander <- function(p, pi0, eps = 1e-8) {
    p <- p[is.finite(p)]
    p <- pmin(pmax(p, 0), 1)
    n <- length(p); stopifnot(n > 0)
    pi0 <- min(max(pi0, 1e-6), 1 - 1e-6)
    
    x <- sort(unique(c(0, sort(p), 1)))
    Fn <- stats::ecdf(p); y <- Fn(x)
    dx <- diff(x); keep <- dx > eps
    xL <- x[-length(x)][keep]; xR <- x[-1][keep]
    yL <- y[-length(y)][keep]; yR <- y[-1][keep]; dx <- xR - xL
    
    s <- (yR - yL) / dx                
    s_hat <- -Iso::pava(-s, w = dx)      
    f1_hat <- pmax((s_hat - pi0) / (1 - pi0), 0)
    area <- sum(f1_hat * dx)
    if (area <= 0) return(function(t) rep(0, length(t)))
    f1_hat <- f1_hat / area
    
    x_knots <- c(0, xR)
    F1_cum  <- c(0, cumsum(f1_hat * dx))
    function(t) {
      t <- pmin(pmax(t, 0), 1)
      v <- stats::approx(x_knots, F1_cum, xout = t, method = "linear",
                         ties = "ordered", rule = 2)$y
      pmin(pmax(v, 0), 1)
    }
  }
  
  F1a <- estimate_F1_grenander(p_a, pi_alpha0)
  F1b <- estimate_F1_grenander(p_b, pi_beta0)
  
  t     <- pmax(p_a, p_b)
  p_mix <- w00 * (t^2) + w10 * (t * F1a(t)) + w01 * (t * F1b(t))
  p_mix <- pmin(pmax(p_mix, 0), 1)
  
  out[keep] <- p_mix
  out
}

p_mediation_hdmt_fdr <- function(p_alpha, p_beta ,exact_p =1) {
  stopifnot(length(p_alpha) == length(p_beta))
  n   <- length(p_alpha)
  out <- rep(NA_real_, n)
  
  keep <- is.finite(p_alpha) & is.finite(p_beta)
  if (!any(keep)) return(out)
  
  pa <- pmin(pmax(p_alpha[keep], 0), 1)
  pb <- pmin(pmax(p_beta [keep], 0), 1)
  input <- cbind(pa, pb)
  
  nullprop <- HDMT::null_estimation(input)  
  fdr <- HDMT::fdr_est(nullprop$alpha00, nullprop$alpha01, nullprop$alpha10,
                       nullprop$alpha1, nullprop$alpha2,
                       input_pvalues = input, exact = exact_p )  
  out[keep] <- fdr
  out
}

# -------------------------------------------------------------------------
# Main CAMRA Function
# -------------------------------------------------------------------------

#' CAMRA: Causal Absolute-abundance Mediation from Relative-Abundance data
#'
#' @description 
#' CAMRA is a method for taxon-level microbiome mediator discovery that delivers 
#' well-calibrated error control. Unlike existing approaches that test mediation effects 
#' directly on the relative-abundance (RA) scale and can suffer from compositionality-induced 
#' false discoveries, CAMRA targets absolute-abundance (AA)–scale mediation effects by 
#' reconstructing AA-level exposure–microbiome and microbiome–outcome associations from RA inputs. 
#' CAMRA then combines the two path-specific signals under a composite-null mediation framework 
#' and reports taxon-wise q-values for mediator discovery.
#'
#' @param mediators A named numeric matrix containing microbiome abundance. Each row is a
#' subject and each column is a taxon. Column name contains the taxon name.
#' @param exposure A numeric vector of the exposure.
#' @param outcome A numeric vector of the outcome variable (currently treated as a continuous variable).
#' @param confounders An optional numeric vector or matrix containing confounders that may affect
#' the treatment, mediators and outcome. Each row is a subject and each column
#' is a specific confounder, e.g., age or sex. Default is NULL.
#' @param pseudo A numeric value for the pseudo-count added before log-ratio transformations. Default is 0.5.
#' @param fdr.alpha An optional numeric value for the desired FDR significance level in identifying
#' mediating nodes on the tree. Default is 0.05.
#' @param hdmt.exact An integer flag passed to HDMT::fdr_est controlling how mixture-null quantiles 
#' are computed: 0 = approximation method, 1 = exact method. Default is 0.
#' @param screen A logical value indicating whether to screen high-dimensional mediators before testing. 
#' Default is FALSE.
#' @param const A numeric multiplier used for the number of filtered bacteria during screening 
#' (const * n/log(n)). Default is 2.
#' @param CClasso A logical value indicating whether to use fastCCLasso to compute correlations 
#' instead of SparCC. Default is FALSE.
#' @param seed An integer seed for reproducibility. Default is 42.
#'
#' @details 
#' CAMRA implements taxon-level mediation testing that targets absolute-abundance (AA)–scale 
#' mediation effects using standard sequencing relative-abundance (RA) count inputs. 
#' For each taxon \eqn{k}, CAMRA first estimates the exposure-to-microbiome association by applying 
#' PALM to obtain a path-specific p-value \eqn{p_{\alpha_k}}. It then estimates the 
#' microbiome-to-outcome association conditional on exposure using a PALAR-style regression 
#' with debiased inference to obtain \eqn{p_{\beta_k}}. 
#' 
#' These two path-specific signals are combined using the max-P joint-significance statistic 
#' \eqn{\max\{p_{\alpha_k}, p_{\beta_k}\}} for testing the composite mediation null 
#' \eqn{H_{0k}: \alpha_k \beta_k = 0}. Finally, CAMRA applies HDMT to obtain taxon-level 
#' q-values by explicitly accounting for the mixture structure of the composite null, enabling 
#' well-calibrated multiple-testing adjustment for mediator discovery.
#'
#' @return A list containing the following components:
#' \item{pval.alpha}{A numeric vector of p-values from the exposure-microbiome association test.}
#' \item{pval.beta}{A numeric vector of p-values from the microbiome-outcome association test (conditional on exposure).}
#' \item{qval.med}{A numeric vector of taxon-level mediation test q-values.}
#' \item{sig.mediators}{A vector of selected microbial mediators significant at the `fdr.alpha` level.}
#' \item{index_detected}{The column indices of the significant mediators in the original matrix.}
#' \item{runtime_sec}{The total execution time in seconds.}
#' 
#' @author Qiyu Wang \email{wang_qy@@mail.ustc.edu.cn}, Yunfei Peng \email{peng228@@wisc.edu}
#'
#' @references 
#' Wang Q, Li Y, Peng Y, Tang, ZZ (2026). Error control in microbiome mediator discovery: benchmark and remedy. Submitted.
#' 
#' Wei Z, Hong Q, Chen G, Hartert TV, Rosas-Salazar C, Das SR, Shilts MH, Levin AM, Tang ZZ (2026). Fast and reliable association discovery in large-scale microbiome studies and meta-analyses using PALM. Genome Biology, accepted.
#' 
#' Li Y, Wang Q, Feng Z, Wang X, Tang ZZ (2025). PALAR: Estimation of absolute abundance effects in regression with relative abundance predictors. Journal of the American Statistical Association (JASA), DOI: 10.1080/01621459.2025.2596250.
#' 
#' Dai JY, Stanford JL, LeBlanc M (2020). A multiple-testing procedure for high-dimensional mediation hypotheses. Journal of the American Statistical Association (JASA), DOI: 10.1080/01621459.2020.1765785.
#'
#' @export
#'
#' @examples
#' \donttest{
#' # Load the real dataset
#' data(data.bmi)
#' 
#' # Run CAMRA mediation screening
#' CAMRA_res <- CAMRA(
#'   mediators = data.bmi$mediators, 
#'   treatment = data.bmi$treatment, 
#'   outcome = data.bmi$outcome,
#'   confounders = NULL,
#'   fdr.alpha = 0.05,
#'   seed = 123
#' )
#' 
#' # View significant microbial mediators
#' print(CAMRA_res$sig.mediators)
#' }
CAMRA <- function(mediators,
                  treatment,
                  outcome,
                  confounders= NULL ,
                  pseudo=0.5,
                  fdr.alpha =0.05,
                  hdmt.exact = 0,
                  screen = FALSE,        
                  const =2,              
                  CClasso = FALSE,       
                  seed=42)
{
  set.seed(seed)
  t0 <- proc.time()[["elapsed"]]
  
  select_otu <- c(1:ncol(mediators))
  
  if (isTRUE(screen)) {
    select_otu <- pre_filter_fun(
      count_m = mediators,
      treat_cov    = treatment,
      y            = outcome,
      const        = const ,
      seed         = seed,
      sudo         = pseudo,
      cov_ad       = confounders
    )
  }
  
  res1 <- recover_l_PALM(count_m=mediators,
                         treat_cov = treatment,
                         cov_ad = confounders)
  
  res2 <- recover_r(count_m=mediators,
                    treat_cov = treatment,
                    y =outcome,
                    cov_ad = confounders,
                    CClasso = CClasso,
                    sudo=pseudo)
  
  p1 <- res1$p
  p2 <- res2$p
  z1 <- res1$z
  z2 <- as.vector(res2$z)
  
  p_matrix <- cbind(p1,p2)
  rownames(p_matrix) <- colnames(mediators)
  rawp.perm <- p_mediation_maxp(p1,p2,pi_method="cp4p",weight_method = "product")   
  p_vec <- stats::p.adjust(rawp.perm, method = "BH")
  
  rawp.perm.rm = stats::na.omit(rawp.perm)
  L = length(rawp.perm.rm)
  rawp.perm.rm[rawp.perm.rm< 1e-10] <- 1e-10
  
  p_vec_f <- NULL
  p_vec_all <- p_vec
  
  if(screen)
  {
    p_vec_f <- stats::p.adjust(rawp.perm[select_otu], method = "BH")
    p_vec_all <- p_vec
    p_vec_all[select_otu] <- p_vec_f 
    p_vec_all[-select_otu] <- 1
  }
  
  stopifnot(hdmt.exact %in% c(0, 1))
  
  exact_first  <- as.integer(hdmt.exact)          
  exact_second <- 1L - exact_first
  
  tmp_locfdr <- try(
    p_mediation_hdmt_fdr(
      p_matrix[select_otu, 1],
      p_matrix[select_otu, 2],
      exact_p = exact_first
    ),
    silent = TRUE
  )
  
  if (inherits(tmp_locfdr, "try-error")) {
    tmp_locfdr <- try(
      p_mediation_hdmt_fdr(
        p_matrix[select_otu, 1],
        p_matrix[select_otu, 2],
        exact_p = exact_second
      ),
      silent = TRUE
    )
  }
  
  if (inherits(tmp_locfdr, "try-error")) {
    selected_values <- p_vec_all
    idx_detected    <- which(selected_values <fdr.alpha)  
    
    locfdr <- rep(NA_real_, length(select_otu))
  } else {
    p_vec_all[select_otu] <- tmp_locfdr
    idx_sub <- which(tmp_locfdr<= fdr.alpha)
    idx_detected <- select_otu[idx_sub]
  }
  
  globalp <- min(p_vec_all, na.rm = TRUE)   
  
  runtime_sec <- as.numeric(proc.time()[["elapsed"]] - t0)
  
  return(list(index_detected=idx_detected,
              qval.med =p_vec_all,
              runtime_sec = runtime_sec,
              sig.mediators = colnames(mediators)[idx_detected],
              pval.alpha = p1,
              pval.beta =p2,
              beta_l = res1$beta_l,
              beta_r = res2$beta_r,
              filter_index = select_otu))
}