PLSIM_test_split_ratio_cross2_ts <- function(
    Y, X, Z, 
    split_ratio = c(4, 4),
    K_si = 4,
    C_lambda_theta = 0.5,
    C_lambda_beta = 0.5,
    df_si_orth = 6,
    do_si_orth = TRUE,
    mu_true = NULL,
    beta_true = NULL,
    oracle_si = NULL,
    gap = 5,
    hac_c_lower = 0.25,
    hac_c_upper = 4,
    seed = 1
) {
  
  set.seed(seed)
  
  Z <- as.matrix(Z)
  X <- as.matrix(X)
  Y <- as.numeric(Y)
  
  n <- nrow(Z)
  p <- ncol(X)
  q <- ncol(Z)
  
  eps_num <- 1e-12
  
  stopifnot(length(Y) == nrow(X), nrow(X) == nrow(Z))
  stopifnot(length(split_ratio) == 2)
  stopifnot(all(split_ratio > 0))
  
  if (K_si < 2L) stop("K_si must be at least 2.")
  
  if (length(gap) != 1L || !is.finite(gap) ||
      gap < 0 || gap != floor(gap)) {
    stop("gap must be a non-negative integer.")
  }
  gap <- as.integer(gap)
  
  if (length(hac_c_lower) != 1L || !is.finite(hac_c_lower) || hac_c_lower <= 0 ||
      length(hac_c_upper) != 1L || !is.finite(hac_c_upper) || hac_c_upper <= hac_c_lower) {
    stop("Require 0 < hac_c_lower < hac_c_upper < Inf.")
  }
  
  np <- reticulate::import("numpy")
  
  
  ############################################################
  ## 0. 时间序列数据：保持时间顺序，并在两个 half 之间留 outer gap
  ############################################################
  
  id_all <- seq_len(n)
  cut0 <- floor(n / 2)
  
  if (cut0 - gap <= 20L || n - cut0 - gap <= 20L) {
    stop("gap is too large relative to the sample size.")
  }
  
  if (gap > 0L) {
    # 分界点左右各丢 gap 个，一共丢弃 2*gap 个数据。
    id_part1 <- seq_len(cut0 - gap)
    id_part2 <- seq(from = cut0 + gap + 1L, to = n)
  } else {
    id_part1 <- seq_len(cut0)
    id_part2 <- seq(from = cut0 + 1L, to = n)
  }
  
  stopifnot(length(id_part1) > 20, length(id_part2) > 20)
  
  
  ############################################################
  ## 1. 先在两个 half 内按时间顺序划分 ce / beta / rf，并在不同阶段之间留 gap
  ############################################################
  
  make_half_split <- function(id_half) {
    
    id_half <- sort(id_half)
    n_half <- length(id_half)
    ratio_use <- split_ratio / sum(split_ratio)
    
    # 为与理论中的 separation 条件一致，ce-beta 和 beta-rf 都使用完整 gap。
    gap_inner <- gap
    n_available <- n_half - 2L * gap_inner
    
    if (n_available <= 30L) {
      stop("Too few observations after removing the internal gaps.")
    }
    
    n_block1 <- floor(n_available * ratio_use[1])
    n_rf <- n_available - n_block1
    
    n_ce <- floor(n_block1 / 2)
    n_beta <- n_block1 - n_ce
    
    stopifnot(n_ce > 5, n_beta > 5, n_rf > 10)
    
    ce_start <- 1L
    ce_end <- n_ce
    
    beta_start <- ce_end + gap_inner + 1L
    beta_end <- beta_start + n_beta - 1L
    
    rf_start <- beta_end + gap_inner + 1L
    rf_end <- n_half
    
    list(
      id_ce = id_half[ce_start:ce_end],
      id_beta = id_half[beta_start:beta_end],
      id_rf = id_half[rf_start:rf_end]
    )
  }
  
  split1 <- make_half_split(id_part1)
  split2 <- make_half_split(id_part2)
  
  # single-index 部分两个 half、所有内部 folds 共用同一套 Z 标准化
  # 只使用 rf 之外的数据，避免 held-out rf fold 进入 scaler
  id_si_scale <- c(
    split1$id_ce, split1$id_beta,
    split2$id_ce, split2$id_beta
  )
  
  Z_scale_ref <- Z[id_si_scale, , drop = FALSE]
  Z_mu_si <- colMeans(Z_scale_ref)
  Z_sd_si <- apply(Z_scale_ref, 2, sd)
  
  if (any(!is.finite(Z_mu_si)) || any(!is.finite(Z_sd_si)) || any(Z_sd_si <= 0)) {
    stop("Invalid common Z scaler for single-index fitting.")
  }
  
  Z_scale_std <- sweep(sweep(Z_scale_ref, 2, Z_mu_si, "-"), 2, Z_sd_si, "/")
  a_si <- as.numeric(
    quantile(sqrt(rowSums(Z_scale_std^2)), probs = 0.95, names = FALSE)
  )
  
  
  ############################################################
  ## 2. 单个 half 内部完整走一遍流程
  ############################################################
  
  run_one_half <- function(split_obj, seed) {
    
    id_ce <- split_obj$id_ce
    id_beta <- split_obj$id_beta
    id_rf <- split_obj$id_rf
    
    ############################################################
    ## 1.1 用 id_ce 估计条件期望；id_beta 估计 beta
    ############################################################
    
    if (!is.null(beta_true)) {
      
      beta_use <- as.numeric(beta_true)
      
    } else {
      
      if (!is.null(mu_true)) {
        
        E_YZ_beta <- as.numeric(mu_true$mu_Y)
        E_XZ_beta <- as.numeric(mu_true$mu_X)
        
        Y_tilde_beta <- Y[id_beta] - E_YZ_beta[id_beta]
        
        X_tilde_beta <- sweep(
          X[id_beta, , drop = FALSE], 1, E_XZ_beta[id_beta], "-"
        )
        
      } else {
        
        fit_ce <- do.call(
          fit_conditional_expectation_models,
          list(
            Z_train = np$array(Z[id_ce, , drop = FALSE]),
            Y_train = np$array(Y[id_ce]),
            X_train = np$array(X[id_ce, , drop = FALSE]),
            verbose = FALSE,
            seed = as.integer(seed + 100)
          )
        )
        
        res_beta <- do.call(
          predict_conditional_expectations,
          list(
            fit_obj = fit_ce,
            Z_new = np$array(Z[id_beta, , drop = FALSE])
          )
        )
        
        E_YZ_beta <- as.numeric(res_beta[[1]])
        E_XZ_beta <- as.matrix(res_beta[[2]])
        
        Y_tilde_beta <- Y[id_beta] - E_YZ_beta
        X_tilde_beta <- X[id_beta, , drop = FALSE] - E_XZ_beta
      }
      
      n_beta_obs <- length(Y_tilde_beta)
      d_beta <- ncol(X_tilde_beta)
      
      if (d_beta > n_beta_obs) {
        
        nfolds_beta <- min(5, n_beta_obs)
        
        # 时间序列下用连续 block folds
        foldid_beta <- cut(
          seq_len(n_beta_obs),
          breaks = nfolds_beta,
          labels = FALSE
        )
        
        fit_cv <- glmnet::cv.glmnet(
          x = X_tilde_beta,
          y = Y_tilde_beta,
          alpha = 1,
          intercept = FALSE,
          standardize = TRUE,
          foldid = foldid_beta,
          nfolds = nfolds_beta
        )
        
        beta_use <- as.numeric(coef(fit_cv, s = "lambda.min")[-1])
        
      } else {
        
        fit_ols <- lm.fit(x = X_tilde_beta, y = Y_tilde_beta)
        
        beta_use <- as.numeric(fit_ols$coefficients)
        beta_use[is.na(beta_use)] <- 0
      }
    }
    
    
    ############################################################
    ## 1.2 在 id_rf 内部做 OOF 单指标 residual
    ############################################################
    
    n_rf_actual <- length(id_rf)
    K_si_use <- min(K_si, n_rf_actual)
    
    # 时间序列下使用连续 block folds。
    fold_si <- cut(
      seq_len(n_rf_actual),
      breaks = K_si_use,
      labels = FALSE
    )
    
    r_rf_oof <- rep(NA_real_, n_rf_actual)
    
    N_si_fold <- rep(NA_integer_, K_si_use)   
    a_si_fold <- rep(NA_real_, K_si_use)   
    N_si_half <- NA_integer_                 
    si_fit_list <- vector("list", K_si_use)
    
    for (k in seq_len(K_si_use)) {
      
      id_si_test_local <- which(fold_si == k)
      
      # hv-block / gap cross-fitting：validation block 保留，
      # 只从 SI training sample 中删除 validation block 及其左右 gap 个时间点。
      test_min <- min(id_si_test_local)
      test_max <- max(id_si_test_local)
      id_gap_local <- which(
        seq_len(n_rf_actual) >= test_min - gap &
          seq_len(n_rf_actual) <= test_max + gap
      )
      id_si_train_local <- setdiff(seq_len(n_rf_actual), id_gap_local)
      
      if (length(id_si_train_local) <= 10L) {
        stop(sprintf(
          paste0(
            "Too few SI training observations after applying gap: ",
            "fold = %d, training size = %d, gap = %d."
          ),
          k, length(id_si_train_local), gap
        ))
      }
      
      id_si_test_global <- id_rf[id_si_test_local]
      id_si_train_global <- id_rf[id_si_train_local]
      
      if (!is.null(oracle_si)) {
        
        # 计算 oracle fitted value
        index_si_test <- as.vector(
          Z[id_si_test_global, , drop = FALSE] %*% oracle_si$theta
        )
        
        Rhat_si_test <- oracle_si$g_fun(index_si_test)
        
      } else {
        
        R_si_train <- Y[id_si_train_global] - 
          as.vector(X[id_si_train_global, , drop = FALSE] %*% beta_use)
        
        fit_si_k <- spline_single_index(
          Z_train = Z[id_si_train_global, , drop = FALSE],
          Y = R_si_train,
          Z_mu = Z_mu_si,
          Z_sd = Z_sd_si,
          a_si = a_si
        )
        
        si_fit_list[[k]] <- fit_si_k
        
        N_si_fold[k] <- as.integer(fit_si_k$N)   # same
        if (is.na(N_si_half)) {
          N_si_half <- as.integer(fit_si_k$N)
        }
        
        a_si_fold[k] <- as.numeric(fit_si_k$a)    
        
        Rhat_si_test <- predict_si(
          theta = fit_si_k$theta_hat,
          g_coef = fit_si_k$g_coef_hat,
          Z_new = Z[id_si_test_global, , drop = FALSE],
          Z_mu = fit_si_k$Z_mu,
          Z_sd = fit_si_k$Z_sd,
          N = fit_si_k$N,
          a = fit_si_k$a,
          d = ncol(Z)
        )
      }
      
      r_rf_oof[id_si_test_local] <- Y[id_si_test_global] - as.vector(
        X[id_si_test_global, , drop = FALSE] %*% beta_use) - Rhat_si_test
    }
    
    
    if (is.null(oracle_si)) {
      
      N_si_unique <- unique(N_si_fold[is.finite(N_si_fold)])
      
      if (length(N_si_unique) > 1L) {
        warning(
          "Different spline dimensions were selected across OOF folds. ",
          "The pooled SI orthogonalization uses N_si = ",
          N_si_half,
          "."
        )
      }
    }
    
    if (any(!is.finite(r_rf_oof))) {
      stop(sprintf(
        "Non-finite OOF residuals detected: %d.",
        sum(!is.finite(r_rf_oof))
      ))
    }
    
    id_rf_valid <- id_rf
    
    X_rf <- X[id_rf_valid, , drop = FALSE]
    Z_rf <- Z[id_rf_valid, , drop = FALSE]
    r_rf <- r_rf_oof
    fold_si_rf <- fold_si
    
    
    ############################################################
    ## 1.3 用这个 half 的 OOF residual 训练 RF
    ############################################################
    
    X_mu <- colMeans(X_rf)
    X_sd <- apply(X_rf, 2, sd)
    X_sd[X_sd == 0 | !is.finite(X_sd)] <- 1
    
    Z_mu <- colMeans(Z_rf)
    Z_sd <- apply(Z_rf, 2, sd)
    Z_sd[Z_sd == 0 | !is.finite(Z_sd)] <- 1
    
    X_rf_std <- sweep(sweep(X_rf, 2, X_mu, "-"), 2, X_sd, "/")
    Z_rf_std <- sweep(sweep(Z_rf, 2, Z_mu, "-"), 2, Z_sd, "/")
    
    rf_final <- randomForest::randomForest(
      x = as.data.frame(cbind(X_rf_std, Z_rf_std)),
      y = r_rf,
      ntree = 500,
      mtry = max(1, floor(sqrt(ncol(X_rf_std) + ncol(Z_rf_std)))),
      nodesize = 5
    )
    
    
    ############################################################
    ## 1.4 返回该 half 的完整结果
    ############################################################
    
    list(
      beta_use = beta_use,
      
      id_ce = id_ce,
      id_beta = id_beta,
      id_rf = id_rf,
      id_rf_valid = id_rf_valid,
      
      r_rf = r_rf,
      rf_final = rf_final,
      
      fold_si_rf = fold_si_rf,
      si_fit_list = si_fit_list,
      Z_mu_si = Z_mu_si,
      Z_sd_si = Z_sd_si,
      
      X_mu = X_mu,
      X_sd = X_sd,
      Z_mu = Z_mu,
      Z_sd = Z_sd,
      
      n_ce = length(id_ce),
      n_beta = length(id_beta),
      n_rf = length(id_rf),
      n_rf_valid = length(id_rf_valid),
      K_si = K_si_use,
      N_si = N_si_half,
      N_si_fold = N_si_fold,
      a_si_fold = a_si_fold
    )
  }
  
  
  ############################################################
  ## 2. 两份数据都完整走一遍
  ############################################################
  
  fit_part1 <- run_one_half(split_obj = split1, seed = seed)
  fit_part2 <- run_one_half(split_obj = split2, seed = seed)
  
  
  ############################################################
  ## 3. 一个 half 的 RF 预测另一个 half
  ############################################################
  
  calc_Zn_between_halves <- function(train_fit, test_fit) {
    
    id_test <- test_fit$id_rf_valid
    
    X_test <- X[id_test, , drop = FALSE]
    Z_test <- Z[id_test, , drop = FALSE]
    
    r_test <- as.numeric(test_fit$r_rf)
    fold_test <- as.integer(test_fit$fold_si_rf)
    
    N_si_fold_test <- test_fit$N_si_fold
    a_si_fold_test <- test_fit$a_si_fold
    
    # 正交化方向必须由 opposite half 的第 k 个 SI fit 在当前 test half 的第 k 折协变量上重新评价
    if (do_si_orth && is.null(oracle_si)) {
      
      if (train_fit$K_si != test_fit$K_si) {
        stop("The two half-samples must use the same number of SI folds.")
      }
      
      u_test <- rep(NA_real_, length(id_test))
      theta_score_test <- matrix(
        NA_real_,
        nrow = length(id_test),
        ncol = q
      )
      
      for (k in seq_len(test_fit$K_si)) {
        
        idx_k <- which(fold_test == k)
        if (length(idx_k) == 0L) next
        
        fit_si_ref <- train_fit$si_fit_list[[k]]
        if (is.null(fit_si_ref)) {
          stop(sprintf("Missing opposite-half SI fit for fold %d.", k))
        }
        
        score_k <- build_oof_theta_score(
          Z_new = Z_test[idx_k, , drop = FALSE],
          fit_si = fit_si_ref,
          eps = 1e-12
        )
        
        u_test[idx_k] <- score_k$u
        theta_score_test[idx_k, ] <- score_k$theta_score
      }
      
    } else {
      
      # oracle / 不做 SI 正交化时，仅用于保持后续长度一致
      u_test <- rep(0, length(id_test))
      theta_score_test <- matrix(0, nrow = length(id_test), ncol = q)
    }
    
    ############################################################
    ## 1. 用 train_fit 的参数标准化 test half
    ############################################################
    
    X_test_std <- sweep(
      sweep(X_test, 2, train_fit$X_mu, "-"), 2, train_fit$X_sd, "/"
    )
    
    Z_test_std <- sweep(
      sweep(Z_test, 2, train_fit$Z_mu, "-"), 2, train_fit$Z_sd, "/"
    )
    
    
    ############################################################
    ## 2. RF 预测 residual predictability
    ############################################################
    
    f_pred <- as.numeric(
      predict(
        train_fit$rf_final,
        newdata = as.data.frame(cbind(X_test_std, Z_test_std))
      )
    )
    
    if (any(!is.finite(r_test))) {
      stop("Non-finite test residual.")
    }
    
    valid <- is.finite(f_pred) & is.finite(u_test) &
      apply(X_test_std, 1, function(v) all(is.finite(v)))
    
    if (do_si_orth && is.null(oracle_si)) {
      valid <- valid & apply(theta_score_test, 1, function(v) all(is.finite(v)))
    }
    
    r_use <- r_test[valid]
    f_use <- f_pred[valid]
    u_use <- u_test[valid]
    
    X_use <- X_test_std[valid, , drop = FALSE]
    theta_score_use <- theta_score_test[valid, , drop = FALSE]
    fold_use <- fold_test[valid]
    
    n_use <- length(r_use)
    p_use <- ncol(X_use)
    
    if (n_use <= 5) {
      stop(sprintf("Effective sample size is too small: n_use = %d.", n_use))
    }
    
    
    ############################################################
    ## 3. 线性正交化辅助函数
    ############################################################
    
    linear_orthogonalize <- function(f_input) {
      
      fit_gam <- glmnet::glmnet(
        x = X_use,
        y = f_input,
        alpha = 1,
        intercept = TRUE,
        standardize = FALSE,
        nlambda = 100
      )
      
      W_mat <- as.matrix(X_use %*% fit_gam$beta)
      W_mat <- sweep(W_mat, 2, fit_gam$a0, "+")
      W_mat <- sweep(W_mat, 1, f_input, function(pred, obs) obs - pred)
      
      scale_vec <- apply(W_mat, 2, function(w) sqrt(mean(w^2)))
      scale_vec[!is.finite(scale_vec) | scale_vec < eps_num] <- eps_num
      
      lambda_theory <- C_lambda_beta * sqrt(log(max(p_use, 2L)) / n_use)
      
      index <- which.min(
        abs(fit_gam$lambda - lambda_theory * scale_vec)
      )
      
      gamma_hat <- as.numeric(fit_gam$beta[, index])
      alpha_hat <- as.numeric(fit_gam$a0[index])
      f_orth <- f_input - alpha_hat - as.vector(X_use %*% gamma_hat)
      
      linear_orth_before <- max(
        abs(as.vector(crossprod(X_use, f_input)) / n_use)
      )
      
      linear_orth_after <- max(
        abs(as.vector(crossprod(X_use, f_orth)) / n_use)
      )
      
      mean_before_linear = mean(f_input)
      mean_after_linear = mean(f_orth)
      
      list(
        f_orth = f_orth,
        gamma_hat = gamma_hat,
        alpha_hat = alpha_hat,
        Orthogonal_beta = sum(abs(gamma_hat) > 1e-8),
        linear_orth_before = linear_orth_before,
        linear_orth_after = linear_orth_after,
        lambda_theory = lambda_theory,
        lambda_selected = fit_gam$lambda[index],
        scale_selected = scale_vec[index],
        mean_before_linear = mean_before_linear,
        mean_after_linear = mean_after_linear
      )
    }
    
    ############################################################
    ## 4. 最终方法：spline + linear + theta-score 联合正交化
    ############################################################
    
    spline_orth_before <- NA_real_
    spline_orth_after <- NA_real_
    theta_orth_before <- NA_real_
    theta_orth_after <- NA_real_
    linear_orth_before_final <- NA_real_
    linear_orth_after_final <- NA_real_
    
    path_index <- NA_real_
    path_ratio <- NA_real_
    path_gap <- NA_real_
    linear_kkt_ratio <- NA_real_
    theta_kkt_ratio <- NA_real_
    theta_kkt_ratio_by_fold <- rep(NA_real_, test_fit$K_si)
    
    
    if (do_si_orth && is.null(oracle_si)) {
      
      joint_fit <- joint_orthogonalize_score(
        f_raw = f_use,
        X = X_use,
        u = u_use,
        theta_score = theta_score_use,
        fold_id = fold_use,
        N_si = test_fit$N_si,
        df_fallback = df_si_orth,
        C_lambda_beta = C_lambda_beta,
        C_lambda_theta = C_lambda_theta,
        eps = eps_num
      )
      
      f_orth_joint <- joint_fit$f_orth   # joint
      
      spline_orth_before <- joint_fit$spline_orth_before
      spline_orth_after <- joint_fit$spline_orth_after
      
      theta_orth_before <- joint_fit$theta_orth_before
      theta_orth_after <- joint_fit$theta_orth_after
      
      linear_orth_before_final <- joint_fit$linear_orth_before
      linear_orth_after_final <- joint_fit$linear_orth_after
      
      path_index <- joint_fit$path_index
      path_ratio <- joint_fit$path_ratio
      path_gap <- joint_fit$path_gap
      linear_kkt_ratio <- joint_fit$linear_kkt_ratio
      theta_kkt_ratio <- joint_fit$theta_kkt_ratio
      theta_kkt_ratio_by_fold <- joint_fit$theta_kkt_ratio_by_fold
      
      
    } else {
      
      # oracle / 关闭 SI 正交化时，保持原来的仅线性正交化
      linear_final_fit <- linear_orthogonalize(f_use)
      f_orth_joint <- linear_final_fit$f_orth
      
      linear_orth_before_final <- linear_final_fit$linear_orth_before
      linear_orth_after_final <- linear_final_fit$linear_orth_after
    }
    
    
    ############################################################
    ## 7. 仅线性正交化结果，用于对比
    ############################################################
    
    linear_only_fit <- linear_orthogonalize(f_use)
    f_orth_linear <- linear_only_fit$f_orth
    
    
    ############################################################
    ## 8. 计算统计
    ############################################################
    
    calc_stat <- function(f_score) {
      
      rw_use <- as.numeric(r_use * f_score)
      n_score <- length(rw_use)
      
      if (any(!is.finite(rw_use))) {
        stop("The score contains non-finite values.")
      }
      
      if (n_score <= 5L) {
        stop("Too few observations for HAC estimation.")
      }
      
      ############################################################
      ## 1. 检验统计量的分子
      ############################################################
      
      S_use <- sum(rw_use) / sqrt(n_score)
      
      ############################################################
      ## 2. 对 score 中心化，用于估计长程方差
      ############################################################
      
      rw_center <- rw_use - mean(rw_use)
      
      ############################################################
      ## 3. Newey-West plug-in bandwidth + n^(1/3) 理论截断
      ############################################################
      
      n_lag <- floor(4 * (n_score / 100)^(2 / 9))
      n_lag <- max(1L, min(as.integer(n_lag), n_score - 1L))
      
      sigma_hat <- numeric(n_lag + 1L)
      sigma_hat[1L] <- sum(rw_center^2) / n_score
      
      for (j in seq_len(n_lag)) {
        sigma_hat[j + 1L] <-
          sum(
            rw_center[(j + 1L):n_score] *
              rw_center[seq_len(n_score - j)]
          ) / n_score
      }
      
      s0_hat <- sigma_hat[1L] + 2 * sum(sigma_hat[-1L])
      s1_hat <- 2 * sum(seq_len(n_lag) * sigma_hat[-1L])
      
      alpha_hat_hac <- NA_real_
      
      if (!is.finite(s0_hat) ||
          !is.finite(s1_hat) ||
          abs(s0_hat) <= eps_num) {
        
        m_raw <- as.numeric(n_lag)
        
      } else {
        
        alpha_hat_hac <- (s1_hat / s0_hat)^2
        m_raw <- 1.1447 * (alpha_hat_hac * n_score)^(1 / 3)
        
        if (!is.finite(m_raw)) {
          m_raw <- as.numeric(n_lag)
        }
      }
      
      m_lower <- min(n_score - 1L, max(1L, floor(hac_c_lower * n_score^(1 / 3))))
      m_upper <- min(n_score - 1L, max(m_lower, ceiling(hac_c_upper * n_score^(1 / 3))))
      
      m <- floor(m_raw)
      m <- max(m_lower, min(as.integer(m), m_upper))
      m <- max(0L, min(as.integer(m), n_score - 1L))
      
      ############################################################
      ## 4. Bartlett HAC long-run variance
      ############################################################
      
      V_use <- sum(rw_center^2) / n_score
      
      if (m >= 1L) {
        for (k in seq_len(m)) {
          gamma_k <- sum(rw_center[(k + 1L):n_score] * rw_center[seq_len(n_score - k)]) / n_score
          weight_k <- 1 - k / (m + 1)
          V_use <- V_use + 2 * weight_k * gamma_k
        }
      }
      
      if (!is.finite(V_use) || V_use <= eps_num) {
        stop(sprintf(
          "The HAC variance is non-positive or non-finite: %.6e.",
          V_use
        ))
      }
      
      ############################################################
      ## 5. Studentized statistic
      ############################################################
      
      Zn_use <- S_use / sqrt(V_use)
      p_use_value <- 1 - pnorm(Zn_use)
      
      list(
        Zn = Zn_use,
        p_value = p_use_value
      )
    }
    
    # 仅线性正交化
    stat_linear <- calc_stat(f_orth_linear)
    
    # 联合正交化
    stat_joint <- calc_stat(f_orth_joint) 
    
    Zn <- stat_joint$Zn
    p_value <- stat_joint$p_value
    
    
    ############################################################
    ## 9. 补回 full-length 的正交化分数
    ############################################################
    
    f_orth_linear_full <- rep(NA_real_, length(f_pred))
    f_orth_joint_full <- rep(NA_real_, length(f_pred))
    
    f_orth_linear_full[valid] <- f_orth_linear
    f_orth_joint_full[valid] <- f_orth_joint
    
    
    ############################################################
    ## 10. 返回结果
    ############################################################
    
    list(
      # 最终结果：联合正交化
      Zn = Zn,
      p_value = p_value,
      n_use = n_use,
      
      # 仅线性正交化
      Zn_linear = stat_linear$Zn,
      p_linear = stat_linear$p_value,
      
      # 线性正交化诊断
      linear_orth_before_final = linear_orth_before_final,
      linear_orth_after_final = linear_orth_after_final,
      
      # KKT
      linear_kkt_ratio = linear_kkt_ratio,
      theta_kkt_ratio = theta_kkt_ratio,
      theta_kkt_ratio_by_fold = theta_kkt_ratio_by_fold
    )
  }
  
  
  ############################################################
  ## 4. 两个方向
  ############################################################
  
  stat_12 <- calc_Zn_between_halves(
    train_fit = fit_part1,
    test_fit = fit_part2
  )
  
  stat_21 <- calc_Zn_between_halves(
    train_fit = fit_part2,
    test_fit = fit_part1
  )
  
  
  ############################################################
  ## 5. CCT 组合两个 p 值
  ############################################################
  
  cct_pvalue <- function(pvals, weights = NULL) {
    
    pvals <- as.numeric(pvals)
    if (any(!is.finite(pvals))) return(NA_real_)
    
    if (length(pvals) == 0) return(NA_real_)
    
    if (is.null(weights)) {
      weights <- rep(1 / length(pvals), length(pvals))
    } else {
      weights <- weights / sum(weights)
    }
    
    pvals <- pmin(pmax(pvals, 1e-15), 1 - 1e-15)
    
    T_cct <- sum(weights * tan((0.5 - pvals) * pi))
    p_cct <- 0.5 - atan(T_cct) / pi
    
    pmin(pmax(p_cct, 0), 1)
  }
  
  # 最终结果：联合正交化
  p_CCT <- cct_pvalue(
    pvals = c(stat_12$p_value, stat_21$p_value)
  )
  
  # 仅线性正交化
  p_CCT_linear <- cct_pvalue(
    pvals = c(stat_12$p_linear, stat_21$p_linear)
  )
  
  
  ############################################################
  ## 6. 输出
  ############################################################
  
  list(
    p_CCT = p_CCT,
    p_CCT_linear = p_CCT_linear,
    
    id_part1 = id_part1,
    id_part2 = id_part2,
    
    gap = gap,
    
    fit_part1 = fit_part1,
    fit_part2 = fit_part2,
    
    stat_12 = stat_12,
    stat_21 = stat_21
  )
}
