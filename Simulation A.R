# Spline single-index
library(splines2)
library(MASS)
library(openxlsx)
library(caret)
library(glmnet)
library(randomForest)
library(reticulate)


source("code/SIM_functions.R")
source("code/Orthogonalization.R")
source("code/PLSIM_test.R")

# Specify the Python interpreter; change this path to your local Python environment.
use_python("D:/Anaconda_envs/envs/myenv/python.exe", required = TRUE)
# Load Python file
source_python("code/python_ts.py")


generate_PLSIM_ts <- function(
    n = 800, p = 20, q = 20, s = 20, 
    rho_z = 0.5, rho_x = 0.5, sigma_u = 1.0,
    scenario = c("H0", "case3"),
    signal_scale = 0, c = 1
){
  scenario <- match.arg(scenario)
  
  # =========================================================
  # 1) Z: q维独立 AR(1)（稳态方差≈1）,每一列都是一个 AR(1) 模型, q 个彼此独立的 AR(1) 时间序列
  # =========================================================
  Z <- matrix(0, n, q)
  Ez <- matrix(rnorm(n * q), n, q)
  Z[1, ] <- Ez[1, ]
  
  for (j in 1:q) {
    for (t in 2:n) {
      Z[t, j] <- rho_z * Z[t - 1, j] + sqrt(1 - rho_z^2) * Ez[t, j]
    }
  }
  
  
  # =========================================================
  # 2) X: 条件均值 muX = Z * Bx + AR(1) 偏差
  # =========================================================
  # Bx  <- build_Bx(q, p, signal_scale = signal_scale, m = m)  # 设置X和Z的相关性，会引入X的自相关性；
  Bx = signal_scale * diag(p)  # 设置X和Z的相关性，同时不引入X的自相关；
  muX <- Z %*% Bx
  
  X  <- matrix(0, n, p)
  Ex <- matrix(rnorm(n*p), n, p)
  X[1,] <- muX[1,] + Ex[1,]
  for (t in 2:n){
    X[t,] <- muX[t,] +
      rho_x * (X[t-1,] - muX[t-1,]) +
      sqrt(1 - rho_x^2) * Ex[t,]
  }
  
  # 2) 单位向量系数
  beta <- numeric(p)               # 全部置 0
  beta[1:s] <- rnorm(s)            # 生成 s 个非零
  beta[1:s] <- beta[1:s] / sqrt(sum(beta[1:s]^2))
  
  theta <- numeric(q)                 # 全部置0
  theta[(q-s+1):q] <- rep(1,s)        # 后s个元素生成非零值
  theta <- theta / sqrt(sum(theta^2)) # 单位化
  
  # 3) 单指标 
  t_val  <- as.vector(Z %*% theta)
  linear_val = as.vector(X %*% beta)
  # Y_base <- linear_val + t_val^2  # data1
  # Y_base <- linear_val + 0.5*exp(-t_val^2)  # data2
  Y_base <- linear_val + cos(2*t_val) + 0.5*t_val^2  # data3
  
  dev <- numeric(n)   # 默认无偏离，满足H0
  
  # 4) H1 场景
  if (scenario == "case3") {
    
    if (p == 10 && q == 10) {
      
      # 低维：Example 3
      dev <- (abs(X[,1]) + abs(Z[,1]))^2 +
        exp(X[,2]) + Z[,2]*X[,3] + abs(Z[,3])^3
      
    } else if (p == q && p %in% c(1000, 2000)) {
      
      # 高维：Model HD3
      dev <- 7*Z[,1]^2 + abs(X[,1]) +
        abs(X[,2]) + abs(Z[,2])
      
    } else {
      stop("未定义的维度设置")
    }
  }
  
  # 5) 生成响应
  eps <- rnorm(n, mean = 0, sd = sigma_u)
  Y   <- Y_base + c*dev + eps   # 通过 c 增加H1偏离程度
  
  list(Y = as.vector(Y), X = X, Z = Z, beta = beta, theta = theta, eps = eps, Bx = Bx)
}



#######################################################################################
start_time <- Sys.time()

n = 800
p = 10
q = 10
s = 10

data <- generate_PLSIM_ts(n = n, p = p, q = q, s = s, scenario = "case3", c = 1, 
                          rho_z = 0.8, rho_x = 0.8, signal_scale = 0)
Y = data[[1]]
X = data[[2]]
Z = data[[3]]

result = PLSIM_test_split_ratio_cross2_ts(
  Y=Y, X=X, Z=Z, 
  split_ratio = c(4, 4),
  K_si = 4,
  C_lambda_beta = 0.3,
  C_lambda_theta = 1,
  gap = 5,
  seed = 1
)


end_time = Sys.time()
end_time - start_time  # 15.87349 secs

# 最终两步正交化结果
cat("joint: \n")
result$stat_12$Zn
result$stat_12$p_value

result$stat_21$Zn
result$stat_21$p_value

result$p_CCT

cat("linear: \n")
result$stat_12$Zn_linear
result$stat_12$p_linear

result$stat_21$Zn_linear
result$stat_21$p_linear

result$p_CCT_linear




# =========================
# 并行计算：PLSIM 拟合优度检验的 size / power
# =========================
library(doParallel)
library(foreach)
library(doRNG)

# ====== Python 解释器与脚本路径 ======
py_bin  <- "D:/Anaconda_envs/envs/myenv/python.exe"  # change this path to your local Python environment.
py_file <- "code/python_ts.py"

.EXPORT_FUNS <- c(
  "generate_PLSIM_ts", 
  "build_si_basis", "build_oof_theta_score", 
  "joint_orthogonalize_score", "PLSIM_test_split_ratio_cross2_ts",
  # 单指标相关函数
  ".safe_inv_sym", "Fd", "dot_Fd", "calc_N", "normalize_theta",
  "build_spline_matrices", "compute_dot_Bp", "compute_R_hat", "compute_hat_S_star",
  ".project_theta_minus_d", ".make_safe_gr", "estimate_theta", "estimate_g",
  "spline_single_index",  "predict_si", 
  "set_R_python_seeds"
)


# -------------------------
# 统一设置 R + Python seed
# -------------------------
set_R_python_seeds <- function(seed, use_python = TRUE) {
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  
  if (use_python) {
    reticulate::py_run_string(sprintf("
import os, random
import numpy as np
try:
    import torch
except Exception:
    torch = None

seed = %d
os.environ['PYTHONHASHSEED'] = str(seed)
random.seed(seed)
np.random.seed(seed)

if torch is not None:
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
    try:
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.benchmark = False
    except Exception:
        pass
", as.integer(seed)))
  }
}


# -------------------------
# 并行功效 / 尺寸评估
# -------------------------
run_power_PLSIM <- function(
    n_sims   = 1000,
    n        = 800,
    p        = 10,
    q        = 10,
    s        = 10,
    rho_x    = 0.5,
    rho_z    = 0.5,
    signal_scale = 0,
    scenario = "H0",
    c        = 1,
    split_ratio = c(4, 4),
    K_si     = 4,
    C_lambda_beta = 0.3,
    C_lambda_theta = 1,
    alpha    = 0.05,
    ncores   = 5,
    gap      = 5,
    seed     = 2026
) {
  
  Sys.setenv(
    OMP_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1",
    NUMEXPR_NUM_THREADS = "1"
  )
  
  suppressPackageStartupMessages(library(reticulate))
  Sys.setenv(RETICULATE_PYTHON = py_bin)
  reticulate::use_python(py_bin, required = TRUE)
  reticulate::source_python(py_file)
  
  set_R_python_seeds(seed)
  
  cl <- parallel::makeCluster(ncores)
  on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
  registerDoParallel(cl)
  
  parallel::clusterExport(
    cl,
    c("py_bin", "py_file", "seed", "set_R_python_seeds"),
    envir = environment()
  )
  
  parallel::clusterEvalQ(cl, {
    Sys.setenv(
      OMP_NUM_THREADS = "1",
      OPENBLAS_NUM_THREADS = "1",
      MKL_NUM_THREADS = "1",
      VECLIB_MAXIMUM_THREADS = "1",
      NUMEXPR_NUM_THREADS = "1"
    )
    
    library(reticulate)
    Sys.setenv(RETICULATE_PYTHON = py_bin)
    reticulate::use_python(py_bin, required = TRUE)
    reticulate::source_python(py_file)
    
    set_R_python_seeds(seed)
    NULL
  })
  
  res_mat <- foreach(
    sim = 1:n_sims,
    .combine = rbind,
    .packages = c(
      "splines2", "MASS", "caret", "randomForest",
      "glmnet", "reticulate"
    ),
    .export = .EXPORT_FUNS,
    .options.RNG = seed
  ) %dorng% {
    
    set_R_python_seeds(seed + sim)
    
    tryCatch({
      
      dat <- generate_PLSIM_ts(
        n = n, p = p, q = q, s = s, scenario = scenario, c = c,
        rho_z = rho_z, rho_x = rho_x, signal_scale = signal_scale
      )
      
      
      test_res <- PLSIM_test_split_ratio_cross2_ts(
        Y = dat$Y, X = dat$X, Z = dat$Z, split_ratio = split_ratio, 
        K_si = K_si, C_lambda_beta = C_lambda_beta, C_lambda_theta = C_lambda_theta, 
        gap = gap, seed = seed + sim
      )
      
      c(
        p_12 = as.numeric(test_res$stat_12$p_value),
        p_21 = as.numeric(test_res$stat_21$p_value),
        p_CCT = as.numeric(test_res$p_CCT)
        
      )
      
    }, error = function(e) {
      
      warning(sprintf("sim %d failed: %s", sim, conditionMessage(e)))
      
      c(
        p_12 = NA_real_,
        p_21 = NA_real_,
        p_CCT = NA_real_
      )
    })
  }
  
  
  list(
    settings = list(
      n_sims = n_sims, n = n, p = p, q = q, s = s, scenario = scenario,
      split_ratio = split_ratio, K_si = K_si, C_lambda_beta = C_lambda_beta,
      C_lambda_theta = C_lambda_theta, alpha = alpha, ncores = ncores, seed = seed
    ),
    
    raw_results = as.data.frame(res_mat)
    
  )
}



# =========================
# 使用示例
# =========================
start_time <- Sys.time()

res_size <- run_power_PLSIM(
  n_sims = 1000, n = 1000, p = 1000, q = 1000, s = 20, scenario = "H0", rho_z = 0.5, rho_x = 0.5, 
  signal_scale = 0.3, split_ratio = c(4, 4), K_si = 4, C_lambda_beta = 0.3, C_lambda_theta = 1,
  alpha = 0.05, ncores = 5, seed = 123
)


save(res_size, file = "results/Simulation_A/cos(2t_val)+0.5t_val^2_n1000_d1000_s20_nsims1000_rho05_sigma03_size.RData")

end_time <- Sys.time()
end_time - start_time  

rr <- res_size$raw_results

# size/power
mean(rr$p_12 < 0.05, na.rm = TRUE)
mean(rr$p_21 < 0.05, na.rm = TRUE)
mean(rr$p_CCT < 0.05, na.rm = TRUE)


