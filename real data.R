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



########################################################################################

library(VGAMextra)

# 1. 读取数据
data("HKdata", package = "VGAMextra")
dat <- as.data.frame(HKdata)

dat$date <- as.Date(dat$date)
dat <- dat[dat$date >= as.Date("1994-01-01") & dat$date <= as.Date("1996-12-31"), ]
dat <- dat[order(dat$date), ]

# 2. 标准化六个环境变量
vars <- c("so2", "no2", "rsp", "o3", "temp", "hum")

for (v in vars) {
  dat[[paste0(v, "_std")]] <- as.numeric(scale(dat[[v]]))
}

# 3. 构造星期虚拟变量
dow <- as.POSIXlt(dat$date)$wday

for (k in 1:6) {
  dat[[paste0("D", k)]] <- as.numeric(dow == k)
}

# 4. 构造滞后变量
lag_n <- function(x, k) {
  c(rep(NA_real_, k), head(as.numeric(x), -k))
}

dat$NO2_lag1  <- lag_n(dat$no2_std, 1)

dat$Temp_lag2 <- lag_n(dat$temp_std, 2)
dat$Temp_lag5 <- lag_n(dat$temp_std, 5)

dat$Hum_lag4  <- lag_n(dat$hum_std, 4)
dat$Hum_lag7  <- lag_n(dat$hum_std, 7)

# 5. 构造响应变量
# 每天的住院总人数
Y_raw <- dat$cardio + dat$resp

t <- seq_along(Y_raw)

trend <- ksmooth(
  x = t,
  y = Y_raw,
  kernel = "normal",
  bandwidth = 30,
  x.points = t
)$y

dat$Y <- Y_raw - trend  # 扣除时间趋势

# 6. 构造 Y、X、Z，对应Xia2006文章p1172页选出来的最终模型
x_names <- c("D1", "D2", "D3", "D4", "D5", "D6", "NO2_lag1")
z_names <- c( "Temp_lag2", "Temp_lag5", "Hum_lag4", "Hum_lag7")

valid <- complete.cases(dat[, c("Y", x_names, z_names)])

Y <- as.numeric(dat$Y[valid])
X <- as.matrix(dat[valid, x_names])
Z <- as.matrix(dat[valid, z_names])

date <- dat$date[valid]

# 7. 检查
cat("n =", length(Y), "| p =", ncol(X), "| q =", ncol(Z), "\n")

stopifnot(
  is.numeric(X), is.numeric(Z),
  all(is.finite(Y)),
  all(is.finite(X)),
  all(is.finite(Z))
)


# =====================================================================
start_time = Sys.time()

result = PLSIM_test_split_ratio_cross2_ts(
  Y=Y, X=X, Z=Z, 
  split_ratio = c(4, 4),
  K_si = 4,
  C_lambda_beta = 0.3,
  C_lambda_theta = 1,
  gap = 10,
  seed = 1
)

end_time = Sys.time()
end_time - start_time  # 12.22344 secs

cat("joint: \n")
result$stat_12$Zn         # -1.33564
result$stat_12$p_value    # 0.9091665

result$stat_21$Zn         # -0.5743284
result$stat_21$p_value    # 0.7171272

result$p_CCT              # 0.8591483





# # 1. R 环境信息
# capture.output(
#   sessionInfo(),
#   file = "R_session_info.txt"
# )
# 
# # 2. Python 包版本
# py <- reticulate::py_config()$python
# 
# system2(
#   py,
#   args = c("-m", "pip", "freeze"),
#   stdout = "requirements.txt"
# )
# 
# # 查看保存位置
# getwd()
# 
# 
# file.exists("R_session_info.txt")
# file.exists("requirements.txt")


