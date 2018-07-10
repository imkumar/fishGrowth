
#' Fit von Bertalanffy growth function
#'
#' @param data   Input data
#'
#' @export
#'

vbgf <- function(data) {

    ndiv <- length(unique(ag_cod$NAFOdiv))
    Lobs <- ag_cod$length
    div <- ag_cod$NAFOdiv
    divs <- as.numeric(as.factor(div)) - 1

    nobs <- length(Lobs)
    age <- as.numeric(ag_cod$age)
    wt <- ag_cod$weight
    yr <- ag_cod$year - min(ag_cod$year) + 1
    ny <- length(unique(yr))
    iy <- yr - 1

    ag_cod$cohort <- ag_cod$survey.year - age
    nc <- length(unique(ag_cod$cohort))
    ic <- as.numeric(as.factor(ag_cod$cohort)) - 1

    fulldate <- as.Date(paste(ag_cod$month,ag_cod$day,ag_cod$year,sep = '/'), "%m/%d/%Y")
    origin <- as.Date(paste("01","01", ag_cod$year, sep = '/'), "%m/%d/%Y")
    ndays <- fulldate - origin
    age_days <- round(as.numeric(age + ndays/365), 2)

    sx <- unclass(as.factor(ag_cod$sex)) - 1
    maturity <- ag_cod$maturity
    maturity <- ifelse(maturity == 100 | maturity == 500, 1, 2)
    mat <- unclass(as.factor(maturity)) - 1

    ############################################################################

    tmb.data <- list(
        ndiv = ndiv,
        nobs = nobs,
        Lobs = Lobs,
        age = age_days,
        Wobs = wt,
        divs = divs,
        nc = nc,
        ic = ic,
        sx = as.numeric(sx),
        mat = as.numeric(mat))

    parameters <- list(
        log_Linf = 5,
        log_k = -1,
        t0 = -0.1,
        log_a = -1,
        log_b = log(3),
        log_std_me1 = log(0.1),
        log_std_me2 = log(0.1),
        log_std_log_Linf = log(0.1),
        log_std_log_k = log(0.1),
        log_std_log_a = log(0.1),
        log_std_log_b = log(0.1),
        log_std_log_Linfc = rep(log(0.1),tmb.data$ndiv),
        log_std_log_kc = rep(log(0.1),tmb.data$ndiv),
        log_std_log_ac = rep(log(0.1),tmb.data$ndiv),
        log_std_log_bc = rep(log(0.1),tmb.data$ndiv),
        log_maty = -1,
        log_sex = -1,
        rlog_Linf = rep(0,tmb.data$ndiv),
        rlog_k = rep(0,tmb.data$ndiv),
        rlog_a = rep(0,tmb.data$ndiv),
        rlog_b = rep(0,tmb.data$ndiv),
        rlog_Linfc = matrix(0, nrow = tmb.data$ndiv, ncol = tmb.data$nc, byrow = TRUE),
        rlog_kc = matrix(0, nrow = tmb.data$ndiv, ncol = tmb.data$nc, byrow = TRUE),
        rlog_ac = matrix(0, nrow = tmb.data$ndiv, ncol = tmb.data$nc, byrow = TRUE),
        rlog_bc = matrix(0, nrow = tmb.data$ndiv, ncol = tmb.data$nc, byrow = TRUE)
    )

    parameters.L <- list(
        log_Linf = 1,
        log_k = -10,
        t0 = -1,
        log_a = -37,
        log_b = log(1),
        log_std_me1 = log(0.01),
        log_std_me2 = log(0.01),
        log_std_log_Linf = log(0.01),
        log_std_log_k = log(0.01),
        log_std_log_a = log(0.01),
        log_std_log_b = log(0.01),
        log_std_log_Linfc = rep(log(0.01),tmb.data$ndiv),
        log_std_log_kc = rep(log(0.01),tmb.data$ndiv),
        log_std_log_ac = rep(log(0.01),tmb.data$ndiv),
        log_std_log_bc = rep(log(0.01),tmb.data$ndiv),
        log_maty = -10,
        log_sex = -10
    )

    parameters.U <- list(
        log_Linf = 10,
        log_k = 1,
        t0 = 0,
        loga = 1,
        log_b = log(10),
        log_std_me1 = log(10),
        log_std_me2 = log(10),
        log_std_log_Linf = log(10),
        log_std_log_k = log(10),
        log_std_log_a = log(10),
        log_std_log_b = log(10),
        log_std_log_Linfc = rep(log(10), tmb.data$ndiv),
        log_std_log_kc = rep(log(10), tmb.data$ndiv),
        log_std_log_ac = rep(log(10), tmb.data$ndiv),
        log_std_log_bc = rep(log(10), tmb.data$ndiv),
        log_maty = 10,
        log_sex = 10
    )

    lower = unlist(parameters.L)
    upper = unlist(parameters.U)

    rname = c("rlog_Linf", "rlog_k", "rlog_a", "rlog_b", "rlog_Linfc",
              "rlog_kc", "rlog_ac", "rlog_bc")

    obj <- MakeADFun(tmb.data, parameters, random = rname, DLL = "fishGrowth",
                     inner.control = list(maxit = 100, trace = TRUE))
    obj$gr(obj$par)

    opt <- nlminb(obj$par, obj$fn, obj$gr, lower = lower, upper = upper,
                  control = list(trace = 0, iter.max = 5000, eval.max = 10000))

    opt$message
    obj$gr(opt$par)

    rep <- obj$report()

    sd.rep <- sdreport(obj)


}

