library(tidyverse)
library(surveillance)
library(googlesheets4)
library(googledrive)
library(DT)
library(purrr)

## Load data GG Drive
load_data <- function(url, sheet) {
  gs4_deauth()
  df <- read_sheet(url, sheet = sheet)
  return(df)
}

calc_serfling <- function(df_h, df_t) {
  tryCatch({
    df_h <- mutate(df_h, time = row_number())
    mod <- lm(cases ~ time + sin(2*pi*week/52) + cos(2*pi*week/52), data = df_h)
    df_t <- mutate(df_t, time = max(df_h$time) + row_number())
    as.numeric(predict(mod, df_t, interval="prediction", level=0.95)[,"upr"])
  }, error = function(e) rep(NA, nrow(df_t)))
}

calc_mem <- function(df_h, df_t) {
  peaks <- df_h %>% group_by(year) %>% summarise(m = max(cases, na.rm=T)) %>% pull(m)
  peaks <- peaks[peaks > 0]
  if(length(peaks) < 2) return(rep(NA, nrow(df_t)))
  thr <- exp(mean(log(peaks)) + 1.96 * sd(log(peaks)))
  rep(thr, nrow(df_t))
}

calc_cdc <- function(df_hist, df_target) {
  stats <- df_hist %>%
    group_by(week) %>%
    summarise(
      mean_val = mean(cases, na.rm = TRUE),
      sd_val = sd(cases, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(threshold = mean_val + 2 * sd_val)
  
  df_merged <- df_target %>%
    left_join(stats, by = "week")
  
  return(df_merged$threshold)
}

## City
run_algos_city <- function(df, target_years) {
  final_list <- list()
  max_w <- 10
  
  for (yr in target_years) {
    if (yr == 2021) next 
    
    hist_years <- setdiff((min(df$year)):(yr - 1), 2021)
    if (length(hist_years) < 2) next 
    
    df_raw <- df %>% filter(year %in% c(hist_years, yr))
    
    df_norm <- df_raw %>%
      filter(week <= 52) %>%
      complete(year = c(hist_years, yr), week = 1:52, fill = list(cases = 0)) %>%
      arrange(year, week)
    
    first_year <- min(hist_years)
    df_pad <- tibble(
      year = first_year - 1,
      week = (52 - max_w + 1):52,
      cases = 0
    )
    
    df_comb <- bind_rows(df_pad, df_norm)
    
    sts_obj <- sts(observed = df_comb$cases, 
                   start = c(first_year - 1, min(df_pad$week)), 
                   frequency = 52)
    dp_obj <- sts2disProg(sts_obj)
    
    idx_start_target <- nrow(df_pad) + length(hist_years)*52 + 1
    idx_end_target   <- nrow(df_comb)
    range_eval <- idx_start_target:idx_end_target
    
    b_curr <- length(hist_years)
    df_target <- df_norm %>% filter(year == yr)
    
    res_rki1 <- algo.rki1(dp_obj, control = list(range = range_eval, b = b_curr, w = 1, actY = FALSE))
    res_rki2 <- algo.rki2(dp_obj, control = list(range = range_eval, b = b_curr, w = 1, actY = TRUE))
    
    res_bayes1 <- algo.bayes1(dp_obj, control = list(range = range_eval, b = b_curr, w = 1, actY = TRUE, alpha = 0.05))
    res_bayes2 <- algo.bayes2(dp_obj, control = list(range = range_eval, b = b_curr, w = 1, actY = TRUE, alpha = 0.05))
    
    res_far_old <- algo.farrington(dp_obj, control = list(range = range_eval, b = b_curr, w = 1, alpha = 0.05))
    res_far_flex <- farringtonFlexible(sts_obj, control = list(range = range_eval, b = b_curr, w = 1, alpha = 0.05))
    
    res_earsC1 <- earsC(sts_obj, control = list(range = range_eval, method = "C1", baseline = 4))
    res_earsC2 <- earsC(sts_obj, control = list(range = range_eval, method = "C2", baseline = 4))
    res_earsC3 <- earsC(sts_obj, control = list(range = range_eval, method = "C3", baseline = 4))
    
    res_glrnb <- algo.glrnb(dp_obj, control = list(range = range_eval, c.ARL = 5, dir = "inc"))
    res_glrpois <- algo.glrpois(dp_obj, control = list(range = range_eval, c.ARL = 5, dir = "inc"))
    
    res_cusum <- algo.cusum(dp_obj, control = list(range = range_eval, k = 1.04, h = 2.26))
    res_cusum_rossi <- algo.cusum(dp_obj, control = list(range = range_eval, k = 1.04, h = 2.26, trans = "rossi"))
    res_cusum_negbin <- algo.cusum(dp_obj, control = list(range = range_eval, k = 1.04, h = 2.26, trans = "pearsonNegBin"))
    res_cusum_rossiGLM <- algo.cusum(dp_obj, control = list(range = range_eval, k = 1.04, h = 2.26, trans = "rossi", m = "glm"))
    res_cusum_negbinGLM <- algo.cusum(dp_obj, control = list(range = range_eval, k = 1.04, h = 2.26, trans = "pearsonNegBin", m = "glm"))
    
    df_res <- df_target %>%
      mutate(
        ub_cdc      = calc_cdc(df_norm %>% filter(year != yr), df_target),
        ub_rki1     = res_rki1$upperbound,
        ub_rki2     = res_rki2$upperbound,
        ub_bayes1   = res_bayes1$upperbound,
        ub_bayes2   = res_bayes2$upperbound,
        ub_far_old  = res_far_old$upperbound,
        ub_far_flex = upperbound(res_far_flex),
        ub_earsC1   = upperbound(res_earsC1),
        ub_earsC2   = upperbound(res_earsC2),
        ub_earsC3   = upperbound(res_earsC3),
        ub_glrnb    = res_glrnb$upperbound,
        ub_glrpois  = res_glrpois$upperbound,
        ub_cusum    = res_cusum$upperbound,
        ub_cusum_rossi = res_cusum_rossi$upperbound,
        ub_cusum_rossiGLM = res_cusum_rossiGLM$upperbound,
        ub_cusum_negbin = res_cusum_negbin$upperbound,
        ub_cusum_negbinGLM = res_cusum_negbinGLM$upperbound,
        ub_serfling = calc_serfling(df_norm %>% filter(year != yr), df_target),
        ub_mem      = calc_mem(df_norm %>% filter(year != yr), df_target)
      )
    
    final_list[[as.character(yr)]] <- df_res
  }
  
  return(bind_rows(final_list))
}

## Districts
run_algos_by_district <- function(df_all_districts, target_years) {
  
  list_districts <- unique(df_all_districts$districts)
  
  final_results <- list()
  
  for (d in list_districts) {
    message(paste("Đang chạy quận:", d, "...")) 
    
    df_sub <- df_all_districts %>% 
      filter(districts == d) %>%
      
      arrange(year, week)

    if (nrow(df_sub) < 52) {
      warning(paste("Quan", d, "khong du du lieu. Bo qua."))
      next
    }
    
    res_sub <- tryCatch({
      run_algos_city(df_sub, target_years)
    }, error = function(e) {
      warning(paste("Loi tai quan:", d, "-", e$message))
      return(NULL)
    })
    
    if (!is.null(res_sub) && nrow(res_sub) > 0) {
      res_sub$districts <- d
      final_results[[d]] <- res_sub
    }
  }
  
  return(bind_rows(final_results))
}
  
  
  