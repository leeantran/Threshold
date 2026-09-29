library(tidyverse)
library(ggrepel)

# ==========================================
# 1. LOAD DATA
# ==========================================
dengue_sim <- readRDS("Dengue_Cases.rds")
dengue_perf <- read.csv("Dengue_Perf.csv") %>% mutate(Disease = "Dengue")

hfmd_sim <- readRDS("HFMD_Cases.rds")
hfmd_perf <- read.csv("HFMD_Perf.csv") %>% mutate(Disease = "HFMD")

# ==========================================
# 2. PLOT 1: SIMULATED CASES (100 Sim Lines per Facet)
# ==========================================
# Function to plot 100 simulation lines (median across k1, k2)
plot_sim_data <- function(df, disease_name, file_name) {
  p_data <- df %>%
    # Group by scenario, simulation ID (sim), and time (t)
    group_by(scenario, sim, t) %>%
    summarise(
      # Calculate median of 'n' across all k1, k2 combinations
      median_n = median(n, na.rm = TRUE),
      mu = first(mu),
      .groups = "drop"
    ) %>%
    mutate(scenario = factor(scenario, levels = paste("Sce", 1:32)))
  
  p <- ggplot(p_data, aes(x = t)) +
    # Draw 100 simulation lines with low alpha (transparency) to show density
    geom_line(aes(y = median_n, group = sim), color = "steelblue", alpha = 0.15, linewidth = 0.4) +
    # Draw the expected baseline (mu) on top with a distinct color (e.g., red dashed line)
    geom_line(aes(y = mu), color = "red", linetype = "dashed", linewidth = 0.6) +
    facet_wrap(~ scenario, ncol = 4, scales = "free_y") +
    theme_bw() +
    theme(
      legend.position = "none", # Hide legend since we only have sim lines and mu
      strip.text = element_text(face = "bold")
    ) +
    labs(
      title = paste("Simulated Cases across 32 Scenarios (100 Simulations) -", disease_name),
      subtitle = "Blue lines: Median cases across k1 & k2 | Red dashed line: Expected baseline (mu)",
      x = "Time (Weeks)",
      y = "Simulated Cases (Median)"
    )
  
  ggsave(file_name, plot = p, width = 16, height = 20, dpi = 300)
  return(p)
}

# Run function to export 2 plots
plot_sim_data(dengue_sim, "Dengue", "Plot1_Dengue_100Sims.png")
plot_sim_data(hfmd_sim, "HFMD", "Plot1_HFMD_100Sims.png")


# ==========================================
# 3. PLOT 2: COMBINED PERFORMANCE (FPR vs Sensitivity - Q1 & Q3)
# ==========================================
# Gom data 2 bệnh
combined_perf <- bind_rows(dengue_perf, hfmd_perf)

# Tính toán Median, Q1 (25%), Q3 (75%)
p2_data <- combined_perf %>%
  mutate(fpr = 1 - spec) %>%
  group_by(Disease, method) %>%
  summarise(
    med_sens = median(sens, na.rm = TRUE),
    q1_sens  = quantile(sens, 0.25, na.rm = TRUE),
    q3_sens  = quantile(sens, 0.75, na.rm = TRUE),
    
    med_fpr  = median(fpr, na.rm = TRUE),
    q1_fpr   = quantile(fpr, 0.25, na.rm = TRUE),
    q3_fpr   = quantile(fpr, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # Phân loại thuật toán bằng str_detect
  mutate(
    Family = case_when(
      str_detect(method, "farrington") ~ "Farrington",
      str_detect(method, "ears") ~ "EARS",
      str_detect(method, "cusum") ~ "CUSUM",
      str_detect(method, "bayes") ~ "Bayesian",
      str_detect(method, "rki") ~ "RKI",
      str_detect(method, "cdc") ~ "CDC",
      str_detect(method, "glr") ~ "GLR",
      TRUE ~ "Other"
    ),
    Category = case_when(
      Family %in% c("CUSUM", "EARS") ~ "Short-term baseline",
      TRUE ~ "Long-term baseline"
    ),
    # Điều kiện để hiện tên (Sens > 0.4 để hiện thêm 1 số thuật toán tạm ổn)
    is_good = med_sens > 0.4 & med_fpr < 0.35
  )

# Vẽ biểu đồ 2
plot2 <- ggplot(p2_data, aes(x = med_fpr, y = med_sens, color = Family, shape = Category)) +
  
  # Thanh ngang và thanh dọc (Map color theo Family)
  geom_errorbarh(aes(xmin = q1_fpr, xmax = q3_fpr), height = 0.02, alpha = 0.5, size = 0.6) +
  geom_errorbar(aes(ymin = q1_sens, ymax = q3_sens), width = 0.02, alpha = 0.5, size = 0.6) +
  
  # Điểm dot (Sẽ tự động nhận color và shape từ aes tổng)
  geom_point(size = 4, alpha = 0.9) +
  
  # Đường chéo y = x
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  
  # # Hiển thị TÊN thuật toán TỐT (cùng màu với Family)
  # geom_text_repel(
  #   data = filter(p2_data, is_good),
  #   aes(label = method),
  #   fontface = "bold",
  #   size = 4.5,
  #   box.padding = 0.8,
  #   point.padding = 0.4,
  #   show.legend = FALSE # Không cần hiện chữ a nhỏ trong legend
  # ) +
  
  # Facet chia 2 bệnh
  facet_wrap(~ Disease) + 
  theme_bw() +
  labs(
    title = "",
    x = "False Positive Rate (FPR)",
    y = "Sensitivity",
    color = "Algorithm Family",
    shape = ""
  ) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  
  # Bảng màu custom cho đẹp mắt (Set 7 màu phân biệt rõ ràng)
  scale_color_brewer(palette = "Dark2") +
  
  # Đổi shape: Tròn cho Short-term, Vuông cho Long-term
  scale_shape_manual(values = c("Long-term baseline" = 15, 
                                "Short-term baseline" = 16)) +
  
  # Làm đẹp theme
  theme(
    strip.text.x.top = element_text(size = 14, face = "bold", color = "black"),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, face = "bold", color = "black"),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_text(face = "bold")
  )

# Lưu ảnh
ggsave("Plot2_Combined_Sens_vs_FPR_Q1Q3 (1).png", plot = plot2, width = 12, height = 6, dpi = 300)


# ==========================================
# 3. TẠO BẢNG KẾT QUẢ CHO HÌNH NÀY (Median, Q1-Q3)
# ==========================================
# Ráp lại thành format: Median (Q1 - Q3)
table_plot2 <- p2_data %>%
  mutate(
    Sensitivity_Med_IQR = paste0(round(med_sens, 3), " (", round(q1_sens, 3), " - ", round(q3_sens, 3), ")"),
    FPR_Med_IQR         = paste0(round(med_fpr, 3), " (", round(q1_fpr, 3), " - ", round(q3_fpr, 3), ")")
  ) %>%
  # Chỉ lấy các cột cần thiết cho bảng báo cáo
  select(Disease, method, Sensitivity_Med_IQR, FPR_Med_IQR) %>%
  arrange(Disease, method)

# Xuất ra CSV
write.csv(table_plot2, "SummaryTable_Plot2_IQR.csv", row.names = FALSE)

# ==========================================
# 4. EXPORT SUMMARY TABLES (Median + Min/Max)
# ==========================================
# Function to format numbers as "Median (Min - Max)"
format_med_range <- function(x) {
  med <- round(median(x, na.rm = TRUE), 3)
  min_val <- round(min(x, na.rm = TRUE), 3)
  max_val <- round(max(x, na.rm = TRUE), 3)
  paste0(med, " (", min_val, " - ", max_val, ")")
}

# Function to generate and save table
generate_summary_table <- function(df, file_name) {
  tbl <- df %>%
    mutate(
      fpr = 1 - spec,
      ppv = ifelse((TP + FP) > 0, TP / (TP + FP), NA_real_),
      npv = ifelse((TN + FN) > 0, TN / (TN + FN), NA_real_)
    ) %>%
    mutate(scenario = factor(scenario, levels = paste("Sce", 1:32))) %>%
    group_by(scenario, method) %>%
    summarise(
      Sensitivity = format_med_range(sens),
      Specificity = format_med_range(spec),
      FPR = format_med_range(fpr),
      PPV = format_med_range(ppv),
      NPV = format_med_range(npv),
      POD1W = format_med_range(pod1w),
      .groups = "drop"
    ) %>%
    arrange(scenario, method)
  
  write.csv(tbl, file_name, row.names = FALSE)
}

# Export 2 report tables
generate_summary_table(dengue_perf, "SummaryTable_Dengue.csv")
generate_summary_table(hfmd_perf, "SummaryTable_HFMD.csv")

################ Re-evaluate performance by using framework ################
############################################################################
run_algos_for_one_sim <- function(df) {
  # 1. Tạo object sts, SAU ĐÓ BẮT BUỘC ÉP SANG disProg để tránh lỗi S4 class
  sts_obj <- sts(observed = df$n, state = df$outbreak, start = c(df$year[1], df$week[1]), frequency = 52)
  disProg_obj <- sts2disProg(sts_obj)
  
  # 2. Khai báo vùng đánh giá (t_current từ 313 đến 364 theo đúng QMD của mày)
  eval_range <- 313:364
  
  # 3. Chạy CUSUM và Farrington bằng disProg_obj (Các tham số này lấy chuẩn từ file của mày)
  ctrl_cusum <- list(range = eval_range, k = 1.04, h = 2.26, m = NULL, trans = "standard")
  res_cusum <- algo.cusum(disProg_obj, control = ctrl_cusum)
  
  ctrl_far <- list(range = eval_range, b = 5, w = 1, reweight = TRUE, verbose = FALSE, alpha = 0.05, trend = TRUE)
  res_far <- algo.farrington(disProg_obj, control = ctrl_far)
  
  # 4. Tạo cột alarm rỗng cho toàn bộ 364 tuần
  df$alarm_cusum <- 0
  df$alarm_farrington <- 0
  
  # 5. Lấp kết quả (52 tuần) vào đúng vị trí t_current
  df$alarm_cusum[eval_range] <- as.numeric(res_cusum$alarm)
  df$alarm_farrington[eval_range] <- as.numeric(res_far$alarm)
  
  # Xóa NA nếu có
  df$alarm_cusum[is.na(df$alarm_cusum)] <- 0
  df$alarm_farrington[is.na(df$alarm_farrington)] <- 0
  
  return(df)
}

# ==========================================
# 2. HÀM XỬ LÝ TOÀN BỘ DATA VÀ TÍNH MEDIAN
# ==========================================
evaluate_framework <- function(file_path, disease_name) {
  print(paste("Đang cày cuốc dữ liệu cho bệnh:", disease_name, "- Chờ xíu nha..."))
  
  sim_data <- readRDS(file_path)
  
  # Chạy thuật toán song song bằng furrr
  sim_with_alarms <- sim_data %>%
    group_by(scenario, k1, k2, sim) %>%
    nest() %>%
    mutate(data_with_alarms = future_map(data, run_algos_for_one_sim, .progress = TRUE)) %>%
    select(-data) %>%
    unnest(cols = c(data_with_alarms)) %>%
    ungroup()
  
  # Áp dụng logic Phễu lọc tuần tự và tính TP, FP cho TỪNG LẦN SIM
  sim_performance <- sim_with_alarms %>%
    # CỰC KỲ QUAN TRỌNG: Chỉ đánh giá trên đúng vùng t_current đã chạy thuật toán
    filter(t >= 313 & t <= 364) %>%
    
    mutate(
      Framework_Signal = case_when(
        alarm_cusum == 0 ~ "Routine",
        alarm_cusum == 1 & alarm_farrington == 0 ~ "Monitor",
        alarm_cusum == 1 & alarm_farrington == 1 ~ "Warning"
      )
    ) %>%
    group_by(scenario, k1, k2, sim) %>%
    summarise(
      Actual_Outbreak_Weeks = sum(outbreak == 1, na.rm = TRUE),
      Actual_Normal_Weeks   = sum(outbreak == 0, na.rm = TRUE),
      
      TP_Warning = sum(Framework_Signal == "Warning" & outbreak == 1, na.rm = TRUE),
      FP_Warning = sum(Framework_Signal == "Warning" & outbreak == 0, na.rm = TRUE),
      
      TP_Early = sum(Framework_Signal %in% c("Warning", "Monitor") & outbreak == 1, na.rm = TRUE),
      FP_Early = sum(Framework_Signal %in% c("Warning", "Monitor") & outbreak == 0, na.rm = TRUE),
      
      .groups = "drop"
    ) %>%
    mutate(
      Sens_Warning_sim = if_else(Actual_Outbreak_Weeks > 0, TP_Warning / Actual_Outbreak_Weeks, NA_real_),
      FPR_Warning_sim  = if_else(Actual_Normal_Weeks > 0, FP_Warning / Actual_Normal_Weeks, NA_real_),
      
      Sens_Early_sim   = if_else(Actual_Outbreak_Weeks > 0, TP_Early / Actual_Outbreak_Weeks, NA_real_),
      FPR_Early_sim    = if_else(Actual_Normal_Weeks > 0, FP_Early / Actual_Normal_Weeks, NA_real_)
    )
  
  # TÍNH MEDIAN (TRUNG VỊ) CHO TỪNG KỊCH BẢN
  final_summary <- sim_performance %>%
    group_by(scenario, k1, k2) %>%
    summarise(
      Disease = disease_name,
      Sens_Warning_Med = median(Sens_Warning_sim, na.rm = TRUE),
      FPR_Warning_Med  = median(FPR_Warning_sim, na.rm = TRUE),
      Sens_Early_Med   = median(Sens_Early_sim, na.rm = TRUE),
      FPR_Early_Med    = median(FPR_Early_sim, na.rm = TRUE),
      .groups = "drop"
    )
  
  return(final_summary)
}

# ==========================================
# 3. BẢNG TỔNG HỢP
# ==========================================
library(furr)
# Bật chế độ chạy song song đa nhân
plan(multisession, workers = availableCores() - 1)

# Chạy hệ thống cho Dengue
dengue_framework_res <- evaluate_framework("Dengue_Cases.rds", "Dengue")
saveRDS(dengue_framework_res, "dengue_framework_res.rds")

# Chạy hệ thống cho HFMD
hfmd_framework_res <- evaluate_framework("HFMD_Cases.rds", "HFMD")
saveRDS(hfmd_framework_res, "hfmd_framework_res.rds")
# Trả CPU về trạng thái bình thường
plan(sequential)

# Gộp kết quả 2 bệnh lại và làm tròn 3 chữ số thập phân cho đẹp báo cáo
final_framework_table <- bind_rows(dengue_framework_res, hfmd_framework_res) %>%
  mutate(
    Sens_Warning_Med = round(Sens_Warning_Med, 3),
    FPR_Warning_Med  = round(FPR_Warning_Med, 3),
    Sens_Early_Med   = round(Sens_Early_Med, 3),
    FPR_Early_Med    = round(FPR_Early_Med, 3)
  ) %>%
  # Sắp xếp lại thứ tự cột cho dễ nhìn
  select(Disease, scenario, k1, k2, Sens_Early_Med, FPR_Early_Med, Sens_Warning_Med, FPR_Warning_Med)


plot_data <- final_framework_table %>%
  pivot_longer(
    cols = c(Sens_Early_Med, FPR_Early_Med, Sens_Warning_Med, FPR_Warning_Med),
    names_to = "variable",
    values_to = "value"
  ) %>%
  mutate(
    # Tách loại chỉ số (Sens vs FPR)
    Metric = if_else(str_detect(variable, "Sens"), "Sensitivity", "False Positive Rate (FPR)"),
    # Tách cấp độ cảnh báo
    Tier = if_else(str_detect(variable, "Early"), "Close Monitoring", "Warning")
  ) %>%
  # Sắp xếp thứ tự cho chuẩn logic
  mutate(
    Tier = factor(Tier, levels = c("Close Monitoring", "Warning")),
    Disease = factor(Disease, levels = c("Dengue", "HFMD"))
  )

p_sens <- plot_data %>%
  filter(Metric == "Sensitivity", !is.na(value)) %>%
  ggplot(aes(x = Tier, y = value, fill = Tier)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.5) +
  geom_jitter(color = "black", alpha = 0.15, width = 0.15, size = 1.2) +
  facet_wrap(~ Disease) +
  scale_fill_manual(values = c("Close Monitoring" = "#f39c12", "Warning" = "red")) +
  theme_bw() +
  labs(
    title = "(A) Sensitivity Distribution",
    subtitle = "Across 288 simulated parameter combinations",
    x = NULL, 
    y = "Sensitivity",
    fill = "" # Thêm tên cho Legend
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  theme(
    title = element_text(face = "bold", size = 14),
    axis.text.x = element_blank(),  # Xóa text trục X
    axis.ticks.x = element_blank(), # Xóa luôn dấu gạch (tick) của trục X
    strip.text = element_text(face = "bold", size = 12)
    # Lưu ý: Đã bỏ lệnh legend.position = "none" ở đây để lát patchwork thu gom
  )

# 3. Vẽ biểu đồ Tỷ lệ hoang báo (FPR)
p_fpr <- plot_data %>%
  filter(Metric == "False Positive Rate (FPR)", !is.na(value)) %>%
  ggplot(aes(x = Tier, y = value, fill = Tier)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.5) +
  geom_jitter(color = "black", alpha = 0.15, width = 0.15, size = 1.2) +
  facet_wrap(~ Disease) +
  scale_fill_manual(values = c("Close Monitoring" = "#f39c12", "Warning" = "red")) +
  theme_bw() +
  labs(
    title = "(B) False Positive Rate Distribution",
    subtitle = "Across 288 simulated parameter combinations",
    x = NULL, 
    y = "False Positive Rate (FPR)",
    fill = "" # Thêm tên cho Legend
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  theme(
    title = element_text(face = "bold", size = 14),
    axis.text.x = element_blank(),  # Xóa text trục X
    axis.ticks.x = element_blank(), # Xóa dấu tick
    strip.text = element_text(face = "bold", size = 12)
  )

# 4. Ghép 2 biểu đồ và gom chung Legend
final_main_plot <- (p_sens | p_fpr) +
  # plot_layout(guides = "collect") sẽ tự động gộp Legend giống nhau của 2 hình làm 1
  plot_layout(guides = "collect") & 
  # Set vị trí Legend chung cho toàn bộ bức hình
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 11)
  )

ggsave("MainText_Framework_Performance_Clean.png", plot = final_main_plot, width = 14, height = 6, dpi = 300)

##### ----------------- ######

p5_data <- combined_perf %>%
  mutate(
    fpr = 1 - spec,
    
    Baseline_State = if_else(k1 == 0, "Baseline: No Outbreak (k1=0)", "Baseline: Outbreak (k1>0)"),
    Current_State  = if_else(k2 == 0, "Current: No Outbreak (k2=0)", "Current: Outbreak (k2>0)"),
    
    Algorithm = case_when(
      str_detect(method, "bayes.*0") ~ "Bayes1",
      str_detect(method, "bayes.*1") ~ "Bayes2",
      str_detect(method, "cdc") ~ "CDC",
      str_detect(method, "pearsonNegBin") ~ "CUSUM NegBin",
      str_detect(method, "rossi") ~ "CUSUM Poisson",
      str_detect(method, "standard") ~ "CUSUM",
      str_detect(method, "earsC1") ~ "EARS C1",
      str_detect(method, "earsC2") ~ "EARS C2",
      str_detect(method, "earsC3") ~ "EARS C3",
      str_detect(method, "farringtonf") ~ "Farrington Flexible",
      str_detect(method, "farrington") ~ "Farrington", 
      str_detect(method, "glrnb") ~ "GLR NegBin",
      str_detect(method, "glrpois") ~ "GLR Poisson",
      str_detect(method, "rki.*0") ~ "RKI1",
      str_detect(method, "rki.*1") ~ "RKI2",
      TRUE ~ method
    )
  )

# Đưa CDC lên đầu tiên trong levels
p5_data$Algorithm <- factor(
  p5_data$Algorithm, 
  levels = c("CDC", "Bayes1", "Bayes2", "CUSUM NegBin", "CUSUM Poisson", "CUSUM", 
             "EARS C1", "EARS C2", "EARS C3", "Farrington", "Farrington Flexible", 
             "GLR NegBin", "GLR Poisson", "RKI1", "RKI2")
)

# Ép kiểu factor cho Baseline và Current
p5_data$Baseline_State <- factor(p5_data$Baseline_State, levels = c("Baseline: No Outbreak (k1=0)", "Baseline: Outbreak (k1>0)"))
p5_data$Current_State  <- factor(p5_data$Current_State, levels = c("Current: No Outbreak (k2=0)", "Current: Outbreak (k2>0)"))

# =================
# Biểu đồ A: Ma trận 2x2 cho False Positive Rate (FPR)
# =================
plot_A_FPR <- ggplot(p5_data, aes(x = Algorithm, y = fpr, fill = Algorithm)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5, outlier.alpha = 0.3) +
  facet_grid(Baseline_State ~ Current_State) + 
  theme_bw() +
  labs(
    title = "(A) False Positive Rate (FPR)",
    subtitle = "Evaluating algorithm noise with and without historical/current outbreaks",
    x = NULL,
    y = "False Positive Rate (FPR)"
  ) +
  coord_cartesian(ylim = c(0, 0.5)) + 
  theme(
    # Điểm mấu chốt: Truyền vector màu (Đỏ cho vị trí 1, Đen cho 14 vị trí còn lại)
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 10, 
                               color = c("red", rep("black", 14))),
    axis.text.y = element_text(color = "black", size = 11),
    axis.title.y = element_text(color = "black", face = "bold", size = 11),
    strip.text = element_text(face = "bold", size = 11, background = element_rect(fill = "gray90")),
    strip.background = element_rect(fill = "lightblue"),
    legend.position = "none" 
  )

# =================
# Biểu đồ B: Sensitivity
# =================
plot_B_Sens <- p5_data %>%
  filter(k2 > 0) %>% 
  ggplot(aes(x = Algorithm, y = sens, fill = Algorithm)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5, outlier.alpha = 0.3) +
  facet_wrap(~ Baseline_State, ncol = 2) + 
  theme_bw() +
  labs(
    title = "(B) Sensitivity",
    subtitle = "When Current Outbreak is Present (k2 > 0)",
    x = "",
    y = "Sensitivity"
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  theme(
    # Truyền vector màu y chang biểu đồ trên
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 9, 
                               color = c("red", rep("black", 14))),
    axis.text.y = element_text(color = "black", size = 11),
    axis.title.y = element_text(color = "black", face = "bold", size = 11),
    strip.text = element_text(face = "bold", size = 11, background = element_rect(fill = "#dff9fb")),
    strip.background = element_rect(fill = "lightblue"),
    legend.position = "none"
  )

# Ghép đồ thị
final_4cases_plot <- plot_A_FPR / plot_B_Sens + 
  plot_layout(heights = c(1.5, 1))

ggsave("Plot5_4_Scenarios_Evaluation_All_Methods.png", plot = final_4cases_plot, width = 14, height = 12, dpi = 300)

