library(tidyverse)
library(ggrepel)

# ==========================================
# 1. LOAD DATA
# ==========================================
dengue_sim <- readRDS("data/Dengue_Cases.rds")
dengue_perf <- read.csv("data/Dengue_Perf.csv") %>% mutate(Disease = "Dengue")

hfmd_sim <- readRDS("data/HFMD_Cases.rds")
hfmd_perf <- read.csv("data/HFMD_Perf.csv") %>% mutate(Disease = "HFMD")

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
# CHUẨN BỊ DỮ LIỆU & ĐỒNG BỘ NHÓM THUẬT TOÁN (FAMILY)
# ==========================================
p5_data <- combined_perf %>%
  mutate(
    fpr = 1 - spec,
    
    Baseline_State = if_else(k1 == 0, "Baseline: No Outbreak (k1=0)", "Baseline: Outbreak (k1>0)"),
    Current_State  = if_else(k2 == 0, "Current: No Outbreak (k2=0)", "Current: Outbreak (k2>0)"),
    
    # SỬA TẠI ĐÂY: Dùng "^" (bắt đầu chuỗi) để phân nhóm chính xác tuyệt đối
    Family = case_when(
      str_detect(method, "^farrington") ~ "Farrington",
      str_detect(method, "^cusum") ~ "CUSUM",
      str_detect(method, "^ears") ~ "EARS",
      str_detect(method, "^bayes") ~ "Bayesian",
      str_detect(method, "^rki") ~ "RKI",
      str_detect(method, "^cdc") ~ "CDC",
      str_detect(method, "^glr") ~ "GLR",
      TRUE ~ "Other"
    ),
    
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

# Ép kiểu factor cho Baseline và Current
p5_data$Baseline_State <- factor(p5_data$Baseline_State, levels = c("Baseline: No Outbreak (k1=0)", "Baseline: Outbreak (k1>0)"))
p5_data$Current_State  <- factor(p5_data$Current_State, levels = c("Current: No Outbreak (k2=0)", "Current: Outbreak (k2>0)"))

# ==========================================
# TÍNH THỨ TỰ SẮP XẾP TRỤC X CHO BIỂU ĐỒ A & B
# ==========================================

# 1. Tính thứ tự FPR riêng biệt cho từng cột (k2=0 riêng, k2>0 riêng)
order_fpr_by_col <- p5_data %>%
  group_by(Current_State, Algorithm) %>%
  summarise(
    med = median(fpr, na.rm = TRUE),
    q3  = quantile(fpr, 0.75, na.rm = TRUE),
    avg = mean(fpr, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # Nếu med và q3 cùng bằng 0 (như CDC và GLR NegBin), ưu tiên giữ CDC ở vị trí đầu tiên
  arrange(Current_State, med, q3, Algorithm != "CDC", avg) %>%
  mutate(Alg_Col = paste(Algorithm, Current_State, sep = "___")) %>%
  pull(Alg_Col)

# Thứ tự cho Plot B (Sensitivity): Cao hơn là tốt hơn -> Xếp từ CAO tới THẤP
order_sens <- p5_data %>%
  filter(k2 > 0) %>%
  group_by(Algorithm) %>%
  summarise(
    med = median(sens, na.rm = TRUE),
    q3  = quantile(sens, 0.75, na.rm = TRUE),
    avg = mean(sens, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(med), desc(q3), desc(avg)) %>%
  pull(Algorithm)

# Tự động định vị chữ "CDC" để tô màu đỏ trên trục X của từng hình
axis_color_A <- if_else(order_fpr %in% c("CDC", "GLR NegBin"), "red", "black")
axis_color_B <- if_else(order_sens %in% c("CDC", "GLR NegBin"), "red", "black")

# ==========================================
# BIỂU ĐỒ A: Ma trận 2x2 cho False Positive Rate (FPR - Thấp đến Cao)
# ==========================================
plot_A_FPR <- p5_data %>%
  mutate(
    Alg_Col = factor(paste(Algorithm, Current_State, sep = "___"), levels = order_fpr_by_col)
  ) %>%
  ggplot(aes(x = Alg_Col, y = fpr, fill = Family)) +
  geom_boxplot(alpha = 0.8, outlier.size = 0.5, outlier.alpha = 0.3) +
  # Bật scales = "free_x" để 2 cột tự xếp thứ tự độc lập
  facet_grid(Baseline_State ~ Current_State, scales = "free_x") + 
  # Cắt bỏ phần hậu tố phụ để hiển thị tên thuật toán sạch đẹp
  scale_x_discrete(labels = function(x) str_remove(x, "___.*")) +
  scale_fill_brewer(palette = "Dark2", name = "Algorithm Family") +
  theme_bw() +
  labs(
    title = "(A) False Positive Rate (FPR)",
    subtitle = "Ordered from lowest to highest FPR within each current-period condition (lower is better)",
    x = NULL,
    y = "False Positive Rate (FPR)"
  ) +
  coord_cartesian(ylim = c(0, 0.5)) + 
  theme(
    # CDC luôn đứng số 1 ở cả 2 cột nên giữ nguyên màu đỏ ở vị trí đầu
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 10, 
                               color = axis_color_A),
    axis.text.y = element_text(color = "black", size = 11),
    axis.title.y = element_text(color = "black", face = "bold", size = 11),
    strip.text = element_text(face = "bold", size = 11),
    strip.background = element_rect(fill = "lightblue")
  )

# ==========================================
# BIỂU ĐỒ B: Sensitivity (Cao đến Thấp)
# ==========================================
plot_B_Sens <- p5_data %>%
  filter(k2 > 0) %>% 
  mutate(Algorithm = factor(Algorithm, levels = order_sens)) %>%
  ggplot(aes(x = Algorithm, y = sens, fill = Family)) +
  geom_boxplot(alpha = 0.8, outlier.size = 0.5, outlier.alpha = 0.3) +
  facet_wrap(~ Baseline_State, ncol = 2) + 
  scale_fill_brewer(palette = "Dark2", name = "Algorithm Family") +
  theme_bw() +
  labs(
    title = "(B) Sensitivity",
    subtitle = "When Current Outbreak is Present (k2 > 0), ordered from highest to lowest (higher is better)",
    x = "",
    y = "Sensitivity"
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 10, 
                               color = axis_color_B),
    axis.text.y = element_text(color = "black", size = 11),
    axis.title.y = element_text(color = "black", face = "bold", size = 11),
    strip.text = element_text(face = "bold", size = 11),
    strip.background = element_rect(fill = "lightblue")
  )

# ==========================================
# GHÉP ĐỒ THỊ VÀ GOM LEGEND CHUNG Ở ĐÁY
# ==========================================
final_4cases_plot <- (plot_A_FPR / plot_B_Sens) + 
  plot_layout(heights = c(1.5, 1), guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10)
  )

ggsave("Plot5_4_Scenarios_Evaluation_All_Methods.png", plot = final_4cases_plot, width = 10.5, height = 12.5, dpi = 300)

