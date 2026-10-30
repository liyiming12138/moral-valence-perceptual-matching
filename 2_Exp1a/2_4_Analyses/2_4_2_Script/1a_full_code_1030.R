# ============================================================
# 1. 环境设置
# ============================================================


# 本文件包含RT和SDT的完整模型、零模型、诊断、后验差值、图和BF。
# 可按小节运行。首次运行brm()且没有相应缓存时，会进行MCMC采样。
# 已有且设置一致的本文件模型缓存，会由brm()直接读取。
# R启动临时目录已通过.Renviron设置在D盘；tempdir()可查看当前目录。

# 加载分析和绘图所需的包。
Sys.setenv(CMDSTAN = "D:/Stan_BF/tools/cmdstan-2.40.0")

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggdist)
library(patchwork)
library(brms)
library(cmdstanr)
library(posterior)
library(bayesplot)
library(bridgesampling)

# 设置项目根目录。
# 先打开项目根目录的moral valence.Rproj，再运行本脚本。
curDir <- getwd()

# 设置数据目录。
data_dir <- paste0(
  "D:/Users/15/Desktop/道德自我/SalientGoodSelf-master/",
  "SalientGoodSelf-master/exp1a"
)

# 设置结果保存路径，并创建文件夹。
result_dir <- file.path(curDir, "2_Exp1a", "2_4_Analyses", "2_4_3_Results", "full_code_1029")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

# 设置模型保存路径，并创建文件夹。
# 新模型缓存保存在本实验分析目录；此目录尚无模型时，运行brm()会拟合模型。
model_dir <- file.path(curDir, "2_Exp1a", "2_4_Analyses", "2_4_5_ModelCache")
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)

# 设置图片保存路径，并创建文件夹。
fig_dir <- file.path(curDir, "2_Exp1a", "2_4_Analyses", "2_4_4_Figures", "full_code_1029")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# 设置Stan源码保存路径；英文路径便于Windows编译。
stan_dir <- "D:/Stan_BF/exp1a_1028/stan_code"
dir.create(stan_dir, recursive = TRUE, showWarnings = FALSE)

# 设置原始采样链CSV保存路径；每个模型使用不同文件名前缀。
chain_dir <- "D:/Stan_BF/exp1a_1028/chains"
dir.create(chain_dir, recursive = TRUE, showWarnings = FALSE)

# 设置随机种子。
set.seed(20260927)

# 设置CmdStan路径及brms采样后端。
cmdstanr::set_cmdstan_path("D:/Stan_BF/tools/cmdstan-2.40.0")

options(
  brms.backend = "cmdstanr",
  brms.chains = 4,
  brms.cores = 4,
  cmdstanr_write_stan_file_dir = stan_dir
)



# ============================================================
# 2. 数据读取与预处理：import data / clear data
# ============================================================


# 2.1 读取两批原始数据。
# 先读为字符，合并后再转换需要的变量。
data_2014 <- read_csv(
  file.path(data_dir, "rawdata_behav_exp1a_201404_export_2019.csv"),
  col_types = cols(.default = col_character()),
  show_col_types = FALSE
) %>%
  mutate(cohort = "2014", source_row = row_number() + 1)

data_2017 <- read_csv(
  file.path(data_dir, "rawdata_behav_exp1a_201704_export_2019.csv"),
  col_types = cols(.default = col_character()),
  show_col_types = FALSE
) %>%
  mutate(cohort = "2017", source_row = row_number() + 1)

raw_data <- bind_rows(data_2014, data_2017)

# 2.2 统一变量名；被试编号前加批次，区分两批同号被试。
data <- raw_data %>%
  transmute(
    cohort,
    source_row,
    subject_uid = paste(cohort, trimws(Subject), sep = "_"),
    block = trimws(`BlockList.Sample`),
    accuracy = as.integer(`Target.ACC`),
    rt_ms = as.numeric(`Target.RT`),
    valence = recode(Shape, "Normal" = "Neutral"),
    matchness = recode(YesNoResp, "Yes" = "Match", "No" = "Mismatch")
  )

# 2.3 排除练习试次；原数据中练习试次的block为空。
formal_data <- data %>%
  filter(!is.na(block), block != "")

# 2.4 先按全部正式试次计算正确率，保留总体正确率不低于60%的被试。
subject_accuracy <- formal_data %>%
  group_by(subject_uid) %>%
  summarise(
    overall_accuracy = mean(accuracy),
    .groups = "drop"
  )

included_subjects <- subject_accuracy %>%
  filter(overall_accuracy >= 0.60)

# 2.5 再保留RT非缺失、200至1500 ms范围内的试次。
valid_data <- formal_data %>%
  semi_join(included_subjects, by = "subject_uid") %>%
  filter(!is.na(rt_ms), rt_ms >= 200, rt_ms <= 1500)

# 2.6 RT分析只使用正确试次；毫秒转换为秒。
rt_data <- valid_data %>%
  filter(accuracy == 1) %>%
  mutate(rt_s = rt_ms / 1000)

# 2.7 SDT分析保留正确和错误试次。
# signal：实际匹配=1，实际不匹配=0。
# response_match：回答匹配=1，回答不匹配=0。
sdt_data <- valid_data %>%
  mutate(
    signal = as.integer(matchness == "Match"),
    response_match = if_else(signal == 1, accuracy, 1L - accuracy)
  )

# 2.8 查看样本量与试次数；年龄缺失不作为行为排除条件。
sample_size <- data.frame(
  recruited = nrow(subject_accuracy),
  included = nrow(included_subjects),
  RT_trials = nrow(rt_data),
  SDT_trials = nrow(sdt_data)
)

print(sample_size)

# 2.9 RT采用中心化对比编码。
rt_data$subject_uid <- factor(rt_data$subject_uid)
rt_data$valence <- factor(rt_data$valence, levels = c("Good", "Neutral", "Bad"))
rt_data$matchness <- factor(rt_data$matchness, levels = c("Match", "Mismatch"))

contrasts(rt_data$valence) <- matrix(
  c(2 / 3, -1 / 3, -1 / 3, -1 / 3, -1 / 3, 2 / 3),
  nrow = 3,
  dimnames = list(c("Good", "Neutral", "Bad"), c("Good", "Bad"))
)

contrasts(rt_data$matchness) <- matrix(
  c(-0.5, 0.5),
  ncol = 1,
  dimnames = list(c("Match", "Mismatch"), "Mismatch")
)

# 2.10 SDT以Neutral为参考水平；signal保持0/1。
sdt_data$subject_uid <- factor(sdt_data$subject_uid)
sdt_data$valence <- factor(sdt_data$valence, levels = c("Neutral", "Good", "Bad"))
contrasts(sdt_data$valence) <- contr.treatment(c("Neutral", "Good", "Bad"), base = 1)

# 2.11 保存筛选记录与建模数据。
write_csv(subject_accuracy, file.path(result_dir, "01_被试总体正确率.csv"))
write_csv(included_subjects, file.path(result_dir, "02_纳入被试.csv"))
write_csv(rt_data, file.path(result_dir, "03_RT建模数据.csv"))
write_csv(sdt_data, file.path(result_dir, "04_SDT建模数据.csv"))
write_csv(sample_size, file.path(result_dir, "05_样本与试次数.csv"))


# 描述统计与雨云图。
rt_subject <- rt_data %>%
  group_by(subject_uid, valence, matchness) %>%
  summarise(value = mean(rt_ms), .groups = "drop")

# 传统dprime只用于描述统计，不作为SDT模型的因变量。
# 极端比例校正沿用原描述统计代码，保持原图数值一致。
# 2.12 按被试、效价汇总命中和虚报。
dprime_subject <- sdt_data %>%
  group_by(subject_uid, valence) %>%
  summarise(
    signal_trials = sum(signal == 1),
    noise_trials = sum(signal == 0),
    hits = sum(signal == 1 & response_match == 1),
    false_alarms = sum(signal == 0 & response_match == 1),
    .groups = "drop"
  ) %>%
  mutate(
    hit_rate = hits / signal_trials,
    false_alarm_rate = false_alarms / noise_trials,
    hit_rate = if_else(hit_rate == 1, 1 - 1 / (2 * signal_trials), hit_rate),
    false_alarm_rate = if_else(false_alarm_rate == 0, 1 / (2 * signal_trials), false_alarm_rate),
    dprime = qnorm(hit_rate) - qnorm(false_alarm_rate)
  )

plot_data <- bind_rows(
  rt_subject %>% transmute(subject_uid, valence, panel = as.character(matchness), value),
  dprime_subject %>% transmute(subject_uid, valence, panel = "dprime", value = dprime)
) %>% mutate(valence = factor(valence, levels = c("Good", "Neutral", "Bad")),
             x = as.integer(valence), x_point = x - 0.23 + runif(n(), -0.06, 0.06))

plot_summary <- plot_data %>% group_by(panel, valence) %>%
  summarise(n = n(), mean = mean(value), sd = sd(value), se = sd / sqrt(n), .groups = "drop")
write_csv(rt_subject, file.path(result_dir, "06_RT被试条件均值.csv"))
write_csv(dprime_subject, file.path(result_dir, "07_传统dprime被试条件值.csv"))
write_csv(plot_summary, file.path(result_dir, "08_描述统计.csv"))

# 各条件的核密度，用于雨云图右侧的半小提琴。
clouds <- plot_data %>% group_by(panel, valence) %>% group_modify(~ {
  den <- density(.x$value, n = 256)
  baseline <- as.integer(.y$valence) + 0.14
  tibble(x = c(rep(baseline, length(den$x)), rev(baseline + 0.42 * den$y / max(den$y))),
         value = c(den$x, rev(den$x)))
}) %>% ungroup()

valence_colours <- c(Good = "#1B9E77", Neutral = "#D95F02", Bad = "#7570B3")
figure_theme <- theme_classic(base_size = 12, base_family = "Times New Roman") +
  theme(axis.text = element_text(colour = "black"), legend.position = "none")

# 2.13 匹配条件下RT雨云图。
p_rt_raincloud <- ggplot(filter(plot_data, panel == "Match"), aes(y = value)) +
  geom_polygon(data = filter(clouds, panel == "Match"),
               aes(x = x, fill = valence, group = valence), alpha = 0.60) +
  geom_line(aes(x = x_point, group = subject_uid), colour = "grey70", alpha = 0.22) +
  geom_point(aes(x = x_point, colour = valence), alpha = 0.48, size = 1.4) +
  geom_boxplot(aes(x = x - 0.04, fill = valence, group = valence),
               width = 0.10, outlier.shape = NA, alpha = 0.65) +
  geom_errorbar(data = filter(plot_summary, panel == "Match"),
    aes(x = as.integer(valence) + 0.08, y = mean, ymin = mean - se, ymax = mean + se,
        colour = valence), width = 0.045) +
  geom_point(data = filter(plot_summary, panel == "Match"),
    aes(x = as.integer(valence) + 0.08, y = mean, colour = valence), size = 2.1) +
  scale_x_continuous(breaks = 1:3, labels = c("Good", "Neutral", "Bad"), limits = c(0.6, 3.65)) +
  scale_colour_manual(values = valence_colours) + scale_fill_manual(values = valence_colours) +
  labs(tag = "A", title = "Match", x = "Moral valence", y = "Reaction time (ms)") + figure_theme

# 2.14 传统dprime雨云图。
p_dprime_raincloud <- ggplot(filter(plot_data, panel == "dprime"), aes(y = value)) +
  geom_polygon(data = filter(clouds, panel == "dprime"),
               aes(x = x, fill = valence, group = valence), alpha = 0.60) +
  geom_line(aes(x = x_point, group = subject_uid), colour = "grey70", alpha = 0.22) +
  geom_point(aes(x = x_point, colour = valence), alpha = 0.48, size = 1.4) +
  geom_boxplot(aes(x = x - 0.04, fill = valence, group = valence),
               width = 0.10, outlier.shape = NA, alpha = 0.65) +
  geom_errorbar(data = filter(plot_summary, panel == "dprime"),
    aes(x = as.integer(valence) + 0.08, y = mean, ymin = mean - se, ymax = mean + se,
        colour = valence), width = 0.045) +
  geom_point(data = filter(plot_summary, panel == "dprime"),
    aes(x = as.integer(valence) + 0.08, y = mean, colour = valence), size = 2.1) +
  scale_x_continuous(breaks = 1:3, labels = c("Good", "Neutral", "Bad"), limits = c(0.6, 3.65)) +
  scale_colour_manual(values = valence_colours) + scale_fill_manual(values = valence_colours) +
  labs(tag = "B", title = "d prime", x = "Moral valence", y = expression(italic(d)*"′")) + figure_theme

# 2.15 将两张图纵向组合，在RStudio显示并保存。
p_raincloud <- p_rt_raincloud / p_dprime_raincloud
print(p_raincloud)
ggsave(file.path(fig_dir, "01_1a_描述统计雨云组合图.png"), p_raincloud, width = 8, height = 9, dpi = 300)
ggsave(file.path(fig_dir, "01_1a_描述统计雨云组合图.jpg"), p_raincloud, width = 8, height = 9, dpi = 300)
ggsave(file.path(fig_dir, "01_1a_描述统计雨云组合图.pdf"), p_raincloud, device = grDevices::cairo_pdf, width = 8, height = 9)


# ============================================================
# 3. RT贝叶斯混合效应模型
# ============================================================


# 3.1 完整模型H1：效价、匹配及交互；被试随机截距、随机斜率及相关。
formula_rt <- brms::bf(
  rt_s ~ valence * matchness + (1 + valence * matchness | subject_uid),
  center = FALSE
)

# 3.2 指定RT先验。
prior_rt <- c(
  prior(normal(0, 0.30), class = "b"),
  prior(normal(-0.35667494, 0.5), class = "b", coef = "Intercept"),
  prior(exponential(5), class = "sd"),
  prior(exponential(5), class = "sigma"),
  prior(lkj(2), class = "cor")
)

# 3.3 拟合RT完整模型；4链，每链4000次，其中500次预热。
# save_pars(all=TRUE)保存全部参数，供后续桥采样使用。
# file保存模型RDS；output_dir保存CmdStan采样链CSV。


model_rt <- brm(
  formula = formula_rt,
  data = rt_data,
  family = lognormal(),
  prior = prior_rt,
  backend = "cmdstanr",
  normalize = TRUE,
  chains = 4,
  iter = 4000,
  warmup = 500,
  cores = 4,
  seed = 20260927,
  refresh = 100,
  sig_figs = 18,
  control = list(adapt_delta = 0.95, max_treedepth = 15),
  save_pars = save_pars(all = TRUE),
  sample_prior = "yes",
  output_dir = chain_dir,
  output_basename = "1a_RT_full",
  file = file.path(model_dir, "1a_RT_full"),
  file_refit = "on_change"
)


# ============================================================
# 4. RT模型诊断
# ============================================================


# 4.1 绘制主要参数的MCMC轨迹图。
p_rt_trace <- mcmc_trace(
  as_draws_array(model_rt),
  pars = c(
    "b_Intercept", "b_valenceGood", "b_valenceBad", "b_matchnessMismatch",
    "b_valenceGood:matchnessMismatch", "b_valenceBad:matchnessMismatch",
    "sd_subject_uid__Intercept", "sigma"
  )
)

print(p_rt_trace)

ggsave(
  file.path(fig_dir, "RT_轨迹图.png"),
  p_rt_trace, width = 11, height = 8, dpi = 200
)

# 4.2 汇总全部模型参数的最大Rhat、最小Bulk ESS。
# 包含固定效应、随机效应、标准差、相关和sigma；排除单独的先验抽样。
rt_parameters <- grep("^(b_|sd_|cor_|r_|sigma$)", variables(model_rt), value = TRUE)
rt_draws <- as_draws_array(model_rt, variable = rt_parameters)

rt_diagnostics <- summarise_draws(
  rt_draws,
  rhat = posterior::rhat,
  ess_bulk = posterior::ess_bulk
)

# 提取发散次数与树深触顶次数。
rt_nuts <- nuts_params(model_rt)

rt_diagnostic_summary <- data.frame(
  max_Rhat = max(rt_diagnostics$rhat),
  min_Bulk_ESS = min(rt_diagnostics$ess_bulk),
  divergences = sum(rt_nuts$Value[rt_nuts$Parameter == "divergent__"]),
  treedepth_hits = sum(rt_nuts$Value[rt_nuts$Parameter == "treedepth__"] >= 15)
)

print(rt_diagnostic_summary)

write_csv(rt_diagnostics, file.path(result_dir, "RT_逐参数诊断.csv"))
write_csv(rt_diagnostic_summary, file.path(result_dir, "RT_诊断汇总.csv"))

# 4.3 后验预测检查：模拟RT分布与实际RT分布进行比较。
p_rt_ppc <- pp_check(model_rt, type = "dens_overlay", ndraws = 50) +
  labs(x = "Reaction time (s)")

print(p_rt_ppc)

ggsave(
  file.path(fig_dir, "RT_后验预测检查.png"),
  p_rt_ppc, width = 9, height = 6, dpi = 200
)


# ============================================================
# 5. RT回归系数、条件预测与后验差值
# ============================================================


# 5.1 回归系数beta：位于log-RT尺度。
# fixef()给出系数的后验均值，以及2.5%和97.5%分位数。
rt_beta <- as.data.frame(fixef(model_rt, probs = c(0.025, 0.975))) %>%
  tibble::rownames_to_column("parameter") %>%
  transmute(
    parameter,
    posterior_mean = Estimate,
    lower95 = Q2.5,
    upper95 = Q97.5,
    CI_contains_zero = lower95 <= 0 & upper95 >= 0
  )

print(rt_beta)
write_csv(rt_beta, file.path(result_dir, "RT_beta_95CI.csv"))

# 5.2 建立六个条件的预测数据。
# 前三列对应匹配的Good/Neutral/Bad，后三列对应不匹配的Good/Neutral/Bad。
rt_grid <- expand.grid(
  valence = c("Good", "Neutral", "Bad"),
  matchness = c("Match", "Mismatch")
)

rt_grid$subject_uid <- rt_data$subject_uid[1]
rt_grid$valence <- factor(rt_grid$valence, levels = levels(rt_data$valence))
rt_grid$matchness <- factor(rt_grid$matchness, levels = levels(rt_data$matchness))
contrasts(rt_grid$valence) <- contrasts(rt_data$valence)
contrasts(rt_grid$matchness) <- contrasts(rt_data$matchness)

# 5.3 由每组后验样本计算模型期望RT；乘1000转换为毫秒。
# re_formula=NA表示预测时把被试随机效应设为0。
# posterior_epred()包含lognormal的sigma校正，得到期望RT。
rt_predicted_ms <- posterior_epred(
  model_rt,
  newdata = rt_grid,
  re_formula = NA
) * 1000

rt_conditions <- cbind(
  rt_grid[, c("valence", "matchness")],
  posterior_mean_ms = colMeans(rt_predicted_ms),
  lower95 = apply(rt_predicted_ms, 2, quantile, 0.025),
  upper95 = apply(rt_predicted_ms, 2, quantile, 0.975)
)

print(rt_conditions)
write_csv(rt_conditions, file.path(result_dir, "RT_六条件后验期望_ms.csv"))

# 5.4 计算六组模型预测差值，单位为ms；这些差值不是回归系数beta。
# GN=积极减中性，GB=积极减消极，NB=中性减消极。
rt_difference_draws <- tibble(
  .draw = seq_len(nrow(rt_predicted_ms)),
  Match_GN = rt_predicted_ms[, 1] - rt_predicted_ms[, 2],
  Match_GB = rt_predicted_ms[, 1] - rt_predicted_ms[, 3],
  Match_NB = rt_predicted_ms[, 2] - rt_predicted_ms[, 3],
  Mismatch_GN = rt_predicted_ms[, 4] - rt_predicted_ms[, 5],
  Mismatch_GB = rt_predicted_ms[, 4] - rt_predicted_ms[, 6],
  Mismatch_NB = rt_predicted_ms[, 5] - rt_predicted_ms[, 6]
) %>%
  pivot_longer(
    cols = -.draw,
    names_to = c("matchness", "comparison"),
    names_sep = "_",
    values_to = "value"
  ) %>%
  mutate(
    comparison = factor(
      comparison,
      levels = c("GN", "GB", "NB"),
      labels = c("Good - Neutral", "Good - Bad", "Neutral - Bad")
    )
  )

# 5.5 汇总后验均值、等尾95%CI、差值为正的后验概率。
# 此处是结果表；下一节才绘制后验差值分布图。
rt_difference_summary <- rt_difference_draws %>%
  group_by(matchness, comparison) %>%
  summarise(
    posterior_mean = mean(value),
    lower95 = quantile(value, 0.025),
    upper95 = quantile(value, 0.975),
    probability_gt_zero = mean(value > 0),
    .groups = "drop"
  ) %>%
  mutate(CI_contains_zero = lower95 <= 0 & upper95 >= 0)

print(rt_difference_summary)

write_csv(
  rt_difference_summary,
  file.path(result_dir, "RT_六组模型差值_ms_95CI.csv")
)

saveRDS(
  list(predicted_ms = rt_predicted_ms, differences = rt_difference_draws),
  file.path(result_dir, "RT_条件与差值后验抽样.rds")
)


# ============================================================
# 6. RT后验差值分布图
# ============================================================


# 6.1 设置RT和SDT共同使用的图形样式。
posterior_theme <- theme_classic(base_size = 12, base_family = "Times New Roman") +
  theme(
    axis.text = element_text(colour = "black"),
    strip.background = element_blank(),
    strip.text = element_text(size = 13),
    legend.position = "bottom"
  )

# 6.2 绘制匹配条件下的三组差值。
# stat_halfeye()同时画密度、后验均值和区间，不再手动生成密度数据。
# mean_qi与.width=0.95指定后验均值和等尾95%CI。
# 正差值为红色，负差值为蓝色；竖虚线表示差值为0。
p_rt_posterior <- rt_difference_draws %>%
  filter(matchness == "Match") %>%
  mutate(panel = paste0(c("A. ", "B. ", "C. ")[as.integer(comparison)], comparison)) %>%
  ggplot(aes(x = value, y = 0)) +
  stat_halfeye(
    aes(fill = after_stat(x > 0)),
    point_interval = mean_qi,
    .width = 0.95,
    normalize = "panels",
    scale = 1,
    justification = 0,
    slab_colour = NA,
    interval_colour = "black",
    point_colour = "black"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey35") +
  # 每个比较单独一栏；横轴按该比较的范围显示。
  facet_wrap(~ panel, ncol = 3, scales = "free_x") +
  scale_fill_manual(
    values = c(`FALSE` = "#3182BD", `TRUE` = "#EF554B"),
    breaks = c(TRUE, FALSE),
    labels = c("> 0", "< 0"),
    name = "Posterior difference"
  ) +
  scale_y_continuous(breaks = seq(0, 1, 0.25)) +
  labs(
    subtitle = "Correct match RT",
    x = "RT difference (ms)",
    y = "Relative density"
  ) +
  posterior_theme

print(p_rt_posterior)

# 6.3 补充分析：绘制不匹配条件下的三组差值，编号D-F。
# 与正文使用相同的后验样本、均值、等尾95%CI及图形样式。
p_rt_mismatch_posterior <- rt_difference_draws %>%
  filter(matchness == "Mismatch") %>%
  mutate(panel = paste0(c("D. ", "E. ", "F. ")[as.integer(comparison)], comparison)) %>%
  ggplot(aes(x = value, y = 0)) +
  stat_halfeye(
    aes(fill = after_stat(x > 0)),
    point_interval = mean_qi,
    .width = 0.95,
    normalize = "panels",
    scale = 1,
    justification = 0,
    slab_colour = NA,
    interval_colour = "black",
    point_colour = "black"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey35") +
  # 每个比较单独一栏；横轴按该比较的范围显示。
  facet_wrap(~ panel, ncol = 3, scales = "free_x") +
  scale_fill_manual(
    values = c(`FALSE` = "#3182BD", `TRUE` = "#EF554B"),
    breaks = c(TRUE, FALSE),
    labels = c("> 0", "< 0"),
    name = "Posterior difference"
  ) +
  scale_y_continuous(breaks = seq(0, 1, 0.25)) +
  labs(
    subtitle = "Correct mismatch RT",
    x = "RT difference (ms)",
    y = "Relative density"
  ) +
  posterior_theme

print(p_rt_mismatch_posterior)

# RT后验差值图在本节只绘制和预览；第11节绘制dprime后统一组合保存。

# 6.4 分别保存正文三组表和补充六组表；数值来自第5.5节。
write_csv(
  filter(rt_difference_summary, matchness == "Match"),
  file.path(result_dir, "RT_正文_匹配三组差值_ms_95CI.csv")
)

write_csv(
  rt_difference_summary,
  file.path(result_dir, "RT_补充_匹配与不匹配六组差值_ms_95CI.csv")
)



# ============================================================
# 7. RT零模型与贝叶斯因子
# ============================================================


# 7.1 取出H1设计矩阵中的两列交互编码，供主效应零模型使用。
# 原样保留交互列，避免R在删除主效应后重新编码并补回相应自由度。
rt_design <- model.matrix(~ valence * matchness, data = rt_data)
rt_data$valenceGood_x_matchness <- rt_design[, "valenceGood:matchnessMismatch"]
rt_data$valenceBad_x_matchness <- rt_design[, "valenceBad:matchnessMismatch"]

# 7.2 效价主效应H0：两个总体效价系数为0。
# 保留匹配主效应、两列固定交互，以及完整被试随机结构。
# 它不要求每个匹配条件下都没有效价差异。
formula_rt_no_valence <- brms::bf(
  rt_s ~ matchness + valenceGood_x_matchness + valenceBad_x_matchness +
    (1 + valence * matchness | subject_uid),
  center = FALSE
)

# 数据、似然、共同参数先验、随机结构与RT完整模型一致。


rt_no_valence <- brm(
  formula = formula_rt_no_valence,
  data = rt_data,
  family = lognormal(),
  prior = prior_rt,
  backend = "cmdstanr",
  normalize = TRUE,
  chains = 4,
  iter = 4000,
  warmup = 500,
  cores = 4,
  seed = 20260927,
  refresh = 100,
  sig_figs = 18,
  control = list(adapt_delta = 0.95, max_treedepth = 15),
  save_pars = save_pars(all = TRUE),
  sample_prior = "yes",
  output_dir = chain_dir,
  output_basename = "1a_RT_no_valence_main",
  file = file.path(model_dir, "1a_RT_no_valence_main"),
  file_refit = "on_change"
)


# 7.3 匹配主效应H0：总体匹配系数为0。
# 保留效价主效应、两列固定交互，以及完整被试随机结构。
formula_rt_no_matchness <- brms::bf(
  rt_s ~ valence + valenceGood_x_matchness + valenceBad_x_matchness +
    (1 + valence * matchness | subject_uid),
  center = FALSE
)


rt_no_matchness <- brm(
  formula = formula_rt_no_matchness,
  data = rt_data,
  family = lognormal(),
  prior = prior_rt,
  backend = "cmdstanr",
  normalize = TRUE,
  chains = 4,
  iter = 4000,
  warmup = 500,
  cores = 4,
  seed = 20260927,
  refresh = 100,
  sig_figs = 18,
  control = list(adapt_delta = 0.95, max_treedepth = 15),
  save_pars = save_pars(all = TRUE),
  sample_prior = "yes",
  output_dir = chain_dir,
  output_basename = "1a_RT_no_matchness_main",
  file = file.path(model_dir, "1a_RT_no_matchness_main"),
  file_refit = "on_change"
)


# 7.4 交互效应H0：两个总体效价×匹配交互系数为0。
# 保留两个固定主效应，以及完整被试随机结构。
formula_rt_no_interaction <- brms::bf(
  rt_s ~ valence + matchness + (1 + valence * matchness | subject_uid),
  center = FALSE
)


rt_no_interaction <- brm(
  formula = formula_rt_no_interaction,
  data = rt_data,
  family = lognormal(),
  prior = prior_rt,
  backend = "cmdstanr",
  normalize = TRUE,
  chains = 4,
  iter = 4000,
  warmup = 500,
  cores = 4,
  seed = 20260927,
  refresh = 100,
  sig_figs = 18,
  control = list(adapt_delta = 0.95, max_treedepth = 15),
  save_pars = save_pars(all = TRUE),
  sample_prior = "yes",
  output_dir = chain_dir,
  output_basename = "1a_RT_no_interaction",
  file = file.path(model_dir, "1a_RT_no_interaction"),
  file_refit = "on_change"
)


# 7.5 汇总每个RT零模型的Rhat和Bulk ESS，供报告及讨论。
rt_null_models <- list(
  no_valence = rt_no_valence,
  no_matchness = rt_no_matchness,
  no_interaction = rt_no_interaction
)

rt_null_diagnostics <- list()

for (model_name in names(rt_null_models)) {
  fit <- rt_null_models[[model_name]]
  parameters <- grep("^(b_|sd_|cor_|r_|sigma$)", variables(fit), value = TRUE)

  diagnostics <- summarise_draws(
    as_draws_array(fit, variable = parameters),
    rhat = posterior::rhat,
    ess_bulk = posterior::ess_bulk
  )

  rt_null_diagnostics[[model_name]] <- data.frame(
    model = model_name,
    max_Rhat = max(diagnostics$rhat),
    min_Bulk_ESS = min(diagnostics$ess_bulk)
  )
}

rt_null_diagnostics <- bind_rows(rt_null_diagnostics)
print(rt_null_diagnostics)
write_csv(rt_null_diagnostics, file.path(result_dir, "RT_零模型诊断.csv"))

# 7.6 为桥采样准备概率密度接口。
# add_rstan_model()内部使用chains=0，不重新进行MCMC。
# 上面完整模型和零模型的实际MCMC后端均为CmdStanR。
rt_full_density <- add_rstan_model(model_rt, overwrite = TRUE)
rt_valence_density <- add_rstan_model(rt_no_valence, overwrite = TRUE)
rt_matchness_density <- add_rstan_model(rt_no_matchness, overwrite = TRUE)
rt_interaction_density <- add_rstan_model(rt_no_interaction, overwrite = TRUE)

# 7.7 分别计算H1和各H0的log边际似然。
# repetitions=1表示单次桥采样；maxiter是桥采样算法的求解上限。


set.seed(20261121)
bridge_rt_full <- bridge_sampler(
  rt_full_density$fit,
  method = "warp3",
  repetitions = 1,
  maxiter = 10000,
  cores = 1
)

# 本模型桥采样完成即保存；中断后可用readRDS()读取这个文件。
saveRDS(bridge_rt_full, file.path(result_dir, "bridge_rt_full.rds"))

# TRUE仅表示log边际似然为有限数值，不表示重复计算稳定。
print(data.frame(
  model = "bridge_rt_full",
  logml = bridge_rt_full$logml,
  logml_is_finite = is.finite(bridge_rt_full$logml)
))


set.seed(20261122)
bridge_rt_no_valence <- bridge_sampler(
  rt_valence_density$fit,
  method = "warp3",
  repetitions = 1,
  maxiter = 10000,
  cores = 1
)

# 本模型桥采样完成即保存；中断后可用readRDS()读取这个文件。
saveRDS(bridge_rt_no_valence, file.path(result_dir, "bridge_rt_no_valence.rds"))

# TRUE仅表示log边际似然为有限数值，不表示重复计算稳定。
print(data.frame(
  model = "bridge_rt_no_valence",
  logml = bridge_rt_no_valence$logml,
  logml_is_finite = is.finite(bridge_rt_no_valence$logml)
))


set.seed(20261123)
bridge_rt_no_matchness <- bridge_sampler(
  rt_matchness_density$fit,
  method = "warp3",
  repetitions = 1,
  maxiter = 10000,
  cores = 1
)

# 本模型桥采样完成即保存；中断后可用readRDS()读取这个文件。
saveRDS(bridge_rt_no_matchness, file.path(result_dir, "bridge_rt_no_matchness.rds"))

# TRUE仅表示log边际似然为有限数值，不表示重复计算稳定。
print(data.frame(
  model = "bridge_rt_no_matchness",
  logml = bridge_rt_no_matchness$logml,
  logml_is_finite = is.finite(bridge_rt_no_matchness$logml)
))


set.seed(20261124)
bridge_rt_no_interaction <- bridge_sampler(
  rt_interaction_density$fit,
  method = "warp3",
  repetitions = 1,
  maxiter = 10000,
  cores = 1
)

# 本模型桥采样完成即保存；中断后可用readRDS()读取这个文件。
saveRDS(bridge_rt_no_interaction, file.path(result_dir, "bridge_rt_no_interaction.rds"))

# TRUE仅表示log边际似然为有限数值，不表示重复计算稳定。
print(data.frame(
  model = "bridge_rt_no_interaction",
  logml = bridge_rt_no_interaction$logml,
  logml_is_finite = is.finite(bridge_rt_no_interaction$logml)
))


# 7.8 BF10 = p(data | H1) / p(data | H0)。
# log_BF10为两个log边际似然之差，exp(log_BF10)即BF10。
# 保留log_BF10；exp()可能溢出为Inf，不能仅凭Inf判定桥采样失败。
# 有限数值检查只随结果输出，不触发停止、重跑或追加MCMC。
rt_bf <- tibble(
  comparison = c(
    "RT效价主效应",
    "RT匹配主效应",
    "RT效价×匹配交互"
  ),
  log_BF10 = bridge_rt_full$logml - c(
    bridge_rt_no_valence$logml,
    bridge_rt_no_matchness$logml,
    bridge_rt_no_interaction$logml
  )
) %>%
  mutate(
    BF10 = exp(log_BF10),
    BF01 = exp(-log_BF10),
    log_BF10_is_finite = is.finite(log_BF10),
    BF10_is_finite = is.finite(BF10),
    note = "单次warp3；未评估重复稳定性"
  )

print(rt_bf)
write_csv(rt_bf, file.path(result_dir, "RT_全部BF.csv"))

saveRDS(
  list(
    full = bridge_rt_full,
    no_valence = bridge_rt_no_valence,
    no_matchness = bridge_rt_no_matchness,
    no_interaction = bridge_rt_no_interaction
  ),
  file.path(result_dir, "RT_桥采样对象.rds")
)


# ============================================================
# 8. SDT贝叶斯混合效应模型
# ============================================================


# 8.1 完整模型H1：对每个试次是否回答“匹配”建模。
# signal×valence对应总体dprime随效价变化；随机结构与原模型一致。
formula_sdt <- brms::bf(
  response_match ~ signal * valence + (1 + signal * valence | subject_uid),
  center = FALSE
)

# 8.2 指定SDT先验；signal系数为中性条件的总体dprime。
prior_sdt <- c(
  prior(normal(0, 1), class = "b"),
  prior(normal(-1, 1.5), class = "b", coef = "Intercept"),
  prior(normal(2, 1.5), class = "b", coef = "signal"),
  prior(exponential(1), class = "sd"),
  prior(lkj(2), class = "cor")
)

# 8.3 拟合SDT完整模型；4链，每链2000次，其中500次预热。
# 因变量是二元回答；使用Bernoulli分布与probit链接。


model_sdt <- brm(
  formula = formula_sdt,
  data = sdt_data,
  family = bernoulli(link = "probit"),
  prior = prior_sdt,
  backend = "cmdstanr",
  normalize = TRUE,
  chains = 4,
  iter = 2000,
  warmup = 500,
  cores = 4,
  seed = 20260927,
  refresh = 100,
  sig_figs = 18,
  control = list(adapt_delta = 0.95, max_treedepth = 15),
  save_pars = save_pars(all = TRUE),
  sample_prior = "yes",
  output_dir = chain_dir,
  output_basename = "1a_SDT_full",
  file = file.path(model_dir, "1a_SDT_full"),
  file_refit = "on_change"
)


# ============================================================
# 9. SDT模型诊断
# ============================================================


# 9.1 绘制主要参数的MCMC轨迹图。
p_sdt_trace <- mcmc_trace(
  as_draws_array(model_sdt),
  pars = c(
    "b_Intercept", "b_signal", "b_valenceGood", "b_valenceBad",
    "b_signal:valenceGood", "b_signal:valenceBad",
    "sd_subject_uid__Intercept"
  )
)

print(p_sdt_trace)

ggsave(
  file.path(fig_dir, "SDT_轨迹图.png"),
  p_sdt_trace, width = 11, height = 8, dpi = 200
)

# 9.2 汇总全部模型参数的最大Rhat和最小Bulk ESS。
sdt_parameters <- grep("^(b_|sd_|cor_|r_)", variables(model_sdt), value = TRUE)
sdt_draws <- as_draws_array(model_sdt, variable = sdt_parameters)

sdt_diagnostics <- summarise_draws(
  sdt_draws,
  rhat = posterior::rhat,
  ess_bulk = posterior::ess_bulk
)

# 提取发散次数与树深触顶次数。
sdt_nuts <- nuts_params(model_sdt)

sdt_diagnostic_summary <- data.frame(
  max_Rhat = max(sdt_diagnostics$rhat),
  min_Bulk_ESS = min(sdt_diagnostics$ess_bulk),
  divergences = sum(sdt_nuts$Value[sdt_nuts$Parameter == "divergent__"]),
  treedepth_hits = sum(sdt_nuts$Value[sdt_nuts$Parameter == "treedepth__"] >= 15)
)

print(sdt_diagnostic_summary)

write_csv(sdt_diagnostics, file.path(result_dir, "SDT_逐参数诊断.csv"))
write_csv(sdt_diagnostic_summary, file.path(result_dir, "SDT_诊断汇总.csv"))

# 9.3 后验预测检查：比较模型模拟的0/1回答与实际回答分布。
p_sdt_ppc <- pp_check(model_sdt, type = "bars", ndraws = 50) +
  labs(x = "Response: 0 = mismatch, 1 = match")

print(p_sdt_ppc)

ggsave(
  file.path(fig_dir, "SDT_后验预测检查.png"),
  p_sdt_ppc, width = 9, height = 6, dpi = 200
)


# ============================================================
# 10. SDT回归系数、dprime与后验差值
# ============================================================


# 10.1 回归系数beta：位于probit线性预测尺度。
sdt_beta <- as.data.frame(fixef(model_sdt, probs = c(0.025, 0.975))) %>%
  tibble::rownames_to_column("parameter") %>%
  transmute(
    parameter,
    posterior_mean = Estimate,
    lower95 = Q2.5,
    upper95 = Q97.5,
    CI_contains_zero = lower95 <= 0 & upper95 >= 0
  )

print(sdt_beta)
write_csv(sdt_beta, file.path(result_dir, "SDT_beta_95CI.csv"))

# 10.2 建立效价×信号状态的六个预测条件。
# 前三列为噪声的Neutral/Good/Bad，后三列为信号的Neutral/Good/Bad。
sdt_grid <- expand.grid(
  valence = c("Neutral", "Good", "Bad"),
  signal = c(0L, 1L)
)

sdt_grid$subject_uid <- sdt_data$subject_uid[1]
sdt_grid$valence <- factor(sdt_grid$valence, levels = levels(sdt_data$valence))
contrasts(sdt_grid$valence) <- contrasts(sdt_data$valence)

# 10.3 计算probit线性预测值；transform=FALSE保留该尺度。
sdt_eta <- posterior_linpred(
  model_sdt,
  newdata = sdt_grid,
  re_formula = NA,
  transform = FALSE
)

# dprime = 同一效价下信号的probit预测值 - 噪声的probit预测值。
# 这里不是直接相减两个回答概率。
dprime_wide <- tibble(
  .draw = seq_len(nrow(sdt_eta)),
  Neutral = sdt_eta[, 4] - sdt_eta[, 1],
  Good = sdt_eta[, 5] - sdt_eta[, 2],
  Bad = sdt_eta[, 6] - sdt_eta[, 3]
)

# 10.4 汇总各效价dprime的后验均值和等尾95%CI。
dprime_draws <- dprime_wide %>%
  pivot_longer(
    cols = -.draw,
    names_to = "valence",
    values_to = "value"
  )

dprime_summary <- dprime_draws %>%
  group_by(valence) %>%
  summarise(
    posterior_mean = mean(value),
    lower95 = quantile(value, 0.025),
    upper95 = quantile(value, 0.975),
    .groups = "drop"
  )

# 10.5 在每组后验样本内计算三个效价间的dprime差值。
sdt_difference_draws <- dprime_wide %>%
  transmute(
    .draw,
    GN = Good - Neutral,
    GB = Good - Bad,
    NB = Neutral - Bad
  ) %>%
  pivot_longer(
    cols = -.draw,
    names_to = "comparison",
    values_to = "value"
  ) %>%
  mutate(
    comparison = factor(
      comparison,
      levels = c("GN", "GB", "NB"),
      labels = c("Good - Neutral", "Good - Bad", "Neutral - Bad")
    )
  )

# 10.6 汇总差值的后验均值、95%CI、差值为正的后验概率。
sdt_difference_summary <- sdt_difference_draws %>%
  group_by(comparison) %>%
  summarise(
    posterior_mean = mean(value),
    lower95 = quantile(value, 0.025),
    upper95 = quantile(value, 0.975),
    probability_gt_zero = mean(value > 0),
    .groups = "drop"
  ) %>%
  mutate(CI_contains_zero = lower95 <= 0 & upper95 >= 0)

print(dprime_summary)
print(sdt_difference_summary)

write_csv(dprime_summary, file.path(result_dir, "SDT_三效价dprime后验.csv"))
write_csv(sdt_difference_summary, file.path(result_dir, "SDT_三组模型差值_95CI.csv"))

saveRDS(
  list(dprime = dprime_wide, differences = sdt_difference_draws),
  file.path(result_dir, "SDT_dprime与差值后验抽样.rds")
)


# ============================================================
# 11. 绘制dprime并统一保存正文六组图与补充九组图
# ============================================================


# 11.1 用与RT相同的stat_halfeye()画三组dprime差值。
# 点表示后验均值，横线表示等尾95%CI；颜色区分差值正负。
p_sdt_posterior <- sdt_difference_draws %>%
  mutate(panel = paste0(c("D. ", "E. ", "F. ")[as.integer(comparison)], comparison)) %>%
  ggplot(aes(x = value, y = 0)) +
  stat_halfeye(
    aes(fill = after_stat(x > 0)),
    point_interval = mean_qi,
    .width = 0.95,
    normalize = "panels",
    scale = 1,
    justification = 0,
    slab_colour = NA,
    interval_colour = "black",
    point_colour = "black"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey35") +
  # 三个比较横向排列，与上方RT对应。
  facet_wrap(~ panel, ncol = 3, scales = "free_x") +
  scale_fill_manual(
    values = c(`FALSE` = "#3182BD", `TRUE` = "#EF554B"),
    breaks = c(TRUE, FALSE),
    labels = c("> 0", "< 0"),
    name = "Posterior difference"
  ) +
  scale_y_continuous(breaks = seq(0, 1, 0.25)) +
  labs(
    subtitle = "d prime",
    x = "Difference in d prime",
    y = "Relative density"
  ) +
  posterior_theme

print(p_sdt_posterior)

# 11.2 正文六组图：匹配RT为A-C，dprime为D-F；共用一个图例。
p_posterior <- (p_rt_posterior / p_sdt_posterior) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

print(p_posterior)

# 11.3 补充分析图的dprime编号改为G-I，后验数值和画图样式保持一致。
sdt_supplement_draws <- sdt_difference_draws %>%
  mutate(panel = paste0(c("G. ", "H. ", "I. ")[as.integer(comparison)], comparison))

# 替换图中数据标签；沿用上面dprime图的密度、均值、95%CI与样式设置。
p_sdt_supplement <- p_sdt_posterior + sdt_supplement_draws

# 11.4 补充九组图：匹配RT为A-C，不匹配RT为D-F，dprime为G-I。
p_posterior_supplement <- (p_rt_posterior / p_rt_mismatch_posterior / p_sdt_supplement) +
  plot_layout(guides = "collect", heights = c(1, 1, 1)) &
  theme(legend.position = "bottom")

print(p_posterior_supplement)

# 11.5 dprime绘制完成后，统一保存正文六组图：PNG、JPG和PDF。
ggsave(
  file.path(fig_dir, "02_1a_RT与dprime后验差值_A-F.png"),
  p_posterior, width = 12.6, height = 7.8, dpi = 300
)

ggsave(
  file.path(fig_dir, "02_1a_RT与dprime后验差值_A-F.jpg"),
  p_posterior, width = 12.6, height = 7.8, dpi = 300
)

ggsave(
  file.path(fig_dir, "02_1a_RT与dprime后验差值_A-F.pdf"),
  p_posterior, device = grDevices::cairo_pdf, width = 12.6, height = 7.8
)


# 11.6 统一保存补充九组图：PNG、JPG和PDF。
# 与正文保持相同宽度；增加一行的高度，保持各行排版一致。
ggsave(
  file.path(fig_dir, "S1_1a_补充_RT匹配_不匹配与dprime后验差值_A-I.png"),
  p_posterior_supplement, width = 12.6, height = 11.4, dpi = 300
)

ggsave(
  file.path(fig_dir, "S1_1a_补充_RT匹配_不匹配与dprime后验差值_A-I.jpg"),
  p_posterior_supplement, width = 12.6, height = 11.4, dpi = 300
)

ggsave(
  file.path(fig_dir, "S1_1a_补充_RT匹配_不匹配与dprime后验差值_A-I.pdf"),
  p_posterior_supplement, device = grDevices::cairo_pdf, width = 12.6, height = 11.4
)


# ============================================================
# 12. SDT效价效应的H0与贝叶斯因子
# ============================================================


# 12.1 dprime效价效应的H0：删除固定signal×valence交互。
# 这是检验效价效应所必需的H0；
# H1：dprime_Neutral=b_signal；其他效价再加各自的signal×valence系数。
# H0：这两个交互系数都为0，因此三个效价的总体dprime均为b_signal。
# 仍保留效价回答基线，以及被试随机截距、随机斜率和相关。
# 这是总体敏感性相等，不要求每名被试的敏感性完全相等。
formula_sdt_equal <- brms::bf(
  response_match ~ signal + valence + (1 + signal * valence | subject_uid),
  center = FALSE
)

# 12.2 拟合SDT零模型；数据、Bernoulli-probit和共同参数先验保持一致。


model_sdt_equal <- brm(
  formula = formula_sdt_equal,
  data = sdt_data,
  family = bernoulli(link = "probit"),
  prior = prior_sdt,
  backend = "cmdstanr",
  normalize = TRUE,
  chains = 4,
  iter = 2000,
  warmup = 500,
  cores = 4,
  seed = 20260927,
  refresh = 100,
  sig_figs = 18,
  control = list(adapt_delta = 0.95, max_treedepth = 15),
  save_pars = save_pars(all = TRUE),
  sample_prior = "yes",
  output_dir = chain_dir,
  output_basename = "1a_SDT_equal_dprime",
  file = file.path(model_dir, "1a_SDT_equal_dprime"),
  file_refit = "on_change"
)


# 12.3 汇总SDT零模型的Rhat和Bulk ESS。
sdt_equal_parameters <- grep(
  "^(b_|sd_|cor_|r_)", variables(model_sdt_equal), value = TRUE
)

sdt_equal_diagnostics <- summarise_draws(
  as_draws_array(model_sdt_equal, variable = sdt_equal_parameters),
  rhat = posterior::rhat,
  ess_bulk = posterior::ess_bulk
)

sdt_equal_summary <- data.frame(
  model = "SDT_equal_dprime",
  max_Rhat = max(sdt_equal_diagnostics$rhat),
  min_Bulk_ESS = min(sdt_equal_diagnostics$ess_bulk)
)

print(sdt_equal_summary)
write_csv(sdt_equal_summary, file.path(result_dir, "SDT_零模型诊断.csv"))

# 12.4 准备概率密度接口；不在此重新进行MCMC。
sdt_full_density <- add_rstan_model(model_sdt, overwrite = TRUE)
sdt_equal_density <- add_rstan_model(model_sdt_equal, overwrite = TRUE)

# 12.5 分别计算完整模型和零模型的log边际似然。
set.seed(20261125)
bridge_sdt_full <- bridge_sampler(
  sdt_full_density$fit,
  method = "warp3",
  repetitions = 1,
  maxiter = 100000,
  cores = 1
)

# 本模型桥采样完成即保存；中断后可用readRDS()读取这个文件。
saveRDS(bridge_sdt_full, file.path(result_dir, "bridge_sdt_full.rds"))

# TRUE仅表示log边际似然为有限数值，不表示重复计算稳定。
print(data.frame(
  model = "bridge_sdt_full",
  logml = bridge_sdt_full$logml,
  logml_is_finite = is.finite(bridge_sdt_full$logml)
))

set.seed(20261126)
bridge_sdt_equal <- bridge_sampler(
  sdt_equal_density$fit,
  method = "warp3",
  repetitions = 1,
  maxiter = 100000,
  cores = 1
)

# 本模型桥采样完成即保存；中断后可用readRDS()读取这个文件。
saveRDS(bridge_sdt_equal, file.path(result_dir, "bridge_sdt_equal.rds"))

# TRUE仅表示log边际似然为有限数值，不表示重复计算稳定。
print(data.frame(
  model = "bridge_sdt_equal",
  logml = bridge_sdt_equal$logml,
  logml_is_finite = is.finite(bridge_sdt_equal$logml)
))

# 12.6 由两个模型的边际似然之比计算BF10。
# 此BF检验整体dprime效价效应，不是三项两两差值各自的BF。
sdt_bf <- tibble(
  comparison = "dprime效价效应",
  log_BF10 = bridge_sdt_full$logml - bridge_sdt_equal$logml
) %>%
  mutate(
    BF10 = exp(log_BF10),
    BF01 = exp(-log_BF10),
    log_BF10_is_finite = is.finite(log_BF10),
    BF10_is_finite = is.finite(BF10),
    note = "单次warp3；未评估重复稳定性"
  )

print(sdt_bf)
write_csv(sdt_bf, file.path(result_dir, "SDT_整体效价BF.csv"))

saveRDS(
  list(full = bridge_sdt_full, equal = bridge_sdt_equal),
  file.path(result_dir, "SDT_桥采样对象.rds")
)

# 保存本次R和各包的版本信息。
writeLines(capture.output(sessionInfo()), file.path(result_dir, "R会话信息.txt"))
