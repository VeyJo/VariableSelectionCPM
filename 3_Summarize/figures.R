##----------------------------------------------------------------------------##
##  create figures
##----------------------------------------------------------------------------##

# Setup
library(tidyverse)
library(scales)
library(kableExtra)
library(cowplot)
library(ggpubr)
library(mfp2)
library(effects)

## color palettes
colors8a <- c("#424949", "#F7D42A", "#F89217", "#CF3E53", "#884EA0", "#91B3D7", "#26897E", "#8DBFA8")
colors7a <- c(           "#F7D42A", "#F89217", "#CF3E53", "#884EA0", "#91B3D7", "#26897E", "#8DBFA8")
colors7b <- c("#424949",            "#F89217", "#CF3E53", "#884EA0", "#91B3D7", "#26897E", "#8DBFA8")
colors6a <- c(                      "#F89217", "#CF3E53", "#884EA0", "#91B3D7", "#26897E", "#8DBFA8")

## oracle
colorsOa <- c("#424949", "#F7D42A", "#F89217", "#CF3E53", "#884EA0", "#91B3D7", "#26897E")
colorsOb <- c("#424949",            "#F89217", "#CF3E53", "#884EA0", "#91B3D7", "#26897E")

## predictor variables
colors_preds <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#91D1C2", "#8491B4", "#DF8F44", "#631879")

# colors8 <- c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628", "#f781bf")

## colors for summary table
color_mapping <- data.frame(
  method = c("LMSS",   "ENET",    "GBM",     "GBT",     "RFB",     "RFH"),
  color = c("#F89217", "#CF3E53", "#884EA0", "#91B3D7", "#26897E", "#8DBFA8"))

## linetype 
ltype8a <- c("solid", "twodash", "solid", "dashed", "dotdash", "solid", "longdash", "dotdash")
ltype7a <- c("solid", "dashed", "longdash", "dotdash", "solid", "longdash", "dotdash")
ltype7b <- c("solid", "dashed", "longdash", "dotdash", "solid", "longdash", "dotdash")
ltype7bb <- c("solid", "solid", "dashed", "dotdash", "solid", "longdash", "dotdash")
ltype6 <- c("solid", "dashed", "dotdash", "solid", "longdash", "dotdash")

## specify predictor variables
preds <- paste0("x", c(1,3,4,5,6,8,10,11))
preds_reg <- ifelse(preds == "x4", "x4_1", preds)
vars <- paste0("x", 1:15)
vars_reg <- paste0("x", c(1:3, "4a", "4b", 5:8, "9a", "9b", 10:15))

## select and specify order of methods for plots
methods <- c("MFP", "LMSS", "ENET", "GBM", "GBT", "RFB", "RFH")
methods8 <- c("TOM", "MFP", "LMSS", "ENET", "GBM", "GBT", "RFB", "RFH")
oracles <- c("TOM", "MFP-O", "LMSS-O", "ENET-O", "GBM-O", "GBT-O", "RF-O")


## themes
theme1 <- theme(
  plot.background = element_rect(fill = "white"),
  panel.background = element_rect(fill = "white"),
  panel.grid.major.y = element_line(color = "#e5e5e5", linewidth = 0.3),
  panel.grid.major.x = element_blank(),
  axis.line = element_line(color = "black"),
  panel.spacing.y = unit(0.5, "cm"),
  axis.text.x = element_text(angle = 90, vjust = 0.25, size = 9),
  axis.title.x = element_text(size = 12),
  axis.title.y = element_text(size = 12),
  legend.text = element_text(size = 12),
  legend.title = element_text(size = 12),
  strip.text = element_text(size = 12))


theme2 <- theme(
  plot.background = element_rect(fill = "white"),
  panel.background = element_rect(fill = "white"),
  panel.grid.major.y = element_line(color = "#e5e5e5", linewidth = 0.3),
  panel.grid.major.x = element_blank(),
  axis.line = element_line(color = "black"),
  panel.spacing.y = unit(0.5, "cm"),
  axis.text.x = element_text(size = 9),
  axis.title.x = element_text(size = 12),
  axis.title.y = element_text(size = 12))

theme3 <- theme(
  plot.background = element_rect(fill = "white"),
  panel.background = element_rect(fill = "white"),
  panel.grid.major.y = element_line(color = "#e5e5e5", linewidth = 0.5),
  panel.grid.major.x = element_blank(),
  axis.line = element_line(color = "black"),
  panel.spacing.y = unit(0.5, "cm"),
  axis.title = element_text(size = 12),
  legend.position = "bottom")



## import summarized and joint performance measures
path_res <- "Results/results_summarized/"


res_A   <- readRDS(paste0(path_res, "results_perf_measures_A.rds"))
res_Aom <- readRDS(paste0(path_res, "results_perf_measures_A_oracles.rds"))
res_Br  <- readRDS(paste0(path_res, "results_perf_measures_B_regression.rds"))
res_Bt  <- readRDS(paste0(path_res, "results_perf_measures_B_trees.rds"))
res_Bo  <- readRDS(paste0(path_res, "results_perf_measures_B_TrueOracle.rds"))




## adapt structure to provide ordering for plots
res_A <- res_A |> 
  mutate(method = factor(method, 
                         levels = c("tom", "mfp", "lmss", "enet", "gbm", "gbt", "rfb", "rfh"),
                         labels = c("TOM", "MFP", "LMSS", "ENET", "GBM", "GBT", "RFB", "RFH"))) |> 
  ## when not converged take minimal model (only x5)
  mutate(RMSPE_min = ifelse(is.na(RMSPE_min), RMSPE_final, RMSPE_min))

res_Aom <- res_Aom |> 
  mutate(method = factor(method,
                         levels = c("mfp_oracle", "lm_oracle", "enet_oracle", "gbm_oracle", "gbt_oracle", "rf_oracle"),
                         labels = c("MFP-O", "LMSS-O", "ENET-O", "GBM-O", "GBT-O", "RF-O")))

res_Br <- res_Br |> 
  mutate(method = factor(method, 
                         levels = c("mfp", "lmss", "enet", "gbm", "gbt", "rfb", "rfh"),
                         labels = c("MFP", "LMSS", "ENET", "GBM", "GBT", "RFB", "RFH")),
         compl = factor(compl, levels = c("A", "B", "C", "D")),
         R2 = as.numeric(R2),
         n = as.numeric(n))
res_Bt <- res_Bt |> 
  mutate(method = factor(method, 
                         levels = c("mfp", "lmss", "enet", "gbm", "gbt", "rfb", "rfh"),
                         labels = c("MFP", "LMSS", "ENET", "GBM", "GBT", "RFB", "RFH")),
         compl = factor(compl, levels = c("A", "B", "C", "D")),
         R2 = as.numeric(R2),
         n = as.numeric(n))


# Results
## Variable selection

## summarize performance measures regarding variable selection - long format
df_VS_sum_l <- res_A |> 
  filter(method != "TOM") |> 
  select(1:6, MS, SV, MSF, TIF, TEF, TCR, tau) |> 
  group_by(compl, R2, n, method) |> 
  summarize(MSN = sum(MSF, na.rm = TRUE),  ## number of correct model selection
            across(MS:tau, \(x) mean(x, na.rm = TRUE)), 
            .groups = "drop") |> 
  ## calculate combined measure of TIF and TEF
  mutate(TIFTEF = (TIF+TEF)/2)


### Model size
ggplot(df_VS_sum_l,
       aes(y = SV, x = n, color = method, linetype = method)) +
  geom_point(size = 1.6) +
  geom_line(linewidth = 0.77) + 
  facet_grid(compl ~ R2, scales = "fixed", 
             labeller = labeller(R2 = function(x) paste0("R² = ", x))) +
  scale_color_manual(values = colors7a) +
  scale_linetype_manual(values = ltype7a) +
  scale_x_continuous(breaks = c(100, 250, 500, 1000),
                     limits = c(90, 1050)) +
  labs(x = "Sample size", y = "Average ratio of selected to available variables",
       color = "Method", linetype = "Method") +
  theme1


### TIF and TEF
p1 <- ggplot(df_VS_sum_l,
             aes(y = TIF, x = n, color = method, linetype = method)) +
  geom_point(size = 1.6) +
  geom_line(linewidth = 0.77) + 
  facet_grid(compl ~ R2, scales = "fixed",
             labeller = labeller(R2 = function(x) paste0("R² = ", x))) +
  scale_color_manual(values = colors7a) +
  scale_linetype_manual(values = ltype7a) +
  scale_x_continuous(breaks = c(100, 250, 500, 1000),
                     limits = c(90, 1050)) +
  labs(x = "Sample size", y = "Average TIF",
       color = "Method", linetype = "Method") +
  theme1

p1 <- p1 + theme(legend.position = "none")

p2 <- ggplot(df_VS_sum_l,
             aes(y = TEF, x = n, color = method, linetype = method)) +
  geom_point(size = 1.6) +
  geom_line(linewidth = 0.77) + 
  facet_grid(compl ~ R2, scales = "fixed", 
             labeller = labeller(R2 = function(x) paste0("R² = ", x))) +
  scale_color_manual(values = colors7a) +
  scale_linetype_manual(values = ltype7a) +
  scale_x_continuous(breaks = c(100, 250, 500, 1000),
                     limits = c(90, 1050)) +
  labs(x = "Sample size", y = "Average TEF",
       color = "Method", linetype = "Method") +
  theme1 + theme(legend.position = "none")


ggarrange(p1, p2, align = "h", 
          common.legend = TRUE, legend = "bottom")


ggplot(df_VS_sum_l,
       aes(y = TIFTEF, x = n, color = method, linetype = method)) +
  geom_point(size = 1.6) +
  geom_line(linewidth = 0.77) + 
  facet_grid(compl ~ R2, scales = "fixed", 
             labeller = labeller(R2 = function(x) paste0("R² = ", x))) +
  scale_color_manual(values = colors7a) +
  scale_linetype_manual(values = ltype7a) +
  scale_x_continuous(breaks = c(100, 250, 500, 1000),
                     limits = c(90, 1050)) +
  labs(x = "Sample size", y = "Mean of averaged TIF and TEF",
       color = "Method", linetype = "Method") +
  theme1


### MSF
ggplot(df_VS_sum_l,
       aes(y = MSF, x = n, color = method, linetype = method)) +
  geom_point(size = 1.6) +
  geom_line(linewidth = 0.77) + 
  facet_grid(compl ~ R2, scales = "fixed", 
             labeller = labeller(R2 = function(x) paste0("R² = ", x))) +
  scale_color_manual(values = colors7a) +
  scale_linetype_manual(values = ltype7a) +
  scale_x_continuous(breaks = c(100, 250, 500, 1000),
                     limits = c(90, 1050)) +
  labs(x = "Sample size", y = "Model selection frequency",
       color = "Method", linetype = "Method") +
  theme1



## VIF
### predictors
## create appropriate data frame containing selection status of variables
df_sel <- res_Br |> 
  select(c(method:converged, ends_with("sel"))) |> 
  ## state variable as selected if one of the dummy variables was selected
  mutate(x4_sel = pmax(x4_1_sel, x4_2_sel),
         x9_sel = pmax(x9_1_sel, x9_2_sel)) |> 
  select(c(method:converged, paste0("x", 1:15, "_sel"))) |> 
  ## merge with results from tree-based methods
  bind_rows(res_Bt |> select(c(method:converged, ends_with("sel")))) |> 
  rename_with(.cols = ends_with("_sel"), ~ gsub("_sel", "", .)) |> 
  ## if not converged, all variables are stated as not selected
  mutate(across(starts_with("x"), ~ifelse(!converged, 0, .)))

## calculate performance measure over repetitions
df_VIF <- df_sel |> 
  group_by(method, compl, R2, n) |> 
  summarize(across(starts_with("x"), \(x) sum(x)/500), .groups = "drop")

x5 <- ggplot(df_VIF,
             aes(y = x5, x = n, color = method, linetype = method)) +
  geom_point(size = 1.6) +
  geom_line(linewidth = 0.77) + 
  facet_grid(compl ~ R2, scales = "fixed", 
             labeller = labeller(R2 = function(x) paste0("R² = ", x))) +
  scale_color_manual(values = colors7a) +
  scale_linetype_manual(values = ltype7a) +
  scale_x_continuous(breaks = c(100, 250, 500, 1000),
                     limits = c(90, 1050)) +
  labs(x = "Sample size", y = "VIF of x5 [%]",
       color = "Method", linetype = "Method") +
  theme1

x3 <- ggplot(df_VIF,
             aes(y = x3, x = n, color = method, linetype = method)) +
  geom_point(size = 1.6) +
  geom_line(linewidth = 0.77) + 
  facet_grid(compl ~ R2, scales = "fixed", 
             labeller = labeller(R2 = function(x) paste0("R² = ", x))) +
  scale_color_manual(values = colors7a) +
  scale_linetype_manual(values = ltype7a) +
  scale_x_continuous(breaks = c(100, 250, 500, 1000),
                     limits = c(90, 1050)) +
  labs(x = "Sample size", y = "VIF of x3 [%]",
       color = "Method", linetype = "Method") +
  theme1

x1 <- ggplot(df_VIF,
             aes(y = x1, x = n, color = method, linetype = method)) +
  geom_point(size = 1.6) +
  geom_line(linewidth = 0.77) + 
  facet_grid(compl ~ R2, scales = "fixed", 
             labeller = labeller(R2 = function(x) paste0("R² = ", x))) +
  scale_color_manual(values = colors7a) +
  scale_linetype_manual(values = ltype7a) +
  scale_x_continuous(breaks = c(100, 250, 500, 1000),
                     limits = c(90, 1050)) +
  labs(x = "Sample size", y = "VIF of x1 [%]",
       color = "Method", linetype = "Method") +
  theme1

x6 <- ggplot(df_VIF,
             aes(y = x6, x = n, color = method, linetype = method)) +
  geom_point(size = 1.6) +
  geom_line(linewidth = 0.77) + 
  facet_grid(compl ~ R2, scales = "fixed", 
             labeller = labeller(R2 = function(x) paste0("R² = ", x))) +
  scale_color_manual(values = colors7a) +
  scale_linetype_manual(values = ltype7a) +
  scale_x_continuous(breaks = c(100, 250, 500, 1000),
                     limits = c(90, 1050)) +
  labs(x = "Sample size", y = "VIF of x6 [%]",
       color = "Method", linetype = "Method") +
  theme1

x11 <- ggplot(df_VIF,
              aes(y = x11, x = n, color = method, linetype = method)) +
  geom_point(size = 1.6) +
  geom_line(linewidth = 0.77) + 
  facet_grid(compl ~ R2, scales = "fixed", 
             labeller = labeller(R2 = function(x) paste0("R² = ", x))) +
  scale_color_manual(values = colors7a) +
  scale_linetype_manual(values = ltype7a) +
  scale_x_continuous(breaks = c(100, 250, 500, 1000),
                     limits = c(90, 1050)) +
  labs(x = "Sample size", y = "VIF of x11 [%]",
       color = "Method", linetype = "Method") +
  theme1

x4 <- ggplot(df_VIF,
             aes(y = x4, x = n, color = method, linetype = method)) +
  geom_point(size = 1.6) +
  geom_line(linewidth = 0.77) + 
  facet_grid(compl ~ R2, scales = "fixed", 
             labeller = labeller(R2 = function(x) paste0("R² = ", x))) +
  scale_color_manual(values = colors7a) +
  scale_linetype_manual(values = ltype7a) +
  scale_x_continuous(breaks = c(100, 250, 500, 1000),
                     limits = c(90, 1050)) +
  labs(x = "Sample size", y = "VIF of x4 [%]",
       color = "Method", linetype = "Method") +
  theme1

x10 <- ggplot(df_VIF,
              aes(y = x10, x = n, color = method, linetype = method)) +
  geom_point(size = 1.6) +
  geom_line(linewidth = 0.77) + 
  facet_grid(compl ~ R2, scales = "fixed", 
             labeller = labeller(R2 = function(x) paste0("R² = ", x))) +
  scale_color_manual(values = colors7a) +
  scale_linetype_manual(values = ltype7a) +
  scale_x_continuous(breaks = c(100, 250, 500, 1000),
                     limits = c(90, 1050)) +
  labs(x = "Sample size", y = "VIF of x10 [%]",
       color = "Method", linetype = "Method") +
  theme1

x8 <- ggplot(df_VIF,
             aes(y = x8, x = n, color = method, linetype = method)) +
  geom_point(size = 1.6) +
  geom_line(linewidth = 0.77) + 
  facet_grid(compl ~ R2, scales = "fixed", 
             labeller = labeller(R2 = function(x) paste0("R² = ", x))) +
  scale_color_manual(values = colors7a) +
  scale_linetype_manual(values = ltype7a) +
  scale_x_continuous(breaks = c(100, 250, 500, 1000),
                     limits = c(90, 1050)) +
  labs(x = "Sample size", y = "VIF of x8 [%]",
       color = "Method", linetype = "Method") +
  theme1 + 
  theme(legend.position = "bottom")

leg <- ggpubr::get_legend(x8)

plots <- lapply(list(x1, x3, x4, x5, x6, x8, x11, x10), function(p) {
  p + 
    theme(legend.position = "none") +
    scale_y_continuous(limits = c(0, 1.05), n.breaks = 4)
})

combined_plot <- plot_grid(
  plotlist = plots,
  ncol = 2,
  align = "hv")

plot_grid(
  combined_plot,
  leg,
  ncol = 1,
  rel_heights = c(1, 0.06))


### non-predictors
x2 <- ggplot(df_VIF,
             aes(y = x2, x = n, color = method, linetype = method)) +
  geom_point(size = 1.6) +
  geom_line(linewidth = 0.77) + 
  facet_grid(compl ~ R2, scales = "fixed", 
             labeller = labeller(R2 = function(x) paste0("R² = ", x))) +
  scale_color_manual(values = colors7a) +
  scale_linetype_manual(values = ltype7a) +
  scale_x_continuous(breaks = c(100, 250, 500, 1000),
                     limits = c(90, 1050)) +
  labs(x = "Sample size", y = "VIF of x2 [%]",
       color = "Method", linetype = "Method") +
  theme1

x7 <- ggplot(df_VIF,
             aes(y = x7, x = n, color = method, linetype = method)) +
  geom_point(size = 1.6) +
  geom_line(linewidth = 0.77) + 
  facet_grid(compl ~ R2, scales = "fixed", 
             labeller = labeller(R2 = function(x) paste0("R² = ", x))) +
  scale_color_manual(values = colors7a) +
  scale_linetype_manual(values = ltype7a) +
  scale_x_continuous(breaks = c(100, 250, 500, 1000),
                     limits = c(90, 1050)) +
  labs(x = "Sample size", y = "VIF of x7 [%]",
       color = "Method", linetype = "Method") +
  theme1

x9 <- ggplot(df_VIF,
             aes(y = x9, x = n, color = method, linetype = method)) +
  geom_point(size = 1.6) +
  geom_line(linewidth = 0.77) + 
  facet_grid(compl ~ R2, scales = "fixed", 
             labeller = labeller(R2 = function(x) paste0("R² = ", x))) +
  scale_color_manual(values = colors7a) +
  scale_linetype_manual(values = ltype7a) +
  scale_x_continuous(breaks = c(100, 250, 500, 1000),
                     limits = c(90, 1050)) +
  labs(x = "Sample size", y = "VIF of x9 [%]",
       color = "Method", linetype = "Method") +
  theme1

x12 <- ggplot(df_VIF,
              aes(y = x12, x = n, color = method, linetype = method)) +
  geom_point(size = 1.6) +
  geom_line(linewidth = 0.77) + 
  facet_grid(compl ~ R2, scales = "fixed", 
             labeller = labeller(R2 = function(x) paste0("R² = ", x))) +
  scale_color_manual(values = colors7a) +
  scale_linetype_manual(values = ltype7a) +
  scale_x_continuous(breaks = c(100, 250, 500, 1000),
                     limits = c(90, 1050)) +
  labs(x = "Sample size", y = "VIF of x12 [%]",
       color = "Method", linetype = "Method") +
  theme1

x13 <- ggplot(df_VIF,
              aes(y = x13, x = n, color = method, linetype = method)) +
  geom_point(size = 1.6) +
  geom_line(linewidth = 0.77) + 
  facet_grid(compl ~ R2, scales = "fixed", 
             labeller = labeller(R2 = function(x) paste0("R² = ", x))) +
  scale_color_manual(values = colors7a) +
  scale_linetype_manual(values = ltype7a) +
  scale_x_continuous(breaks = c(100, 250, 500, 1000),
                     limits = c(90, 1050)) +
  labs(x = "Sample size", y = "VIF of x13 [%]",
       color = "Method", linetype = "Method") +
  theme1

x14 <- ggplot(df_VIF,
              aes(y = x14, x = n, color = method, linetype = method)) +
  geom_point(size = 1.6) +
  geom_line(linewidth = 0.77) + 
  facet_grid(compl ~ R2, scales = "fixed", 
             labeller = labeller(R2 = function(x) paste0("R² = ", x))) +
  scale_color_manual(values = colors7a) +
  scale_linetype_manual(values = ltype7a) +
  scale_x_continuous(breaks = c(100, 250, 500, 1000),
                     limits = c(90, 1050)) +
  labs(x = "Sample size", y = "VIF of x14 [%]",
       color = "Method", linetype = "Method") +
  theme1

x15 <- ggplot(df_VIF,
              aes(y = x15, x = n, color = method, linetype = method)) +
  geom_point(size = 1.6) +
  geom_line(linewidth = 0.77) + 
  facet_grid(compl ~ R2, scales = "fixed", 
             labeller = labeller(R2 = function(x) paste0("R² = ", x))) +
  scale_color_manual(values = colors7a) +
  scale_linetype_manual(values = ltype7a) +
  scale_x_continuous(breaks = c(100, 250, 500, 1000),
                     limits = c(90, 1050)) +
  labs(x = "Sample size", y = "VIF of x15 [%]",
       color = "Method", linetype = "Method") +
  theme1

plots <- lapply(list(x2, x7, x9, x12, x13, x14, x15), function(p) {
  p + 
    theme(legend.position = "none") +
    scale_y_continuous(limits = c(0, 1.05), n.breaks = 4)
})
plots[[8]] <- leg

cowplot::plot_grid(
  plotlist = plots,
  ncol = 2)




### Kendall's tau
ggplot(df_VS_sum_l,
       aes(y = tau, x = n, color = method, linetype = method)) +
  geom_point(size = 1.6) +
  geom_line(linewidth = 0.77) + 
  facet_grid(compl ~ R2, scales = "fixed", 
             labeller = labeller(R2 = function(x) paste0("R² = ", x))) +
  scale_color_manual(values = colors7a) +
  scale_linetype_manual(values = ltype7a) +
  scale_x_continuous(breaks = c(100, 250, 500, 1000),
                     limits = c(90, 1050)) +
  labs(x = "Sample size", y = bquote("Average Kendall's  " * tau),
       color = "Method", linetype = "Method") +
  theme1


### Rankings
## create appropriate data frame containing selection status of variables
df_rank <- res_Br |> 
  select(c(method:converged, ends_with("rank"))) |> 
  ## state variable as selected if one of the dummy variables was selected
  mutate(x4_rank = pmax(x4_1_rank, x4_2_rank),
         x9_rank = pmax(x9_1_rank, x9_2_rank)) |> 
  select(c(method:converged, paste0("x", 1:15, "_rank"))) |> 
  ## merge with results from tree-based methods
  bind_rows(res_Bt |> select(c(method:converged, ends_with("rank")))) |> 
  rename_with(.cols = ends_with("_rank"), ~ gsub("_rank", "", .)) |> 
  filter(converged == TRUE)

r3 <- df_rank_sum |> 
  select(method:n, x3, x13) |> 
  filter(compl %in% c("B", "C"), R2 == 0.5, n == 500)

## calculate performance measure over repetitions
## median rank
df_rank_sum <- df_rank |> 
  group_by(method, compl, R2, n) |> 
  summarize(across(starts_with("x"), \(x) median(x, na.rm = TRUE)), .groups = "drop")


p3 <- ggplot(df_rank_sum, aes(y = x3, x = factor(n), fill = method)) + 
  geom_bar(position = "dodge", stat = "identity") + 
  labs(y = "Median rank of x3", x = "Sample size",
       fill = "Method") +
  facet_grid(compl ~ R2) + 
  scale_y_continuous(limits = c(0,15)) +
  scale_fill_manual(values = colors7a) +
  theme2 

legend_ranks <- cowplot::get_legend(p3)
p3 <- p3 + theme(legend.position = "none")

p13 <- ggplot(df_rank_sum, aes(y = x13, x = factor(n), fill = method)) + 
  geom_bar(position = "dodge", stat = "identity") + 
  labs(y = "Median rank of x13", x = "Sample size",
       fill = "Method") +
  facet_grid(compl ~ R2) + 
  scale_y_continuous(limits = c(0,15)) +
  scale_fill_manual(values = colors7a) +
  theme2 +
  theme(legend.position = "none")

combined <- cowplot::plot_grid(p3, p13, ncol = 1)
cowplot::plot_grid(combined, legend_ranks, rel_widths = c(1, 0.2))




## Prediction
### RMSPE

## long format for ggplots
df_pred_final_l <- res_A |> 
  select(1:6, RMSPE_final, calA_final, calB_final) |> 
  mutate(calB_wins = ifelse(calB_final < 0.01, 0.01, calB_final),
         AD_slope = abs(log(calB_wins)),       # absolute deviation of log slope from target 0
         SD_slope = log(calB_wins)^2) |>       # squared distance of log slope from target 0
  group_by(compl, R2, n, method) |> 
  summarize(RMSPE = mean(RMSPE_final, na.rm = TRUE),
            MCSE = sd(RMSPE_final, na.rm = TRUE)/sqrt(500),
            calA = median(calA_final, na.rm = TRUE),
            calB = median(calB_final, na.rm = TRUE),
            MAD  = median(AD_slope, na.rm = TRUE),
            RMSD = sqrt(mean(SD_slope, na.rm = TRUE)),
            .groups = "drop")



ggplot(res_A |> filter(method != "MFP"),
       aes(y = RMSPE_final, x = factor(n), fill = method)) +
  geom_jitter(
    alpha = 0.2, size = 0.5, aes(color = method),
    position = position_jitterdodge(jitter.width = 0.45, dodge.width = 0.8)) +
  stat_summary(
    fun = median,
    geom = "errorbar",
    color = "black",
    aes(ymax = after_stat(y), ymin = after_stat(y)),
    position = position_dodge(width = 0.8),
    width = 0.6, linewidth = 0.25) +
  facet_grid(compl ~ R2, scales = "free_y",
             labeller = labeller(R2 = function(x) paste0("R² = ", x))) +
  scale_fill_manual(values = colors7b) +
  scale_color_manual(values = colors7b) +
  labs(x = "Sample size", y = "RMSPE", color = "Method", fill = "Method") +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 3))) +
  theme2




### Calibration 
#### Intercept
ggplot(df_pred_final_l,
       aes(y = calA, x = n, color = method, linetype = method)) +
  geom_point(size = 1.6) +
  geom_line(linewidth = 0.77) + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_grid(compl ~ R2, scales = "fixed",
             labeller = labeller(R2 = function(x) paste0("R² = ", x))) +
  scale_color_manual(values = colors8a) +
  scale_linetype_manual(values = ltype8a) +
  scale_x_continuous(breaks = c(100, 250, 500, 1000),
                     limits = c(90, 1050)) +
  labs(x = "Sample size", y = "Median calibration intercept",
       color = "Method", linetype = "Method") +
  theme1


#### Slope
ggplot(df_pred_final_l,
       aes(y = calB, x = n, color = method, linetype = method)) +
  geom_point(size = 1.6) +
  geom_line(linewidth = 0.77) + 
  geom_hline(yintercept = 1, linetype = "dashed") +
  facet_grid(compl ~ R2, scales = "fixed",
             labeller = labeller(R2 = function(x) paste0("R² = ", x))) +
  scale_color_manual(values = colors8a) +
  scale_linetype_manual(values = ltype8a) +
  scale_x_continuous(breaks = c(100, 250, 500, 1000),
                     limits = c(90, 1050)) +
  labs(x = "Sample size", y = "Median calibration slope",
       color = "Method", linetype = "Method") +
  theme1


#### MAD 
ggplot(df_pred_final_l |> filter(method != "MFP"),
       aes(y = MAD, x = n, color = method, linetype = method)) +
  geom_point(size = 1.6) +
  geom_line(linewidth = 0.77) + 
  facet_grid(compl ~ R2, scales = "fixed",
             labeller = labeller(R2 = function(x) paste0("R² = ", x))) +
  scale_color_manual(values = colors7b) +
  scale_linetype_manual(values = ltype7b) +
  scale_x_continuous(breaks = c(100, 250, 500, 1000),
                     limits = c(90, 1050)) +
  labs(x = "Sample size", y = "MAD",
       color = "Method", linetype = "Method") +
  theme1


#### RMSD
ggplot(df_pred_final_l |> filter(method != "MFP"),
       aes(y = RMSD, x = n, color = method, linetype = method)) +
  geom_point(size = 1.6) +
  geom_line(linewidth = 0.77) + 
  facet_grid(compl ~ R2, scales = "fixed",
             labeller = labeller(R2 = function(x) paste0("R² = ", x))) +
  scale_color_manual(values = colors7b) +
  scale_linetype_manual(values = ltype7b) +
  scale_x_continuous(breaks = c(100, 250, 500, 1000),
                     limits = c(90, 1050)) +
  labs(x = "Sample size", y = "RMSD",
       color = "Method", linetype = "Method") +
  theme1

### Correlation
df_shrinkage <- res_A |>
  select(method:rep, calB_final, calB_full) |> 
  ## winserize at 0.01 to avoid problems of negative slopes
  mutate(final = ifelse(calB_final < 0.01, 0.01, calB_final),
         full = ifelse(calB_full < 0.01, 0.01, calB_full)) |> 
  select(-calB_final, -calB_full) |> 
  pivot_wider(names_from = method, 
              values_from = c(final, full)) |> 
  select(-full_MFP, -full_ENET, -full_GBM) |> 
  mutate(LMSS_e = log(final_LMSS) - log(full_LMSS),
         ENET_e = log(final_ENET) - log(full_LMSS),
         GBM_e = log(final_GBM) - log(full_LMSS),
         GBT_e = log(final_GBT) - log(full_GBT),
         RFB_e = log(final_RFB) - log(full_RFB),
         RFH_e = log(final_RFH) - log(full_RFH),
         TOM_e = log(final_TOM) - log(full_TOM))

df_shrinkage_cor <- df_shrinkage |>
  group_by(compl, R2, n) |> 
  summarise(LMSS_cor = cor(-log(full_LMSS), LMSS_e, method = "spearman"),
            ENET_cor = cor(-log(full_LMSS), ENET_e, method = "spearman"),
            GBM_cor = cor(-log(full_LMSS), GBM_e, method = "spearman", use = "complete.obs"),
            GBT_cor = cor(-log(full_GBT), GBT_e, method = "spearman", use = "complete.obs"),
            RFB_cor = cor(-log(full_RFB), RFB_e, method = "spearman", use = "complete.obs"),
            RFH_cor = cor(-log(full_RFH), RFH_e, method = "spearman", use = "complete.obs"),
            TOM_cor = cor(-log(full_TOM), TOM_e, method = "spearman", use = "complete.obs"),)




## Sparse models
### TIF

## first, transfer to long format and filter 4 best ranked variables
### regression-based methods
df_ranks_Br_4 <- res_Br |> 
  select(c(method:converged, ends_with("rank"))) |> 
  filter(converged == TRUE) |> 
  pivot_longer(cols = x1_rank:x15_rank, 
               names_to = "variable", values_to = "rank",
               ## remove everything from the first underscore (removes "rank" and dummy-coding)
               names_transform = list(variable = ~ sub("_.*$", "", .))
  ) |> 
  filter(rank <= 4)

df_ranks_Bo_4 <- res_Bo |> 
  select(c(method:converged, ends_with("rank"))) |> 
  filter(converged == TRUE) |> 
  pivot_longer(cols = x1_rank:x11_rank, 
               names_to = "variable", values_to = "rank",
               ## remove everything from the first underscore (removes "rank" and dummy-coding)
               names_transform = list(variable = ~ sub("_.*$", "", .))
  ) |> 
  filter(rank <= 4)

### tree-based methods
df_ranks_Bt_4 <- res_Bt |> 
  select(c(method:converged, ends_with("rank"))) |> 
  filter(converged == TRUE) |> 
  pivot_longer(cols = x1_rank:x15_rank, 
               names_to = "variable", values_to = "rank",
               names_transform = list(variable = ~ sub("_.*$", "", .))) |> 
  filter(rank <= 4)

## join together
df_ranks_4 <- bind_rows(df_ranks_Br_4, df_ranks_Bt_4)

## record the model size to correctly calculate TIF when <4 vars selected
df_ranks_4_size <- df_ranks_4 |> 
  group_by(method, compl, R2, n, rep) |> 
  summarise(size = n(), .groups = "drop") |> 
  mutate(size = ifelse(size>4, 4, size)) ## one rep where MFP 6 variables on same rank

## filter predictor variables and count rows per repetition - that is number of predictors in top 4
df_ranks_4_TI <- df_ranks_4 |> 
  filter(variable %in% preds) |> 
  group_by(method, compl, R2, n, rep) |> 
  summarise(TI = n(), .groups = "drop") 

## join together and calculate TIF across rep
df_TIF_sparse <- left_join(df_ranks_4_size, df_ranks_4_TI, 
                           by = c("method", "compl", "R2", "n", "rep")) |> 
  mutate(TIF = TI/size) |> 
  group_by(compl, R2, n, method) |> 
  summarize(TIF = mean(TIF, na.rm = TRUE), .groups = "drop")


TIF_4_table <- df_ranks_4 |> 
  select(method, variable) |> 
  group_by(method, variable) |> 
  summarise(count = n()/24000)


ggplot(df_TIF_sparse,
       aes(y = TIF, x = n, color = method, linetype = method)) +
  geom_point(size = 1.6) +
  geom_line(linewidth = 0.77) + 
  facet_grid(compl ~ R2, scales = "fixed", 
             labeller = labeller(R2 = function(x) paste0("R² = ", x))) +
  scale_color_manual(values = colors7a) +
  scale_linetype_manual(values = ltype7a) +
  scale_x_continuous(breaks = c(100, 250, 500, 1000),
                     limits = c(90, 1050)) +
  labs(x = "Sample size", y = "Average TIF of top 4 selected variables",
       color = "Method", linetype = "Method") +
  theme1


### RMSPE
## long format for ggplots
df_pred_sparse_l <- res_A |> 
  select(1:6, RMSPE_min, calA_min, calB_min) |> 
  group_by(compl, R2, n, method) |> 
  summarize(RMSPE = mean(RMSPE_min, na.rm = TRUE),
            calA = median(calA_min, na.rm = TRUE),
            calB = median(calB_min, na.rm = TRUE),
            .groups = "drop")

ggplot(res_A |> filter(method != "MFP"),
       aes(y = RMSPE_min, x = factor(n), fill = method)) +
  geom_jitter(
    alpha = 0.2, size = 0.5, aes(color = method),
    position = position_jitterdodge(jitter.width = 0.45, dodge.width = 0.8)) +
  stat_summary(
    fun = median,
    geom = "errorbar",
    color = "black",
    aes(ymax = after_stat(y), ymin = after_stat(y)),
    position = position_dodge(width = 0.8),
    width = 0.6, linewidth = 0.25) +
  facet_grid(compl ~ R2, scales = "free_y",
             labeller = labeller(R2 = function(x) paste0("R² = ", x))) +
  scale_fill_manual(values = colors7b) +
  scale_color_manual(values = colors7b) +
  labs(x = "Sample size", y = "RMSPE", color = "Method", fill = "Method") +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 3))) +
  theme2


### Slope
ggplot(df_pred_sparse_l,
       aes(y = calB, x = n, color = method, linetype = method)) +
  geom_point(size = 1.6) +
  geom_line(linewidth = 0.77) + 
  geom_hline(yintercept = 1, linetype = "dashed") +
  facet_grid(compl ~ R2, scales = "fixed",
             labeller = labeller(R2 = function(x) paste0("R² = ", x))) +
  scale_color_manual(values = colors8a) +
  scale_linetype_manual(values = ltype8a) +
  scale_x_continuous(breaks = c(100, 250, 500, 1000),
                     limits = c(90, 1050)) +
  scale_y_continuous(breaks = c(0.5, 0.75, 1, 1.25),
                     limits = c(0.5, 1.25)) +
  labs(x = "Sample size", y = "Median calibration slope",
       color = "Method", linetype = "Method") +
  theme1

### Comparison
df_pred_sparse_comparison_l <- left_join(
  df_pred_final_l |> rename(RMSPE_final = RMSPE), 
  df_pred_sparse_l |> rename(RMSPE_sparse = RMSPE),
  by = c("compl", "R2", "n", "method")) |> 
  mutate(RMSPE_diff = RMSPE_final - RMSPE_sparse)

ggplot(df_pred_sparse_comparison_l |> filter(method != "MFP"),
       aes(y = RMSPE_diff, x = n, color = method, linetype = method)) +
  geom_point(size = 1.6) +
  geom_line(linewidth = 0.77) + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_grid(compl ~ R2, scales = "free_y",
             labeller = labeller(R2 = function(x) paste0("R² = ", x))) +
  scale_color_manual(values = colors7b) +
  scale_linetype_manual(values = ltype7b) +
  scale_x_continuous(breaks = c(100, 250, 500, 1000),
                     limits = c(90, 1050)) +
  labs(x = "Sample size", y = "Average difference in RMSPE",
       color = "Method", linetype = "Method") +
  theme1


## Oracle models
## add true oracle model stored in res_A to other oracles in res_Aom
res_A_oms <- res_A |> 
  filter(method == "TOM") |> 
  select(method:rep, RMSPE_final, calA_final, calB_final) |> 
  rename(RMSPE = RMSPE_final,
         calA = calA_final,
         calB = calB_final) |> 
  bind_rows(res_Aom) |> 
  mutate(method = factor(method, levels = oracles))

df_pred_oms_l <- res_A_oms |> 
  group_by(compl, R2, n, method) |> 
  summarize(RMSPE = mean(RMSPE, na.rm = TRUE),
            calA = median(calA, na.rm = TRUE),
            calB = median(calB, na.rm = TRUE),
            .groups = "drop")

### RMSPE
ggplot(res_A_oms |> filter(method != "MFP-O"),
       aes(y = RMSPE, x = factor(n), fill = method)) +
  geom_jitter(
    alpha = 0.2, size = 0.5, aes(color = method),
    position = position_jitterdodge(jitter.width = 0.45, dodge.width = 0.8)) +
  stat_summary(
    fun = median,
    geom = "errorbar",
    color = "black",
    aes(ymax = after_stat(y), ymin = after_stat(y)),
    position = position_dodge(width = 0.8),
    width = 0.6, linewidth = 0.25) +
  facet_grid(compl ~ R2, scales = "free_y",
             labeller = labeller(R2 = function(x) paste0("R² = ", x))) +
  scale_fill_manual(values = colors7b) +
  scale_color_manual(values = colors7b) +
  labs(x = "Sample size", y = "RMSPE", color = "Method", fill = "Method") +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 3))) +
  theme2

### Calibration slope
ggplot(df_pred_oms_l,
       aes(y = calB, x = n, color = method, linetype = method)) +
  geom_point(size = 1.6) +
  geom_line(linewidth = 0.77) + 
  geom_hline(yintercept = 1, linetype = "dashed") +
  facet_grid(compl ~ R2, scales = "fixed",
             labeller = labeller(R2 = function(x) paste0("R² = ", x))) +
  scale_color_manual(values = colors8a) +
  scale_linetype_manual(values = ltype8a) +
  scale_x_continuous(breaks = c(100, 250, 500, 1000),
                     limits = c(90, 1050)) +
  scale_y_continuous(breaks = c(0.5, 0.75, 1, 1.25),
                     limits = c(0.5, 1.28)) +
  labs(x = "Sample size", y = "Median calibration slope",
       color = "Method", linetype = "Method") +
  theme1


### Comparison RMPSE
mapping <- data.frame(
  m = methods8,
  method_O = c(oracles, "RF-O")
)

df_pred_final_l_m <- df_pred_final_l |> 
  rename(RMSPE_final = RMSPE,
         calA_final = calA,
         calB_final = calB) |> 
  left_join(mapping, by = c("method" = "m"))

df_pred_comparison_l <- left_join(
  df_pred_final_l_m, 
  df_pred_oms_l |> rename(method_O = method,
                          RMSPE_O = RMSPE,
                          calA_O = calA,
                          calB_O = calB),
  by = c("compl", "R2", "n", "method_O")) |> 
  mutate(RMSPE_diff = RMSPE_O - RMSPE_final,
         calA_diff  = calA_O  - calA_final,
         calB_diff  = calB_O  - calB_final,
         RMSPE_diff_p = RMSPE_diff/RMSPE_O*100,
         calA_diff_p  = calA_diff/calA_O*100,
         calB_diff_p  = calB_diff/calB_O*100,
         method = factor(method, levels = methods8))

ggplot(df_pred_comparison_l |> filter(method_O != "MFP-O" & method_O != "TOM"),
       aes(y = RMSPE_diff, x = n, color = method, linetype = method)) +
  geom_point(size = 1.6) +
  geom_line(linewidth = 0.77) + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_grid(compl ~ R2, scales = "free_y",
             labeller = labeller(R2 = function(x) paste0("R² = ", x))) +
  scale_color_manual(values = colors6a) +
  scale_linetype_manual(values = ltype6) +
  scale_x_continuous(breaks = c(100, 250, 500, 1000),
                     limits = c(90, 1050)) +
  labs(x = "Sample size", y = "Average difference in RMSPE",
       color = "Method", linetype = "Method") +
  theme1



## Duration
## transform results data to wide format
df_duration <- res_A |>
  select(1:6, duration) |> 
  pivot_wider(names_from = method, values_from = duration, id_cols = c(compl, R2, n, rep))

df_duration_long <- res_A |> 
  select(1:6, duration) |> 
  filter(method != "TOM", !is.na(duration))


ggplot(df_duration_long, aes(y = as.numeric(duration), x = factor(n), fill = method)) + 
  geom_boxplot(outliers = FALSE) + 
  labs(y = "Computational time [minutes]", x = "Sample size") +
  facet_grid(compl ~ R2) + 
  scale_fill_manual(values = colors7a) +
  scale_y_log10(breaks = c(0.1, 1, 10, 100), labels = c("0.1", "1", "10", "100")) +
  theme(plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"),
        panel.grid.major.y = element_line(color = "#e5e5e5", linewidth = 0.5),
        panel.grid.major.x = element_blank(),
        axis.line = element_line(color = "black"),
        panel.spacing.y = unit(0.5, "cm")) 
