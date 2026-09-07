# load library
library(tidyverse)
library(ggplot2)
library(reshape2)
library(svglite)

# load dataset
prvcont.logfc <- read.delim('result/fig4c/LBR/LBR_v_Cont_logFC.tsv', header = TRUE, sep = '\t')
prvcont.pvalue <- read.delim('result/fig4c/LBR/LBR_v_Cont_pval.tsv', header = TRUE, sep = '\t')
prvlbr.logfc <- read.delim('result/fig4c/LBR/LBR_v_PR_logFC.tsv', header = TRUE, sep = '\t')
prvlbr.pvalue <- read.delim('result/fig4c/LBR/LBR_v_PR_pval.tsv', header = TRUE, sep = '\t')

prvcont.logfc
# ---- inputs ----
logfc <- as.matrix(prvcont.logfc)
pval  <- as.matrix(prvcont.pvalue)

colnames(logfc) <- colnames(pval) <- c("24", "48", "72", "96")
rownames(logfc) <- rownames(pval) <- paste0("T", 1:8)

# ---- long format, merged so fill and stars stay aligned ----
lf <- melt(logfc); colnames(lf) <- c("Template", "Time", "logFC")
pv <- melt(pval);  colnames(pv) <- c("Template", "Time", "pvalue")
dat <- merge(lf, pv, by = c("Template", "Time"))

dat$Time     <- factor(dat$Time, levels = c("24", "48", "72", "96"))
dat$Template <- factor(dat$Template, levels = paste0("T", 1:8))

# stars from p-values
dat$stars <- ifelse(dat$pvalue < 0.001, "***",
                    ifelse(dat$pvalue < 0.01,  "**",
                           ifelse(dat$pvalue < 0.05,  "*", "")))

# symmetric limits so white sits exactly at 0
lim <- max(abs(dat$logFC), na.rm = TRUE)

lim_min <- min(dat$logFC, na.rm = TRUE)
lim_max <- max(dat$logFC, na.rm = TRUE)

lim_min
lim_max

p <- ggplot(dat, aes(x = Template,
                     y = factor(Time, levels = rev(c("24", "48", "72", "96"))),
                     fill = logFC)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = stars), color = "black", size = 7, vjust = 0.78) +
  scale_fill_gradient2(low = "#2E5F8C", mid = "white", high = "#E8704B",
                       midpoint = 0, limits = c(-lim, lim),
                       name = expression(log[2]*FC)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  labs(x = "Template", y = "Time (hpi)") +
  theme_void() +
  theme(
    legend.position   = "top",
    legend.direction  = "horizontal",
    legend.text       = element_text(size = 11),
    legend.title      = element_text(size = 13),
    legend.key.width  = unit(1.5, "cm"),
    legend.key.height = unit(0.5, "cm"),
    axis.text.x = element_text(color = "black", size = 15),
    axis.text.y = element_text(color = "black", size = 15),
    panel.grid  = element_blank()
  ) +
  coord_fixed(ratio = 1)
p

# p <- ggplot(dat, aes(x = Template,
#                      y = factor(Time, levels = rev(c("24", "48", "72", "96"))),
#                      fill = logFC)) +
#   geom_tile(color = "white", linewidth = 1) +
#   geom_text(aes(label = stars), color = "black", size = 7, vjust = 0.78) +
#   scale_fill_gradient2(low = "#2E5F8C", mid = "white", high = "#E8704B",
#                        midpoint = mean(range(dat$logFC, na.rm = TRUE)),
#                        limits   = range(dat$logFC, na.rm = TRUE),
#                        name     = expression(log[2]*FC)) +
#   scale_y_discrete(expand = c(0, 0)) +
#   scale_x_discrete(expand = c(0, 0)) +
#   labs(x = "Template", y = "Time (hpi)") +
#   theme_void() +
#   theme(
#     legend.position   = "top",
#     legend.direction  = "horizontal",
#     legend.text       = element_text(size = 11),
#     legend.title      = element_text(size = 13),
#     legend.key.width  = unit(1.5, "cm"),
#     legend.key.height = unit(0.5, "cm"),
#     axis.text.x = element_text(color = "black", size = 15),
#     axis.text.y = element_text(color = "black", size = 15),
#     panel.grid  = element_blank()
#   ) +
#   coord_fixed(ratio = 1)
# p
# 
# p <- ggplot(dat, aes(x = Template,
#                      y = factor(Time, levels = rev(c("24", "48", "72", "96"))),
#                      fill = logFC)) +
#   geom_tile(color = "white", linewidth = 1) +
#   geom_text(aes(label = stars), color = "black", size = 7, vjust = 0.78) +
#   scale_fill_gradient2(low = "#2E5F8C", mid = "white", high = "#E8704B",
#                        midpoint = 0, limits = c(lim_min, lim_max),
#                        name = expression(log[2]*FC)) +
#   scale_y_discrete(expand = c(0, 0)) +
#   scale_x_discrete(expand = c(0, 0)) +
#   labs(x = "Template", y = "Time (hpi)") +
#   theme_void() +
#   theme(
#     legend.position   = "top",
#     legend.direction  = "horizontal",
#     legend.text       = element_text(size = 11),
#     legend.title      = element_text(size = 13),
#     legend.key.width  = unit(1.5, "cm"),
#     legend.key.height = unit(0.5, "cm"),
#     axis.text.x = element_text(color = "black", size = 15),
#     axis.text.y = element_text(color = "black", size = 15),
#     panel.grid  = element_blank()
#   ) +
#   coord_fixed(ratio = 1)
# 
# p
ggsave(p, file = 'result/svg/lb10_v_cont.png', height = 6, width = 6, dpi = 300)
