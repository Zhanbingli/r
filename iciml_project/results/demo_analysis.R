# demo_pipeline.R  -----------------------------------------------------------
library(readxl); library(tidyverse); library(stringr)
install.packages("UpSetR")
library(pheatmap); library(ggplot2); library(ggrepel); library(UpSetR)
dir.create("results", showWarnings = FALSE)

## 1. read data --------------------------------------------------------------
fn     <- "data/NormalizedData.xlsx"
sheet  <- "VSN_NormalyzedData_Mapped"
df     <- read_excel(fn, sheet = sheet)

ann_cols <- c("Input Name","Standardized.name","Formula","Exact.mass",
              "Super.class","Main.class","Sub.class",
              "PubChem_CID","ChEBI_ID","HMDB_ID","LM_ID","KEGG_ID",
              "INCHI_KEY","RefMet_ID")
sample_cols <- setdiff(names(df), ann_cols)

# parse group from column name (…_sampleP_pos…)
grp_raw <- str_match(sample_cols, "sample([A-Za-z0-9]+)_")[,2]
stopifnot(!any(is.na(grp_raw)))
group <- factor(grp_raw, levels = unique(grp_raw))          # P, PP, V

## 2. detection matrix (≥2/3 rule) ------------------------------------------
mat      <- as.matrix(df[, sample_cols])
mode(mat) <- "numeric"

# per-group replicate counts that are non-NA
det_by_rep <- !is.na(mat)
det_by_grp <- sapply(levels(group), function(g)
  rowSums(det_by_rep[, group == g, drop = FALSE]) >= 2)

det_mat <- as.data.frame(det_by_grp * 1)   # 1/0
write.csv(det_mat, "results/detection_matrix.csv", row.names = FALSE)

## 3. counts & percentages ---------------------------------------------------
tot_feat   <- nrow(det_mat)
per_group  <- colSums(det_mat)
common_all <- rowSums(det_mat) == nlevels(group)
unique_tbl <- apply(det_mat, 1, function(x){
  if(sum(x)==1) names(which(x==1)) else NA_character_
}) |> na.omit()

counts_tbl <- tibble(
  Group        = names(per_group),
  Detected_n   = per_group,
  Percent      = round(100*per_group/tot_feat,1)
)
counts_tbl <- add_row(counts_tbl,
                      Group      = "Common_3groups",
                      Detected_n = sum(common_all),
                      Percent    = round(100*sum(common_all)/tot_feat,1)
)
write.csv(counts_tbl, "results/counts_overview.csv", row.names=FALSE)

## 4. UpSet / overlap figure -------------------------------------------------
up_input <- as.data.frame(det_mat*1)
png("results/Upset_overlap.png", 1200, 700)
upset(up_input, nsets = 3, mb.ratio = c(0.6,0.4),
      sets = colnames(up_input),
      order.by = "freq")
dev.off()

## 5. PCA of samples (for QC) ------------------------------------------------
library(stats)
mat_imp <- replace(mat, is.na(mat), NA)             # keep NAs
mat_sc  <- t(scale(t(mat_imp), center = TRUE, scale = TRUE))
pc      <- prcomp(t(mat_sc), na.action = na.omit)
## --- PCA (samples) : keep only complete rows -------------------------------
mat_pca <- mat_sc[complete.cases(mat_sc), ]

pc <- prcomp(t(mat_pca), center = FALSE, scale. = FALSE)  # 已经手动 scale 过
pc_df <- data.frame(pc$x[,1:2],
                    group  = group,
                    sample = sample_cols)

library(ggplot2); library(ggrepel)
p <- ggplot(pc_df, aes(PC1, PC2, colour = group, label = sample)) +
  geom_point(size = 4) +
  geom_text_repel() +
  theme_bw() +
  labs(title = "PCA – samples (complete-case features)")
ggsave("results/PCA_samples.png", p, width = 6, height = 5, dpi = 300)


## 6. One differential contrast demo (PP vs P) -------------------------------
keep_feat <- (det_mat$PP + det_mat$P) >= 2          # present in ≥2 farms
expr      <- mat[keep_feat, c(group %in% c("P","PP"))]
expr      <- expr[, match(sample_cols[group %in% c("P","PP")],
                          colnames(expr))]

# design & limma
library(limma)
g_sub <- factor(group[group %in% c("P","PP")], levels = c("P","PP"))
design <- model.matrix(~0 + g_sub); colnames(design) <- levels(g_sub)
fit    <- eBayes(lmFit(expr, design))
cont   <- makeContrasts(PPvsP = PP - P, levels = design)
fit2   <- eBayes(contrasts.fit(fit, cont))
tt     <- topTable(fit2, coef = 1, number = Inf, adjust.method = "BH")

res <- bind_cols(df[keep_feat, ann_cols], tt)
write.csv(res, "results/DE_PPvsP.csv", row.names = FALSE)

# Volcano
vol <- ggplot(tt, aes(logFC, -log10(adj.P.Val))) +
  geom_point(alpha=.6) +
  geom_vline(xintercept = c(-1,1), linetype="dashed") +
  geom_hline(yintercept = -log10(0.05), linetype="dashed") +
  theme_bw() + labs(title="Volcano – PP vs P")
ggsave("results/Volcano_PPvsP.png", vol, width=6, height=5, dpi=300)

# Heatmap (Top 30)
top30 <- tt |> arrange(adj.P.Val) |> slice(1:30) |> rownames()
pheatmap(expr[top30,], scale="row",
         annotation_col = data.frame(group = g_sub),
         show_rownames = FALSE,
         main = "Top 30 – PP vs P",
         filename = "results/Heatmap_PPvsP.png",
         width = 6, height = 8)
# 取 adj.P.Val 最小的前 30 行对应的“行号”：
top_idx <- order(tt$adj.P.Val)[1:30]   # 这里得到的是整数向量

pheatmap(expr[top_idx, ],
         scale           = "row",
         annotation_col  = data.frame(group = g_sub, row.names = colnames(expr)),
         show_rownames   = FALSE,
         main            = "Top 30 – PP vs P",
         filename        = "results/Heatmap_PPvsP.png",
         width = 6, height = 8)

# Save top20 table
top20 <- res %>%                     # res 已经包含注释列
  arrange(adj.P.Val) %>% 
  slice_head(n = 20) %>%
  select(Standardized.name, Main.class, Sub.class,
         KEGG_ID, HMDB_ID, logFC, adj.P.Val)

write.csv(top20, "results/DE_PPvsP_top20.csv", row.names = FALSE)


message("Demo run complete → see results/")
# ---------------------------------------------------------------------------#
