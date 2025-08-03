# ----------------- 0 载入包 -----------------
if(!requireNamespace("pacman")) install.packages("pacman")
pacman::p_load(readxl, dplyr, stringr, tidyr, httr, jsonlite)

# ----------------- 1 读数据 -----------------
fn    <- "data/NormalizedData.xlsx"
sheet <- "VSN_NormalyzedData_Mapped"
raw   <- read_excel(fn, sheet = sheet)

ann_cols <- c("Input Name","Standardized.name","Formula","Exact.mass",
              "Super.class","Main.class","Sub.class",
              "PubChem_CID","ChEBI_ID","HMDB_ID","LM_ID","KEGG_ID")
anno <- raw[, ann_cols]

# ----------------- 2 Level-1：已有 ID -----------------
l1 <- anno %>% 
  filter(!is.na(LM_ID) | !is.na(HMDB_ID) | !is.na(KEGG_ID)) %>% 
  mutate(confidence = "L1")

# ----------------- 3 Level-4：仅类信息 -----------------
l4 <- anno %>% 
  filter(row_number() %in% setdiff(row_number(), row_number(l1))) %>% 
  mutate(confidence = "L4")

annot <- bind_rows(l1, l4)

# ----------------- 4 缺省检查 -----------------
missing_stats <- annot %>% 
  mutate(has_KEGG = !is.na(KEGG_ID),
         has_HMDB = !is.na(HMDB_ID),
         has_LM   = !is.na(LM_ID)) %>% 
  summarise(across(starts_with("has_"), ~sum(!.x)),
            .groups = "drop")

print(missing_stats)

# ----------------- 5 输出 -----------------
dir.create("results", showWarnings = FALSE)
write.csv(annot, "results/annotated_lipids.csv", row.names = FALSE)
write.csv(missing_stats, "results/ID_missing_summary.csv", row.names = FALSE)

message("annotated_lipids.csv + 缺省统计已写入 results/")
