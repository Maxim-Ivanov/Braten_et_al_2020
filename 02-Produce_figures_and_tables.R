library(tidyverse)
library(scales)

# Load the haplotype data:
data_dir <- "." # change to the directory with HAP/LEGEND/SAMPLE files which were returned by "01-Preprocess_1000G_haplotypes.sh"
hap <- file.path(data_dir, "CYP2C19_10Mb_ALL.hap.gz") %>% read_delim(delim = " ", col_names = FALSE, col_types = cols(.default = col_integer())) %>% as.matrix() # 285864 x 5008
legend <- file.path(data_dir, "CYP2C19_10Mb_ALL.legend.gz") %>% read_delim(delim = " ")
legend_id <- legend %>% unite(id, id, position, a0, a1, sep = ":") %>% .$id
rownames(hap) <- legend_id

# Skip *2-*n haplotypes:
blacklist <- read_tsv("CYP2C19_blacklist.txt", col_names = FALSE) # this file is available from the same GitHub page
bad_rows <- which(legend$id %in% blacklist$X1)
bad_cols <- which(colSums(hap[bad_rows, ]) > 0) # n = 2076
hap_star1 <- hap[-bad_rows, -bad_cols] # 285847 x 2932

# Classify *1 haplotypes into TG, TA, CA and CG:
m1 <- legend_id[legend$id == "rs2860840"] # ref = C, alt = T
m2 <- legend_id[legend$id == "rs11188059"] # ref = G, alt = A
hap_m1 <- hap_star1[rownames(hap_star1) == m1, ]
hap_m2 <- hap_star1[rownames(hap_star1) == m2, ]
cg <- hap_m1 == 0 & hap_m2 == 0
ca <- hap_m1 == 0 & hap_m2 == 1
tg <- hap_m1 == 1 & hap_m2 == 0
grp <- ifelse(cg, "CG", ifelse(ca, "CA", ifelse(tg, "TG", "TA")))

# Skip both mutations from the haplotype table:
hap_star1 <- hap_star1[!rownames(hap_star1) %in% c(m1, m2), ] # 285845 x 2932

# Skip non-variable positions:
rmn <- rowMeans(hap_star1)
hap_star1 <- hap_star1[rmn != 0 & rmn != 1, ] # 220112 x 2932

### Define custom functions ----------------------------------------------------------------------------------

calculate_padj <- function(hap, grp, g1, g2, method = "fisher") {
  stopifnot(method %in% c("fisher", "chi"))
  # Extract two relevant haplotypes:
  h1 <- hap[, grp == g1]
  h2 <- hap[, grp == g2]
  message(ncol(h1), " ", g1, " vs ", ncol(h2), " ", g2, " haplotypes;")
  # Count frequency of each SNP in these two groups:
  k1 <- rowSums(h1) %>% as.integer()
  n1 <- ncol(h1)
  freq1 <- round(k1 / n1, 3)
  k2 <- rowSums(h2) %>% as.integer()
  n2 <- ncol(h2)
  freq2 <- round(k2 / n2, 3)
  # Calculate p-values:
  if (method == "fisher") {
    pvals <- mapply(function(k1, n1, k2, n2) { fisher.test(x = rbind(c(k1, n1 - k1), c(k2, n2 - k2)))$p.value }, k1, n1, k2, n2, SIMPLIFY = FALSE) %>% unlist()
  } else if (method == "chi") {
    pvals <- mapply(function(k1, n1, k2, n2) { chisq.test(x = c(k1, n1 - k1), y = c(k2, n2 - k2))$p.value }, k1, n1, k2, n2, SIMPLIFY = FALSE) %>% unlist()
  }
  padj <- p.adjust(pvals, method = "BH")
  coord <- rownames(hap) %>% str_split(pattern = ":") %>% lapply(`[`, 2) %>% unlist() %>% as.numeric()
  out <- tibble("id" = rownames(hap), "coord" = coord, "pval" = pvals, "padj" = padj, "mllk" = -log10(padj), "n1" = n1, "k1" = k1, "freq1" = freq1, 
                "n2" = n2, "k2" = k2, "freq2" = freq2, "dir" = paste0("Enriched in ", ifelse(freq1 > freq2, g1, g2)))
  return(out)
}

get_results <- function(hap, grp, genotypes = list(c("TA", "TG"), c("TA", "CG"), c("TG", "CG")), method = "fisher") {
  tbl_list <- vector("list", length(genotypes))
  for (i in seq_along(genotypes)) {
    gt <- genotypes[[i]]
    tbl <- calculate_padj(hap, grp, g1 = gt[[1]], g2 = gt[[2]], method = method)
    tbl$gt <- paste(gt[[1]], "vs", gt[[2]], sep = " ")
    tbl_list[[i]] <- tbl
  }
  out <- bind_rows(tbl_list)
  return(out)
}

### Draw 3-panel Manhattan plot (Supplementary Figure S1) -----------------------------------------------------

ymax <- 200
xmin <- 94000000
xmax <- 99000000
cbp <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") # colorblind-friendly palette


tbl <- get_results(hap_star1, grp)
tbl$mllk <- ifelse(tbl$mllk > ymax, ymax, tbl$mllk)
tbl$gt <- tbl$gt %>% as.factor() %>% relevel(ref = "TA vs TG")

p <- ggplot(tbl, aes(x = coord, y = mllk, colour = dir)) + geom_point(size = 1) + 
  scale_colour_manual(values = cbp[c(1, 6, 7)]) + 
  scale_x_continuous(labels = format_format(big.mark = ".", decimal.mark = ",", scientific = FALSE), 
                     limits = c(xmin, xmax)) + 
  xlab("chr10") + 
  ylab("-log10(padj)") + 
  ylim(0, ymax) + 
  facet_grid(gt ~ .) + 
  theme_bw() + 
  theme(legend.position = "top", 
        legend.title = element_blank(), 
        text = element_text(size = 12), 
        axis.text.x = element_text(size = 8), 
        axis.text.y = element_text(size = 8))

ggsave("Fig_S2.tiff", plot = p, width = 180, height = 150, units = "mm", dpi = 300)


##### Save the tables ------------------------------------------------------------------------------------------

# 1) Big table of all SNPs in all comparisons:
t1 <- calculate_padj(hap_star1, grp, g1 = "TA", g2 = "TG") %>% select(id, freq1, freq2, mllk) %>% mutate(mllk = abs(round(mllk, 3))) %>% rename(freq_TA = freq1, freq_TG = freq2, mllk_TA_TG = mllk)
t2 <- calculate_padj(hap_star1, grp, g1 = "TA", g2 = "CG") %>% select(freq2, mllk) %>% mutate(mllk = abs(round(mllk, 3))) %>% rename(freq_CG = freq2, mllk_TA_CG = mllk)
t3 <- calculate_padj(hap_star1, grp, g1 = "TG", g2 = "CG") %>% select(mllk) %>% mutate(mllk = abs(round(mllk, 3))) %>% rename(mllk_TG_CG = mllk)
big <- bind_cols(t1, t2, t3) %>% select(id, contains("freq"), everything())
write_tsv(big, "All_variable_SNPs_in_10Mb_window.txt.gz") # for convenience, this table was uploaded to the same GitHub page

# 2) Small table of the most significant SNPs in TA vs TG comparison (Supplementary File S2):
small <- big %>% filter(mllk_TA_TG >= 50) %>% select(id, freq_TA, freq_TG, mllk_TA_TG) %>% mutate(mllk_TA_TG = round(mllk_TA_TG, 3), enrichment = round(freq_TA / freq_TG, 3))
write_tsv(small, "File_S2.txt")
