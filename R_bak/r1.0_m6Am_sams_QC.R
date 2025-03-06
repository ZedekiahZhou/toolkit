# Date: Feb 11, 2025
# Author: Zhou Zhe
# motif QC for sams

rm(list=ls())
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(scales)
source("~/toolkit/funs2.R")
mydir = "04_m6Am/"
bca = paste0(c("G", "C", "T"), "CA")

# I. load data ===========================
file = "03_Sites/merged_m6Am.m6Am.passed"
m6A = data.table::fread(file, na.strings = ".") %>%
  data.frame() %>%
  fill_na_as_0() %>%
  dplyr::rename(Gene = geneID) %>%
  mutate(ID = paste(Chr, Pos, sep = "."))
summary(m6A[, grep("Passed_", colnames(m6A))])
sams = sub("AGcov_", "", grep("AGcov_", colnames(m6A), value = T))

## motif
motif = data.table::fread(paste0(file, "_len5_motif.bed")) %>%
  data.frame() %>%
  mutate(ID = paste(V1, V2+3, sep = "."), 
         motif = toupper(V7)) %>%
  dplyr::select(ID, motif)
m6A = left_join(m6A, motif, by = "ID")
m6A$BCA = ifelse(stringr::str_sub(m6A$motif, 1, 3) %in% bca, "BCA", "Non-BCA")
table(m6A$BCA)

## loc info
loc = data.table::fread(paste0(file, "_dist.measures.txt")) %>%
  data.frame() %>%
  mutate(rel_location2 = rel_location)
utr5.SF <- median(loc$utr5_size, na.rm = T)/median(loc$cds_size, na.rm = T)
utr3.SF <- median(loc$utr3_size, na.rm = T)/median(loc$cds_size, na.rm = T)
loc$rel_location2[loc$rel_location < 1] = rescale(loc$rel_location[loc$rel_location < 1], 
                                                  to = c(1-utr5.SF, 1), from = c(0,1))
loc$rel_location2[loc$rel_location >= 2] = rescale(loc$rel_location[loc$rel_location >= 2], 
                                                   to = c(2, 2+utr3.SF), from = c(2,3))


# II. QC =====================================================
pdf(paste0(mydir, "QC_sams_motif.pdf"), width = 6, height = 6)
for (sam in sams) {
  df = m6A[m6A[[paste0("Passed_", sam)]], ]
  
  title = sam
  freq_BCA = table(df$BCA)/nrow(df)*100
  freq_BCA = signif(freq_BCA["BCA"], digits = 3)
  title = paste0(title, " #", nrow(df), " BCA_", freq_BCA, "%")
  
  motif_freq = table(df$motif) %>% 
    sort() %>%
    data.frame() %>%
    mutate(Prop = Freq/sum(Freq)*100) %>%
    slice_tail(n = 15) %>%
    mutate(color = ifelse(stringr::str_sub(Var1, 1, 3) %in% bca, "black", "red"))
  
  motif_color = colorRampPalette(RColorBrewer::brewer.pal(7, "YlGnBu"))(15)
  p1 = ggplot(motif_freq, aes(y = Var1, x = Prop, fill = Var1)) + 
    geom_bar(stat = "identity", show.legend = F) + 
    theme_classic() + 
    scale_fill_manual(values = motif_color) + 
    scale_x_continuous(expand = c(0, 0), breaks = seq(0, 20, 2), position = "top") + 
    ylab("Context 5-mers") + xlab("Percentage (%)") + 
    theme(axis.text.y = element_text(color = motif_freq$color))
  
  tmpdf = df %>% filter(motif %in% motif_freq$Var1)
  tmpdf$motif = factor(tmpdf$motif, levels = motif_freq$Var1)
  p2 = ggplot(tmpdf, 
              aes(y = motif, x = .data[[paste0("Ratio_", sam)]])) + 
    geom_violin(aes(fill = motif), show.legend = F) + 
    scale_fill_manual(values = motif_color) + 
    scale_x_continuous(expand = c(0,2), position = "top") + 
    geom_boxplot(width = 0.1, outlier.shape = NA) + 
    ylab("Context 5-mers") + xlab("m6Am level (%)") + 
    theme_classic() + 
    theme(axis.text.y = element_text(color = motif_freq$color))
  
  print(p1 + p2 + plot_annotation(title))
  
  # prepare input for meme
  write(paste0(">", df$ID, "\n", df$motif), 
        file = paste0(mydir, "motif/Motif_", sam, ".fa"))
  
  # location 
  df = left_join(df, loc, join_by("Chr"=="chr", "Pos"=="coord"))
  minloc = min(df$rel_location2, na.rm = T)
  maxloc = max(df$rel_location2, na.rm = T)
  p3 = ggplot(df, aes(rel_location)) + 
    geom_density(adjust = 0.6) + 
    theme_classic() + 
    xlab(NULL) + ylab("Density") + 
    scale_x_continuous(breaks = c(0.5, 1.5, 2.5), 
                       labels = c("UTR5", "CDS", "UTR3")) + 
    geom_vline(xintercept = 1, linetype = "dashed") + 
    geom_vline(xintercept = 2, linetype = "dashed")
  p4 = ggplot(df, aes(rel_location2)) + 
    geom_density(adjust = 0.6) + 
    theme_classic() + 
    xlab(NULL) + ylab("Density") + 
    scale_x_continuous(breaks = c(0.5+minloc/2, 1.5, maxloc/2+1), 
                       labels = c("UTR5", "CDS", "UTR3")) + 
    geom_vline(xintercept = 1, linetype = "dashed") + 
    geom_vline(xintercept = 2, linetype = "dashed")
  
  onecolor = ggsci::pal_npg()(10)[6]
  p5 = ggplot(df, aes(.data[[paste0("Ratio_", sam)]])) + 
    geom_density(color = onecolor) + 
    xlab(NULL) + ylab("Density") + 
    scale_x_continuous(limits = c(0, 105), breaks = seq(0, 100, 10)) + 
    theme_classic()
  p6 = ggplot(df, aes(y = "", x = .data[[paste0("Ratio_", sam)]])) + 
    geom_boxplot(outlier.shape = NA, fill = onecolor) + 
    scale_x_continuous(limits = c(0, 105), breaks = seq(0, 100, 10), position = "top") + 
    theme_classic() + 
    xlab("m6Am level (%)") + ylab("") + 
    theme(axis.line.y = element_blank(), 
          axis.ticks.y = element_blank(), 
          axis.text.y = element_blank()) 
  p7 = wrap_plots(list(p6, p5), nrow = 2, heights = c(1, 5))
  print((p7 | plot_spacer()) / (p3 | p4) + plot_annotation(title))
  
}
dev.off()
