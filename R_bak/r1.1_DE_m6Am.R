# 筛选条件：在**任意一组**2/3个样本中被鉴定为m6Am的点；在**两组**的2/3个样本中都有AGcov>=15

rm(list=ls())
library(dplyr)
library(ggplot2)
library(ggpointdensity)
library(methylSig)
library(RColorBrewer)
source("~/toolkit/funs2.R")
mydir = "Duo/04_m6Am/"


# I. load data ===========================
file = "Duo/03_Sites/merged_m6Am.m6Am.passed"
m6A = data.table::fread(file, na.strings = ".") %>%
  data.frame() %>%
  fill_na_as_0() %>%
  dplyr::rename(Gene = geneID)
summary(m6A[, grep("Passed_", colnames(m6A))])
sams = sub("AGcov_", "", grep("AGcov_", colnames(m6A), value = T))
write2BisCoverage(m6A, x_sams = sams, outdir = mydir)




# II. compare reps =======================
sams
qc_reps(m6A, sams[1:3], fout = paste0(mydir, "QC_reps_AsPC.1.pdf"))
qc_reps(m6A, sams[4:6], fout = paste0(mydir, "QC_reps_HPNE.pdf"))


# III. select samples to perform DE =======================
sams
cutoff_fdr = 0.05
cutoff_diff = 15
for (case in c("AsPC.1")) {
  # change index of ctrl samples
  ctrl = "HPNE"
  x_sams = sams[c(grep(case, sams), 4:6)]
  treatment = c(rep(case, 3), rep(ctrl, 3))
  
  required_sams = c(2, 2)
  names(required_sams) = c(case, ctrl)
  
  
  pdf(paste0(mydir, "DEplot_", case, ".pdf"), width = 4, height = 4)
  ## pre-filter sites passed in 2 samples of 3 (EITHER group)
  cand = select_sams(m6A, x_sams, addtional_col = T, 
                     treatment = treatment, min_samples_per_group = required_sams) %>%
    filter(Passed.Either)
  
  ## euler of passed sites
  lsites = lapply(c(case, ctrl), function(x) {
    c(1:nrow(cand))[cand[[paste0("PassedSam_", x)]] >= 2]
  })
  names(lsites) = c(case, ctrl)
  p = plot(euler(lsites), quantities = T, 
           fill = RColorBrewer::brewer.pal(3, "Set2")[2:3], 
           main = "Sites passed in 2 of 3 reps\n(Either Group)")
  p$vp$width = unit(0.7, "npc")
  p$vp$height = unit(0.7, "npc")
  print(p)
  
  G_cand = GenomicRanges::GRanges(
    seqnames = cand$Chr, 
    ranges = IRanges::IRanges(
      start = cand$Pos, end = cand$Pos
    )
  )
  
  ## filter sites covered in 2 samples of 3 (BOTH group)
  bs_raw = bsseq::read.bismark(
    files = paste0(mydir, "/", x_sams, ".biscov"),
    loci = G_cand, 
    colData = data.frame(Type= treatment,row.names = x_sams),
    rmZeroCov = FALSE,
    strandCollapse = FALSE
  )
  bs = filter_loci_by_coverage(bs_raw, min_count = 15, max_count = Inf) # only set counts below 15 to 0
  bs = filter_loci_by_group_coverage(
    bs = bs,
    group_column = 'Type', 
    min_samples_per_group = required_sams)
  
  # methylSig chi-square
  res_diff_sig = diff_methylsig(
    bs = bs,
    group_column = 'Type',
    comparison_groups = c('case' = case, 'control' = ctrl),
    disp_groups = c('case' = TRUE, 'control' = TRUE),
    local_window_size = 0,
    t_approx = FALSE,
    n_cores = 4)
  df_diff_sig = data.frame(res_diff_sig) %>% mutate(
    Alteration = ifelse(fdr < cutoff_fdr,
                        ifelse(
                          meth_diff > cutoff_diff, "Up",
                          ifelse(meth_diff < -cutoff_diff, "Down", "NonSig")
                        ), "NonSig"),
    .after = fdr
  )
  table(df_diff_sig$Alteration)
  df_diff_sig = inner_join(df_diff_sig, cand, 
                           by = join_by("seqnames"=="Chr", "end"=="Pos"))
  write.table(df_diff_sig, paste0(mydir, "DEtable_", case, ".tsv"), 
              quote = F, sep = "\t", row.names = F)
  
  ## Euler of used 
  lsites = lapply(c(case, ctrl), function(x) {
    c(1:nrow(df_diff_sig))[df_diff_sig[[paste0("PassedSam_", x)]] >= 2]
  })
  names(lsites) = c(case, ctrl)
  p = plot(euler(lsites), quantities = T, 
           fill = RColorBrewer::brewer.pal(3, "Set2")[2:3], 
           main = "Sites covered in 2 of 3 reps\n(Both Group, for DA)")
  p$vp$width = unit(0.7, "npc")
  p$vp$height = unit(0.7, "npc")
  print(p)
  
  
  ## Cor plot
  pden_color = rev(brewer.pal(11, "RdYlBu"))
  df_short = df_diff_sig
  p1 = ggplot(df_short, aes(.data[[paste0("Ratio_", ctrl)]], .data[[paste0("Ratio_", case)]])) + 
    geom_pointdensity(adjust = 3, size = 0.3) +
    geom_abline(slope = 1, linetype = "dashed") + 
    geom_abline(slope = 1, intercept = 15, linetype = "dashed", color = "red") + 
    geom_abline(slope = 1, intercept = -15, linetype = "dashed", color = "red") + 
    scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 20), expand = c(0.02, 0.02)) + 
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20), expand = c(0.02, 0.02)) + 
    coord_fixed() + 
    xlab(paste("m6Am level in", ctrl, "(%)")) + 
    ylab(paste("m6Am level in", case, "(%)")) + 
    scale_color_gradientn(colours = pden_color) + 
    theme_bw() + 
    theme(panel.grid = element_blank()) + 
    ggtitle(paste0(case, " vs ", ctrl, "\n", 
                   "Up_", sum(df_short$Alteration == "Up"), " | ", 
                   "NonSig_", sum(df_short$Alteration == "NonSig"), " | ", 
                   "Down_", sum(df_short$Alteration == "Down")))
  
  valcano_color = c(ggsci::pal_npg()(10)[c(1,4)], "gray")
  names(valcano_color) = c("Up", "Down", "NonSig")
  df_short = df_short %>% mutate(mlog10FDR = ifelse(-log10(fdr) > 20, 20, -log10(fdr)))
  p2 = ggplot(df_short, aes(meth_diff, mlog10FDR, color = Alteration)) + 
    geom_point(size = 0.3) + 
    coord_cartesian(ylim = c(0, 20), xlim = c(-60, 60)) + 
    theme_classic() + 
    geom_vline(xintercept = cutoff_diff, linetype = "dashed") + 
    geom_vline(xintercept = -cutoff_diff, linetype = "dashed") + 
    geom_hline(yintercept = -log10(cutoff_fdr), linetype = "dashed") + 
    scale_color_manual(values = valcano_color) + 
    xlab("Differential m6Am level (%)") + 
    ylab("-log10(FDR)") + 
    ggtitle(paste0(case, " vs ", ctrl, "\n", 
                   "Up_", sum(df_short$Alteration == "Up"), " | ", 
                   "NonSig_", sum(df_short$Alteration == "NonSig"), " | ", 
                   "Down_", sum(df_short$Alteration == "Down")))
  
  ecdf_color = c(ggsci::pal_npg()(10)[c(1,4)])
  names(ecdf_color) = c(case, ctrl)
  df_long = reshape2::melt(df_short, 
                           id.vars = c("seqnames", "end"), 
                           measure.vars = paste0("Ratio_", c(case, ctrl)), 
                           value.name = "Ratio") %>%
    mutate(Group = factor(sub("Ratio_", "", variable), levels = c(ctrl, case)))
  p3 = ggplot(df_long, aes(Ratio, color = Group)) + 
    stat_ecdf() + 
    xlab("m6Am level (%)") + ylab("Cumulative fraction") + 
    scale_color_manual(values = ecdf_color) + 
    theme_classic()
  
  box_color = c(ggsci::pal_npg(alpha = 0.8)(10)[c(1,4)])
  names(box_color) = c(case, ctrl)
  p4 = ggplot(df_long, aes(x = Group, y = Ratio, fill = Group)) + 
    geom_boxplot() + 
    scale_fill_manual(values = box_color) + 
    theme_classic() + 
    ylab("m6Am level (%)") + 
    ggpubr::stat_compare_means()
  
  
  print(p1)
  print(p2)
  print(p3)
  print(p4)
  dev.off()
}

