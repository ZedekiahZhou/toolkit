require(dplyr)
require(eulerr)
require(ggplot2)
require(ggpointdensity)
require(patchwork)
require(methylKit)
require(MASS)

# pden_color = rev(RColorBrewer::brewer.pal(11, "RdBu"))
pden_color = rev(RColorBrewer::brewer.pal(11, "RdYlBu"))

# I. preprocessing ==========
## fill na as 0
fill_na_as_0 <- function(df, x_sams = NULL) {
  if (is.null(x_sams)) {  
    x_sams = colnames(df)[grepl("AGcov_", colnames(df))]
    x_sams = sub("AGcov_", "", x_sams)
  }
  for (sam in x_sams) {
    df[, paste0("AGcov_", sam)] = ifelse(is.na(df[, paste0("AGcov_", sam)]), 
                                         0, df[, paste0("AGcov_", sam)])
    df[, paste0("Acov_", sam)] = ifelse(is.na(df[, paste0("Acov_", sam)]), 
                                        0, df[, paste0("Acov_", sam)])
    df[, paste0("Ratio_", sam)] = ifelse(df[, paste0("AGcov_", sam)] == 0, 0, 
                                         df[, paste0("Acov_", sam)]/df[, paste0("AGcov_", sam)]*100)
  }
  return(df)
}


## select samples, return a clean data frame
select_sams <- function(
    df, x_sams,
    info_col = c("Chr", "Pos", "Strand", "Gene", "txBiotype"), 
    addtional_col = FALSE, 
    treatment = NULL,
    min_samples_per_group = NULL
) {
  info_col = info_col[info_col %in% colnames(df)]
  data_col = paste0(c("AGcov_", "Acov_", "Ratio_", "Passed_"), 
                    rep(x_sams, each = 4))
  info_df = df[, colnames(df) %in% info_col]
  data_df = df[, colnames(df) %in% data_col]
  

  # add additional columns of group means
  if (addtional_col) {
    if (is.null(treatment)) {
      stop("treatment must be specified if addtional_col is TRUE")
    }
    if (is.null(min_samples_per_group)) {
      message("Require sites passed in all samples of group1 or group0 if min_samples_per_group is NULL!")
      min_samples_per_group = table(treatment)
    }
    group1 = unique(treatment)[1]
    group0 = unique(treatment)[2]
    info_df[[paste0("Ratio_", group1)]] = rowMeans(data_df[, paste0("Ratio_", x_sams[treatment == group1]),
                                                           drop = F])
    info_df[[paste0("Ratio_", group0)]] = rowMeans(data_df[, paste0("Ratio_", x_sams[treatment == group0]), 
                                                           drop = F])
    info_df[[paste0("PassedSam_", group1)]] = apply(data_df[, paste0("Passed_", x_sams[treatment == group1]), 
                                                         drop = F], 1, sum)
    info_df[[paste0("PassedSam_", group0)]] = apply(data_df[, paste0("Passed_", x_sams[treatment == group0]), 
                                                         drop = F], 1, sum)
    info_df$Passed.Either = info_df[, paste0("PassedSam_", group1)] >= min_samples_per_group[group1] | 
      info_df[, paste0("PassedSam_", group0)] >= min_samples_per_group[group0]
    info_df$Passed.Both = info_df[, paste0("PassedSam_", group1)] >= min_samples_per_group[group1] & 
      info_df[, paste0("PassedSam_", group0)] >= min_samples_per_group[group0]
  }

  info_df$Passed.Any = apply(data_df[, paste0("Passed_", x_sams)], 1, any)
  df = cbind(info_df, data_df)
  return(df %>% filter(Passed.Any))
}


## write out bismark coverage format (no strand info)
## chr  start end ratio(%)  methyled  unmethyled
write2BisCoverage <- function(df, x_sams, outdir = "../03_Sites") {
  for (sam in x_sams) {
    used_col = c("Chr", "Start", "End", "Ratio", "Acov", "Gcov")
    # all positions are 1-based
    bis = df %>% mutate(
      Start = Pos,
      End = Pos,
      AGcov = .data[[paste0("AGcov_", sam)]],
      Acov = .data[[paste0("Acov_", sam)]],
      Gcov = AGcov - Acov,
      Ratio = Acov/AGcov*100
    )
    bis = bis[, used_col]
    write.table(bis, file = paste0(outdir, sam, ".biscov"), sep = "\t",
                row.names = F, col.names = F, quote = F)
  }
}


## QC between replicates
qc_reps = function(df, x_sams, fout) {
  pdf(fout, height = 5, width = 5)
  df = select_sams(m6A, x_sams)
  lsites = lapply(x_sams, function(x) {
    c(1:nrow(df))[df[[paste0("Passed_", x)]]]
  })
  names(lsites) = x_sams
  
  # plot(euler(lsites), quantities = T, fill = ggsci::pal_npg(alpha = 0.6)(length(x_sams)))
  p = plot(euler(lsites), quantities = T, fill = RColorBrewer::brewer.pal(length(x_sams), "Set2"))
  p$vp$width = unit(0.7, "npc")
  p$vp$height = unit(0.7, "npc")
  print(p)
  
  for (i in combn(1:length(x_sams), 2, simplify = FALSE)) {
    x = x_sams[i[1]]
    y = x_sams[i[2]]
    df_plot = df[df[[paste0("Passed_", x)]] & df[[paste0("Passed_", y)]], ]
    p = ggplot(df_plot, aes(.data[[paste0("Ratio_", x)]], .data[[paste0("Ratio_", y)]])) + 
      geom_pointdensity(adjust = 5, size = 0.6) + 
      scale_color_gradientn(colors = pden_color) + 
      geom_abline(slope = 1, intercept = 0, linetype = "dashed") + 
      scale_x_continuous(breaks = seq(10, 100, 30)) + 
      scale_y_continuous(breaks = seq(10, 100, 30)) + 
      coord_fixed() + 
      xlab(paste0("m6A level in ", x)) + 
      ylab(paste0("m6A level in ", y)) + 
      ggpubr::stat_cor() + 
      theme_bw() + 
      theme(panel.grid = element_blank())
    print(p)
  }
  
  dev.off()
}