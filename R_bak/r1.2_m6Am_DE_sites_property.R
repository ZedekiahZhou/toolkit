rm(list=ls())
library(ggplot2)
library(dplyr)
library(patchwork)
library(scales)

mydir = "04_m6Am/"
bca = paste0(c("G", "C", "T"), "CA")

# motif info
motif = data.table::fread(paste0("03_Sites/merged_m6Am.m6Am.passed_len5_motif.bed")) %>%
  data.frame() %>%
  mutate(ID = paste(V1, V2+3, sep = "."), 
         motif = toupper(V7)) %>%
  dplyr::select(ID, motif)

## loc info
loc = data.table::fread(paste0("03_Sites/merged_m6A.m6A.passed_dist.measures.txt")) %>%
  data.frame() %>%
  mutate(rel_location2 = rel_location)
utr5.SF <- median(loc$utr5_size, na.rm = T)/median(loc$cds_size, na.rm = T)
utr3.SF <- median(loc$utr3_size, na.rm = T)/median(loc$cds_size, na.rm = T)
loc$rel_location2[loc$rel_location < 1] = rescale(loc$rel_location[loc$rel_location < 1], 
                                                  to = c(1-utr5.SF, 1), from = c(0,1))
loc$rel_location2[loc$rel_location >= 2] = rescale(loc$rel_location[loc$rel_location >= 2], 
                                                   to = c(2, 2+utr3.SF), from = c(2,3))

for (case in c("FTO.KO", "CS1", "FB23.2")) {
  ctrl = "CTRL"
  dfshort = read.delim(paste0(mydir, "DEtable_", case, ".tsv")) %>%
    mutate(ID = paste(seqnames, end, sep = "."))
  dfshort = left_join(dfshort, motif, by = "ID") %>% 
    mutate(BCA = ifelse(stringr::str_sub(motif, 1, 3) %in% bca, "BCA", "Non-BCA"))
  table(dfshort$BCA)
  
  ## ratio
  dflong = reshape2::melt(dfshort, measure.vars = paste0("Ratio_", c(case, ctrl)), 
                          id.vars = c("seqnames", "end", "Alteration"), 
                          value.name = "Ratio") %>%
    mutate(Group = factor(sub("Ratio_", "", variable), levels = c(ctrl, case)))
  
  alter_color = c(ggsci::pal_npg()(10)[c(1,4)], "gray")
  names(alter_color) = c("Up", "Down", "NonSig")
  
  pdf(paste0(mydir, "DEsites_", case, "_vs_", ctrl, ".pdf"), width = 9, height = 4)
  p1 = ggplot(dflong, aes(Ratio, color = Alteration)) + 
          geom_density(aes(linetype = Group)) + 
          scale_color_manual(values = alter_color) + 
          scale_x_continuous(breaks = seq(0, 100, 10)) + 
          ylab("Density") + 
          theme_classic()
  
  group_color = c(ggsci::pal_npg(alpha = 0.8)(10)[c(1,4)])
  names(group_color) = c(case, ctrl)
  p2 = ggplot(dflong, aes(x = Alteration, Ratio, fill = Group)) + 
          geom_boxplot(outlier.shape = NA) + 
          scale_fill_manual(values = group_color) + 
          scale_y_continuous(breaks = seq(0, 100, 10)) + 
          ylab("m6A level (%)") + 
          theme_classic()
  
  tbl_alter = table(dfshort$Alteration)
  title = paste0(case, " vs CTRL: Up #", tbl_alter["Up"], 
                 " / NonSig #", tbl_alter["NonSig"], " / Down #", tbl_alter["Down"])
  print(p1 + p2 + plot_annotation(title))
  
  lp = lapply(c("Up", "NonSig", "Down"), function(group) {
    df = dfshort %>% filter(Alteration == group)
    
    title = group
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
    p = ggplot(motif_freq, aes(y = Var1, x = Prop, fill = Var1)) + 
      geom_bar(stat = "identity", show.legend = F) + 
      theme_classic() + 
      scale_fill_manual(values = motif_color) + 
      scale_x_continuous(expand = c(0, 0), breaks = seq(0, 20, 2), position = "top") + 
      ylab("Context 5-mers") + xlab("Percentage (%)") + 
      theme(axis.text.y = element_text(color = motif_freq$color)) + 
      ggtitle(title)
    
    # prepare input for meme
    write(paste0(">", df$ID, "\n", df$motif), 
          file = paste0(mydir, "motif/Motif_", case, "_vs_", ctrl, "_", group, ".fa"))
    return(p)
  })
  print(wrap_plots(lp) + plot_annotation(paste0(case, " vs ", ctrl)))
  
  # location
  dfshort = left_join(dfshort, loc, join_by("seqnames"=="chr", "end"=="coord"))
  minloc = min(dfshort$rel_location2, na.rm = T)
  maxloc = max(dfshort$rel_location2, na.rm = T)
  alter_color = c(ggsci::pal_npg()(10)[c(1,4)], "gray")
  names(alter_color) = c("Up", "Down", "NonSig")
  p3 = ggplot(dfshort, aes(rel_location, color = Alteration)) + 
    geom_density(adjust = 0.6) + 
    theme_classic() + 
    scale_x_continuous(breaks = c(0.5, 1.5, 2.5), 
                       labels = c("UTR5", "CDS", "UTR3")) + 
    scale_color_manual(values = alter_color) + 
    geom_vline(xintercept = 1, linetype = "dashed") + 
    geom_vline(xintercept = 2, linetype = "dashed")
  p4 = ggplot(dfshort, aes(rel_location2, color = Alteration)) + 
    geom_density(adjust = 0.6) + 
    theme_classic() + 
    scale_x_continuous(breaks = c(0.5+minloc/2, 1.5, maxloc/2+1), 
                       labels = c("UTR5", "CDS", "UTR3")) + 
    scale_color_manual(values = alter_color) + 
    geom_vline(xintercept = 1, linetype = "dashed") + 
    geom_vline(xintercept = 2, linetype = "dashed")
  print(p3 + p4) + plot_annotation(paste0(case, " vs ", ctrl))
  
  dev.off()
}

