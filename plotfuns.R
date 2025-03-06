require(dplyr)
require(eulerr)
require(ggplot2)
require(ggpointdensity)
require(patchwork)
require(methylKit)
require(MASS)

pden_color = rev(RColorBrewer::brewer.pal(11, "RdBu"))

get_density <- function(x, y, ...) {
    dens <- MASS::kde2d(x, y, ...)
    ix <- findInterval(x, dens$x)
    iy <- findInterval(y, dens$y)
    ii <- cbind(ix, iy)
    return(dens$z[ii])
}

# fill na as 0
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


# write out as Bismark format ----
write2Bis <- function(df, x_sams, outdir = "../03_Sites") {
    for (sam in x_sams) {                                                                                               
        used_col = c("Chr", "Start", "End", "Strand", "Ratio", "AGcov", "Acov")
        bis = df %>% mutate(
            Start = Pos - 1,
            End = Pos,
            AGcov = .data[[paste0("AGcov_", sam)]],
            Acov = .data[[paste0("Acov_", sam)]],
            Ratio = Acov/AGcov*100
        )
        bis = bis[, used_col]
        write.table(bis, file = paste0(outdir, sam, ".cov"), sep = "\t",
                    row.names = F, col.names = F, quote = F)
    }
}

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


# differential
fun_diff <- function(x_sams, treatment, name1 = "Treat", name0 = "Ctrl",
                     outdir = "../03_Sites", 
                     mincov = 10,
                     info_col = c("Chr", "Pos", "Strand", "Gene"),
                     df_meta, 
                     q_cutoff = 0.05, diffRatio_cutoff = 10
) {
    files = file.path(outdir, paste0(x_sams, ".cov"))
    myobj = methRead(as.list(files), sample.id = as.list(x_sams), assembly = "hg38", 
                     treatment = treatment, mincov = mincov, header = F, 
                     pipeline = list("fraction"=FALSE, "chr.col"=1, "start.col"=2, "end.col"=3,
                                     "coverage.col"=6, "strand.col"=4, "freqC.col"=5))
    meth=unite(myobj, destrand=FALSE)
    myDiff=calculateDiffMeth(meth)
    myDiff = left_join(data.frame(myDiff), 
                       df_meta[, c(info_col, unlist(lapply(x_sams, grep, x = colnames(df_meta), value = T)))], 
                       by = join_by("chr"=="Chr", "end"=="Pos", "strand"=="Strand"))
    myDiff = myDiff %>% mutate(
        Alteration = ifelse(qvalue >= q_cutoff, "NonSig", 
                            ifelse(meth.diff >= diffRatio_cutoff, "Up", 
                                   ifelse(meth.diff <= -diffRatio_cutoff, "Down", "NonSig"))), 
        label = ifelse(Alteration == "NonSig", "", Gene), 
    )
    
    myDiff[[paste0("Ratio_", name1)]] = rowMeans(myDiff[, paste0("Ratio_", x_sams[treatment == 1]), drop = F])
    myDiff[[paste0("Ratio_", name0)]] = rowMeans(myDiff[, paste0("Ratio_", x_sams[treatment == 0]), drop = F])
    myDiff[[paste0("Passed_", name1)]] = apply(myDiff[, paste0("Passed_", x_sams[treatment == 1]), drop = F], 1, all)
    myDiff[[paste0("Passed_", name0)]] = apply(myDiff[, paste0("Passed_", x_sams[treatment == 0]), drop = F], 1, all)
    myDiff[[paste0("Passed_either")]] = myDiff[, paste0("Passed_", name1)] | myDiff[, paste0("Passed_", name0)]
    myDiff[[paste0("Passed_both")]] = myDiff[, paste0("Passed_", name1)] & myDiff[, paste0("Passed_", name0)]
    return(myDiff %>% filter(Passed_either))
}


# euler plot of sites of two samples
fun_euler <- function(df, sam1, sam0) {
    ltmp = list(rownames(df)[df[[paste0("Passed_", sam1)]]], rownames(df)[df[[paste0("Passed_", sam0)]]])
    names(ltmp) = c(sam1, sam0)
    p = plot(euler(ltmp), quantities = TRUE)
    p$vp$height = unit(0.4, "npc")
    p$vp$height = unit(0.4, "npc")
    print(p)
}


# point 2d density 
fun_pden2d <- function(df, sam1, sam0, dotsize = 0.3,
                     label = FALSE, adjust = 0.03, plot_count = FALSE, main = "", 
                     use_common = FALSE, split_coding = FALSE, useMASS = FALSE, MASS_para_n = 100
) {
    
    
    get_plot = function(dfplot) {
        if (use_common) dfplot = dfplot[dfplot$Passed_both, ]
        p1 = ggplot(dfplot, 
                    aes(y = .data[[paste0("Ratio_", sam1)]], 
                        x = .data[[paste0("Ratio_", sam0)]]))
        if (useMASS) {
            dfplot$density = get_density(dfplot[[paste0("Ratio_", sam1)]], 
                                         dfplot[[paste0("Ratio_", sam0)]], n = MASS_para_n)
            p1 = p1 + geom_point(data = dfplot, aes(color = density), size = dotsize) 
        } else {
            p1 = p1 + geom_pointdensity(size = dotsize, adjust = adjust)
        }
        p1 = p1 + 
            geom_abline(slope = 1, linetype = "dashed", color = "black") + 
            geom_abline(slope = 1, intercept = 10, linetype = "dashed", color = "red4") + 
            geom_abline(slope = 1, intercept = -10, linetype = "dashed", color = "red4") +
            scale_color_gradientn(colors = pden_color) + 
            theme_classic() + xlim(0, 100) + ylim(0, 100)
        if (label) {
            p1 = p1 + ggrepel::geom_text_repel(mapping = aes(label = label), min.segment.length = 0)
        }
        if (plot_count) {
            main = paste0(main, " -- ", "Up_", sum(dfplot$Alteration == "Up"), " | ", 
                              "Down_", sum(dfplot$Alteration == "Down"))
        }
        return(p1 + ggtitle(main))
    }
        
    if (split_coding) {
        print(get_plot(df %>% filter(txBiotype %in% c("protein_coding"))) + 
                  get_plot(df %>% filter(txBiotype %in% c("snRNA", "snoRNA"))))
    } else {
        print(get_plot(df))
    }
}


# plot 1d density
fun_pden1d <- function(df, sam1, sam0, 
    use_respective = FALSE, use_common = FALSE, width = 0.001,
    main = ""
) {
    mycolor = rev(RColorBrewer::brewer.pal(11, "RdBu"))
    if (use_common) df = df[df$Passed_both, ]
    dfplot = rbind(
        data.frame(Ratio = df[[paste0("Ratio_", sam1)]], 
                   Passed = df[[paste0("Passed_", sam1)]], 
                   group = sam1),
        data.frame(Ratio = df[[paste0("Ratio_", sam0)]], 
                   Passed = df[[paste0("Passed_", sam0)]], 
                   group = sam0)
    )
    if (use_respective) dfplot = dfplot[dfplot$Passed, ]
    dfplot$group = factor(dfplot$group, levels = c(sam1, sam0))
    
    p1 = ggplot(dfplot, aes(x = Ratio, color = group)) + 
        geom_density() + 
        geom_boxplot(width = width, outlier.shape = NA) + 
        scale_color_brewer(palette = "Set1") + 
        theme_classic() + xlim(0, 100) + 
        ggtitle(main)
    print(p1)
}


# plot volcano
fun_volcano <- function(df, ylim = c(0,30)) {
    tmp_main = paste0("Up_", sum(df$Alteration == "Up"), " | ", 
                      "Down_", sum(df$Alteration == "Down"))
    p1 = ggplot(df, aes(meth.diff, -log10(qvalue), color = Alteration)) + 
        geom_point() + theme_classic() + coord_cartesian(ylim= ylim) + 
        geom_vline(xintercept = 10, linetype = "dashed", color = "red4") + 
        geom_vline(xintercept = -10, linetype = "dashed", color = "red4") + 
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red4") +
        scale_color_manual(values = c("Up" = "red4", "Down" = "blue4", "NonSig" = "gray")) + 
        ggtitle(tmp_main)
    print(p1)
}