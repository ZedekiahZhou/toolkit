# for 2D density plot
get_density <- function(x, y, n = 1000, ...) {
        dens <- MASS::kde2d(x, y, n = n, ...)
        ix <- findInterval(x, dens$x)
        iy <- findInterval(y, dens$y)
        ii <- cbind(ix, iy)
        return(dens$z[ii])
}

plot_density_2d <- function(df, ax, ay, nbreaks = 1000) {
    require(ggplot2)
    mycolor = rev(RColorBrewer::brewer.pal(11, "RdYlBu"))
    
    df$density = get_density(df[[ax]], df[[ay]], n = nbreaks)
    ggplot(df, mapping = aes(x = .data[[ax]], y = .data[[ay]])) + 
        geom_point(aes(color = density), size = 0.3) + 
        scale_color_gradientn(colors = mycolor, limits = c(0, quantile(df$density, 0.98)), na.value = mycolor[11]) + 
        ggpubr::stat_cor() + theme_classic()
}


# plot for feature of overlapped elements in a list
corplot_lofd = function(x, # list of data.frame
                      idx1 = NULL, idx2 = NULL, # index or name of element to plot correlation
                      col_merge = "ID", # columns to merge
                      col_plot, # columns to plot
                      title = "", 
                      log = FALSE, return_df = FALSE) {
    require(ggplot2)
    if (is.null(idx1)) {
        idx1 = names(x)[1]
        idx2 = names(x)[2]
    }
    
    df = dplyr::inner_join(x[[idx1]], x[[idx2]], suffix = paste0(".", c(idx1, idx2)), by = col_merge)
    p = ggplot(df, mapping = aes(x = .data[[paste0(col_plot, ".", idx1)]], y = .data[[paste0(col_plot, ".", idx2)]])) + 
        geom_point(size = 1) + 
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red4") + 
        xlim(0,1) + ylim(0, 1) + 
        xlab(idx1) + ylab(idx2) + ggtitle(title) + 
        ggpubr::stat_cor() + theme_classic() + theme(legend.position = "None")
    if (log) {
        p = p + scale_x_log10() + scale_y_log10()
    }
    print(p)
    
    if (return_df) {
        return(df)
    }
}

density_2d_lofd = function(x, # list of data.frame
                        idx1 = NULL, idx2 = NULL, # index or name of element to plot correlation
                        col_merge = "ID", # columns to merge
                        col_plot, # columns to plot
                        nbreaks = 1000, title = "", 
                        axis_lim = NULL,
                        log = FALSE, return_df = FALSE) {
    require(ggplot2)
    if (is.null(idx1)) {
        idx1 = names(x)[1]
        idx2 = names(x)[2]
    }
    
    df = dplyr::inner_join(x[[idx1]], x[[idx2]], suffix = paste0(".", c(idx1, idx2)), by = col_merge)
    if (nrow(df) == 0) {
        print(paste("No overlap between", idx1, "and", idx2, "!"))
        return(NULL)
    }
    
    mycolor = rev(RColorBrewer::brewer.pal(11, "RdYlBu"))
    df$density = get_density(df[[paste0(col_plot, ".", idx1)]], df[[paste0(col_plot, ".", idx2)]], n = nbreaks)
    p = ggplot(df, mapping = aes(x = .data[[paste0(col_plot, ".", idx1)]], y = .data[[paste0(col_plot, ".", idx2)]])) + 
        geom_point(aes(color = density), size = 0.3) + 
        scale_color_gradientn(colors = mycolor, limits = c(0, quantile(df$density, 0.98)), na.value = mycolor[11]) + 
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red4") + 
        xlab(idx1) + ylab(idx2) + ggtitle(title) + 
        ggpubr::stat_cor() + theme_classic() + theme(legend.position = "None")
    if (!is.null(axis_lim)) {
        p = p + xlim(axis_lim[1], axis_lim[2]) + ylim(axis_lim[1], axis_lim[2])
    }
    if (log) {
        p = p + scale_x_log10() + scale_y_log10()
    }
    print(p)
    
    if (return_df) {
        return(df)
    }
}

euler_lofd = function(x, idx1 = NULL, idx2 = NULL, col_ID = "ID") {
    require(eulerr)
    if (is.null(idx1)) {
        idx1 = names(x)[1]
        idx2 = names(x)[2]
    }
    
    ltmp = list(x[[idx1]][[col_ID]], x[[idx2]][[col_ID]])
    names(ltmp) = c(idx1, idx2)
    p = plot(euler(ltmp), quantities = TRUE)
    p$vp$height = unit(0.7, "npc")
    p$vp$height = unit(0.7, "npc")
    print(p)
}

density_lofd = function(x, # list of data.frame
                        idx1 = NULL, idx2 = NULL, # index or name of element to plot correlation
                        col_merge = "ID", # columns to merge
                        col_plot, # columns to plot
                        log = FALSE, return_df = FALSE, overlap = FALSE) {
    require(ggplot2)
    if (is.null(idx1)) {
        idx1 = names(x)[1]
        idx2 = names(x)[2]
    }
    
    if (overlap) {
        used = intersect(x[[idx1]][[col_merge]], x[[idx2]][[col_merge]])
        x[[idx1]] = x[[idx1]][match(used, x[[idx1]][[col_merge]]), ]
        x[[idx2]] = x[[idx2]][match(used, x[[idx1]][[col_merge]]), ]
    }
    
    df = rbind(data.frame(value = x[[idx1]][[col_plot]], group = idx1), data.frame(value = x[[idx2]][[col_plot]], group = idx2))
    p = ggplot(df, mapping = aes(x = value, color = group)) + geom_density() + theme_classic()
    if (log) {
        p = p + scale_x_log10()
    }
    print(p)
    if (return_df) {
        return(df)
    }
}
