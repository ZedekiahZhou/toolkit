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
        theme_classic()
}


# plot for feature of overlapped elements in a list
corplot_lofd = function(x, # list of data.frame
                          idx1 = NULL, idx2 = NULL, # index or name of element to plot correlation
                          col_merge = "ID", # columns to merge
                          col_plot, # columns to plot
                          log = FALSE, return_df = FALSE) {
    require(ggplot2)
    if (is.null(idx1)) {
        idx1 = names(x)[1]
        idx2 = names(x)[2]
    }
    
    df = dplyr::inner_join(x[[idx1]], x[[idx2]], suffix = c(idx1, idx2), by = col_merge)
    if (log) {
        print(ggplot(df, mapping = aes(x = .data[[paste0(col_plot, idx1)]], y = .data[[paste0(col_plot, idx2)]])) + 
                  geom_point(size = 1) + ggpubr::stat_cor()) + scale_x_log10() + scale_y_log10()
    } else {
        print(ggplot(df, mapping = aes(x = .data[[paste0(col_plot, idx1)]], y = .data[[paste0(col_plot, idx2)]])) + 
                  geom_point(size = 1) + ggpubr::stat_cor())
    }
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
    plot(euler(ltmp), quantities = TRUE)
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
        x[[idx1]] == x[[idx1]][used, ]
        x[[idx2]] == x[[idx2]][used, ]
    }
    
    df = rbind(data.frame(value = x[[idx1]][[col_plot]], group = idx1), data.frame(value = x[[idx2]][[col_plot]], group = idx2))
    if (log) {
        print(ggplot(df, mapping = aes(x = value, color = group)) + geom_density() + scale_x_log10())
    } else {
        print(ggplot(df, mapping = aes(x = value, color = group)) + geom_density())
    }
    if (return_df) {
        return(df)
    }
}
