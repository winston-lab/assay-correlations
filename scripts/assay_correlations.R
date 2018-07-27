library(tidyverse)
library(forcats)
library(viridis)
library(GGally)

scatter = function(intable='coding-genes_WT-37C_allassays.tsv',
                   pcount=0.01,
                   anno_label='coding genes',
                   condition = 'WT-37C',
                   cutoff_high = 0.999,
                   cutoff_low = 0.01,
                   outpath = 'test.svg'){
    df = read_tsv(intable,
                  col_names=c('chrom', 'start', 'end', 'name', 'score', 'strand', 'assay')) %>%
        mutate(score=score/(end-start),
               assay = fct_inorder(assay, ordered=TRUE)) %>%
        group_by(assay) %>%
        mutate(score=log10(score+pcount)) %>%
        mutate(score=scale(score))

    maxsignal = df %>% pull(score) %>% quantile(probs=cutoff_high)
    minsignal = df %>% pull(score) %>% quantile(probs=cutoff_low)
    n_assays = n_distinct(df[["assay"]])

    df = df %>%
        select(-c(chrom, start, end, strand)) %>%
        spread(key=assay, value=score) %>%
        select(-name)

    mincor = min(cor(df, use="complete.obs")) * .99
    plots = list()

    #for each row
    for (i in 1:ncol(df)){
        #for each column
        for (j in 1:ncol(df)){
            idx = ncol(df)*(i-1)+j
            #upper right (correlation)
            if (i < j){
                c = cor(df[,i], df[,j], use = "complete.obs") %>% as.numeric()
                plot = ggplot(data = tibble(x=c(0,1), y=c(0,1), corr=c)) +
                        geom_rect(aes(fill=corr), xmin=0, ymin=0, xmax=1, ymax=1) +
                        annotate("text", x=0.5, y=0.5, label=sprintf("%.2f",round(c,2)),
                                 size=16/72*25.4*c, fontface="plain") +
                        scale_x_continuous(breaks=NULL) +
                        scale_y_continuous(breaks=NULL) +
                        scale_fill_distiller(palette="Blues", limits = c(mincor,1), direction=1)
                plots[[idx]] = plot
            }
            #top left to bot right diag (density)
            else if (i == j){
                subdf = df %>% select(i) %>% gather(sample, value)
                plot = ggplot(data = subdf, aes(x=(value+pcount))) +
                        geom_density(aes(y=..scaled..), fill="#114477", size=0.8) +
                        scale_y_continuous(breaks=c(0,.5,1)) +
                        scale_x_continuous(limits=c(minsignal, maxsignal)) +
                        annotate("text", x=maxsignal, y=0.9, hjust=1,
                                 label=unique(subdf$sample), size=10/72*25.4, fontface="plain")
                plots[[idx]] = plot
            }
            #bottom left (scatter)
            else {
                subdf = df %>% select(i,j) %>% gather(xsample, xvalue, -1) %>%
                            gather(ysample, yvalue, -c(2:3))
                plot = ggplot(data = subdf, aes(x=xvalue+pcount, y=yvalue+pcount)) +
                    geom_abline(intercept = 0, slope=1, color="grey80", size=1) +
                    stat_bin_hex(geom="point", aes(color=log10(..count..)),
                    # stat_bin_hex(geom="point", aes(color=..count..),
                                 binwidth=c(.075,.075), size=0.5, shape=16, stroke=0) +
                    scale_fill_viridis(option="inferno") +
                    scale_color_viridis(option="inferno") +
                    scale_x_continuous(limits = c(minsignal, maxsignal)) +
                    scale_y_continuous(limits = c(minsignal, maxsignal))
                plots[[idx]] = plot
            }
        }
    }

    mat = ggmatrix(plots, nrow=ncol(df), ncol=ncol(df),
                   title = paste0(nrow(df), " ", anno_label, ", ", condition),
                   xAxisLabels = names(df), yAxisLabels = names(df), switch="both") +
        theme_light() +
        theme(plot.title = element_text(size=12, color="black", face="bold"),
              axis.text = element_text(size=10, margin=margin(0,0,0,0,"pt")),
              strip.background = element_blank(),
              strip.text = element_text(size=12, color="black", face="bold"),
              strip.text.y = element_text(angle=180, hjust=1),
              strip.placement="outside")

    w = 3+n_assays*4
    h = 9/16*w+0.5
    ggsave(outpath, plot=mat, width=w, height=h, units="cm", limitsize=FALSE)
    return(mat)
}

plot = scatter(intable = snakemake@input[["table"]],
               pcount = snakemake@params[["pcount"]],
               anno_label = snakemake@params[["anno_label"]],
               condition = snakemake@wildcards[["condition"]],
               cutoff_high = snakemake@params[["cutoff_high"]],
               cutoff_low = snakemake@params[["cutoff_low"]],
               outpath = snakemake@output[["scatter"]])
