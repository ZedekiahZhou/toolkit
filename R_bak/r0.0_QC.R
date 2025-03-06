rm(list=ls())
library(dplyr)
library(ggplot2)
library(patchwork)
library(ggrepel)

# I. Count the number of multi-tags -----
mydir = "05_m6A/"
dir.create(mydir)
group = "m6A"
allsam = readLines(paste0("src/config_", group))
df <- plyr::ldply(allsam, function(sam) {
  if (group == "m6A") {
    res <- read.delim(paste0("src/log/pre_", sam, ".out"), skip = 37, nrows = 133)
  } else {
    res <- read.delim(paste0("src/log/pre_", sam, ".out"), skip = 38, nrows = 137)
  }
  
  res$Group = sam
  return(res)
})

df <- df %>% 
  group_by(Group) %>%
  mutate(Density=count/sum(count), 
         label = ifelse(Density == max(Density), Group, NA))
tag_stats <- df %>% group_by(Group) %>% 
  summarize(total_tag = sum(count), 
            one_tag = sum(count[length <= 27]),
            percent_1tag = one_tag/total_tag*100, 
            multi_tag = sum(count[length >27]), 
            percent_mtag = multi_tag/total_tag*100)
tag_stats <- tag_stats[order(tag_stats$Group), ]
write.table(tag_stats, file = paste0(mydir, "tag_stats.tsv"), row.names = F, sep = "\t")

ggplot(df, mapping = aes(x = length, y = Density, color = Group)) + 
  geom_line() + 
  geom_text_repel(aes(label = label), min.segment.length = 0, show.legend=FALSE) + 
  scale_x_continuous(breaks = c(seq(6, 150, 17))) + 
  geom_vline(xintercept = 27, color = "red4", linetype = "dashed") + 
  theme_classic() 
ggsave(paste0(mydir, "tag_length.pdf"), width = 7, height = 5)


# II. Reads distribution ------

## read mapping length
len_tf <- plyr::ldply(allsam, function(sam) {
  res = read.csv(paste0("03_Sites/QC/", sam, "_align_len.csv"))
  res$mapping_len_Prop <- res$mapping_len/sum(res$mapping_len)
  res$read_len_Prop <- res$read_len/sum(res$read_len)
  res$Group = sam
  return(res)
})

len_tf = len_tf %>% group_by(Group) %>%
  mutate(label = ifelse(mapping_len_Prop == max(mapping_len_Prop), Group, NA))

ggplot(len_tf, aes(X, mapping_len_Prop, color = Group)) +
  geom_line() + 
  scale_x_continuous(breaks = c(seq(0, 150, 30))) +
  # geom_vline(xintercept = 10, linetype = "dashed", color = "red4") + 
  xlab("Mapping length") + ylab("Proportion") + 
  geom_text_repel(aes(label = label), min.segment.length = 0, show.legend=FALSE)+
  theme_classic() + geom_vline(xintercept = 26, linetype = "dashed")
ggsave(paste0(mydir, "Distribution_of_mapping_length.pdf"), width = 7, height = 5)


## read distribution in transcriptome
distr_tf <- plyr::ldply(allsam, function(sam) {
  res = read.delim(paste0("03_Sites/QC/qualimap/", sam, "/raw_data_qualimapReport/coverage_profile_along_genes_(total).txt"))
  res$Density <- res$Transcript.coverage.profile/sum(res$Transcript.coverage.profile)
  res$Group = sam
  return(res)
})

distr_stats <- distr_tf %>% group_by(Group) %>%
  summarize(Prop05 = sum(Density[X.Transcript.position < 5]), 
            Prop10 = sum(Density[X.Transcript.position < 10]))
distr_tf$label = paste0(distr_tf$Group, " ",
                        signif(distr_stats$Prop10[match(distr_tf$Group, distr_stats$Group)]*100, 3), "%")
distr_tf = distr_tf %>% group_by(Group) %>%
  mutate(label2 = ifelse(Density == max(Density), label, NA))

label_data = subset(distr_tf, X.Transcript.position == 2)
ggplot(distr_tf, aes(X.Transcript.position, Density, color = label)) +
  geom_line() + 
  scale_x_continuous(breaks = c(seq(0, 100, 25), 10)) +
  geom_vline(xintercept = 10, linetype = "dashed", color = "red4") + 
  xlab("Quantile of Gene Body") + ylab("Density") + 
  geom_text_repel(aes(label = label2), min.segment.length = 0, show.legend=FALSE, max.overlaps = 100)+
  theme_classic() 
ggsave(paste0(mydir, "Distribution_of_reads.pdf"), width = 7, height = 5)
