#!/usr/bin/env Rscript
library(ggplot2)
library(dplyr) 



################
#     MAIN     #
################

args <- commandArgs(TRUE)
p_dataset = args[1]
pr_out = args[2]


#pr_out = "./../../results/PFAS/"
#p_dataset = "./../../data/PFAS.csv"
d = read.csv(p_dataset, sep = ",")

d_count = table(d$Group)
d_count = cbind(names(d_count), d_count)
colnames(d_count) = c("Group", "Count")

d_count = as.data.frame(d_count)
d_count$Count = as.double(as.character(d_count$Count))

d_count <- d_count %>%
  arrange(desc(Group)) %>%
  mutate(lab.ypos = cumsum(Count) - 0.5*Count)


ggplot(d_count, aes(x = "", y = Count, fill = Group)) +
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar("y", start = 0)+
  geom_text(aes(y = lab.ypos, label = Count), color = "white")+
  theme_void()

ggsave(paste(pr_out, "pie_group.png", sep = ""), dpi=300, height = 10, width = 16)

