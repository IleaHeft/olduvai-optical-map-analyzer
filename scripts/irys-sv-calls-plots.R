library(tidyverse)
library(ggbeeswarm)
library(ggsci)



svs <- read.delim("~/LabProjects/UCSF-visit/sv-calls-HLS-1000.txt", header = FALSE)
svs <- as_tibble(svs)
colnames(svs) <- c("gene","sample","type","size","support","min","max","range","nick_type","median.conf")

svs <- svs %>% mutate(size2 = ifelse(type == "del", size*-1,size),
                      type_group = ifelse(type == "nd", "nd",
                                          ifelse(type == "ref", "ref","sv")),
                      num_molecules = ifelse(support == 1,"1",">1"))

num_type_per_gene <- svs %>% group_by(gene,type,nick_type) %>% summarise(count = n()) %>%
  mutate(x_pos = 100000)
colnames(num_type_per_gene) <- c("gene","type","nick_type","samples_without_data","x_pos")

ggplot(data = svs, mapping = aes(y=gene, x = size2, color = type)) +
  geom_point(data = filter(num_type_per_gene, type == "nd",nick_type == "std"), mapping = aes(x = x_pos, y = gene, size = samples_without_data), color = "burlywood4") + 
  geom_quasirandom(data = filter(svs, type_group == "ref",nick_type == "std"),size = 2, color = "gray27", alpha = 0.65, varwidth = TRUE) +
  geom_quasirandom(data = filter(svs, type_group == "sv",nick_type == "std"),mapping = aes(shape = num_molecules), size = 2,alpha = 0.65,varwidth = TRUE) +
  scale_color_aaas() +
  labs(x = "size change relative to reference") + 
  geom_vline(xintercept = 95000, color = "gray27") +
  scale_x_continuous(breaks = seq(-75000,95000,10000)) +
  scale_y_discrete(limits = c("NBPF20",
                              "NBPF19",
                              "NBPF10",
                              "NBPF14",
                              "NBPF12",
                              "NBPF26",
                              "NBPF9",
                              "NBPF1",
                              "NBPF25P",
                              "NBPF15",
                              "NBPF8",
                              "NBPF3")) +
  theme(axis.text.y = element_text(size = 13,
                                   face = "bold")) + 
  theme(axis.text.x = element_text(size = 13,
                                   face = "bold")) +
  theme(axis.title.y = element_text(size = 13,
                                    face = "bold")) + 
  theme(axis.title.x = element_text(size = 13,
                                    face = "bold"))


# Make a plot of the HLS region varition for the genes that we had to use an alterate nick for 
ggplot(data = svs, mapping = aes(y=gene, x = size2, color = type)) +
  geom_point(data = filter(num_type_per_gene, type == "nd",nick_type == "alt"), mapping = aes(x = x_pos, y = gene, size = samples_without_data), color = "burlywood4") + 
  geom_quasirandom(data = filter(svs, type_group == "ref",nick_type == "alt"),size = 2, color = "gray27", alpha = 0.65, varwidth = TRUE) +
  geom_quasirandom(data = filter(svs, type_group == "sv",nick_type == "alt"),size = 2,alpha = 0.65,varwidth = TRUE) +
  scale_color_aaas() +
  labs(x = "size change relative to reference") + 
  geom_vline(xintercept = 95000, color = "gray27") +
  scale_x_continuous(breaks = seq(-75000,95000,10000)) +
  scale_y_discrete(limits = c("NBPF2P",
                              "NBPF11")) +
  theme(axis.text.y = element_text(size = 13,
                                   face = "bold")) + 
  theme(axis.text.x = element_text(size = 13,
                                   face = "bold")) +
  theme(axis.title.y = element_text(size = 13,
                                    face = "bold")) + 
  theme(axis.title.x = element_text(size = 13,
                                    face = "bold"))

###### Do same plots as above, but for CON1 region

svs.con1 <- read.delim("~/LabProjects/UCSF-visit/sv-calls-CON1-1000.txt", header = FALSE)
svs.con1 <- as_tibble(svs.con1)
colnames(svs.con1) <- c("gene","sample","type","size","support","min","max","range","nick_type","median.conf")

svs.con1 <- svs.con1 %>% mutate(size2 = ifelse(type == "del", size*-1,size),
                      type_group = ifelse(type == "nd", "nd",
                                          ifelse(type == "ref", "ref","sv")),
                      num_molecules = ifelse(support == 1,"1",">1"))

num_type_per_gene.con1 <- svs.con1 %>% group_by(gene,type,nick_type) %>% summarise(count = n()) %>%
  mutate(x_pos = 100000)
colnames(num_type_per_gene.con1) <- c("gene","type","nick_type","samples_without_data","x_pos")

ggplot(data = svs.con1, mapping = aes(y=gene, x = size2, color = type)) +
  geom_point(data = filter(num_type_per_gene.con1, type == "nd"), mapping = aes(x = x_pos, y = gene, size = samples_without_data), color = "burlywood4") + 
  geom_quasirandom(data = filter(svs.con1, type_group == "ref"), size = 2, color = "gray27", alpha = 0.65, varwidth = TRUE) +
  geom_quasirandom(data = filter(svs.con1, type_group == "sv"), size = 2,alpha = 0.65,varwidth = TRUE) +
  scale_color_aaas() +
  labs(x = "size change relative to reference") + 
  geom_vline(xintercept = 95000, color = "gray27") +
  scale_x_continuous(breaks = seq(-75000,95000,10000)) +
  scale_y_discrete(limits = c("NBPF1",
                              "NBPF12",
                              "NBPF8",
                              "NBPF10",
                              "NBPF11",
                              "NBPF14",
                              "NBPF26",
                              "NBPF9",
                              "NBPF15",
                              "NBPF17P",
                              "NBPF19",
                              "NBPF20",
                              "NBPF25P",
                              "NBPF3")) +
  theme(axis.text.y = element_text(size = 13,
                                   face = "bold")) + 
  theme(axis.text.x = element_text(size = 13,
                                   face = "bold")) +
  theme(axis.title.y = element_text(size = 13,
                                    face = "bold")) + 
  theme(axis.title.x = element_text(size = 13,
                                    face = "bold"))







# Plot data for every sample 
ggplot(data = svs, mapping = aes(y=gene, x = size, color = type)) + geom_point() + facet_wrap(~ sample) + theme(axis.text.y = element_text(size = 5)) + scale_color_aaas()


############ Code below here not likely to be necessary -- generates the histogram plots ########
# Plot all samples with facets (for a given gene)
ggplot(data = lengths, mapping = aes(x = length, fill = conf.group)) +
  geom_histogram(data = filter(lengths, gene == "NBPF10"),binwidth = 1000) + 
  geom_vline(xintercept = 63134) + 
  facet_wrap(~ sample)

geom_dotplot(stackgroups = TRUE, binwidth = 1, method = "histodot")

#geom_vline(xintercept = nbpf10.ref)

samples <- lengths %>% select(sample) %>% distinct(sample)
genes <- lengths %>% select(gene) %>% distinct(gene)

# Make a separate plot for every sample for a given gene
for (i in samples$sample) {
  each_sample <- i
  print(each_sample)
  print(ggplot(data = lengths, mapping = aes(x = length)) +
          geom_histogram(data = filter(lengths, gene == "NBPF10", sample == each_sample), binwidth = 1000) +
          labs(title = each_sample))}

# Plot one sample, one gene (or all genes if you take that out from the filter)
ggplot(data = lengths, mapping = aes(x = length, y = gene, color = conf.group)) +
  geom_point(data = filter(lengths, gene == "NBPF10", sample == "HG02283"))

geom_histogram(data = filter(lengths, gene == "NBPF19", sample == "HG01139"),binwidth = 1000) 


#dotplot line
#geom_dotplot(data = filter(lengths,sample =="HG02283", gene == "NBPF10"),
# stackgroups = TRUE,binwidth = 1000,method = "histodot")


geom_histogram(data = filter(lengths, sample == "HG02283", gene == "NBPF9"), binwidth = 1000) +
  facet_wrap(~ gene)

