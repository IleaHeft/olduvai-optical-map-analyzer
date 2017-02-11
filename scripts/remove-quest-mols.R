library(dplyr)


quest <- read.delim("test5.txt", header = FALSE)
colnames(quest) <- c("mol","conf")
quest <- unique(quest$mol)

all <- read.delim("GM11994-all-mols.txt", header = FALSE)


stripped <- all %>% filter(!V3 %in% quest) 

quest
all %>% arrange(V3)

stripped %>% arrange(V3)


all %>% group_by(V3) %>% summarise(count = n()) %>% filter(count != 1) 


length(all$V2)
length(quest)
length(stripped$V2)


length(unique(all$V3))
length(quest)
length(unique(stripped$V3))
