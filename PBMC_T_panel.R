# Load packages

library(tidyverse)
library(ggpubr)
library(FSA)
library(factoextra)

# READ IN TABLE FROM FLOWJO ----

setwd("~/Desktop/T_STATS_211118")
PBMC <- read_csv("PBMC_stat.csv")
summarise(PBMC)
PBMC <- PBMC[-c(28:29),] # delete summary rows from original

# CD4 EXPLORATORY ----

# 1. Calculate ratios

CD4_PBMC <- mutate(PBMC,
                   ID = PBMC$...1,
                   CD4_ratio = PBMC$`CD4+ | Count` / PBMC$Count * 100,
                   Th2_in_memory_ratio = PBMC$`CD4+/T_memory/CD45RA-CCR6-/Th2 | Count` / PBMC$`CD4+/T_memory | Count` * 100,
                   Th17_in_memory_ratio = PBMC$`CD4+/T_memory/CD45RA-CCR6+/Th17 | Count` / PBMC$`CD4+/T_memory | Count` * 100,
                   Th1_in_memory_ratio = PBMC$`CD4+/T_memory/CD45RA-CCR6-/Th1 | Count` / PBMC$`CD4+/T_memory | Count` * 100,
                   Th1Th17_in_memory_ratio = PBMC$`CD4+/T_memory/CD45RA-CCR6+/Th17+1 | Count` / PBMC$`CD4+/T_memory | Count` * 100,
                   TCM_in_CD4_ratio = PBMC$`CD4+/T_CM | Count` / PBMC$`CD4+ | Count` * 100,
                   TEM_in_CD4_ratio = PBMC$`CD4+/T_EM | Count` / PBMC$`CD4+ | Count` * 100,
                   Teff_in_CD4_ratio = PBMC$`CD4+/T_effector | Count` / PBMC$`CD4+ | Count` * 100,
                   Tnaiv_in_CD4_ratio = PBMC$`CD4+/T_naive | Count` / PBMC$`CD4+ | Count` * 100,
                   Tmemory_in_CD4_ratio = PBMC$`CD4+/T_memory | Count` / PBMC$`CD4+ | Count`*100)
       
ratios <- select(CD4_PBMC, ID, group, CD4_ratio:Tmemory_in_CD4_ratio)

# 2. Group by participant status

ratios %>%
  group_by(group) %>%
  summarise_at(vars(CD4_ratio:Tmemory_in_CD4_ratio), list(name = mean))

ratios$group <- as.factor(ratios$group)
head(ratios)

# 3. Exploratory plot

proba_melt <- pivot_longer(ratios, CD4_ratio:Tmemory_in_CD4_ratio, names_to = "variables", values_to = "values")
proba_melt

q <- ggplot(proba_melt, aes(x = variables, y = values, fill = group)) +  # ggplot function
  geom_boxplot()+
  labs(y="Percentages (%)", x = "Populations", colour = "Groups") + 
  scale_fill_manual(breaks = c("Gabon_neg", "Gabon_pos", "German"), 
                    values=c("royalblue", "red2", "seagreen4")) +
  scale_x_discrete(labels=c("CD4_ratio" = "CD4", "Th2_in_memory_ratio" = "Th2/memory", "Th17_in_memory_ratio" = "Th17/memory", 
                            "Th1Th17_in_memory_ratio" = "Th1Th17/memory",
                            "Th1_in_memory_ratio" = "Th1/memory",
                            "TEM_in_CD4_ratio" = "TEM",
                            "TCM_in_CD4_ratio" = "TCM",
                            "Teff_in_CD4_ratio" = "T effector",
                            "Tnaiv_in_CD4_ratio" = "T naiv",
                            "Tmemory_in_CD4_ratio" = "T memory"
  )) +
  coord_flip() +
  stat_compare_means(method = "kruskal.test", label = "p.signif", hide.ns = TRUE) +
  theme_bw()

# CD4 STATS ----

cols <- names(ratios)[3:ncol(ratios)]
cols

all_test <- lapply(cols, function(x) 
  kruskal.test(reformulate("group", x), data = ratios))

all_test

# post-hoc test
dunnTest(ratios$TEM_in_CD4_ratio~ratios$group,data=ratios,method="bonferroni")



# PRETTY PLOT for statistically significant differences ----


p <- ggplot(ratios, aes(x=group, y=TEM_in_CD4_ratio, fill = ratios$group)) + 
  geom_boxplot() +
  labs(title = "Effector memory CD4 cells", y="Percentages (%)", x = "Group", fill = "group") + 
  scale_fill_manual(breaks = c("Gabon_neg", "Gabon_pos", "German"), 
                    values=c("royalblue", "red2", "seagreen4")) +
  scale_x_discrete(labels=c("Gabon_neg" = "Gabonese non-infected", 
                            "Gabon_pos" = "Gabonese infected",
                            "German" = "German"
  )) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  
my_comparisons <- list( c("Gabon_neg", "Gabon_pos"), c("Gabon_pos", "German"), c("Gabon_neg", "German") )
p + stat_compare_means(comparisons = my_comparisons, label = "p.signif")+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 25)
 




# TESTS ONE BY ONE -----

group_by(ratios, group) %>%
  summarise(
    count = n(),
    mean = mean(CD4_ratio, na.rm = TRUE),
    sd = sd(CD4_ratio, na.rm = TRUE),
    median = median(CD4_ratio, na.rm = TRUE),
    IQR = IQR(CD4_ratio, na.rm = TRUE)
  )

kruskal.test(CD4_ratio ~ group, data = ratios)

# TEM

group_by(ratios, group) %>%
  summarise(
    count = n(),
    mean = mean(TEM_in_CD4_ratio, na.rm = TRUE),
    sd = sd(TEM_in_CD4_ratio, na.rm = TRUE),
    median = median(TEM_in_CD4_ratio, na.rm = TRUE),
    IQR = IQR(TEM_in_CD4_ratio, na.rm = TRUE)
  )

kruskal.test(TEM_in_CD4_ratio ~ group, data = ratios)

pairwise.wilcox.test(ratios$TEM_in_CD4_ratio, ratios$group,
                     p.adjust.method = "BH")

my_comparisons <- list( c("Gabon_neg", "Gabon_pos"), c("Gabon_pos", "German"), c("Gabon_neg", "German") )
ggboxplot(ratios, x = "group", y = "TEM_in_CD4_ratio",
          color = "group", palette = "jco") +
  stat_compare_means(comparisons = my_comparisons, label =  "p.signif") +# Add pairwise comparisons p-value
  stat_compare_means(label.y = 30)  # Add global p-value
       

# NCSM B 

group_by(ratios, Group) %>%
  summarise(
    count = n(),
    mean = mean(NCSM_ratio, na.rm = TRUE),
    sd = sd(NCSM_ratio, na.rm = TRUE),
    median = median(NCSM_ratio, na.rm = TRUE),
    IQR = IQR(NCSM_ratio, na.rm = TRUE)
  )

kruskal.test(NCSM_ratio ~ Group, data = ratios)

# NCSM B / memory cells

group_by(ratios, Group) %>%
  summarise(
    count = n(),
    mean = mean(NCSM_in_memory, na.rm = TRUE),
    sd = sd(NCSM_in_memory, na.rm = TRUE),
    median = median(NCSM_in_memory, na.rm = TRUE),
    IQR = IQR(NCSM_in_memory, na.rm = TRUE)
  )

kruskal.test(NCSM_in_memory ~ Group, data = ratios)

pairwise.wilcox.test(ratios$NCSM_in_memory, ratios$Group,
                     p.adjust.method = "BH")

# switch / memory cells

group_by(ratios, Group) %>%
  summarise(
    count = n(),
    mean = mean(switch_in_memory, na.rm = TRUE),
    sd = sd(switch_in_memory, na.rm = TRUE),
    median = median(switch_in_memory, na.rm = TRUE),
    IQR = IQR(switch_in_memory, na.rm = TRUE)
  )

kruskal.test(switch_in_memory ~ Group, data = ratios)

pairwise.wilcox.test(ratios$switch_in_memory, ratios$Group,
                     p.adjust.method = "BH")

# switch B 

group_by(ratios, Group) %>%
  summarise(
    count = n(),
    mean = mean(switch_ratio, na.rm = TRUE),
    sd = sd(switch_ratio, na.rm = TRUE),
    median = median(switch_ratio, na.rm = TRUE),
    IQR = IQR(switch_ratio, na.rm = TRUE)
  )

kruskal.test(switch_ratio ~ Group, data = ratios)

pairwise.wilcox.test(ratios$switch_ratio, ratios$Group,
                     p.adjust.method = "BH")

# transB 

group_by(ratios, Group) %>%
  summarise(
    count = n(),
    mean = mean(transB_ratio, na.rm = TRUE),
    sd = sd(transB_ratio, na.rm = TRUE),
    median = median(transB_ratio, na.rm = TRUE),
    IQR = IQR(transB_ratio, na.rm = TRUE)
  )

kruskal.test(transB_ratio ~ Group, data = ratios)



# CD8 EXPLORATORY ----

# 1. Calculate ratios
CD8_PBMC <- mutate(PBMC,
                   ID = PBMC$...1,
                   CD8_ratio = PBMC$`CD8+ | Count` / PBMC$Count * 100,
                   TCM_in_CD8_ratio = PBMC$`CD8+/CD8T_CM | Count` / PBMC$`CD8+ | Count` * 100,
                   TEM_in_CD8_ratio = PBMC$`CD8+/CD8T_EM | Count`/ PBMC$`CD8+ | Count` * 100,
                   Teff_in_CD8_ratio = PBMC$`CD8+/CD8T_eff | Count`/ PBMC$`CD8+ | Count` * 100,
                   Tnaiv_in_CD8_ratio = PBMC$`CD8+/CD8T_naive | Count` / PBMC$`CD8+ | Count` * 100,
                   Tmemory_in_CD8_ratio = PBMC$`CD8+/CD8T_memory | Count` / PBMC$`CD8+ | Count`*100
)

CD8_ratios <- select(CD8_PBMC, ID, group, CD8_ratio:Tmemory_in_CD8_ratio)

# 2. Group by participant status

CD8_ratios %>%
  group_by(group) %>%
  summarise_at(vars(CD8_ratio:Tmemory_in_CD8_ratio), list(name = mean))

CD8_ratios$group <- as.factor(CD8_ratios$group)
head(CD8_ratios)





CD8_melt <- pivot_longer(CD8_ratios, CD8_ratio:Tmemory_in_CD8_ratio, names_to = "variables", values_to = "values")
CD8_melt

ggplot(CD8_melt, aes(x = variables, y = values, color = group)) +  # ggplot function
  geom_boxplot()


q <- ggplot(CD8_melt, aes(x = variables, y = values, fill = group)) +  # ggplot function
  geom_boxplot()+
  labs(y="Percentages (%)", x = "Populations", colour = "Groups") + 
  scale_fill_manual(breaks = c("Gabon_neg", "Gabon_pos", "German"), 
                    values=c("royalblue", "red2", "seagreen4")) +
  scale_x_discrete(labels=c("CD8_ratio" = "CD8", 
                            "TEM_in_CD8_ratio" = "TEM",
                            "TCM_in_CD8_ratio" = "TCM",
                            "Teff_in_CD8_ratio" = "T effector",
                            "Tnaiv_in_CD8_ratio" = "T naiv",
                            "Tmemory_in_CD8_ratio" = "T memory"
  )) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  stat_compare_means(method = "kruskal.test", label = "p.signif", hide.ns = TRUE, label.y = 75) #+
  

q


# CD8 STATS ----

CD8_cols <- names(CD8_ratios)[3:ncol(CD8_ratios)]
CD8_cols

CD8_all_test <- lapply(CD8_cols, function(x) 
  kruskal.test(reformulate("group", x), data = CD8_ratios))

CD8_all_test

# Post-hoc tests

dunnTest(CD8_ratios$TCM_in_CD8_ratio~ratios$group,data=CD8_ratios,method="bonferroni")
dunnTest(CD8_ratios$Teff_in_CD8_ratio~ratios$group,data=CD8_ratios,method="bonferroni")
dunnTest(CD8_ratios$Tnaiv_in_CD8_ratio~ratios$group,data=CD8_ratios,method="bonferroni")


# PRETTY PLOT - CD8 ----

TCM_plot <- ggplot(CD8_ratios, aes(x=group, y=TCM_in_CD8_ratio, fill = CD8_ratios$group)) + 
  geom_boxplot() +
  labs(title = "Central memory CD8 cells", y="Percentages (%)", x = "Group", fill = "group") + 
  scale_fill_manual(breaks = c("Gabon_neg", "Gabon_pos", "German"), 
                    values=c("royalblue", "red2", "seagreen4")) +
  scale_x_discrete(labels=c("Gabon_neg" = "Gabonese non-infected", 
                            "Gabon_pos" = "Gabonese infected",
                            "German" = "German"
  )) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

my_comparisons <- list( c("Gabon_neg", "Gabon_pos"), c("Gabon_pos", "German"), c("Gabon_neg", "German") )
TCM_plot + stat_compare_means(comparisons = my_comparisons, label = "p.signif")+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 35)

Teff_plot <- ggplot(CD8_ratios, aes(x=group, y=Teff_in_CD8_ratio, fill = CD8_ratios$group)) + 
  geom_boxplot() +
  labs(title = "Effector CD8 cells", y="Percentages (%)", x = "Group", fill = "group") + 
  scale_fill_manual(breaks = c("Gabon_neg", "Gabon_pos", "German"), 
                    values=c("royalblue", "red2", "seagreen4")) +
  scale_x_discrete(labels=c("Gabon_neg" = "Gabonese non-infected", 
                            "Gabon_pos" = "Gabonese infected",
                            "German" = "German"
  )) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

my_comparisons <- list( c("Gabon_neg", "Gabon_pos"), c("Gabon_pos", "German"), c("Gabon_neg", "German") )
Teff_plot + stat_compare_means(comparisons = my_comparisons, label = "p.signif")+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 60)

Tnaiv_plot <- ggplot(CD8_ratios, aes(x=group, y=Tnaiv_in_CD8_ratio, fill = CD8_ratios$group)) + 
  geom_boxplot() +
  labs(title = "Naive CD8 cells", y="Percentages (%)", x = "Group", fill = "group") + 
  scale_fill_manual(breaks = c("Gabon_neg", "Gabon_pos", "German"), 
                    values=c("royalblue", "red2", "seagreen4")) +
  scale_x_discrete(labels=c("Gabon_neg" = "Gabonese non-infected", 
                            "Gabon_pos" = "Gabonese infected",
                            "German" = "German"
  )) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

my_comparisons <- list( c("Gabon_neg", "Gabon_pos"), c("Gabon_pos", "German"), c("Gabon_neg", "German") )
Tnaiv_plot + stat_compare_means(comparisons = my_comparisons, label = "p.signif")+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 110)




# PCA ----

# Combine T cell tables
total_T <- full_join(ratios, CD8_ratios, by = "ID")
combined_T <- select(total_T, -c(ID, group.x, group.y))

# Calculate PCA
res.pca <- prcomp(combined_T, scale = TRUE)
fviz_eig(res.pca)

# Plots 

fviz_pca_ind(res.pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

groups <- as.factor(total_T$group.x)
fviz_pca_ind(res.pca,
             col.ind = groups, # color by groups
             palette = c("royalblue", "red2", "seagreen4"),
             addEllipses = TRUE, # Concentration ellipses
             ellipse.type = "confidence",
             legend.title = "Groups",
             title = "PCA of T cell subpopulations",
             repel = TRUE
)

fviz_pca_var(res.pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             select.var = list(contrib = 10),
             repel = TRUE     # Avoid text overlapping
)

fviz_pca_biplot(res.pca, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969"  # Individuals color
)

