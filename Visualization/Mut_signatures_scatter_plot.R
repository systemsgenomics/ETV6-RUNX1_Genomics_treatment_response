#Scatter plots to show relationship between MRD on APOBEC

# Load necessary library
library(ggplot2)
library(hrbrthemes)
library(ggpmisc)
library(ggpubr)
library(dplyr)

# Mid-induction (day15) #

# Read data (ID, signature contributions, MRDs, response groups)
sig_da15 <- read.table("scatter_data_day15.txt", header=T)
sig_da15

# Create a scatter plot using ggplot2

#SBS1 (Deamination of 5-methylcytosine)
sbs1 <- ggplot(sig_da15_mod, aes(x = SBS1, y = day15_flow)) +
  geom_point(aes(color = day15_flow < 0.1), size = 3) + # Plot points, set point size, set cut off for point colors
  labs(title = "SBS1 vs mrd_day15",
       x = "SBS1",
       y = "day15_flow") +
  geom_smooth(method=lm,se = FALSE, colour="deepskyblue4") +
  scale_color_manual(values = c("TRUE" = "#A1D76A", "FALSE" = "#B2182B")) + # Map colors for points
  theme_minimal() 

#calculate Spearman's rank correlation coefficient
sbs1+stat_cor(method = "spearman")


#SBS2 (APOBEC activity)
sbs2 <- ggplot(sig_da15_mod, aes(x = SBS2, y = day15_flow)) +
  geom_point(aes(color = day15_flow < 0.1), size = 3) + 
  labs(title = "SBS2 (APOBEC activity) vs Mid induction MRD",
       x = "SBS2 ",
       y = "Mid indcution MRD") +
  geom_smooth(method=lm,se = FALSE, colour="deepskyblue4") +
  scale_color_manual(values = c("TRUE" = "#A1D76A", "FALSE" = "#B2182B")) + 
  theme_minimal() 

#calculate Spearman's rank correlation coefficient
sbs2+stat_cor(method = "spearman")


#SBS8 (Unknown)
sbs8 <- ggplot(sig_da15_mod, aes(x = SBS8, y = day15_flow)) +
  geom_point(aes(color = day15_flow < 0.1), size = 3) + 
  labs(title = "SBS8 vs mrd_day15",
       x = "SBS8",
       y = "day15_flow") +
  geom_smooth(method=lm,se = FALSE, colour="deepskyblue4") +
  scale_color_manual(values = c("TRUE" = "#A1D76A", "FALSE" = "#B2182B")) + 
  theme_minimal() 

#calculate Spearman's rank correlation coefficient
sbs8+stat_cor(method = "spearman")


#SBS13 (APOBEC activity)
sbs13 <- ggplot(sig_da15_mod, aes(x = SBS13, y = day15_flow)) +
  geom_point(aes(color = day15_flow < 0.1), size = 3) + 
  labs(title = "SBS13 vs mrd_day15",
       x = "SBS13",
       y = "day15_flow") +
  geom_smooth(method=lm,se = FALSE, colour="deepskyblue4") +
  scale_color_manual(values = c("TRUE" = "#A1D76A", "FALSE" = "#B2182B")) + 
  theme_minimal() 

#calculate Spearman's rank correlation coefficient
sbs13+stat_cor(method = "spearman")


# EOI (day29)

# Read data (ID, signature contributions, MRDs, response groups)
sig_day29 <- read.table("scatter_data_day29.txt", header=T)
sig_day29

#SBS2 (APOBEC activity)
sbs2 <- ggplot(sig_day29_mod, aes(x = SBS2, y = day29_flow)) +
  geom_point(aes(color = case_when(
    day29_flow == 0 ~ "fast",
    day29_flow < 0.001 & day29_flow > 0 ~ "intermediate",
    day29_flow > 0.001 ~ "slow" )), size = 3) +
    labs(title = "SBS2 vs mrd_day29",
         x = "SBS2",
         y = "day29_flow") +
    geom_smooth(method = lm, se = FALSE, colour = "deepskyblue4") +
    scale_color_manual(values = c("fast" = "#A1D76A", "intermediate" = "#E69F00", "slow" = "#B2182B")) +
    theme_minimal()

#calculate Spearman's rank correlation coefficient 
  sbs2 + stat_cor(method = "spearman")

