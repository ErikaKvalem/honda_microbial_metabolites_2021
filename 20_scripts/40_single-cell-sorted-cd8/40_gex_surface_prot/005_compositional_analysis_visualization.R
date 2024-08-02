############ Compositional analysis for TIL Honda project 

plot_df_inf = read_csv("plot_df_INF.csv")

frac_by_condition_inf <- read_csv("frac_by_condition_INF.csv")

plot_df= read_csv("plot_df.csv")
colnames(plot_df)[colnames(plot_df) == 'Cell type'] <- 'cell_type'

frac_by_condition <- read_csv("frac_by_condition.csv")


###################### SCCODA Cell proportion for 10mix, 11mix, GF 

desired_order <- c("MPEC_Progenitor", "MPEC_Intermediate", "MPEC_Effector",
                   "SLEC_Plastic", "SLEC_Progenitor", "SLEC_Intermediate",
                   "SLEC_Ifn", "SLEC_Effector", "SLEC_Terminal")

# Set the cell_type column to a factor with the specified order
plot_df$"Cell type" <- factor(plot_df$"cell_type", levels = desired_order)
# Customize facet labels
facet_labels <- c("10mix", "11mix","GF")


p1 <- ggboxplot(plot_df, x = "cell_type", y = "Proportion",
                color = "cell_type", palette = "jco",
                add = "jitter",
                facet.by = "condition", short.panel.labs = FALSE) +
  rotate_x_text(angle = 45) +
  facet_wrap(~ condition, labeller = as_labeller(setNames(facet_labels, unique(plot_df$condition)))) +
  theme(legend.position = "none", # Remove legend
        axis.text.x = element_text(size = 8)) + # Set x-axis label size
  ylab("Cell proportion") + # Set y-axis label
  xlab("") # Set x-axis label
p1

ggsave("cell_proportion.png", p1, width =6, height = 4)

######################### Number of cells 

desired_order <- c("MPEC_Progenitor", "MPEC_Intermediate", "MPEC_Effector",
                   "SLEC_Plastic", "SLEC_Progenitor", "SLEC_Intermediate",
                   "SLEC_Ifn", "SLEC_Effector", "SLEC_Terminal")

# Set the cell_type column to a factor with the specified order
frac_by_condition$cell_type <- factor(frac_by_condition$cell_type, levels = desired_order)


p2 <- ggbarplot(frac_by_condition, "cell_type", "n_cells",
                fill = "cell_type", color = "cell_type",   palette = "jco") + rotate_x_text(angle = 45) +  xlab("") +   ylab("Number of cells") + # Set y-axis label
  theme(legend.position = "none",  axis.title.x = element_text(size = 10))

p2 

ggsave("number_of_cells.png", p2, width =6, height = 4)

######################### # Is there any difference between groups in cell_count 
p3_results = compare_means(n_cells ~ cell_type,  data = frac_by_condition,
                           ref.group = ".all.", method = "t.test")
# Calculate the maximum value for cell_fraction
max_cell_num <- max(frac_by_condition$n_cells, na.rm = TRUE)

p3= ggboxplot(frac_by_condition, x = "cell_type", y = "n_cells", color = "cell_type", 
              add = "jitter", legend = "none",  palette = "jco") +
  rotate_x_text(angle = 45)+
  geom_hline(yintercept = mean(frac_by_condition$n_cells), linetype = 2)+ # Add horizontal line at base mean
  # stat_compare_means(method = "anova", label.y = 1600)+        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE)   +
  ylim(0, max_cell_num * 1.1)+   ylab("Number of cells") + xlab("") +  theme(legend.position = "none", # Remove legend
                                                                             axis.title.x = element_text(size = 8)) # Set x-axis label size # Set y-axis limits

p3

ggsave("number_of_cells_baseline.png", p3, width =6, height = 4)

##########################  Cell proportion Treated vs Untreated GF all cell types 

filtered_df <- subset(plot_df, condition == "GF")

desired_order <- c("MPEC_Progenitor", "MPEC_Intermediate", "MPEC_Effector",
                   "SLEC_Plastic", "SLEC_Progenitor", "SLEC_Intermediate",
                   "SLEC_Ifn", "SLEC_Effector", "SLEC_Terminal")

# Set the cell_type column to a factor with the specified order
filtered_df$cell_type <- factor(filtered_df$cell_type, levels = desired_order)

p4_results = compare_means(Proportion ~ treatment, data = filtered_df, method = "wilcox.test")

my_comparisons <- list( c("treated", "naive") )

# Is there any difference between the conditions for cell_fraction 
p4 = ggboxplot(filtered_df, x = "cell_type", y = "Proportion",color="treatment",
               palette = "jco")  + 
  stat_compare_means(comparisons = my_comparisons, label="p.signif", method="wilcox.test", ref.group = "naive")+
  rotate_x_text(angle = 45) +   ylab("Cell proportion") +   xlab("")  +  theme(
        axis.title.x = element_text(size = 7)) # Set x-axis label size # Set y-axis label


p4
ggsave("cell_proportion_treated_vs_untreated_GF.png", p4, width =6, height = 4)

############# Cell type: SLEC_Ifn TREATED only 
conflicts_prefer(rstatix::filter)

filtered_df <- plot_df
#filtered_df <- plot_df %>% filter(cell_type == "SLEC_Ifn")
#filtered_df <- filtered_df %>% filter(treatment == "naive")
p5_results = compare_means(Proportion ~ condition, data = filtered_df, method = "wilcox.test")

my_comparisons <- list( c("10mix", "11mix"), c("11mix", "GF"), c("10mix", "GF") )
#my_comparisons <- list( c("treated", "naive") )
# Is there any difference between the conditions for cell_fraction 
p5 = ggboxplot(filtered_df, x = "condition", y = "Proportion",color="condition",
               palette = "jco", add = "jitter")  + 
  stat_compare_means(comparisons = my_comparisons, label="p.signif", method="wilcox.test.test", ref.group = "GF")+
  theme(legend.position = "none")+  # Remove legend
  ylab("Cell proportion") + # Set y-axis label
  xlab("Treatment") # Set x-axis label

p5

ggsave("SLEC_Ifn_treated.png", p5, width =4, height = 4)















