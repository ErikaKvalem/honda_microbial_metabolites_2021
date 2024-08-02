



df_num = read_csv("df_num.csv")
df_num  = as.data.frame(df_num)
df_num$log_cf = log(df_num$cell_fraction)


## ALL 
p0_results = compare_means(cell_fraction ~ condition, data = df_num, method = "wilcox.test")

my_comparisons <- list( c("10mix", "11mix"), c("11mix", "GF"), c("10mix", "GF") )

# Is there any difference between the conditions for cell_fraction 
p0 = ggboxplot(df_num, x = "condition", y = "cell_fraction",color="condition",
                palette = "jco", add = "jitter")  + 
  stat_compare_means(comparisons = my_comparisons, label="p.signif", method="wilcox.test", ref.group = "GF")

p0



# Is there any difference between groups in cell_fraction/ cell_count 
p1_results = compare_means(cell_fraction ~ cell_type,  data = df_num,
              ref.group = ".all.", method = "t.test")
# Calculate the maximum value for cell_fraction
max_cell_fraction <- max(df_num$cell_fraction, na.rm = TRUE)

p1= ggboxplot(df_num, x = "cell_type", y = "cell_fraction", color = "cell_type", 
          add = "jitter", legend = "none") +
  rotate_x_text(angle = 45)+
  geom_hline(yintercept = mean(df_num$cell_fraction), linetype = 2)+ # Add horizontal line at base mean
 # stat_compare_means(method = "anova", label.y = 1600)+        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE)   +
  ylim(0, max_cell_fraction * 1.1)  # Set y-axis limits

p1
#Can conclude that cell fraction is significantly increased in SLCE_Effector group and
#it’s significantly decreased in MPEC_Effector, SLEC_Terminal, MPEC_intermediate & SLEC_Plastic

## Grouping variables  

desired_order <- c("MPEC_Progenitor", "MPEC_Intermediate", "MPEC_Effector",
                   "SLEC_Plastic", "SLEC_Progenitor", "SLEC_Intermediate",
                   "SLEC_Ifn", "SLEC_Effector", "SLEC_Terminal")

# Set the cell_type column to a factor with the specified order
df_num$cell_type <- factor(df_num$cell_type, levels = desired_order)

# Remove rows where condition is "GF"
#df_num<- df_num[df_num$condition != "GF", ]


p2_results = compare_means(cell_fraction ~ condition, data = df_num, 
              group.by = "cell_type")

# Customize facet labels
facet_labels <- c("10mix", "11mix")

# Box plot facetted by "dose"
p2 <- ggboxplot(df_num, x = "cell_type", y = "cell_fraction",
               color = "cell_type", palette = "jco",
               add = "jitter",
               facet.by = "condition", short.panel.labs = FALSE)+  rotate_x_text(angle = 45) + 
  facet_wrap(~ condition, labeller = as_labeller(setNames(facet_labels, unique(df_num$condition)))) +
  theme(legend.position = "none") + # Remove legend
  ylab("Cell fraction")  + # Set y-axis label
xlab("Cell type")

p2 
