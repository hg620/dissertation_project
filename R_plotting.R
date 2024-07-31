#1. ANALYSING CORE GENES IN THE ALIGNMENT OF 500 SEQUENCES 
gene_presence_counts$percentage <- (gene_presence_counts$X2 / 499) * 100
ggplot(gene_presence_counts, aes(x = percentage, y = reorder(X1, percentage))) +
  geom_bar(stat = "identity", fill = "skyblue") +
  theme_minimal() +
  labs(
    title = "Percentage of Genomes with Gene Present",
    x = "Percentage of Genomes",
    y = "Gene Name"
  )
#2. ANALYSING GAP COUNTS IN THE ALIGNMENT 
ggplot(gap_counts, aes(x = Gap_Count, y = Frequency)) +
  geom_bar(stat = "identity", fill = "blue", alpha = 0.7) +
  theme_minimal() +
  labs(title = "Gap Count Frequencies",
       x = "Number of Gaps",
       y = "Frequency") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#3. ANALYSIS OF PHASTCONS SCORES
#Plotting scores in wig format 
library(rphast)
library(ggplot2)
library(dplyr)
wig_file <- read.wig("C:****/scores.wig")
genes_boundaries <- read_table("C:****/genes_boundaries.txt", 
                               +     col_names = FALSE)
wig_df <- data.frame(
  start = wig_file$start,
  end = wig_file$end,
  score = wig_file$score
)
wig_df$position <- seq(from = min(wig_df$start), to = max(wig_df$end), length.out = nrow(wig_df))
ggplot(wig_df, aes(x = position, y = score)) +
  geom_line(color = "blue") +
  labs(title = "WIG Scores", x = "Position", y = "Score") +
  theme_minimal()

# Aggregating data 
breaks <- seq(0, 1, by = 0.05)
labels <- paste0(seq(0, 0.95, by = 0.05), "-", seq(0.05, 1, by = 0.05))
wig_file$category <- cut(wig_file$score, breaks = breaks, labels = labels, include.lowest = TRUE)
score_agg <- wig_file %>%
  group_by(category) %>%
  summarise(count = n())

# Plotting agreggated data 
ggplot(score_agg, aes(x = category, y = count)) +
  geom_bar(stat = "identity", fill = "blue") +
  theme_minimal() +
  labs(title = "Distribution of PhastCons scores across 50 E. coli genomes",
       x = "PhastCons score",
       y = "Count") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


#plotting scores across whole genome
ggplot(wig_file, aes(x = start, y = score)) +
  geom_line(color = "blue", alpha = 0.5) +
  labs(title = "PhastCons Scores across Genome Length", x = "Genome Position", y = "PhastCons Score") +
  theme_minimal()


# Filtering the dataframe to include only the first 45,000 positions
filtered_wig_file <- wig_file[wig_file$start <= 45000, ]

# Plotting wig scores across the first 45,000 genome positions
ggplot(filtered_wig_file, aes(x = start, y = score)) +
  geom_line(color = "blue") +
  labs(title = "PhastCons Scores across the First 45,000 Genome Positions", x = "Genome Position", y = "PhastCons Score") +
  theme_minimal()

# Plotting annotated genes in K12 reference at first 45 000 positions 
genes_filtered <- genes_boundaries[genes_boundaries$X2 <= 45000 | genes_boundaries$X3 <= 45000, ]
positions <- unique(c(genes_filtered$X2, genes_filtered$X3))
ggplot(data.frame(position = positions), aes(x = position)) +
  geom_bar(stat = "identity", width = 1) +
  labs(title = "Annotated Genes across Genome Positions",
       x = "Genome Position",
       y = "Annotated Genes") +
  theme_minimal() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

#plotting combined phastcons scores and annotated genes 
ggplot() +
  geom_line(data = filtered_wig_file, aes(x = start, y = score), color = "blue") +
  geom_segment(data = genes_filtered, aes(x = X2, xend = X3, y = 0.8, yend = 0.8), color = "red", size = 1) +
  labs(title = "Comparison of PhastCons Scores and Annotated Genes across the First 45,000 Genome Positions",
       x = "Genome Position", y = "PhastCons Score / Annotated Genes") +
  theme_minimal() +
  theme(axis.text.y = element_text(), axis.ticks.y = element_line())
