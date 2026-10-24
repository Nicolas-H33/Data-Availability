library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(multcompView)
library(ggplot2)

setwd("~/Library/CloudStorage/OneDrive-IowaStateUniversity/Nick_MS_Thesis/lab presentations, posters, writing, pictures/My writing/Objective 1 Manuscript")

selected <- "fusarium"

df <- read_excel("Supplementary_Tables_Manuscript.xlsx", sheet = "Table S8") %>%
  filter(str_squish(as.character(Study)) == selected) %>%
  mutate(Treatment = factor(
    case_when(
      str_detect(str_to_lower(str_squish(`mock/inoculated`)), "^mock") ~ "Mock",
      str_detect(str_to_lower(str_squish(`mock/inoculated`)), "^inoc") ~ "Inoculated"),
    levels = c("Mock", "Inoculated"))) %>%
  rename(`GmROS2-1` = `Glyma.10G065900.Wm82.a6.v1`, `GmROS2-2` = `Glyma.13G151000.Wm82.a6.v1`,
         `GmROS1`   = `Glyma.03G190800.Wm82.a6.v1`, `GmDME1-1` = `Glyma.20G188300.Wm82.a6.v1`,
         `GmDME1-2` = `Glyma.10G202150.Wm82.a6.v1`)

tags <- c("GmROS2-1", "GmROS2-2", "GmROS1", "GmDME1-1", "GmDME1-2")

#Statistics
anova_results <- list()
tukey_results <- list()
letters_results <- list()

for (g in tags) {
  
  fit <- aov(as.formula(paste0("`", g, "` ~ Treatment")), data = df)
  
  anova_results[[g]] <- summary(fit)
  tukey_results[[g]] <- TukeyHSD(fit)
  
  tuk_group <- tukey_results[[g]][["Treatment"]]
  
  p_vals <- tuk_group[, "p adj"]
  names(p_vals) <- rownames(tuk_group)
  
  letters_results[[g]] <- multcompLetters(p_vals)$Letters
  print(letters_results[[g]])
}

#Mean, standard error, and significance letters
plot_df <- df %>%
  pivot_longer(cols = all_of(tags), names_to = "Gene", values_to = "TPM") %>%
  group_by(Gene, Treatment) %>%
  summarise(
    mean_TPM = mean(TPM, na.rm = TRUE),
    se_TPM = sd(TPM, na.rm = TRUE) / sqrt(sum(!is.na(TPM))),
    .groups = "drop") %>%
  mutate(
    Letters = mapply(function(g, trt) {
      letters_results[[g]][[as.character(trt)]]
    }, Gene, Treatment),
    Letters = ifelse(
      ave(Letters, Gene, FUN = function(x) length(unique(x))) == 1,
      NA_character_, Letters),
    Gene = factor(Gene, levels = tags))

#ggplot
y <- 5
pos_dodge <- position_dodge(width = 0.8)

ggplot(plot_df, aes(x = Gene, y = mean_TPM, fill = Treatment)) +
  geom_col(position = pos_dodge, width = 0.78,
           color = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = mean_TPM - se_TPM,
                    ymax = mean_TPM + se_TPM),
                position = pos_dodge, width = 0.2) +
  geom_text(aes(label = Letters, y = mean_TPM + se_TPM + 1),
            position = pos_dodge, size = 4, na.rm = TRUE) +
  scale_fill_manual(
    values = c(
      "Mock"       = "#9bbfe0",
      "Inoculated" = "#f4a96c"),
    name = "Treatment") +
  scale_y_continuous(breaks = function(x) seq(0, max(x), by = y)) +
  labs(
    title = "Study: Fusarium virguliforme",
    y = "Mean TPM ± SE",
    x = "Gene") +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15))

#Results
list(ANOVA = anova_results, Tukey = tukey_results, Letters = letters_results,
     Mean_SE = plot_df %>%
       select(Gene, Treatment, mean_TPM, se_TPM))
