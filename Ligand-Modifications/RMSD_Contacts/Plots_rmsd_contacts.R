### Visualizations of Molecular Simulation Results of Ligand Modifications ####
###############################################################################

# Load libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)

# Load and label the data
mod1 <- read_table("MOD1_final_contacts.tsv", col_names = FALSE) %>%
  mutate(ligand = "MOD1")
mod2 <- read_table("MOD2_final_contacts.tsv", col_names = FALSE) %>%
  mutate(ligand = "MOD2")
adn  <- read_table("ADN_final_contacts.tsv", col_names = FALSE) %>%
  mutate(ligand = "ADN")

# Assign column names
colnames(mod1) <- colnames(mod2) <- colnames(adn) <- c("residue", "residue_number", "frame", "contacts", "ligand")

# Combine data
all_data <- bind_rows(mod1, mod2, adn) %>%
  mutate(residue_label = paste(residue, residue_number, sep = " "))

# Plot without forcing residue factor levels
ggplot(all_data, aes(x = frame, y = residue_label, fill = contacts)) +
  geom_tile() +
  facet_wrap(~ ligand, nrow = 1, scales = "free_y") +  # allow each facet to have its own residues
  scale_fill_viridis_c(option = "plasma", limits = c(0, 60)) +
  labs(
    title = "Contacts per Residue and Frame (within 3Å) ",
    x = "Simulation Frame",
    y = "Residue",
    fill = "Contacts"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 16),
    plot.title = element_text(face = "bold", size = 18)
  )

### ADN Ligand ################################################################


# Reload data with correct column names
rmsd_data <- read.table("ADN_ligand_rmsd.dat", header = FALSE)
colnames(rmsd_data) <- c("Frame", "RMSD")  # Rename columns
mean_rmsd <- mean(rmsd_data$RMSD)
print(mean_rmsd)


# Verify column names
print(colnames(rmsd_data))  # Should show "Frame" and "RMSD"

# Check min/max values
min(rmsd_data$Frame)
max(rmsd_data$Frame)

plot(rmsd_data$Frame, rmsd_data$RMSD, type = "l", col = "blue", lwd = 2,
     xlab = "Frame", ylab = "Ligand RMSD (Å)", main = "Ligand RMSD Over Time (ADN)",
     xlim = c(0, 501))  # Ensures the plot stops at the last frame


### ADN protein ##############################################################

# Read data and assign to rmsd_data2 directly
rmsd_data2 <- read.table("protein_rmsd_ADN.dat", header = FALSE)
colnames(rmsd_data2) <- c("Frame", "RMSD")

# Convert columns to numeric (in case they are characters)
rmsd_data2$Frame <- as.numeric(rmsd_data2$Frame)
rmsd_data2$RMSD <- as.numeric(rmsd_data2$RMSD)

# Plot
plot(rmsd_data2$Frame, rmsd_data2$RMSD, type = "l", col = "blue", lwd = 2,
     xlab = "Frame", ylab = "Protein RMSD (Å)", main = "Protein RMSD Over Time (ADN)")



### MOD1 Ligand ################################################################


# Reload data with correct column names
rmsd_data3 <- read.table("MOD1_ligand_rmsd.dat", header = FALSE)
colnames(rmsd_data3) <- c("Frame", "RMSD")  # Rename columns
mean_rmsd <- mean(rmsd_data3$RMSD)
print(mean_rmsd)


# Verify column names
print(colnames(rmsd_data3))  # Should show "Frame" and "RMSD"

# Check min/max values
min(rmsd_data3$Frame)
max(rmsd_data3$Frame)

plot(rmsd_data3$Frame, rmsd_data3$RMSD, type = "l", col = "blue", lwd = 2,
     xlab = "Frame", ylab = "Ligand RMSD (Å)", main = "Ligand RMSD Over Time (MOD1)",
     xlim = c(0, 501))  # Ensures the plot stops at the last frame


### MOD1 protein ##############################################################

rmsd_data4 <- read.table("protein_rmsd_mod1.dat", header = FALSE)
colnames(rmsd_data4) <- c("Frame", "RMSD")

# Convert columns to numeric (in case they are character)
rmsd_data4$Frame <- as.numeric(rmsd_data4$Frame)
rmsd_data4$RMSD <- as.numeric(rmsd_data4$RMSD)

# Now plot
plot(rmsd_data4$Frame, rmsd_data4$RMSD, type = "l", col = "blue", lwd = 2,
     xlab = "Frame", ylab = "Protein RMSD (Å)", main = "Protein RMSD Over Time (MOD1)")


### MOD2 Ligand ################################################################


# Reload data with correct column names
rmsd_data5 <- read.table("MOD2_ligand_rmsd.dat", header = FALSE)
colnames(rmsd_data5) <- c("Frame", "RMSD")  # Rename columns
mean_rmsd <- mean(rmsd_data5$RMSD)
print(mean_rmsd)


# Verify column names
print(colnames(rmsd_data5))  # Should show "Frame" and "RMSD"

# Check min/max values
min(rmsd_data5$Frame)
max(rmsd_data5$Frame)

plot(rmsd_data5$Frame, rmsd_data5$RMSD, type = "l", col = "blue", lwd = 2,
     xlab = "Frame", ylab = "Ligand RMSD (Å)", main = "Ligand RMSD Over Time (MOD2)",
     xlim = c(0, 501))  # Ensures the plot stops at the last frame


### MOD2 protein ##############################################################

rmsd_data6 <- read.table("protein_rmsd_mod2.dat", header = FALSE)
colnames(rmsd_data6) <- c("Frame", "RMSD")

# Convert columns to numeric (in case they are character)
rmsd_data6$Frame <- as.numeric(rmsd_data6$Frame)
rmsd_data6$RMSD <- as.numeric(rmsd_data6$RMSD)

# Now plot
plot(rmsd_data6$Frame, rmsd_data6$RMSD, type = "l", col = "blue", lwd = 2,
     xlab = "Frame", ylab = "Protein RMSD (Å)", main = "Protein RMSD Over Time (MOD2)")


## Smoothed contacts  ######################################################
library(dplyr)
library(zoo)
library(ggplot2)

### Summarized RMSD ligand ################################################

# Fix ADN protein (rmsd_data2)
rmsd_data2_clean <- rmsd_data2[-1, ]
colnames(rmsd_data2_clean) <- c("Frame", "RMSD")
rmsd_data2_clean$Frame <- as.numeric(rmsd_data2_clean$Frame)
rmsd_data2_clean$RMSD <- as.numeric(rmsd_data2_clean$RMSD)
rmsd_data2_clean <- mutate(rmsd_data2_clean, System = "ADN")


# --- Combine ligand RMSD data ---
ligand_rmsd <- bind_rows(
  mutate(rmsd_data, System = "ADN"),
  mutate(rmsd_data3, System = "MOD1"),
  mutate(rmsd_data5, System = "MOD2")
)

# --- Plot Ligand RMSD ---
p_lig <- ggplot(ligand_rmsd, aes(x = Frame, y = RMSD, color = System)) +
  geom_line(linewidth = 1.2) +
  labs(
    title = "Ligand RMSD Over Time",
    x = "Frame",
    y = "Ligand RMSD (Å)"
  ) +
  scale_color_manual(values = c("ADN" = "#b2182b", "MOD1" = "#762a83", "MOD2" = "#2166ac")) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold"),
    legend.title = element_blank()
  )
### Summarized RMSD ligand ################################################

# Fix ADN protein (rmsd_data2)
rmsd_data2_clean <- rmsd_data2[-1, ]
colnames(rmsd_data2_clean) <- c("Frame", "RMSD")
rmsd_data2_clean$Frame <- as.numeric(rmsd_data2_clean$Frame)
rmsd_data2_clean$RMSD <- as.numeric(rmsd_data2_clean$RMSD)
rmsd_data2_clean <- mutate(rmsd_data2_clean, System = "ADN")

# --- Combine protein RMSD data ---
protein_rmsd <- bind_rows(
  mutate(rmsd_data2_clean, System = "ADN"),
  mutate(rmsd_data4, System = "MOD1"),
  mutate(rmsd_data6, System = "MOD2")
)

# --- Plot Protein RMSD ---
p_prot <- ggplot(protein_rmsd, aes(x = Frame, y = RMSD, color = System)) +
  geom_line(linewidth = 1.2) +
  labs(
    title = "Protein RMSD Over Time",
    x = "Frame",
    y = "Protein RMSD (Å)"
  ) +
  scale_color_manual(values = c("ADN" = "#b2182b", "MOD1" = "#762a83", "MOD2" = "#2166ac")) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold"),
    legend.title = element_blank()
  )

# --- Display both plots ---
print(p_lig)
print(p_prot)

