
# Reload data with correct column names
rmsd_data <- read.table("lig_rmsd.dat", header = FALSE)
colnames(rmsd_data) <- c("Frame", "RMSD")  # Rename columns
mean_rmsd <- mean(rmsd_data$RMSD)
print(mean_rmsd)


# Verify column names
print(colnames(rmsd_data))  # Should show "Frame" and "RMSD"

# Check min/max values
min(rmsd_data$Frame)
max(rmsd_data$Frame)

plot(rmsd_data$Frame, rmsd_data$RMSD, type = "l", col = "blue", lwd = 2,
     xlab = "Frame", ylab = "Ligand RMSD (Å)", main = "Ligand RMSD Over Time",
     xlim = c(0, 1021))  # Ensures the plot stops at the last frame

rmsd_data <- read.table("protein_rmsd.dat", header = FALSE)
colnames(rmsd_data) <- c("Frame", "RMSD")

plot(rmsd_data$Frame, rmsd_data$RMSD, type = "l", col = "blue", lwd = 2,
     xlab = "Frame", ylab = "Protein RMSD (Å)", main = "Protein RMSD Over Time")

### CONTACTS ###################################################################
###############################################################################

# Load required libraries
library(readr)
library(zoo)
library(ggplot2)
library(reshape2)
library(dplyr)

contacts_data <- read_tsv("final_contacts.tsv", col_names = FALSE)
colnames(contacts_data) <- c("Residue", "Frame", "Contacts")

contacts_data$Frame <- as.numeric(contacts_data$Frame)
contacts_data$Contacts <- as.numeric(contacts_data$Contacts)

# Summarize total contacts per frame
contact_summary <- contacts_data %>%
  group_by(Frame) %>%
  summarize(TotalContacts = sum(Contacts))

# Apply smoothing
contact_summary$Smoothed <- rollmean(contact_summary$TotalContacts, k = 10, fill = NA)

# Plot raw + smoothed
plot(contact_summary$Frame, contact_summary$TotalContacts, type = "l", col = "gray", lwd = 1,
     xlab = "Frame", ylab = "Total Protein-Ligand Contacts", main = "Contacts Over Time")
lines(contact_summary$Frame, contact_summary$Smoothed, col = "red", lwd = 2)
legend("topright", legend = c("Raw Contacts", "Smoothed Contacts"), col = c("gray", "red"), lty = 1, lwd = 2)

# === Heatmap ===
# Pivot to wide format (Residue = columns)
contacts_wide <- dcast(contacts_data, Frame ~ Residue, value.var = "Contacts", fill = 0)

# Convert to long format for ggplot and filter out zero-contact entries
contacts_long <- melt(contacts_wide, id.vars = "Frame", variable.name = "AA", value.name = "Contacts") %>%
  filter(Contacts > 0)

# Plot heatmap
ggplot(contacts_long, aes(x = Frame, y = AA, fill = Contacts)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red", name = "Contacts") +
  labs(title = "Residues Involved in ADN Binding",
       x = "Frame", 
       y = "Amino Acid",
       fill = "Contacts") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, size = 8),
    axis.text.y = element_text(size = 10, face = "bold"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

