start<-Sys.time()
rep_time <- 100
n<-5
RMSE<-array(NA,dim = c(rep_time,6))

source('MGP_stoch.R')
source('MGP_stoch(unpenalized).R')
source('Minimal Transfer.R')
source('Multioutput with identity.R')
source('Single SK.R')
source('FGP_19.R')

colMeans(RMSE) 
end<-Sys.time()
end-start



#Code below are just for generating boxplot, the simulation has been done
par(mfrow = c(1,1))
par(mar = c(8, 4.5, 2, 2)) 
colnames(RMSE)<-c("Proposed","Unpenalized","Minimal Transfer","i.i.d MGP","Individualized SK","i.i.d Enforced")
bp<-boxplot(RMSE, xlab = "", ylab = "RMSE", col = "lightblue",xaxt = "n",ylim = c(0.2,2))
tick <- seq_along(bp$names)
axis(1, at = tick, labels = FALSE)
text(tick, par("usr")[3]-0.3 , bp$names, srt = 45, xpd = TRUE)


# Load the ggplot2 library
library(ggplot2)

#Create sample data frames

data1<-data.frame(RMSE)
# Combine the two data frames
combined_data <- cbind(data1)
colnames(combined_data)<-c("A","B","C","D","E","F")

# Convert the data to long format
combined_data_long <- tidyr::pivot_longer(combined_data, cols = c(A,B,C,D,E,F), names_to = "Variable", values_to = "RMSE")

# Create the boxplot with parallel boxplots
ggplot(combined_data_long, aes(x = Variable, y = RMSE)) +
  geom_boxplot(fill = "skyblue", width = 0.8) + # Adjust width parameter
  theme_minimal()+
  xlab(NULL)+
  coord_cartesian(ylim = c(0.2, 2.0))+ # Edit legend name here
  scale_x_discrete(labels = c("Proposed", "Unpenalized", "Minimal Transfer", "i.i.d. MGP", "Individualized SK", "i.i.d Enforced"))+ 
  theme(axis.title.x = element_text(size = 12),  # Increase x-axis label size
        axis.title.y = element_text(size = 12),  # Increase y-axis label size
        axis.text.x = element_text(size = 12, angle = 40, hjust = 1),   # Increase x-axis tick size
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        text = element_text(family = "Times"))



