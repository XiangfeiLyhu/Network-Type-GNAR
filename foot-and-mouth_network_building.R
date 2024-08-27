#############################################################################
############## UK foot-and-mouth disease outbreak network ##############
#############################################################################
library(maps)
library(mapdata)
library(parallel)
library(GNAR)

options(mc.cores=6)
load("/Users/lvxiangfei/Desktop/Xiangfei.RData")
load("/Users/lvxiangfei/Desktop/FMD_MVTS.RData")


####################################################################################
################# Create 1967 Network and simulate the data  #######################
####################################################################################
library(sp)
library(geodist)
library(igraph)
library(GNAR)

# Double the number of farms by simulating new locations around existing ones
set.seed(123)
locations <- as.data.frame(Xiangfei.coords)
new_locations <- data.frame(Longitude = numeric(), Latitude = numeric())
for (i in 1:nrow(locations)) {
  new_locations <- rbind(new_locations, locations[i, ])
  new_lon <- rnorm(1, mean = locations$Longitude[i], sd = 0.1)
  new_lat <- rnorm(1, mean = locations$Latitude[i], sd = 0.1)
  new_locations <- rbind(new_locations, data.frame(Longitude = new_lon, Latitude = new_lat))
}


# Plot original and new locations for verification
plot(new_locations$Longitude, new_locations$Latitude, col = "blue", pch = 16, cex = 0.6, main = "Farm Locations")
points(locations$Longitude, locations$Latitude, col = "red", pch = 16, cex = 0.6)
legend("topright", legend = c("2001 Farm Location", "1967 Farm Location"), col = c("red", "blue"), pch = 16, cex = 1, pt.cex = 0.7, text.width = 1)

# Calculate the distance matrix
####### Way 1
distance_matrix <- geodist(new_locations, measure = 'haversine' )/1000 #converting it to km

# Create edges for any two farms less than 30 km apart
threshold <- 30 # 30 km 

# Initialize adjacency matrix
n_farms <- nrow(new_locations)
adj_matrix <- matrix(0, nrow = n_farms, ncol = n_farms)

# Randomly pick 2 columns for each row where the distance is less than 30 km
set.seed(123) # For reproducibility

for (i in 1:n_farms) {
  # Get the indices of farms within 30 km distance
  nearby_farms <- which(distance_matrix[i, ] < threshold & distance_matrix[i, ] > 0)
  
  # If there are fewer than 2 farms within 30 km, select all
  if (length(nearby_farms) > 1) {
    # Randomly pick 2 farms from the nearby farms
    selected_farms <- sample(nearby_farms, size = min(2, length(nearby_farms)), replace = FALSE)
    adj_matrix[i, selected_farms] <-  distance_matrix[i, selected_farms]
    adj_matrix[selected_farms, i] <- distance_matrix[selected_farms, i] # Ensure the adjacency is symmetric
  } else if (length(nearby_farms) == 1) {
    adj_matrix[i, nearby_farms] <- distance_matrix[i, nearby_farms]
    adj_matrix[nearby_farms, i] <- distance_matrix[nearby_farms, i]
  }
}
###### CHECK
dim(adj_matrix) # should be 4036 * 4036
net_1967 <- matrixtoGNAR(adj_matrix)
save(net_1967, file = "net_1967.RData")


net_1967 <- graph.adjacency(adj_matrix, mode = "undirected")
net_1967 <- igraphtoGNAR(net_1967)

library(maps)
library(mapdata)
library(parallel)
options(mc.cores=6)

Xiangfei.map(net = net_1967, loc.coords = new_locations)
Xiangfei.map(net = Xiangfei.net, loc.coords = Xiangfei.coords)

#################### Generate a new farm size distribution. ###############
# New farm size distribution with mean half of the 2001 distribution
mean_2001 <- mean(Xiangfei.farmsize)
sd(Xiangfei.farmsize)

# Create the density plot using ggplot2
ggplot(as.data.frame(Xiangfei.farmsize), aes(x = farm_size)) +
  geom_density(fill = "blue", alpha = 0.5) +
  labs(title = "Density Plot of farm size in 2001",
       x = "Farm Size",
       y = "Density") +
  theme_minimal()



# Estimate the density of the original data
density_est <- density(Xiangfei.farmsize, na.rm = TRUE)
# Calculate the mean of the original data
original_mean <- mean(Xiangfei.farmsize, na.rm = TRUE)
# Desired mean
desired_mean <- original_mean / 2

# Scale factor to adjust the mean
scale_factor <- desired_mean / original_mean

# Simulate new data points by scaling the original data
simulated_data <- Xiangfei.farmsize * scale_factor

if (length(simulated_data) < 4036) {
  simulated_data <- sample(simulated_data, 4036, replace = TRUE)
} else if (length(simulated_data) > 4036) {
  simulated_data <- sample(simulated_data, 4036)
}

new_farm_1967 <- round(simulated_data, 0)
# Calculate the density of the simulated data
simulated_density <- density(new_farm_1967, na.rm = TRUE)
mean(new_farm_1967)
length(new_farm_1967)
# Create the density plot using ggplot2
ggplot() +
  geom_density(aes(x = Xiangfei.farmsize, fill = "2001 farm size"), alpha = 0.5) +
  geom_density(aes(x = new_farm_1967, fill = "1967 farm size"), alpha = 0.5) +
  labs(title = "Density Plot of 2001 and Simulated 1967 farm size",
       x = "Farm Size",
       y = "Density") +
  scale_fill_manual(name = "Legend", values = c("2001 farm size" = "blue", "1967 farm size" = "red")) +
  theme_minimal() +
  scale_x_continuous(limits = c(0, max(Xiangfei.farmsize, na.rm = TRUE))) +
  theme(plot.title = element_text(hjust = 0.5))

#save(simulated_data, file="/Users/lvxiangfei/Desktop/farmsize_1967.csv")
new_farm_1967 <- round(simulated_data, 0)
#write.csv(new_farm_1967, file = "/Users/lvxiangfei/Desktop/farmsize_1967.csv", row.names = FALSE)


####################################################################################
################# True Time series GNAR(1,[1]) Fitting #############################
####################################################################################

#########################  Create a new Network for 2002 nodes ########################
# Function to find the missing number
find_missing_nodes <- function(x) {
  full_sequence <- 1:2018
  missing_number <- setdiff(full_sequence, x)
  return(missing_number)
}

true_ts_nodes <- dimnames(FMD.MVTS)[[2]]
# Find the missing number
missing_nodes <- find_missing_nodes(true_ts_nodes)
print(missing_nodes)


head(Xiangfei.coords)
tail(Xiangfei.coords)

length(Xiangfei.coords)/2 # 2018
Xiangfei.coords_2002 <- Xiangfei.coords
for (i in indices_to_delete){
  Xiangfei.coords_2002 <- Xiangfei.coords_2002[-i, ]
}
print(length(Xiangfei.coords_2002)/2) # 2002

##### convert coords to km distance
library(geodist)
# Calculate the distance matrix for Xiangfei.coords_2002
distance_matrix1 <-geodist(Xiangfei.coords_2002, measure = 'haversine' )/1000 

# Display the distance matrix
print(dim(distance_matrix1))

# Initialize adjacency matrix
n_farms <- nrow(Xiangfei.coords_2002)
adj_matrix1 <- matrix(0, nrow = n_farms, ncol = n_farms)

# Randomly pick 2 columns for each row where the distance is less than 30 km
set.seed(123) # For reproducibility
threshold<-30
for (i in 1:n_farms) {
  # Get the indices of farms within 30 km distance
  nearby_farms <- which(distance_matrix1[i, ] < threshold & distance_matrix1[i, ] > 0)
  
  # If there are fewer than 2 farms within 30 km, select all
  if (length(nearby_farms) > 1) {
    # Randomly pick 2 farms from the nearby farms
    selected_farms <- sample(nearby_farms, size = min(2, length(nearby_farms)), replace = FALSE)
    adj_matrix1[i, selected_farms] <- distance_matrix1[i, selected_farms]
    adj_matrix1[selected_farms, i] <- distance_matrix1[selected_farms, i] # Ensure the adjacency is symmetric
  } else if (length(nearby_farms) == 1) {
    adj_matrix1[i, nearby_farms] <- distance_matrix1[i, nearby_farms]
    adj_matrix1[nearby_farms, i] <- distance_matrix1[nearby_farms, i]
  }
}


dim(adj_matrix1)
Xiangfei.net_2002 <- matrixtoGNAR(adj_matrix1)
save(Xiangfei.net_2002, file= "Xiangfei.net_2002.RData")

Xiangfei.net_2002 <- graph.adjacency(adj_matrix1, mode = "undirected", weighted = TRUE)
Xiangfei.net_2002 <- igraphtoGNAR(Xiangfei.net_2002)
is.GNARnet(Xiangfei.net_2002)
Xiangfei.map(net = Xiangfei.net_2002, loc.coords = Xiangfei.coords_2002)




