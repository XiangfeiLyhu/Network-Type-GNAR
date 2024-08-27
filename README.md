# Network-Type-GNAR
# Overview
This repository contains the code and data used to reproduce the results of the master's thesis titled "Exploring the Impact of Network Type and Characteristics on the Accuracy of Network Time Series Forecasting". The study investigates how different network types and their characteristics affect the accuracy of time series forecasting, specifically in the context of foot-and-mouth disease outbreaks.

## Data
The repository includes several important datasets:

Xiangfei.net_2002.RData: This file contains a time series dataset representing the number of infected farms during a foot-and-mouth disease outbreak. The dataset consists of 2002 nodes (farms) and 230 time stamps.

Xiangfei.RData: This file includes several key components:

Xiangfei.net: The network structure associated with the 2002 outbreak.
Xiangfei.coords: A set of coordinates for plotting the network.
Xiangfei.map: A function to plot the edges and locations of the network.
getSegET: A subsidiary function used within the plotting process.
Xiangfei.farmsize: Additional data related to farm sizes, used in the analysis.
net_1967.RData: This file includes the network structure (GNARnet) for the 1967 foot-and-mouth disease outbreak.

## Code
### R Scripts
pred_perf.R: This script contains the code for all analyses performed in the predictive performance section of the thesis.

GNARsim_t_noise.R: This script modifies the GNARsim function from the GNAR package by changing the innovation distribution from normal to t-distribution. It is used to simulate network time series data under different assumptions.

foot-and-mouth_network_building.R: This script includes the code necessary to create the network structures for the 2001 and 1967 foot-and-mouth disease outbreaks.

### Python Scripts
FMD_data_simulation.py: This Python script uses the NGNAR model to simulate time series data for the 2001 and 1967 foot-and-mouth disease outbreaks. It also includes model selection and statistical analysis for the NGNAR model.

NGNAR Net 2: This file contains all the code related to the NGNAR model, which is used for forecasting in the study.

## Usage Instructions
Loading Data
Before running any analyses, you must load the necessary RData files. This can be done by running the following command in your R session:

load("Xiangfei.RData")


Required Libraries
Make sure you have the following R libraries installed. You can install any missing libraries using 

install.packages("package_name").


library(maps)
library(mapdata)
library(parallel)
options(mc.cores=6)


Running the Code

After loading the required data and libraries, you can generate the network plots using the following command:

Xiangfei.map(lapplyfn=mclapply)

This command uses the loaded network data and plotting functions to visualize the foot-and-mouth disease outbreak network.

Predictive Performance Analysis
To run the predictive performance analysis, execute the pred_perf.R script in your R environment.

GNAR Simulation
If you want to run the GNAR simulations with a t-distribution, use the GNARsim_t_noise.R script.

Network Building
For recreating the 2001 and 1967 network structures, run the foot-and-mouth_network_building.R script.

Python Simulations
For time series simulations and model selection using the NGNAR model, run the FMD_data_simulation.py script in your Python environment.

Notes
Please do not distribute this information further.
The code is intended for academic and research purposes only.
Contact
If you have any questions or need further assistance, please contact the author of the thesis.

