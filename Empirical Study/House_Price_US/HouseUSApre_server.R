
####### 1. House Index #########
HouseUSA <- read.csv("HouseUSA.csv")

# Extract the relevant columns
Year <- HouseUSA$Year
Month <- as.character(HouseUSA$Month)

GEO_Name <- unique(HouseUSA$GEO_Name)
date <- as.Date(paste0(Year, "-", ifelse(nchar(Month) == 1, paste0("0", Month), Month), "-01"))
date_all <- unique(date)

HouseUSA[,1] = date

NSA <- matrix(NA, nrow = length(GEO_Name), ncol = length(date_all ) )
SA <- matrix(NA, nrow = length(GEO_Name), ncol = length(date_all ) )


# Create the NSA matrix
for (i in 1:length(GEO_Name)) {
  geo_rows = HouseUSA[HouseUSA$GEO_Name == GEO_Name[i],]
  NSA[i,] = geo_rows$Index_NSA
  SA[i,] = geo_rows$Index_SA
}

rownames(NSA) = GEO_Name
colnames(NSA) = as.character(date_all)
rownames(SA) = GEO_Name
colnames(SA) = as.character(date_all)

logdiff = function(matt,dolog = TRUE,dopercent=FALSE){
  matt_diff = matt[,-1]
  if(dolog == TRUE && dopercent==FALSE){
    for (i in 1:nrow(matt)) {
      matt_diff[i,] = diff(log(matt[i,]))
    }
  }else if(dopercent==FALSE && dolog == FALSE){
    for (i in 1:nrow(matt)) {
      matt_diff[i,] = diff(matt[i,])
    }
  }else if(dolog == FALSE && dopercent==TRUE){
    for (i in 1:nrow(matt)) {
      matt_diff[i,] = diff(matt[i,]) / matt[i,1:ncol(matt_diff)]
    }
  }
  return(matt_diff)
}

NSA_diff = logdiff(NSA)

SA_diff = logdiff(SA)

is_na = function(input){
  any_na_values <- any(is.na(input))
  
  # print the result
  if (any_na_values) {
    print("There are NA values in the new_output matrix.")
  } else {
    print("There are no NA values in the new_output matrix.")
  }
}

is_na(NSA_diff)

####### 2. CPI###########
##### 2.1 Extract CPI #####
library(readxl)

# set directory path
dir_path <- "/home/zhpu/Panel CH Test/empricial/Data/USHouse/CPI_USA/"

# get list of file names in directory
file_names <- list.files(path = dir_path, pattern = "*.xlsx")

# create empty list to store data frames
df_list <- list()

# loop through file names and read XLSX files
for (file_name in file_names) {
  file_path <- paste0(dir_path, "/", file_name)
  df <- read_excel(file_path)
  area_name <- gsub(",", "", df[6, 2])
  df <- read_excel(file_path, skip = 11)
  df_list[[area_name]] <- df
}


# loop through data frames and extract quarterly data
for (i in seq_along(df_list)) {
  df_out <- data.frame()
  df <- df_list[[i]]
  year_col <- which(colnames(df) == "Year")
  month_cols <- which(colnames(df) %in% month.abb)
  half1_col <- which(colnames(df) == "HALF1")
  half2_col <- which(colnames(df) == "HALF2")
  annual_col <- which(colnames(df) == "Annual")
  
  # create datetime columns for quarters
  col_indices <- month_cols
  col_names <- colnames(df)[col_indices]
  col_indices <- (col_indices - 1) * 4 + 1
  col_names <- format(seq(as.Date(paste0(df[1, year_col], "-01-01")), as.Date(paste0(df[nrow(df), year_col], "-12-31")), by = "quarter"), "%Y-%m-%d")
  #col_names <- col_names[col_indices]
  
  for (j in 1:nrow(df)) {
    if (sum(!is.na(df[j, month_cols])) >=4) {
      for (k in 1:4) {
        col_current_name <- col_names[(j*4-3):(j*4)]
        start_month <- (k - 1) * 3 + 1
        end_month <- k * 3 
        quarter_values <- as.numeric(df[j, month_cols[start_month:end_month]])
        quarter_avg <- mean(quarter_values, na.rm = TRUE)
        df_out[1, col_current_name[k]] <- quarter_avg
      }
    } else if (all(is.na(df[j, month_cols]))) {
      col_current_name <- col_names[(j*4-3):(j*4)]
      df_out[1, col_current_name[1]] <- df[j, half1_col]
      df_out[1, col_current_name[2]] <- df[j, half1_col]
      df_out[1, col_current_name[3]] <- df[j, half2_col]
      df_out[1, col_current_name[4]] <- df[j, half2_col]
    } else {
      annual_value <- df[j, annual_col]
      col_current_name <- col_names[(j*4-3):(j*4)]
      for (k in 1:4) {
        df_out[1, col_current_name[k]] <- annual_value 
      }
    }
  }
  
  df_list[[i]] <- df_out
}

##### 2.2 CPI MSA #####
# create CPI_MSA matrix
start_year <- 1975
end_year <- 2014
num_quarters <- (end_year - start_year + 1) * 4
CPI_MSA <- matrix(NA, nrow = length(GEO_Name), ncol = num_quarters)

# extract column names for quarters only
col_names <- format(seq(as.Date(paste0(start_year, "-01-01")), as.Date(paste0(end_year, "-12-31")), by = "quarter"), "%Y-%m-%d")
col_indices <- which(substr(col_names, 6, 7) %in% c("01", "04", "07", "10"))
col_names <- col_names[col_indices]
col_indices <- (col_indices - 1) * 4 + 1

# assign column names to CPI_MSA matrix
colnames(CPI_MSA) <- col_names
rownames(CPI_MSA) <- GEO_Name

# loop through data frames and fit data into CPI_MSA matrix
start_col <- which(colnames(df_list$`U.S. city average`) == col_names[1])
end_col <- which(colnames(df_list$`U.S. city average`) == col_names[length(col_names)])

for (i in 1:nrow(CPI_MSA)) {
  CPI_MSA[i, ] <- as.matrix(df_list$`U.S. city average`[1,c(start_col:end_col)])
}

ss = 0
for (i in seq_along(df_list)) {
  df <- df_list[[i]]
  area_name <- names(df_list)[i]
  
  if (area_name %in% GEO_Name) {
    ss = ss+ 1
    row_index <- which(GEO_Name == area_name)
  
    intersect_date <- intersect(colnames(df),col_names)
    
    CPI_MSA[row_index, intersect_date] = as.matrix(df[1,intersect_date])
  }
}


####### 3. real House Index###########
Real_NSA <- NSA[,col_names]
Real_SA <- SA[,col_names]

Real_NSA <- logdiff(Real_NSA / CPI_MSA)
Real_SA <- logdiff(Real_SA / CPI_MSA)

####### 4. Population ###########
Population <- read.csv("Population.csv", header = TRUE, stringsAsFactors = FALSE)
Population <- Population[1:385,]

# remove parentheses and their contents from GeoName column
Population$GeoName <- gsub(" \\(.*\\)", "", Population$GeoName)

# remove commas from GeoName column
Population$GeoName <- gsub(",", "", Population$GeoName)

Population$GeoName <- gsub("\\s+\\*$", "", Population$GeoName)

# set row names to GeoName column
rownames(Population) <- Population$GeoName

# remove GeoName column from data frame
Population <- Population[, -1]

# set column names to the first row starting from the third column
Years <- colnames(Population)[-1]
Years <- gsub("X", "", Years)

Population <- Population[,-1]
colnames(Population) <- Years



# convert data frame to matrix
Population <- data.matrix(Population)

######4.1 interpolation######
library(zoo)

# Convert the character column names to quarterly dates
quarterly_dates <- as.yearqtr(paste0(colnames(Population), "-01-01"))

# Set the column names of the Population matrix to the quarterly dates
colnames(Population) <- as.character(quarterly_dates)

# Create a sequence of quarterly dates from 1975Q1 to 2014Q4
quarterly_dates_interp <- seq(as.yearqtr("1975Q1"), as.yearqtr("2014Q4"), by = 1/4)

# Convert the quarterly dates to character strings in the desired format
colnames_interp <- format(as.Date(as.yearqtr(quarterly_dates_interp)), "%Y-%m-%d")

# Interpolate the population numbers using the quarterly dates
Population_interpolation <- matrix(NA, nrow = nrow(Population), ncol =length(colnames_interp))

for (i in 1:nrow(Population)) {
  Population_zoo <- zoo(as.vector(t(Population[i,])), order.by = quarterly_dates)
  Population_interpolation[i,] <- na.approx(Population_zoo, xout = quarterly_dates_interp)
}

rownames(Population_interpolation) = rownames(Population)
colnames(Population_interpolation) = colnames_interp

# Convert the interpolated population matrix back to a regular matrix and set the new column names
Population_interpolation <- as.matrix(Population_interpolation)

MSA_Population <- intersect(rownames(Population_interpolation) , rownames(Real_NSA))

percentchange_Population <- logdiff(Population_interpolation[MSA_Population,],dolog = FALSE, dopercent = TRUE)


#######5. Income ###########
Income <- read.csv("Income.csv", header = TRUE, stringsAsFactors = FALSE)
Income <- Income[1:385,]

# remove parentheses and their contents from GeoName column
Income$GeoName <- gsub(" \\(.*\\)", "", Income$GeoName)

# remove commas from GeoName column
Income$GeoName <- gsub(",", "", Income$GeoName)

Income$GeoName <- gsub("\\s+\\*$", "", Income$GeoName)

# set row names to GeoName column
rownames(Income) <- Income$GeoName

# remove GeoName column from data frame
Income <- Income[, -1]

# set column names to the first row starting from the third column
Years <- colnames(Income)[-1]
Years <- gsub("X", "", Years)

Income <- Income[,-1]
colnames(Income) <- Years

# convert data frame to matrix
Income <- data.matrix(Income)

######5.1 interopolation######
quarterly_dates <- as.yearqtr(paste0(colnames(Income), "-01-01"))

# Set the column names of the Population matrix to the quarterly dates
colnames(Income) <- as.character(quarterly_dates)

# Create a sequence of quarterly dates from 1975Q1 to 2014Q4
quarterly_dates_interp <- seq(as.yearqtr("1975Q1"), as.yearqtr("2014Q4"), by = 1/4)

# Convert the quarterly dates to character strings in the desired format
colnames_interp <- format(as.Date(as.yearqtr(quarterly_dates_interp)), "%Y-%m-%d")

# Interpolate the population numbers using the quarterly dates
Income_interpolation <- matrix(NA, nrow = nrow(Income), ncol =length(colnames_interp))

for (i in 1:nrow(Income)) {
  Income_zoo <- zoo(as.vector(t(Income[i,])), order.by = quarterly_dates)
  Income_interpolation[i,] <- na.approx(Income_zoo, xout = quarterly_dates_interp)
}

rownames(Income_interpolation) = rownames(Income)
colnames(Income_interpolation) = colnames_interp

# Convert the interpolated population matrix back to a regular matrix and set the new column names
Income_interpolation <- as.matrix(Income_interpolation)

MSA_Income <- intersect(rownames(Income_interpolation) , rownames(Real_NSA))

percentchange_Income <- logdiff(Income_interpolation[MSA_Income,],dolog = FALSE, dopercent = TRUE)

##########Real House Index############
Real_NSA <- Real_NSA[MSA_Income,]
Real_SA <- Real_SA[MSA_Income,]

ymat <- Real_NSA
xmat1 <- percentchange_Population
xmat2 <- percentchange_Income

seasondiff <- function(mat) {
  # Get the number of rows and columns in the input matrix
  n_rows <- nrow(mat)
  n_cols <- ncol(mat)
  
  # Create an output matrix of the correct dimensions
  out_mat <- mat[,5:n_cols]
  
  # Loop over each row of the input matrix
  for (i in 1:n_rows) {
    # Loop over each column of the output matrix
    for (j in 1:(n_cols-4)) {
      # Compute the seasonal difference
      out_mat[i, j] <- mat[i, (j+4)] - mat[i, (j)]
    }
  }
  
  # Return the output matrix
  return(out_mat)
}

ymat <- seasondiff(Real_NSA)
xmat2 <- seasondiff(percentchange_Income)
xmat1 <- seasondiff(percentchange_Population)

#### 6. W ######
USCity <- read.csv("uscities.csv")
USCity <- as.data.frame(USCity)

MSA_info <- matrix(0, nrow = length(MSA_Income), ncol = 3, dimnames = list(MSA_Income, c("lat", "lng", "num_cities")))

# loop through each MSA name and compute the latitude and longitude


for (z in 1:length(MSA_Income)) {
  latitude <- c()
  longitude <- c()
  
  # split the MSA name into its components
  split_string <- strsplit(MSA_Income[z], " ")[[1]]
  split_index <- length(split_string)
  left_side <- paste(split_string[1:(split_index-1)], collapse = " ")
  right_side <- split_string[split_index]
  
  # state names
  state_names <- strsplit(right_side,  "--|-" )[[1]]
  
  # check if there are any county names in the split
  county_index <- grep("County", left_side, ignore.case = TRUE)
  if (length(county_index) > 0) {
    # if there is a county name, split it by the "/"'

    county.name <- strsplit(left_side[county_index], "-|/" )[[1]][2]
    left_side <-  strsplit(left_side[county_index], "-|/" )[[1]][1]
    
    county.name <- strsplit(county.name, " ")[[1]][1]
    
    county.state <- state_names[length(state_names)]
    
    state_names <- state_names[1] # the state names for city but not for county
    
    county_data <- subset(USCity, county_name == county.name & state_id == county.state)
    
    # extract the latitude and longitude values from the subset
    if(nrow(county_data) ==1){
      latitude[1] <- county_data$lat
      longitude[1] <- county_data$lng
    }else if (nrow(county_data) >=1){
      county_data <- subset(USCity, city == left_side & county_name == county.name & state_id == county.state)
      latitude[1] <- county_data$lat
      longitude[1] <- county_data$lng
    }

    county_num <- 1
  } else {
    county.name <- ""
    county_num <- 0
  }
  
  # loop through each split and determine if it's a city or state
  city_names <- strsplit(left_side, "-" )[[1]]
  
  if(length(state_names) == length(city_names)){ # for city = states
    for (i in 1:length(state_names)) {
      city_data <-subset(USCity, city == city_names[i] & state_id == state_names[i])
      if(nrow(city_data) != 0){
        latitude = c(latitude,city_data$lat)
        longitude = c(longitude,city_data$lng)
        county_num <- county_num + 1
      }
    }
  }else if(length(state_names) > length(city_names)){  # for city < states
    if (length(city_names) == 1){  # for 1 city, many states
      for (i in 1:length(state_names)) {
        city_data <-subset(USCity, city == city_names[1] & state_id == state_names[i])
        if(nrow(city_data) != 0){
          latitude = c(latitude,city_data$lat)
          longitude = c(longitude,city_data$lng)
          county_num <- county_num + 1
        }
      }
    }else{ # for some cities, more states
      for (i in 1:length(state_names)) {
        if(i == 1 ){
          city_data <-subset(USCity, city == city_names[1] & state_id == state_names[1])
          if(nrow(city_data) != 0){
            latitude = c(latitude,city_data$lat)
            longitude = c(longitude,city_data$lng)
            county_num <- county_num + 1
          }
        }else if( i == length(state_names)){
          city_data <-subset(USCity, city == city_names[length(city_names)] & state_id == state_names[i])
          if(nrow(city_data) != 0){
            latitude = c(latitude,city_data$lat)
            longitude = c(longitude,city_data$lng)
            county_num <- county_num + 1
          }
        }else{
          if(length(city_names) >=3){
            city_data <-subset(USCity, city == city_names[2] & state_id == state_names[i])
            if(nrow(city_data) != 0){
              latitude = c(latitude,city_data$lat)
              longitude = c(longitude,city_data$lng)
              county_num <- county_num + 1
            }              
          }  
        }
      }
      
    }
  }else if(length(state_names) < length(city_names)){ # city > state
    if (length(state_names) == 1){  # for many city, 1 states
      for (i in 1:length(city_names)) {
        city_data <-subset(USCity, city == city_names[i] & state_id == state_names[1])
        if(nrow(city_data) != 0){
          latitude = c(latitude,city_data$lat)
          longitude = c(longitude,city_data$lng)
          county_num <- county_num + 1
        }
      }
    }else{ # for more cities, some states
      for (i in 1:length(city_names)) {
        if(i == 1 ){
          city_data <-subset(USCity, city == city_names[1] & state_id == state_names[1])
          if(nrow(city_data) != 0){
            latitude = c(latitude,city_data$lat)
            longitude = c(longitude,city_data$lng)
            county_num <- county_num + 1
          }
        }else if( i == length(city_names)){
          city_data <-subset(USCity, city == city_names[i] & state_id == state_names[length(state_names)])
          if(nrow(city_data) != 0){
            latitude = c(latitude,city_data$lat)
            longitude = c(longitude,city_data$lng)
            county_num <- county_num + 1
          }
        }else{
          if(length(state_names) >=3){
            city_data <-subset(USCity, city == city_names[i] & state_id == state_names[2])
            if(nrow(city_data) != 0){
              latitude = c(latitude,city_data$lat)
              longitude = c(longitude,city_data$lng)
              county_num <- county_num + 1
            }              
          }  
        }
      }
      
    }
  }
  
  MSA_info[MSA_Income[z],"lat"] = mean(latitude)
  MSA_info[MSA_Income[z],"lng"] = mean(longitude)
  MSA_info[MSA_Income[z],"num_cities"] = county_num
}

which(MSA_info[,3] == 0)
MSA_info[26,1] = 41.6923
MSA_info[26,2] = -70.0694

MSA_info[42,1] = 43.6005
MSA_info[42,2] = -116.2308


MSA_info[355,1] =21.3294
MSA_info[355,2] = -157.846


haversine <- function(lat1, lon1, lat2, lon2) {
  # convert latitudes and longitudes to radians
  lat1 <- lat1 * pi / 180
  lon1 <- lon1 * pi / 180
  lat2 <- lat2 * pi / 180
  lon2 <- lon2 * pi / 180
  
  # haversine formula
  dlon <- lon2 - lon1
  dlat <- lat2 - lat1
  a <- sin(dlat/2)^2 + cos(lat1) * cos(lat2) * sin(dlon/2)^2
  c <- 2 * asin(sqrt(a))
  r <- 3959 # radius of earth in miles
  distance <- r * c
  
  return(distance)
}

W = matrix(0,nrow = length(MSA_Income), ncol = length(MSA_Income))
rownames(W) <- colnames(W) <- MSA_Income
d = 100

for (i in 1:length(MSA_Income)) {
  for (j in (i+1):length(MSA_Income)) {
    distance <- haversine(MSA_info[i,1],MSA_info[i,2],MSA_info[j,1],MSA_info[j,2])
    if(distance <= d){
      W[i,j] <- W[j,i] <- 1
    }
  }
}

# reweight W
for (i in 1:nrow(W)) {
  sumrow = sum(W[i,])
  if(sumrow !=0){
    W[i,] = W[i,] / sumrow
  }
}

# ymat = ymat[1:50,1:50]
# xmat1 = xmat1[1:50,1:50]
# xmat2 = xmat2[1:50,1:50]
# W = W[1:50,1:50]
# 
# rm(list = setdiff(ls(), c("ymat", "xmat1", "xmat2", "W")))


#####7. State#########
all_state_names <- data.frame()
for (z in 1:length(MSA_Income)) {
  split_string <- strsplit(MSA_Income[z], " ")[[1]]
  split_index <- length(split_string)
  left_side <- paste(split_string[1:(split_index-1)], collapse = " ")
  right_side <- split_string[split_index]
  
  # state names
  state_names <- strsplit(right_side,  "--|-" )[[1]]
  
  for (name_i in state_names) {
    if(name_i %in% rownames(all_state_names)){
      all_state_names[name_i,1] = all_state_names[name_i,1] + 1
    }else{
      all_state_names[name_i,1] = 1
    }
  }
}

east_coast = c("ME", "NH", "MA", "RI", "CT", "NY", "NJ", "PA", "DE", "MD", "VA", "NC", "SC", "GA", "FL")
east_coast = c("WA", "OR", "CA")
MSA_eastcoast = c()
for (z in 1:length(MSA_Income)) {
  split_string <- strsplit(MSA_Income[z], " ")[[1]]
  split_index <- length(split_string)
  left_side <- paste(split_string[1:(split_index-1)], collapse = " ")
  right_side <- split_string[split_index]
  
  # state names
  state_names <- strsplit(right_side,  "--|-" )[[1]]
  
  if(length(intersect(state_names , east_coast)) != 0){
    MSA_eastcoast = c(MSA_eastcoast, MSA_Income[z])
  }
}

ymat = ymat[MSA_eastcoast,]
xmat1 = xmat1[MSA_eastcoast,]
xmat2 = xmat2[MSA_eastcoast,]
W = W[MSA_eastcoast,MSA_eastcoast]

rm(list = setdiff(ls(), c("ymat", "xmat1", "xmat2", "W")))
