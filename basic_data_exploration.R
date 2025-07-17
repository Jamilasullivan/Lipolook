## Data given in .xlsx format ##################################################

setwd("C:/Users/jamsu/OneDrive - Cardiff University/University/Masters/Big Data Biology/Modules/Dissertation/R Programme/Dissertation/test_data")

## data provided to the client

raw_data <- read.csv("raw_data_1.csv", header = T)
names(raw_data)[c(1,2,3,4,5)] <- c("Filename", "CMaLL Sample Name", "Group", "code", "Cell Count")
raw_data <- raw_data[-1,]
raw_data$code <- NULL
raw_data$`Cell Count` <- NULL
View(raw_data)
summary(raw_data)
str(raw_data) # they were all characters here
# Convert columns 4 to end (since the first 3 are metadata)
raw_data[, 4:ncol(raw_data)] <- lapply(raw_data[, 4:ncol(raw_data)], function(x) as.numeric(as.character(x)))
sum(is.na(raw_data[, 4:ncol(raw_data)])) # no NA values introduced by coercion
str(raw_data)
summary(raw_data)

sapply(raw_data[, 4:ncol(raw_data)], class)
colSums(is.na(raw_data[, 4:ncol(raw_data)]))
colSums(is.infinite(as.matrix(raw_data[, 4:ncol(raw_data)])))

summary(raw_data[, 4:10])

## lipids tested for the client

lipids_tested <- read.csv("lipids_tested_1.csv", header = T)
View(lipids_tested)


## subset the data by groups ###################################################

group_names <- unique(raw_data$Group)
View(group_names)
