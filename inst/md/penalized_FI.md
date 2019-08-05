

url <- "https://data.mendeley.com/datasets/dmxs3gpd5k/5/files/329c8022-0a37-4e6f-96e6-9c25b260b9e9/Dataset.RData?dl=1"
download.file(url,"./hybrid_data.RData",mode="wb")

load("./hybrid_data.RData")