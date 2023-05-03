library(GEOquery)
library(reshape2)
library(limma)
library(vsn)
library(dplyr)
library(readr)


source("scripts/support_functions.R")

Sys.setenv(VROOM_CONNECTION_SIZE=500072) #lol what

#Tried this dataset because the ERBB2 isn't predicted to be active by moon.
#The data actually looks alright, but it's weird that ERBB2 activity isn't up, even though ERBB2 score goes significantly down when an ERBB2 inhibtor is applied. 
#Possibly, the experiement didn't work well (I think moon actually give an acurate transcriptional activity score downstream of ERBB2)

gse <- getGEO("GSE50760", GSEMatrix = T)

df_list <- list()

i <- 1
for(file_name in list.files(pattern = "FPKM", recursive = TRUE))
{
  print(file_name)
  df <- as.data.frame(
    read_delim(file_name, 
               delim = "\t", escape_double = FALSE, 
               trim_ws = TRUE))
  dubs <- df[duplicated(df$symbol),"genes"]
  if(length(dubs) > 0)
  {
    df <- df[-which(df$genes %in% dubs),]
  }
  df <- df[!grepl("[a-z]",df$genes),]
  df_list[[i]] <- df
  i <- i+1
}

df <- merge(df_list[[1]], df_list[[2]], by = "genes")

for(i in 3:length(df_list))
{
  print(i)
  print(dim(df)[1])
  df <- merge(df,df_list[[i]], by = "genes")
}

rm(df_list)

targets <- pData(gse[[1]])
targets <- targets[,c(1,40)]
names(targets) <- c("sample","condition")
targets$sample <- gsub(".* AMC","AMC",targets$sample)

names(df) <- gsub("[.]","-",names(df))
names(df) <- gsub("_FPKM","",names(df))

row.names(df) <- df$genes
eset <- df[,-1]

eset <- eset[,targets$sample]

plot(density(log2(as.numeric(unlist(eset)))))

eset[log2(eset) <= 0] <- NA 
eset <- eset[complete.cases(eset),]

fit <- vsnMatrix(as.matrix(eset))
meanSdPlot(fit)
eset <- as.data.frame(vsn::predict(fit,as.matrix(eset)))

unique(targets$condition)

comparisons <- list("PMC_vs_ADJHEALTHY" = c(1,-2))



limmaRes <- runLimma(eset, targets, comparisons = comparisons)

ttop <- limma_res_to_ttop_list(limma_res = limmaRes,
                       comp_names = names(comparisons),
                       number = length(eset[,1]),
                       adjust.method = "fdr")[[1]]

write_csv(ttop[,c(1,4)],file = "example_inputs/differential_stats.csv")
write_csv(ttop[,c(1,5)],file = "example_inputs/differential_pvals.csv")
write_csv(ttop[ttop[,6] < 0.05,c(1),drop = F],file = "example_inputs/differential_genes.csv", col_names = F)
write_csv(df[,1,drop = F],file = "example_inputs/background_genes.csv", col_names = F)

