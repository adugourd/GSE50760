library(GEOquery)
library(reshape2)
library(limma)
library(vsn)
library(dplyr)
library(GSEABase)
library(org.Hs.eg.db)
library(readr)
library(decoupleR)
library(ggplot2)

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

comparisons <- list("PMC_vs_ADJHEALTHY" = c(1,-2),
                    "METASTATIC_vs_ADJHEALTHY" = c(3,-2),
                    "METASTATIC_vs_PMC" = c(3,-1))



limmaRes <- runLimma(eset, targets, comparisons = comparisons)


t_table <- ttop_list_to_t_table(
  limma_res_to_ttop_list(limma_res = limmaRes,
                         comp_names = names(comparisons),
                         number = length(eset[,1]),
                         adjust.method = "fdr"))

row.names(t_table) <- t_table$ID
t_table <- t_table[,-1]

cp_pathways <- import_gmt(gmtfile = 'support/c2.cp.v7.2.symbols.gmt')

names(cp_pathways) <- c("target","tf")
cp_pathways$mor <- 1
cp_pathways$likelihood <- 1

cp_pathways <- cp_pathways[grepl("NABA_",cp_pathways$tf) |  grepl("KEGG_",cp_pathways$tf),]

cp_pathways <- intersect_regulons(mat = as.matrix(t_table), network = cp_pathways, minsize = 10, .source = tf, .target = target)
cp_pathway_mean <- run_wmean(mat = as.matrix(t_table), network = cp_pathways, .source = tf, times = 1000)
cp_pathway_mean <- cp_pathway_mean[cp_pathway_mean$statistic == "norm_wmean",]
cp_pathway_mean_df <- dcast(cp_pathway_mean, formula = source~condition, value.var = "score")
row.names(cp_pathway_mean_df) <- cp_pathway_mean_df$source
cp_pathway_mean_df <- cp_pathway_mean_df[,-1]

cp_pathway_mean_df <- cp_pathway_mean_df[apply(abs(cp_pathway_mean_df), 1, max) > 3.5,]
cp_pathway_mean_df <- cp_pathway_mean_df[order(cp_pathway_mean_df$METASTATIC_vs_ADJHEALTHY,decreasing = F),]

to_plot <- melt(cbind(cp_pathway_mean_df, row.names(cp_pathway_mean_df)))
names(to_plot)[1] <- "pathway"
to_plot$pathway <- gsub("KEGG_","",to_plot$pathway)
to_plot$pathway <- gsub("NABA_","",to_plot$pathway)
to_plot$pathway <- gsub("_"," ",to_plot$pathway)

to_plot$pathway <- factor(to_plot$pathway, levels = unique(to_plot$pathway))

ggplot(to_plot, aes(x = value, y = pathway, color = variable, size = abs(value))) + 
  geom_point(alpha = 0.5) + 
  geom_vline(xintercept = 0) + 
  theme_minimal()

pathways_df <- as.data.frame(cp_pathways)
BCCA <- pathways_df[pathways_df$tf == "KEGG_VALINE_LEUCINE_AND_ISOLEUCINE_DEGRADATION","target"]

ttop <- ttopFormatter(topTable(limmaRes[[1]], coef = 3, number = 9180, adjust.method = "fdr")) #metastasis vs PMC
ttop_BCCA <- ttop[ttop$ID %in% BCCA,]

volcano_nice(ttop_BCCA, FCIndex = 2, pValIndex = 5, IDIndex = 1, )
