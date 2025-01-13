library(CCNV)

# path to capper description
path_csv <- "/beegfs/home/s/schumany/AntoniaCCNV/Data/Capper_Methylation/capper_et_al_clinical_2.csv"
output_dir <- "/beegfs/home/s/schumany/AntoniaCCNV/Scripts/Simulation_Output"
# description
descr <- read.table(path_csv, header=TRUE, sep=";")

CapperTypes <- unique(descr$methylation_class)
CapperCtrl <- CapperTypes[grep('contr', CapperTypes)]
controls <- descr[which(descr$methylation_class %in% CapperCtrl),]

args <- commandArgs(trailingOnly=TRUE)
my_index <- as.integer(args[1])
my_job <- as.integer(args[2])
my_gamma <- as.numeric(args[3])

output_list = cumul.CNV(controls, conumee.version = 2, segmentationMode = "multi", out="data", gamma = my_gamma)
df <- output_list[[1]]
gene_name <- output_list[[2]]
fraction <- output_list[[3]]
strength <- output_list[[4]]
write.table(df, paste(output_dir, "/", "output_",my_job,"_df_", my_index, "_", my_gamma, "_", gene_name, "_", fraction, "_", strength, ".csv", sep=""), sep=",")

