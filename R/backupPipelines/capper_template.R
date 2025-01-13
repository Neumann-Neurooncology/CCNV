library(CCNV)

# path to capper description
path_csv <- "/beegfs/home/s/schumany/AntoniaCCNV/Data/Capper_Methylation/capper_et_al_clinical_2.csv"
output_dir <- "/beegfs/home/s/schumany/AntoniaCCNV/Scripts/Capper_Output"
# description
descr <- read.table(path_csv, header=TRUE, sep=";")

# subgroups

groups <- sort(unique(descr$methylation_class))

# select subset

args <- commandArgs(trailingOnly=TRUE)
my_group <- groups[as.integer(args[1])]

descr <- descr[which(descr$methylation_class==my_group),]

print(paste("Processing group", my_group, "with", nrow(descr), "IDATs"))

# create output directory

my_output_dir <- paste(output_dir, my_group, sep="/")
dir.create(my_output_dir, showWarnings = FALSE)

# apply CCNV

v1_list <- cumul.CNV(descr, conumee.version=1, out="data")
v2_list <- cumul.CNV(descr, conumee.version=2, out="data")

df1_v1 <- v1_list[[1]]
df2_v1 <- v1_list[[2]]

df1_v2 <- v2_list[[1]]
df2_v2 <- v2_list[[2]]

# store

write.table(df1_v1, paste(my_output_dir, "multiseg_v1.csv", sep="/"), sep=",")
write.table(df2_v1, paste(my_output_dir, "singleseg_v1.csv", sep="/"), sep=",")
write.table(df1_v2, paste(my_output_dir, "multiseg_v2.csv", sep="/"), sep=",")
write.table(df2_v2, paste(my_output_dir, "singleseg_v2.csv", sep="/"), sep=",")

