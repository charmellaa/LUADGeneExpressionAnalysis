## Just some visualization of the clinical data (information about the patients)

rm(list = ls())

tumor = "LUAD"
path <- "./TCGA_datasets/clinical/"
filename_in <- paste0(path, "/clinical_", tumor, ".txt")

tmp <- read.table(filename_in, header = T, sep = "\t", check.names = F, row.names = 1, quote = "", nrows = 10)
classes <- sapply(tmp, class)

tmp <- read.table(filename_in, header = T, sep = "\t", check.names = F, row.names = 1, quote = "", colClasses = classes)

genes <- rownames(tmp)
col <- colnames(tmp)
dev.off()
barplot(table(tmp$ajcc_pathologic_stage))
barplot(table(tmp$tissue_or_organ_of_origin), col = "brown")

cpd <- tmp$cigarettes_per_day[!is.na(tmp$cigarettes_per_day)]
hist(cpd, main = "", xlab = "cigarettes a day", xlim = c(0,10), ylim = c(0,120), col="orange")
#plot(table(tmp$years_smoked))
#plot(tmp$year_of_death)

age_at_diag <- tmp$age_at_diagnosis[!is.na(tmp$age_at_diagnosis)] 
age_at_diag <- sapply(age_at_diag, function(x) {x/365.24})
hist(age_at_diag, col = "green4", main = "", xlab = "Age at diagnosis")
#hist(tmp$age_at_index, xlim = c(20,100), ylim = c(0, 120))

count<- table((tmp$primary_diagnosis))
pie(count, labels = paste0(names(count), " ", round(100 * count/sum(count), 2), "%"), cex = 0.45)

count<- table((tmp$tissue_or_organ_of_origin))
pie(count, labels = paste0(names(count), " ", round(100 * count/sum(count), 2), "%"), cex = 0.45)

count<- table((tmp$vital_status))
pie(count, labels = paste0(names(count), " ", round(100 * count/sum(count), 2), "%"), col = c("red", "blue4"))

count<- table((tmp$ethnicity))
pie(count, labels = paste0(names(count), " ", round(100 * count/sum(count), 2), "%"))

count<- table((tmp$gender))
pie(count, labels = paste0(names(count), " ", round(100 * count/sum(count), 2), "%"), col = c("cyan", "pink"))

count<- table((tmp$race))
pie(count, labels = paste0(names(count), " ", round(100 * count/sum(count), 2), "%"))

count<- table(tmp$tissue_or_organ_of_origin)
pie(count, labels = paste0(names(count), " ", round(100 * count/sum(count), 2), "%"))



# i want to see: those who underwent treatment or therapy, how many are alive

c <- tmp[tmp$treatments_pharmaceutical_treatment_or_therapy == "yes", "vital_status"]
count <- table(c)
pie(count, labels = paste0(names(count), " ", round(100 * count/sum(count), 2), "%"))

c <- tmp[tmp$treatments_pharmaceutical_treatment_or_therapy == "no", "vital_status"]
count <- table(c)
pie(count, labels = paste0(names(count), " ", round(100 * count/sum(count), 2), "%"))

