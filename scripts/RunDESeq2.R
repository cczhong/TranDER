args=commandArgs()
#args
#stop()

file=args[6]
group=args[7]
name=args[8]
folder=args[9]

print(file)
print(group)
print(name)

library(DESeq2)

# load in read count and condition table
countData = read.table(file, row.names=1, header=T)
head(countData)
colData = read.table(group, row.names=1, header=T)
head(colData)

dds = DESeqDataSetFromMatrix(countData=countData, colData=colData, design=~Condition)
dds = DESeq(dds)

# wirte differential expressed genes
diffexpr = results(dds, contrast=c("Condition", "treatment", "control"))
diffexpr = diffexpr[complete.cases(diffexpr),]
de_file_name = paste("diffexpr", name, "tab", sep=".")
write.table(diffexpr, paste(folder, de_file_name, sep="/"))

# print out MA-plot
ma_file_name = paste("MAplot", name, "pdf", sep=".")
pdf(paste(folder, ma_file_name, sep="/"))
plotMA(diffexpr, main=name, ylim=c(-10,10));
dev.off()

