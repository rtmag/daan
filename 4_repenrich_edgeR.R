library('edgeR')

# In the case of seperate outputs, load the RepEnrich results - fraction counts
ob1 <- read.delim('OB1_fraction_counts.txt', header=FALSE)
we1 <- read.delim('WE1_fraction_counts.txt', header=FALSE)
yb1 <- read.delim('YB1_fraction_counts.txt', header=FALSE)
ob2 <- read.delim('OB2_fraction_counts.txt', header=FALSE)
we2 <- read.delim('WE2_fraction_counts.txt', header=FALSE)
yb2 <- read.delim('YB2_fraction_counts.txt', header=FALSE)
##################
#' Build a counts table
counts <- data.frame(
  row.names = ob1[,1],  
  ob1 = ob1[,4], we1 = we1[,4], 
  yb1 = yb1[,4], ob2 = ob2[,4], 
  we2 = we2[,4], yb2 = yb2[,4]
)

# Build a meta data object. I am comparing young, old, and veryold mice.
# I manually input the total mapping reads for each sample.
# The total mapping reads are calculated using the bowtie logs:
# # of reads processed - # reads that failed to align
meta <- data.frame(
	row.names=colnames(counts),
	condition=c("OB","WE","YB","OB","WE","YB"),
	libsize=c(88017941,95590198,105286740,98478413,98608413,104906147)
)

# Define the library size and conditions for the GLM
libsize <- meta$libsize
condition <- factor(meta$condition)
design <- model.matrix(~0+condition)
colnames(design) <- levels(meta$condition)

# Build a DGE object for the GLM
y <- DGEList(counts=counts, lib.size=libsize)

# Normalize the data
y <- calcNormFactors(y)
y$samples
plotMDS(y)

# Estimate the variance
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTrendedDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
plotBCV(y)

# Build an object to contain the normalized read abundance
logcpm <- cpm(y, log=TRUE, lib.size=libsize)
logcpm <- as.data.frame(logcpm)
colnames(logcpm) <- factor(meta$condition)

# Conduct fitting of the GLM
yfit <- glmFit(y, design)

# Initialize result matrices to contain the results of the GLM
results <- matrix(nrow=dim(counts)[1],ncol=0)
logfc <- matrix(nrow=dim(counts)[1],ncol=0)

# Make the comparisons for the GLM
# Young VS Old
# Young VS WEHI
my.contrasts <- makeContrasts( y_o = YB - OB,
			       y_w = YB - WE,
			      levels = design )  

# Define the contrasts used in the comparisons
allcontrasts = c("y_o","y_w")


# Conduct a for loop that will do the fitting of the GLM for each comparison
# Put the results into the results objects
for(current_contrast in allcontrasts) {
	lrt <- glmLRT(yfit, contrast=my.contrasts[,current_contrast])
	plotSmear(lrt, de.tags=rownames(y))
	title(current_contrast)
	res <- topTags(lrt,n=dim(c)[1],sort.by="none")$table
	colnames(res) <- paste(colnames(res),current_contrast,sep=".")
	results <- cbind(results,res[,c(1,5)])
	logfc <- cbind(logfc,res[c(1)])
}

# Add the repeat types back into the results.
# We should still have the same order as the input data
results$class <- ob1[,2]
results$type <- ob1[,3]

# Sort the results table by the logFC
results <- results[with(results, order(-abs(logFC.y_o))), ]
results <- cbind(results,rownames(results))
colnames(results)[7]="repeats"
# Save the results
#write.table(results, 'results.txt', quote=FALSE, sep="\t")
results_ori=results

#FDR 5%
results=results_ori

for(current_contrast in allcontrasts) {
  logFC <- results[ results[, paste0("FDR.", current_contrast)]<0.05, 
		   paste0("logFC.", current_contrast)]
  # Plot the repeat classes
	
pdf(paste(current_contrast,"_class_replicates_repenrich_FDR_5_noFilter.pdf",sep=""))
  classes <- with(results[results[, paste0("FDR.", current_contrast)]<0.05,], reorder(class, -logFC, median))
  par(mar=c(6,10,4,1))
  boxplot(logFC ~ as.vector(classes), data=results, outline=FALSE, horizontal=TRUE,
          las=2, xlab="log2(Fold Change)", main=paste("Class",current_contrast,"FDR 5%") )
  abline(v=0)
dev.off()
	
pdf(paste(current_contrast,"_type_replicates_repenrich_FDR_5_noFilter.pdf",sep=""))
  # Plot the repeat types
	
	  par(mar=c(6,10,4,1))
  types <- with(results[results[, paste0("FDR.", current_contrast)]<0.05,], reorder(type, -logFC, median))
    boxplot(logFC ~ as.vector(types), data=results, outline=FALSE, horizontal=TRUE,
          las=2, xlab="log2(Fold Change)", main=paste("Type",current_contrast,"FDR 5%") )
  abline(v=0)
dev.off()

pdf(paste(current_contrast,"_repeats_replicates_repenrich_FDR_5_noFilter.pdf",sep=""))
	#plot repeats
		  par(mar=c(6,10,4,1),cex.axis=.4)
  repe <- with(results[results[, paste0("FDR.", current_contrast)]<0.05,], reorder(repeats, -logFC, median))
boxplot(logFC ~ as.vector(repe), data=results, outline=FALSE, horizontal=TRUE,
          las=2, xlab="log2(Fold Change)", main=paste("Repeat",current_contrast,"FDR 5%"))
  abline(v=0)
		  par(mar=c(6,10,4,1),cex.axis=1)
	dev.off()

}
