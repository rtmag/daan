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
	condition=c("OB1","WE1","YB1","OB2","WE2","YB2"),
	libsize=c(88017941,95590198,105286740,98478413,98608413,104906147)
)

# Define the library size and conditions for the GLM
libsize <- meta$libsize
condition <- factor(meta$condition)
design <- model.matrix(~0+condition)
colnames(design) <- levels(meta$condition)


