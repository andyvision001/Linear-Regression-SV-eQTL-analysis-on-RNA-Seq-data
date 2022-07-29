#Loading the required libraries
library(devtools)
library(Biobase)
library(VariantAnnotation)
install.packages("tidymodels")
library("GenomicFeatures")
library(tidyverse)
library(broom)
library(R.utils)
library(ggplot2)
library(readr)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("SNPlocs.Hsapiens.dbSNP.20101109")

#Loading both the SNP and the gene expression

Gene_Expr <- read_csv("C:/Users/HP/Desktop/MSC FILES/Genomics/Andra_Umoru_Project_2/Expression.txt")
#Gene_Expr <- read_csv("C:/Users/HP/Desktop/MSC FILES/Genomics/Andra_Umoru_Project_2/Expression.txt")
SV <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
my_vcf <- readVcf(SV, "/C:/Users/HP/Desktop/MSC FILES/Genomics/Andra_Umoru_Project_2/Genotypes")
COV <- read_csv("C:/Users/HP/Desktop/MSC FILES/Genomics/Andra_Umoru_Project_2/Covariates.txt")
my_output_file = tempfile()

#Getting the SVs information
my_vcf
header(my_vcf)
samples(header(my_vcf))
geno(header(my_vcf))
head(rowRanges(my_vcf), 3)
ref(my_vcf)[1:10]
qual(my_vcf)[1:466]
alt(my_vcf)[1:10]
geno(my_vcf)
sapply(geno(my_vcf), class)

#Plotting histogram showing Genotype Dosage
geno(header(my_vcf))["DS",]
DS <-geno(my_vcf)$DS
dim(DS)
DS[1:10,]
fivenum(DS)
length(which(DS==0))/length(DS)
hist(DS[DS !=0], breaks=seq(0, 2, by=0.05), main="Genotype Dosage (DS) non-zero values", xlab="Genotype Dosage - DS")

#Fetching info data frame
info(my_vcf)[1:20, 1:5]
rd <- rowRanges(my_vcf)
seqlevels(rd) <- "ch22"
ch22snps <- getSNPlocs("ch22")
dbsnpchr22 <- sub("rs", "", names(rd)) %in% ch22snps$RefSNP_id
table(dbsnpchr22)
info(header(my_vcf))[c("VT", "LDAF", "RSQ"),]
metrics <- data.frame(QUAL=qual(my_vcf), inDbSNP=dbsnpchr22, VT=info(my_vcf)$VT, LDAF=info(my_vcf)$LDAF, RSQ=info(my_vcf)$RSQ)


#Histogram for Chromosome Dosage
My_Hist<-(LDAF=info(my_vcf)$LDAF)
dim(My_Hist)
My_Hist[1:1000]
my_data<-My_Hist[1:1000]
hist(my_data[my_data !=0], breaks=seq(0, 2, by=0.017), main="Chromosome/Average Dosage", xlab="Average Dosage", ylab="Chromosome")


#Getting the Gene TSS
header(my_vcf)
my_pos<-(AVGPOST=info(my_vcf)$AVGPOST)
head(my_pos)
my_var<-(LDAF=info(my_vcf)$LDAF)
head(my_var)
lm1 = lm(my_pos ~ my_var)
tidy(lm1)
hist(my_pos, breaks=seq(0, 2, by=0.03), main="Gene TSS > 0.05", xlab="P-Value")


#Visualizing the distribution of quality
ggplot(metrics, aes(x=RSQ, fill=inDbSNP)) +
  geom_density(alpha=0.25) +
  scale_x_continuous(name="Input Quality") +
  scale_y_continuous(name="Density") +
  theme(legend.position="right") +
  ggtitle("GG Plot Showing Input Quality")

SV2<- DS <-geno(my_vcf)$DS
head(Gene_Expr, sep="\t", header=True, row.names=1)
The_Exp<-head(Gene_Expr, sep="\t", header=True, row.names=1)
head(SNP, sep="\t", header=True, row.names=1)
my_cov = read.table(COV, sep="\t", header=T,row.names=1)
My_expression = as.numeric(The_Exp, row.names=1)
My_SNPS = as.numeric(SNP2[1:466])
lm1 = lm(My_expression ~ My_SNPS)
tidy(lm1, conf.int = TRUE, conf.level = 0.95)
hist(My_Hist, breaks=seq(0, 2, by=0.05), main="Genotype Dosage (DS) non-zero values", xlab="Genotype Dosage - DS")
hist(My_Hist, breaks=0.2, xlim=c(10,120), ylim=c(0,5000))
summary(lm1)


#Plot for the Genotype
plot(My_expression ~ jitter(My_SNPS), col=(My_SNPS+1),xaxt="n",xlab="The Genotype",ylab="The Gene Expression")
axis(1,at=c(0:2),labels=c("AA","Aa","aa"))
lines(lm1$fitted ~ My_SNPS,type="b",pch=15,col="grey")


#Performing the Linear Regression and getting all the P-Value
n = 1000;

# Number of variables
no_col <- max(count.fields(SV, "/C:/Users/HP/Desktop/MSC FILES/Genomics/Andra_Umoru_Project_2/Genotypes", sep = "\t"))
no_rows <- sapply(SV, "/C:/Users/HP/Desktop/MSC FILES/Genomics/Andra_Umoru_Project_2/Genotypes",countLines)

# Common signal in all variables (population stratification)
pop = 0.2 * rnorm(n);

# Setting the Gene-SNPs Columns and rows matrices
snps.mat = matrix(rnorm(n*no_col), ncol = no_col) + pop;
gene.mat = matrix(rnorm(n*no_col), ncol = no_col) + pop + snps.mat*((1:no_col)/no_col)^9/2;


#Getting the data from the files (The Expression and SNPs)
snps.mat = matrix(SV, "/C:/Users/HP/Desktop/MSC FILES/Genomics/Andra_Umoru_Project_2/Genotypes", ncol = no_col) + pop;
gene.mat = matrix(Gene_Expr, ncol = no_col) + pop + snps.mat*((1:no_col)/no_col)^9/2;
##


# data objects slice
snps1 = SlicedData$new( t( snps.mat ) );
gene1 = SlicedData$new( t( gene.mat ) );
cvrt1 = SlicedData$new( );
rm(snps.mat, gene.mat)

# Slice data in blocks of 500 variables
snps1$ResliceCombined(500);
gene1$ResliceCombined(500);

#Setting name output file
my_filename = tempfile();

#significant level
pvOutputThreshold = 1e-2;

# Perform analysis recording information for
# a histogram
snpspos = read.table(my_vcf, header = TRUE, stringsAsFactors = FALSE);
genepos = read.table(Gene_Expr, header = TRUE, stringsAsFactors = FALSE);

eQTL = My_Linear_Regression(
  snps = snps1,
  gene = gene1,
  cvrt = cvrt1,
  snpspos = snpspos,
  genepos = genepos,
  output_file_name = my_filename,
  pvOutputThreshold = 1e-2,
  useModel = modelLINEAR,
  errorCovariance = numeric(),
  verbose = TRUE,
  pvalue.hist = 100);
  unlink( my_filename );

#Getting the analysis result
  cat('The Analysis done in: ', eQTL$time.in.sec, ' seconds', '\n');
  cat('Detected eQTLs:', '\n');
  show(eQTL$all$eqtls)

# Plotting the Histogram for all P-Values
hist(eQTL, breaks=seq(0, 2, by=0.05), main="Histogram for all ", eqtl, " p-values", xlab="P-values")
plot(eQTL, col="darkgrey")

eQTL = My_Linear_Regression(
  snps = snps1,
  gene = gene1,
  cvrt = cvrt1,
  snpspos = snpspos,
  genepos = genepos,
  output_file_name = my_filename,
  pvOutputThreshold = 1e-6,
  useModel = modelLINEAR,
  errorCovariance = numeric(),
  verbose = TRUE,
  pvalue.hist = "qqplot");
unlink( my_filename );

# Q-Q Plot for all P-Values
ggplot(eQTL) +
  geom_density(alpha=0.25) +
  scale_x_continuous(name="-log$$10$$(p-value), expected") +
  scale_y_continuous(name="-log$$10$$(p-value), observed") +
  theme(legend.position="top") +
  ggtitle("QQ-plot for all ", eQTL, "p-values")
plot(eQTL, pch = 16, cex = 0.5)