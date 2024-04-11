
setwd("C:/Users/lvisinti/OneDrive - UGent/Projects/TeA OPLS-DA/Orbitrap - Veterinary medicine/OPLSDA/Github")

suppressPackageStartupMessages({ 

	# Bioconductor packages
	library(structToolbox)
	library(ropls)
  
 	 # CRAN libraries
  	library(ggplot2)
  	library(gridExtra)
  	library(cowplot)
  	library(openxlsx)
  	library(caret)
    library(lattice)
  })


# import the data
data = vroom::vroom("C:/Users/lvisinti/OneDrive - UGent/Projects/TeA OPLS-DA/Orbitrap - Veterinary medicine/OPLSDA/Github/data_neg.csv")

#edit the dataset and create a dataframe
variables = (data[,1])
samples = (data[1,-1])
Class = data[30509,]
data = as.data.frame(t(data))
colnames(data) = data[1,]
data <- data.frame(sapply(data[-1,-30509], function(x) as.numeric(as.character(x))))
row.names(data) = colnames(samples)
data[,30509] = t(Class[,-1])
names(data)[names(data) == "V30509"] <- "Class"

data$Class = as.factor(data$Class)
data.class(data)
Class = data$Class
data = data[,-c(30509)] # no class

#delete outliers 
data = data[-c(36,38,55,72),]

#remove empty columns, infinite values, NAs, and centre the data
X <- data[,colSums(data != 0) != 0] 
is.na(X)<-sapply(X, is.infinite)
X[is.na(X)]<-0
X = as.data.frame(scale(X,center = TRUE, scale = FALSE))

# define dataset and metadata
dataMatrix = data.frame(X)
dataMatrix <- dataMatrix[,colSums(dataMatrix != 0) != 0] 
sampleMetadata = data.frame(Class)
variable_meta = data.frame(colnames(dataMatrix))

# prepare the OPLS-DA model
oplsda <- opls(dataMatrix, Class, permI = 1000,
                        predI = 1, orthoI = NA)

# assess the predictive performance
oplsda <- opls(dataMatrix, Class,
                        predI = 1, orthoI = NA,
                        subset = "odd")

# prediction on the training subset
trainVi <- getSubsetVi(oplsda)
confusion_train.tb <- table(Class[trainVi], fitted(oplsda))
confusion_train.tb


#create a confusion matrix for the training set
conf_matrix <- as.matrix(confusion_train.tb)

image(1:ncol(conf_matrix), 1:nrow(conf_matrix), conf_matrix, col = c("white", "gray"),
      xlab = "", ylab = "", axes=FALSE)
for(i in 1:nrow(conf_matrix)) {
  for(j in 1:ncol(conf_matrix)) {
    text(j, i, conf_matrix[i, j], cex = 0.7, col = "black")
  }
}
axis(1, at = 1:ncol(conf_matrix), labels = c("Control", "TeA"))
axis(2, at = 1:nrow(conf_matrix), labels = c("Control", "TeA"), las = 2)

# prediction on the test subset
confusion_test.tb <- table(Class[-trainVi],
                           predict(oplsda, dataMatrix[-trainVi, ]))
confusion_test.tb

#create a confusion matrix for the test set
conf_matrix <- as.matrix(confusion_test.tb)

image(1:ncol(conf_matrix), 1:nrow(conf_matrix), conf_matrix, col = c("white", "gray"),
      xlab = "", ylab = "", axes=FALSE)
for(i in 1:nrow(conf_matrix)) {
  for(j in 1:ncol(conf_matrix)) {
    text(j, i, conf_matrix[i, j], cex = 0.7, col = "black")
  }
}
axis(1, at = 1:ncol(conf_matrix), labels = c("Control", "TeA"))
axis(2, at = 1:nrow(conf_matrix), labels = c("Control", "TeA"), las = 2)

# calculate the p-value of each feature using a t-test
ttest <- function(df, grp1, grp2) {
  x = df[grp1]
  y = df[grp2]
  x = as.numeric(x)
  y = as.numeric(y)  
  results = t.test(x, y)
  results$p.value
}
rawpvalue = apply(dataMatrix, 2, ttest, grp1 = c(1:35), grp2 = c(36:68))
hist(rawpvalue)

# False Discovery Rate correction
pvalue = p.adjust(rawpvalue, method = "fdr", n = length(rawpvalue))
hist(pvalue)

#transform our data into log2 base and replace NAs
logdataMatrix = log2(dataMatrix)
logdataMatrix[is.na(logdataMatrix)]<-0

#calculate the mean of each feature for control group
control = apply(logdataMatrix[1:35,], 2, mean)

#calculate the mean of each gene for TeA group
test = apply(logdataMatrix[36:68,], 2, mean) 

#Calculate the log2 Fold Change or log2 Ratio == log2(control / test)
foldchange <- test - control
class(foldchange)
hist(foldchange, xlab ="log2 Fold Change (Control vs TeA)")

#volcano plot
results = cbind(foldchange, pvalue)
results = as.data.frame(results)
results$probename <- rownames(results)

# color based on excretion
results$diffexpressed[results$pvalue > 0.05] <- "NS"
results$diffexpressed[abs(results$foldchange) < 1] <- "NS"
results$diffexpressed[results$foldchange > 1 & results$pvalue < 0.05] <- "UP"
results$diffexpressed[results$foldchange < -1 & results$pvalue < 0.05] <- "DOWN"

p <- ggplot(data=results, aes(x = foldchange, y = -log10(pvalue), col=diffexpressed)) + geom_point(size = 0.1) + theme_minimal()
p2 <- p + geom_vline(xintercept=c(-1, 1), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")
p2

#create reduced dataset with the 100 features most influenced by the expsure to TeA
logpvalue = -log10(pvalue)
dataMatrix_reduced = cbind(t(dataMatrix), pvalue, foldchange)
dataMatrix_reduced = as.data.frame(dataMatrix_reduced)

# label the significant variables
dataMatrix_reduced$drop[dataMatrix_reduced$foldchange < -1 & dataMatrix_reduced$pvalue < 0.05] <- 1
dataMatrix_reduced$drop[dataMatrix_reduced$foldchange > 1 & dataMatrix_reduced$pvalue < 0.05] <- 1

#normalization of pvalue and FC
process <- preProcess(as.data.frame(logpvalue), method=c("range"))
pvalue_norm <- predict(process, as.data.frame(logpvalue))
absFC <- abs(foldchange)
process <- preProcess(as.data.frame(absFC), method=c("range"))
foldchange_norm <- predict(process, as.data.frame(absFC))
score <- pvalue_norm + foldchange_norm
dataMatrix_reduced = cbind(dataMatrix_reduced, score)

# delete not labelled variables
reduced = dataMatrix_reduced[which(dataMatrix_reduced$drop >= 1), ]

Top100features = reduced[order(reduced$logpvalue, decreasing = TRUE),]
Top100features = Top100features[1:100,]

#save the dataset as .xlsx file
openxlsx::write.xlsx(Top100features, "C:/Users/lvisinti/OneDrive - UGent/Projects/TeA OPLS-DA/Orbitrap - Veterinary medicine/OPLSDA//Github/Top100_neg.xlsx", rowNames=TRUE)

reduced = as.data.frame(t(Top100features[, 1:68]))
p3 <-  ggplot(data=Top100features, aes(x = foldchange, y = -log10(pvalue))) + geom_point(size = 0.1) + theme_minimal()
p3

#prepare the OPLS-DA model using the reduced dataset
oplsda_red <- opls(reduced, Class, permI = 2000,
                   predI = 1, orthoI = NA)

# assess the predictive performance
oplsda_red <- opls(reduced, Class,
                   predI = 1, orthoI = NA,
                   subset = "odd")

# prediction on the training subset
trainVi <- getSubsetVi(oplsda_red)
confusion_train.tb <- table(Class[trainVi], fitted(oplsda_red))
confusion_train.tb

#create a confusion matrix for the test set
conf_matrix <- as.matrix(confusion_train.tb)

image(1:ncol(conf_matrix), 1:nrow(conf_matrix), conf_matrix, col = c("white", "gray"),
      xlab = "", ylab = "", axes=FALSE)
for(i in 1:nrow(conf_matrix)) {
  for(j in 1:ncol(conf_matrix)) {
    text(j, i, conf_matrix[i, j], cex = 0.7, col = "black")
  }
}
axis(1, at = 1:ncol(conf_matrix), labels = c("Control", "TeA"))
axis(2, at = 1:nrow(conf_matrix), labels = c("Control", "TeA"), las = 2)

# prediction on the test subset
confusion_test.tb <- table(Class[-trainVi],
                           predict(oplsda_red, reduced[-trainVi, ]))
confusion_test.tb

conf_matrix <- as.matrix(confusion_test.tb)

image(1:ncol(conf_matrix), 1:nrow(conf_matrix), conf_matrix, col = c("white", "gray"),
      xlab = "", ylab = "", axes=FALSE)
for(i in 1:nrow(conf_matrix)) {
  for(j in 1:ncol(conf_matrix)) {
    text(j, i, conf_matrix[i, j], cex = 0.7, col = "black")
  }
}
axis(1, at = 1:ncol(conf_matrix), labels = c("Control", "TeA"))
axis(2, at = 1:nrow(conf_matrix), labels = c("Control", "TeA"), las = 2)
