# Here we will revisit the blood pressure data and illustrate regression with principal components

setwd("~/Desktop/2020 Winter/STA437/project")
#########################################
# Principal Components Regression
#########################################
# Recall that
# We had data on blood pressure, age, weight, body surface area, duration of hypertension,  
# basal  pulse, stress index for 20 individuals with high blood pressure
#########################################
BP.data= read.table(file="bp.txt", header=TRUE)[,-1]
BP.data
n=20

colnames(BP.data)
X= as.matrix(BP.data[2:7])      # exclude the response variable!!

# Obtain the sample mean vector and sample covariance matrix of X:

x.bar = apply(X,2,mean)
x.bar

S= cov(X)
S
round(S, 3)
# the scales are quite different, suggesting that we may be better off using R to get the PCs
# nevertheless, let's continue with S. We will redo the analysis with R later!


# Principal Components:
#######################################

# to get the PCs, first find the eigenvalues and eigenvectors

Val = eigen(S)$values
Vec = eigen(S)$vectors

# we have 6 PCs which are linear combinations of X1, ..., X6 with weights given by each eigenvector
# we can obtain the data for PCs using:

W = X  # just to create a data matrix of the same size of X

# now fill in the entries by calculating sample PCs

for(i in 1:6){
 	for(j in 1:20){
	W[j,i] = Vec[,i] %*% ( X[j,] -x.bar)  # centered PCs
 }}

colnames(W) = paste("W", 1:6, sep="")
# Principal Components have zero correlation:

plot(data.frame(W)) 
round( cor(W),3)

# Centered Principal Components have zero mean:

round(apply(W, 2, mean),3)



# How many components should we keep? 
#######################################
# screeplot:

plot(Val, type="b")  # suggests keeping the first PC only!

# Proportion of variation explained by each PC:

round( Val/sum(Val),3)  # 97.3 % of the sample variation in X is explained by the first PC.

# If you like, you can use built-in functions in R for a summary:

summary( prcomp(X,cov.mat=TRUE))
screeplot( prcomp(X,cov.mat=TRUE))


# Now let's run regression with the first PC as the explanatory variable:
##############################################################################

PC.model = lm(BP.data$BP ~ W[,1])
summary(PC.model)  # W1 is not found significant! adj.R^2 is too low (negative!)

# however, if we add W2 (or even more PCs) we would find them significant:

PC.model.2 = lm(BP.data$BP ~  W[,1]+ W[,2])
summary(PC.model.2)  # W1 is not significant!

# check the importance of each variable in the first two PCs:

round(Vec[,1],3)   # Stress level dominates the first PC
round(Vec[,2],3)   # weight and pulse contribute to the second PC more than other variables. 





#######################################
# Principal Components using R 
#######################################

# obtain the standardized variables:

Z=X
for(i in 1:6){
Z[,i] = (X[,i]-x.bar[i])/sqrt(diag(S)[i])
}

# obtain correlation matrix
R = cor(X)
R
cov(Z)  # they should be the same!

# obtain eigenvalues and eigenvectors of R
Val.new = eigen(R)$values
round(Val.new ,2)

Vec.new = eigen(R)$vectors
rownames(Vec.new) = colnames(X)
colnames(Vec.new) = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6")
round(Vec.new ,2)


# Obtain sample PC values:

W.new = X  # just to create a data matrix of the same size of X
colnames(W.new) = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6")

# now fill in the entries by calculating sample PCs

for(i in 1:6){ # PC's
 	for(j in 1:20){
	W.new[j,i] = Vec.new[,i] %*% Z[j,]   # no need to center when using normalized PCCs 
 }}


# Principal Components have zero correlation:

plot(data.frame(W.new)) 
round( cor(W.new),3)




# How many components should we keep? 
#######################################
# screeplot:

plot(Val.new, type="b", pch=19, xlab="",ylab="Variances")  # suggests keeping the first 4 or 5 PCs.

# Proportion of variation explained by each PC:

round( Val.new/sum(Val.new),3)  

# If you like, you can use built-in functions in R for a summary:

summary( prcomp(Z))
screeplot(prcomp(Z),npcs = 6, type = "lines")



# Regression with all standardized PCs as the explanatory variables
##############################################################################

PC.model.new.1 = lm(BP.data$BP ~  W.new[,1] + W.new[,2] + W.new[,3] + W.new[,4] + W.new[,5]+ W.new[,6])
summary(PC.model.new.1)  # W2 and W5 seem not to be significant

# Let's remove W3.new and W5.new

PC.model.new.2 = lm(BP.data$BP ~  W.new[,1] + W.new[,2] + W.new[,4] + W.new[,6])
summary(PC.model.new.2) 

# Note that removing PCs from the model does not change the coefficients!

vif(PC.model.new.2)



# check the importance of each variable in standardized PCs:

round(Vec.new[,1],3)   # all variables contribute to the first PC

round(Vec.new[,2],3)   # stress vs (weight and BSA)  

round(Vec.new[,3],3)   # Dur dominates the third PC  (perhaps vs stress)

round(Vec.new[,4],3)   # Age dominates the fourth PC (perhaps vs stress)

round(Vec.new[,5],3)   # Pulse vs (BSA and stress)  

round(Vec.new[,6],3)   # Weight vs (BSA and pulse)


# the correlations between PCs and variables also indicate importance:

cor.WX = function(mat){
	vals= mat
	rownames(vals) = paste("X", 1:nrow(mat), sep="")
	colnames(vals) = paste("W", 1:nrow(mat), sep="")
	
	for(i in 1: ncol(vals)){
		val= t(eigen(mat)$vectors[,i])* sqrt( eigen(mat)$values[i])
		vals[,i] = val/sqrt(diag(mat))
		}
	return(vals)
	
}


cor.WX(R)







