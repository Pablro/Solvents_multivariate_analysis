#Solvents projects
#install packages
install.packages("dplyr")
install.packages("car")
install.packages("ggplot2")
install.packages("gridExtra")
install.packages("psych")
install.packages("factoextra")
install.packages("dendextend")

#libraries
library("dplyr")
library("car")
library(ggplot2)
library("gridExtra")
library("MASS")
library(psych)
library("factoextra")
library("dendextend")
#This is the R script project
#Locating and loading my input file with respect to my directory.
getwd()
setwd("C:/Users/toros/Desktop/multivariate-project")
#Reading input file
solvents<-read.csv("solvents.csv",header = T)
head(solvents)
dim(solvents)
names(solvents)
View(solvents)

#pairs plot- Inter relation exploration-Exploratory
pairs(solvents[-1])
# Clearly logP and solubility show a strong linear relation
# Dielectric point seems to show poor relation with other variables
#Melting point and boiling point also seems shows a noisy linear relation

#Checking for normality assumptions by exploratory methods
attach(solvents)
names(solvents)
hist(MeltingPoint)
summary(MeltingPoint)
boxplot(MeltingPoint)
hist(ET30)
summary(ET30)
boxplot(ET30)
hist(Solubility)
summary(Solubility)
boxplot(Solubility)
hist(DipoleMoment)
summary(DipoleMoment)
boxplot(DipoleMoment)
hist(Density)
summary(Density)
boxplot(Density)
hist(BoilingPoint)
summary(BoilingPoint)
boxplot(BoilingPoint)
hist(RefractiveIndex)
summary(RefractiveIndex)
boxplot(RefractiveIndex)
hist(logP)
summary(logP)
boxplot(logP)
hist(Dielectric)#This is definetely not normal
summary(Dielectric)
boxplot(Dielectric)
# #Detecting outliers according to Johson & Wichern (step 1)
par(mfcol=c(2,5),cex=0.4,mar=c(4,4,1,1))
dotchart(MeltingPoint,ylab = "solvents",xlab = "Melting Point",main = "Dotchart")
dotchart(ET30,ylab = "solvents",xlab = "ET30",main = "Dotchart")
dotchart(Solubility,ylab = "solvents",xlab = "Solubility",main = "Dotchart")
dotchart(DipoleMoment,ylab = "solvents",xlab = "Dipole Moment",main = "Dotchart")
dotchart(Density,ylab = "solvents",xlab = "Density",main = "Dotchart")
dotchart(BoilingPoint,ylab = "solvents",xlab = "Boiling Point",main = "Dotchart")
dotchart(RefractiveIndex,ylab = "solvents",xlab = "Refractive index",main = "Dotchart")
dotchart(logP,ylab = "solvents",xlab = "logP",main = "Dotchart")
dotchart(Dielectric,ylab = "solvents",xlab = "Dielectric",main = "Dotchart")
#I find no clear evidence of outliers
par(mfcol=c(2,5))
qqPlot(MeltingPoint)
qqPlot(ET30)
qqPlot(Solubility)
qqPlot(DipoleMoment)
qqPlot(Density)
qqPlot(BoilingPoint)
qqPlot(RefractiveIndex)
qqPlot(logP)
qqPlot(Dielectric)
#Maybe some outliers. Most variables can be consider normal except for dielectric
#no possible transformation
#Detecting outliers according to Johson & Wichern (step 3)
#Standarized values
stdsolvents <- solvents %>% mutate_at(c('MeltingPoint', 'ET30','Solubility','DipoleMoment'
                                        ,'Density','BoilingPoint','RefractiveIndex','logP',
                                        'Dielectric'), ~(scale(.) %>% as.vector))
View(stdsolvents)
#My impression is that there is not clear indication of outliers
#like values range from -1 to 3. Maybe most frequent values are 
#between 0 and 1 but it is not clear indication that a value of 3 is certainly an outlier
#However because of the univariate distribution saw, I will exclude Dielectric variable
shapiro.test(MeltingPoint)
shapiro.test(ET30)
shapiro.test(Solubility)
shapiro.test(DipoleMoment)
shapiro.test(Density)
shapiro.test(BoilingPoint)
shapiro.test(RefractiveIndex)
shapiro.test(logP) # strictly speaking only normal variable 
shapiro.test(Dielectric)

#I will remove the dielectric variable as it has a clear not normality nor linear 
#relationship to ther variables.
#I will avoid presuming about outliers because there is not strong evidence and I 
#do not want to remove chemical components from the analysis
#Removal of dielectric variable
solvents2<-solvents[-4]
#centered solvents
solventscenter <- solvents2 %>% mutate_at(c('MeltingPoint', 'ET30','Solubility',
                                            'DipoleMoment','Density','BoilingPoint',
                                            'RefractiveIndex','logP'),
                                          ~(scale(.,scale = F) %>% as.vector))
attach(solventscenter)
head(solventscenter)
#Detecting outliers according to Johson & Wichern (step 2 and 4)
#pairs plot
pairs(solventscenter[-1])
#Mahalanobis distance for detection of outliers
#covariance matrix of center data
names(solventscenter)
s2=cov(solventscenter[-1])
s2
#Data already center
solventscenter$dist=mahalanobis(solventscenter[-1],center = F, cov = s2)
par(mfcol=c(1,1))
plot(dchisq(solventscenter$dist,df=7),ylab = "mahalanobis distance",xlab = "solvents",
     main = " outliers Dotchart")
View(solventscenter)
#PCA
pairs(solventscenter[-1])
#PCA analysis
#I will work with correlation matrix as I would expected important variance 
#due to scaling in different properties measurments
pca_results=princomp(solvents2[-1],cor = T)
pca_results
summary(pca_results)
pca_results
#Screeplot
#calculate total variance explained by each principal component
names(pca_results)
var_explained = pca_results$sdev^2 / sum(pca_results$sdev^2)

#create scree plot
components=c(1:8)
screeplotdata=data.frame(components,var_explained,row.names = NULL)
ggplot(data = screeplotdata,mapping = aes(components,var_explained)) + 
  geom_line() + geom_point()+
  xlab("Principal Component") + 
  ylab("Variance Explained") +
  ggtitle("Scree Plot") +
  ylim(0, 1)
#The elbow will be in 3.However most of the variance is capture by the first 3 Components. 
#However the total variance capture is still quite low around 80%
#eigen values
pca_results$sdev^2
pca_results$loadings
#This is the cor(z,x) which will yield the observation about the share variance in PC1 and PC2
importance_var_comp=cor(pca_results$scores,solvents2[-1])
#write.csv(importance_var_comp,file="score_x_cor.csv")#for the pdf report
View(importance_var_comp)
#Basically we see that from the 70% of the variance capture
#the PC1 has a strong linear correlation on Dipole moment, Et30,logP and solubility
#while the PC2 has a strong linear correlation on the remaining variables
#as the score are linear combinations. it would be logical that variables with strong linear correlation
#would also have higher weights.
#The screeplot shows that PC1 and PC2 yield arround 70% of the variance of the data
colnames(pca_results$scores)
pca_scores=as.data.frame(pca_results$scores)
pca_scores
#5PC visualization (captures around 90% variance)
p1<-ggplot(data = pca_scores,mapping = aes(Comp.1,Comp.2)) + 
  geom_point()+
  xlab("Pc1") + 
  ylab("PC2") +
  ggtitle("PCA: PC1 vs PC2")+
  theme(text=element_text(size=8)) #this is the plot of the report (PC1 vs PC2)
#p1
p2<-ggplot(data = pca_scores,mapping = aes(Comp.1,Comp.3)) + 
  geom_point()+
  xlab("Pc1") + 
  ylab("PC3") +
  ggtitle("PCA: PC1 vs PC3")+
  theme(text=element_text(size=8))
p3<-ggplot(data = pca_scores,mapping = aes(Comp.1,Comp.4)) + 
  geom_point()+
  xlab("Pc1") + 
  ylab("PC4") +
  ggtitle("PCA: PC1 vs PC4")+
  theme(text=element_text(size=8))
p4<-ggplot(data = pca_scores,mapping = aes(Comp.1,Comp.5)) + 
  geom_point()+
  xlab("Pc1") + 
  ylab("PC5") +
  ggtitle("PCA: PC1 vs PC5")+
  theme(text=element_text(size=8))
p5<-ggplot(data = pca_scores,mapping = aes(Comp.2,Comp.3)) + 
  geom_point()+
  xlab("Pc2") + 
  ylab("PC3") +
  ggtitle("PCA: PC2 vs PC3")+
  theme(text=element_text(size=8))
p6<-ggplot(data = pca_scores,mapping = aes(Comp.2,Comp.4)) + 
  geom_point()+
  xlab("Pc2") + 
  ylab("PC4") +
  ggtitle("PCA: PC2 vs PC4")+
  theme(text=element_text(size=8))
p7<-ggplot(data = pca_scores,mapping = aes(Comp.2,Comp.5)) + 
  geom_point()+
  xlab("Pc2") + 
  ylab("PC5") +
  ggtitle("PCA: PC2 vs PC5")+
  theme(text=element_text(size=8))
p8<-ggplot(data = pca_scores,mapping = aes(Comp.3,Comp.4)) + 
  geom_point()+
  xlab("Pc3") + 
  ylab("PC4") +
  ggtitle("PCA: PC3 vs PC4")+
  theme(text=element_text(size=8))
p9<-ggplot(data = pca_scores,mapping = aes(Comp.3,Comp.5)) + 
  geom_point()+
  xlab("Pc3") + 
  ylab("PC5") +
  ggtitle("PCA: PC3 vs PC5")+
  theme(text=element_text(size=8))
p10<-ggplot(data = pca_scores,mapping = aes(Comp.4,Comp.5)) + 
  geom_point()+
  xlab("Pc4") + 
  ylab("PC5") +
  ggtitle("PCA: PC4 vs PC5")+
  theme(text=element_text(size=8))
grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10, ncol = 2) 
#There is no clear evidence of grouping structure on PCA.
#However with all these pcalots 90% of the variance is explained

#Factor analysis
#It seems pointless to employ principal component extraction with the predefined components 
# As the number of components to explain at least 90% of the variance 
#will be 5 components which will be hard to interpret given the physical and chemical 
#properties. However give what we observed between the correlation of the scores and the variables
#(importance_var_comp)
# I will choose to do a factor analysis on two factor variables as 70% of the variance is capture in
#two components (summary(pca_results))
# and the respective variables are strongly located in each of these components suggesting strong 
#shared variance.
#I will first attempt factor extraction with PCA method  with two PC and further look for 
#improvement with another
#factor extraction technique.
#PCA 
pca_results$loadings
solvents_cor<-cor(solvents2[-1])
solvents_cor_eigen<-eigen(solvents_cor)
#Same eigen values
solvents_cor_eigen
pca_results$sdev^2
pca_comp_fa<-solvents_cor_eigen$vectors[,1:2]
pca_results$loadings[,1:2]
#Similar results with different signs for the PC both should be valid
pca_diag_sd<-diag(sqrt(solvents_cor_eigen$values[1:2]))
pca_diag_sd^2 #verified
#Scaled principal components for factor loading estimation
factor_loadings<- pca_comp_fa%*%pca_diag_sd
rownames(factor_loadings)<-names(solvents2[-1])
colnames(factor_loadings)<-c("F1","F2")
factor_loadings
#As expected, there is an interesting strong underlying weight (factor loadings) 
#for this two factor variable suggesting this common variance among the association of
# observed variables in this two groups which coincides with the correlation of the 
#variables with the scores
#seen in PCA analysis and the great proportion of total variance explained in this this
#two first principal components

#According to what observed and investigated (simple google search IA output). 
#The factor variables names could be first attempt to be describe as follows:
#F2 which strongly holds Melting Point, Boiling Point, Density and refractive index can be further
# name as phase transition and strong density dependence measures.
#F1 which strongly holds Dipole Moment, ET30, logP and Solubility which then can be name:
#polarity and solubility measures.
colnames(factor_loadings)<-c("polarity/solubility","pha_tra/den_dep")
factor_loadings
#However this naming convention base on this grouping structure would be better assess and analyzed
#by chemist (intepretation improvements of factor variables)

#Estimated correlate matrix
est_cor_mat_PCA<-factor_loadings%*%t(factor_loadings)
est_cor_mat_PCA
#communalities:
diag(est_cor_mat_PCA)
#residual calculation between estimated and observed cor- matrix
mat_res<-solvents_cor-est_cor_mat_PCA
mat_res
#Specific variance:
diag(mat_res)
#In this sense density will represent the variable which is less explain by the common factor.
#However for the rest fo variable the is a an important of total variance share among variables
#and capture by
#the common factors.
#RMSE calculation
RMSE_res <- (mat_res - diag(mat_res))^2 #Make the head diagonal equal to 0, and square the 
#off-diagonal residuals       
solvents_faPCA_RMSE_overall <- sqrt(sum(RMSE_res)/length(RMSE_res))  #Take the square root of 
#the average residual
solvents_faPCA_RMSE_overall
#Seems good but there is doubt if it can be improve by other technique given the limitation 
#of the given RMSE value and the 70% of total variance explain with just 
# 2 PC

#Alternative Factor extraction method:
#Given the strong loading factor observed with PCA I will not rotate the loadings factors
#default factor score estimation ("Thurdstone"<- W = R^{-1} F)-regression based weights
#SMC use squared multiple correlations as initial communality
#Factor extraction: weighted least squares method
WLSmethod<-fa(r=solvents_cor,nfactors = 2,rotate = "none",SMC=TRUE,fm="wls")
#factor loadings and variance explanation
WLSmethod$Structure
#Communalities
WLSmethod$communality
#specific variation
1-WLSmethod$communality
#Comparable to RMSE
WLSmethod$rms
#clearly improvement
#Third approach:
#Given the strong loading factor observed with PCA I will not rotate the loadings factors
#default factor score estimation ("Thurdstone"<- W = R^{-1} F)-regression based weights
#SMC use squared multiple correlations as initial communality
#Factor extraction method: Principal factor method
PFmethod<-fa(r=solvents_cor,nfactors = 2,rotate = "none",SMC=TRUE,fm="pa")
#factor loadings and variance explanation
PFmethod$Structure
#Communalities
PFmethod$communality
#specific variation
1-PFmethod$communality
#Comparable to RMSE
PFmethod$rms
#No further improvement.
#It seems that using two factors will never explain more than 70% of total variances.
#Furthermore regarding the loading factors all methods yield similar results with similar 
#communalities and specific variance.
#yet iterative methods yields better root mean square of the residual correlation matrix

###############
#Biplot
#Vector of factor
fA_var=c("pha_tra/den_dep","pha_tra/den_dep","polarity/solubility","pha_tra/den_dep",
         "polarity/solubility","pha_tra/den_dep","polarity/solubility","polarity/solubility")
fviz_pca_biplot(pca_results,addEllipses = T,col.var =fA_var,
                title = "PCA Biplot and intituive relation to FA inference")
#Hierarchical clustering- looking for structure in solvents
#standarized observations-standarized euclidean distance to be used
solvents3<-scale(solvents2[-1],center=T,scale=T)
solvents3
dim(solvents3)
#Euclidean distance as default
#However because of previous step I will be performing 
#standarized euclidean distance to make result not sensible to scaling
dist1<-dist(solvents3,method = "euclidean")#distance between rows
#hierachical cluster for finding structure in data (solvents)
#I will use different method approaches
cl1<-hclust(dist1,method="single")
plot(cl1)
#does not seem conclusive
cl2<-hclust(dist1,method="complete")
plot(cl2)
#good enough
cl3<-hclust(dist1,method="average")
plot(cl3)
#good enough
#Illustrating clusters
#complete linkage
comp_dend_clust <- as.dendrogram(cl2)
#height cutoff 6
comp_col_dend <- color_branches(comp_dend_clust,h=6)
plot(comp_col_dend)
abline(h=6,lty=2)
#good enough
#average linkage
avg_dend_clust <- as.dendrogram(cl3)
#height cutoff 4
avg_col_dend <- color_branches(avg_dend_clust,h=3.8)
plot(avg_col_dend)
abline(h=3.8, lty=2)
#almost good but for this case one data point seems like an outlier.
#Conclusion: Complete linkage seems the best approach for hierarchical clustering
#of this data
# 5 clusters seems the solution
#Cutting tree of complete linkage
solv_cut_comp_tree<-cutree(cl2,h=6)
solv_cut_comp_tree
#members count
table(solv_cut_comp_tree)
#solvents and their respective cluster
solvents2["cluster"]=solv_cut_comp_tree
View(solvents2)
#some pairwise plots,
par(mfcol=c(2,2))
plot(solvents2$MeltingPoint,solvents2$BoilingPoint,col=solvents2$cluster,
     xlab = "Melting Point",ylab = "Boiling Point",main = "Pairwise plot")
plot(solvents2$RefractiveIndex,solvents2$Density,col=solvents2$cluster,
     xlab = "Refractive Index",ylab = "Density",main = "Pairwise plot")
plot(solvents2$DipoleMoment,solvents2$BoilingPoint,col=solvents2$cluster,
     xlab = "Dipole Moment",ylab = "Boiling Point",main = "Pairwise plot")
plot(solvents2$Solubility,solvents2$Density,col=solvents2$cluster,
     xlab = "Solubility",ylab = "Density",main = "Pairwise plot")
# 2 cluster groups seems confident enough across variables. 
#The others clusters may lack more points to distinguish them clearly
# or better physical and chemical properties to distinguish such solvents
# or another technique to distinguish them
# looking again at the biplot
par(mfcol=c(1,1))
labels=as.factor(solv_cut_comp_tree)
#Final output for chemical intepretation (in the report)
fviz_pca_biplot(pca_results,addEllipses = T,col.ind = labels,
                col.var =fA_var,title = "Biplot and intituive relation to FA  and  cluster groups")
