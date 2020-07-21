# utils

df <- read.csv("/Users/Frankie/Desktop/syfrankie/20ss/542/project/county_data_apr22.csv", header=T)
id = sort(df[,2])

pred <- read.csv("/Users/Frankie/Desktop/syfrankie/20ss/542/project/usafacts_infections.csv",header = T)
pred = pred[which(pred[,1]%in%id),]
true_y = cbind(pred[,1],pred[,94:100],pred[,200:206])
true_y$tot_cases = apply(true_y[,2:8],1,sum)
true_y$tot_deaths = apply(true_y[,9:15],1,sum)

write.table(true_y, "/Users/Frankie/Desktop/syfrankie/20ss/542/project/cases-deaths-23-29.csv", row.names=F, sep=",")

info_dat = df[,3:88]
covid_dat = df[,89:(ncol(df)-4)]

library(VIM)
info_dat <- kNN(info_dat, k=10, imp_var=F)
which(apply(info_dat[,-(1:12)], 2, sd)==0)
info_dat = info_dat[,-which(names(info_dat) %in% c("federal.guidelines", "foreign.travel.ban"))]


# unsupervise

library(ggplot2)
library(usmap)
#clus1 = kmeans(log1p(covid_dat[,1:92]), centers=4, iter.max=1e4, nstart=10)
clus1 = kmeans(log1p(covid_dat[,93:184]), centers=4, iter.max=1e4, nstart=10)
us_covid = as.data.frame(cbind(df[,11], df[,10], df[,2], clus1$cluster, df[,1]))
colnames(us_covid) = c("lon", "lat", "countyFIPS", "cluster", "group")
us_covid = usmap_transform(us_covid)
plot_usmap(regions="counties", color="grey") + 
  geom_point(data=us_covid, aes(x=lon.1, y=lat.1, group=group, color=factor(cluster)), 
             shape=18, size=2) + theme(legend.position="none")

clus2 = kmeans(scale(cbind(info_dat[,29:33])), centers=4, iter.max=1e4)
us_info = as.data.frame(cbind(df[,11], df[,10], df[,2], clus2$cluster, df[,1]))
colnames(us_info) = c("lon", "lat", "countyFIPS", "cluster", "group")
us_info = usmap_transform(us_info)
plot_usmap(regions="counties", color="grey") + 
  geom_point(data=us_info, aes(x=lon.1, y=lat.1, group=group, color=factor(cluster)), 
             shape=18, size=2) + theme(legend.position="none")


library(anocva) #spectral
W = as.matrix(exp(-dist(as.matrix(scale(cbind(info_dat[,20:33], info_dat[,81:84]))))^2) / 4)
clus2 = spectralClustering(W, k=10)
us_info = as.data.frame(cbind(df[,11], df[,10], df[,2], clus2, df[,1]))
colnames(us_info) = c("lon", "lat", "countyFIPS", "cluster", "group")
us_info = usmap_transform(us_info)
plot_usmap(regions="counties", color="grey") + 
  geom_point(data=us_info, aes(x=lon.1, y=lat.1, group=group, color=factor(cluster)), 
             shape=18, size=2) + theme(legend.position="none")

library(FactoMineR)
#info.pc = prcomp(scale(info_dat[,23:25]))
#info.pc = prcomp(scale(info_dat[,29:33]))
info.pc = prcomp(scale(info_dat[,55:62]))
info.pc = prcomp(scale(info_dat[,49:54]))
plot(info.pc, type = "l", pch = 19)
clus2pc = PCA(scale(info_dat[,49:54]), ncp = 2, graph = FALSE)
clus2 = HCPC(clus2pc, graph = FALSE, min=4)
us_info = as.data.frame(cbind(df[,11], df[,10], df[,2], clus2$data.clust[,7], df[,1]))
colnames(us_info) = c("lon", "lat", "countyFIPS", "cluster", "group")
us_info = usmap_transform(us_info)
plot_usmap(regions="counties", color="grey") + 
  geom_point(data=us_info, aes(x=lon.1, y=lat.1, group=group, color=factor(cluster)), 
             shape=18, size=2) + theme(legend.position="none")

library(kohonen)
#clus3 = som(as.matrix(cbind(scale(info_dat[,23:27]), log1p(covid_dat[,93:184]))), somgrid(xdim=2, ydim=2))
clus3 = som(as.matrix(cbind(scale(info_dat[,81:84]), log1p(covid_dat[,93:184]))), somgrid(xdim=2, ydim=2))
us_covid = as.data.frame(cbind(df[,11], df[,10], df[,2], clus3[["unit.classif"]], df[,1]))
colnames(us_covid) = c("lon", "lat", "countyFIPS", "cluster", "group")
us_covid = usmap_transform(us_covid)
plot_usmap(regions="counties", color="grey") + 
  geom_point(data=us_covid, aes(x=lon.1, y=lat.1, group=group, color=factor(cluster)), 
             shape=18, size=2) + theme(legend.position="none")

#new_data = data.frame(scale(info_dat[,23:27]), death=log1p(covid_dat[,184]), cluster=clus3[["unit.classif"]])
new_data = data.frame(scale(info_dat[,81:84]), death=scale(covid_dat[,184]), cluster=clus3[["unit.classif"]])
new_dat = aggregate(.~cluster, new_data, mean)
#diabetes
plot(new_dat$cluster, new_dat$death, col="orange", pch=20, ylim=c(-1,6), xaxt="n", 
     xlab="Cluster Label", ylab="")
axis(1, at=seq(1,4,by=1), las=1)
points(new_dat$X3.YrDiabetes2015.17, col="skyblue", pch=18)
legend("topleft",c("Mean Diabetes", "Mean Log Death"), col=c("skyblue", "orange"), pch=c(18, 20))
#hpsa shortage
plot(new_dat$cluster, new_dat$death, col="orange", pch=20, ylim=c(-1,4), xaxt="n",
     xlab="Cluster Label", ylab="")
axis(1, at=seq(1,4,by=1), las=1)
points(new_dat$HPSAShortage, col="skyblue", pch=18)
legend("topright",c("Mean HPSA Shortage", "Mean Log Death"), col=c("skyblue", "orange"), pch=c(18, 20))

# supervise
## classification

info_dat$death.rate = covid_dat$tot_deaths/info_dat$PopulationEstimate2018*1e5
info_dat$death.rate = as.factor(ifelse(info_dat$death.rate>1,1,0))
trn_idx = sample(1:3141, 2500, replace=F)
info_trn = info_dat[trn_idx,]
info_tst = info_dat[-trn_idx,]

library(e1071)
cost = c(seq(1,9,by=1), seq(10,90,by=10), seq(100,1000,by=100))
folds <- rep_len(1:10,2500)
acc1 = matrix(NA, length(cost), 10)
for (i in 1:length(cost)){
  for (j in 1:10){
    val_idx = which(folds==j)
    info_est = info_trn[-val_idx,]
    info_val = info_trn[val_idx,]
    class1 = svm(death.rate~., data=info_est, type='C-classification', 
                  kernel = 'radial', scale=T, cost=cost[i])
    pred_c1 = predict(class1, info_val)
    acc1[i,j] = mean(pred_c1==info_val$death.rate)
  }
}
cost[which.max(apply(acc1, 1, mean))]
class1 = svm(death.rate~., data=info_trn, type='C-classification', 
             kernel = 'radial', scale=T, cost=30)
pred_c1 = predict(class1, info_tst)
mean(pred_c1==info_tst$death.rate)
table(pred_value=pred_c1, true_value=info_tst$death.rate)

library(caret)
all_k = c(1:3, seq(5,60,by=5))
acc2 = matrix(NA, length(all_k), 10)
for (i in 1:length(all_k)){
  for (j in 1:10){
    val_idx = which(folds==j)
    info_est = info_trn[-val_idx,]
    info_val = info_trn[val_idx,]
    class2 = knn3(death.rate~., data=info_est, k=all_k[i])
    pred_c2 = predict(class2, info_val, type="class")
    acc2[i,j] = mean(pred_c2==info_val$death.rate)
  }
}
all_k[which.max(apply(acc2, 1, mean))]
class2 = knn3(death.rate~., data=info_trn, k=30)
pred_c2 = predict(class2, info_tst, type="class")
mean(pred_c2==info_tst$death.rate)
table(pred_value=pred_c2, true_value=info_tst$death.rate)

##regression
x = matrix(c(rep(1,15), c(1:15)), 15, 2)
newx = matrix(c(rep(1,7), c(16:22)), 7, 2)
y = matrix(as.numeric(t(covid_dat[16, 163:177])), 15, 1)
newy = matrix(as.numeric(t(covid_dat[16, 178:184])), 7, 1)
library(glmnet)
x = c(1:15)
y = as.numeric(t(covid_dat[16, 163:177]))
newx = c(16:22)
reg1 = lm(y~exp(x))
reg2 = glmnet(x, y, lambda=0.5, alpha=0)
predy = predict(reg1, data.frame(x=exp(newx)))
