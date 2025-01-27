#PROTEIN
library(caret)
library(MASS)
library(car)
library(randomForest)
library(glmnet)
library(gbm)
load("/home/jasiek/Dokumenty/Grupa1/protein.RData")
X = as.matrix(data.train[,-ncol(data.train)])
czy.numeryczne=sapply(X[1,],is.numeric)
length(X[1,czy.numeryczne])==(ncol(data.train)-1)
Y = as.matrix(data.train[,ncol(data.train)])
colnames(data.train)
sum(is.na(data.train))
dane.min = apply(data.train, 2, min)
dane.max = apply(data.train, 2, max)
plot(rep(1:ncol(data.train),2),c(dane.min,dane.max))
ncol(X)
#dane.protein=as.data.frame(cbind(X,Y))
dane.protein=as.data.frame(data.train)
liczba.predyktorow=c(1:20)
cv.lasy = rfcv(dane.protein[,-ncol(dane.protein)], dane.protein$Y, cv.fold =5 , type="regression")
cv.lasy$error.cv
las = randomForest(Y~. ,dane.protein, mtry=8, importance=TRUE)
plot(las)
varImpPlot(las)
varUsed(las)

data.split = createFolds(dane.protein$Y, k=5, returnTrain = FALSE)
train.set = lapply(data.split, function(x) as.data.frame(dane.protein[-x,]) )

full.formula1=terms(Y~.,data=train.set$Fold1)
full.formula1=formula(full.formula)
null.model1=lm(Y~1, train.set$Fold1)
full.formula2=terms(Y~.,data=train.set$Fold2)
full.formula2=formula(full.formula)
null.model2=lm(Y~1, train.set$Fold2)
full.formula3=terms(Y~.,data=train.set$Fold3)
full.formula3=formula(full.formula)
null.model3=lm(Y~1, train.set$Fold3)
full.formula4=terms(Y~.,data=train.set$Fold4)
full.formula4=formula(full.formula)
null.model4=lm(Y~1, train.set$Fold4)
full.formula5=terms(Y~.,data=train.set$Fold5)
full.formula5=formula(full.formula)
null.model5=lm(Y~1, train.set$Fold5)


full.formulas=lapply(train.set, function (x) terms(Y~.,data=x) )
full.formulas=lapply(full.formulas, function(x) formula(x))

#full.formulas=lapply(train.set, function(x) as.formula(x))

full.formula=terms(Y~.,data=dane.protein)
full.formula=formula(full.formula)


my.cv = function() {
  l=numeric(48)
  i = 1
  forward.model1 = stepAIC(null.model1, full.formula1, direction = "forward", steps=1)
  forward.model2 = stepAIC(null.model1, full.formula2, direction = "forward", steps=1)
  forward.model3 = stepAIC(null.model1, full.formula3, direction = "forward", steps=1)
  forward.model4 = stepAIC(null.model1, full.formula4, direction = "forward", steps=1)
  forward.model5 = stepAIC(null.model1, full.formula5, direction = "forward", steps=1)
  mse_old = 1000000000000000
  mse_new = 900000000000
  while (i < 49){
    mse_old = mse_new
    forward.model1 = stepAIC(forward.model1, full.formula1, direction = "forward", steps=1)
    forward.model2 = stepAIC(forward.model2, full.formula2, direction = "forward", steps=1)
    forward.model3 = stepAIC(forward.model3, full.formula3, direction = "forward", steps=1)
    forward.model4 = stepAIC(forward.model4, full.formula4, direction = "forward", steps=1)
    forward.model5 = stepAIC(forward.model5, full.formula5, direction = "forward", steps=1)
    mse1 = sum((predict(forward.model1, newdata = dane.protein[data.split$Fold1,-ncol(dane.protein)])-dane.protein$Y[data.split$Fold1])^2)/length(data.split$Fold1)
    mse2 = sum((predict(forward.model2, newdata = dane.protein[data.split$Fold2,-ncol(dane.protein)])-dane.protein$Y[data.split$Fold2])^2)/length(data.split$Fold2)
    mse3 = sum((predict(forward.model3, newdata = dane.protein[data.split$Fold3,-ncol(dane.protein)])-dane.protein$Y[data.split$Fold3])^2)/length(data.split$Fold3)
    mse4 = sum((predict(forward.model4, newdata = dane.protein[data.split$Fold4,-ncol(dane.protein)])-dane.protein$Y[data.split$Fold4])^2)/length(data.split$Fold4)
    mse5 = sum((predict(forward.model5, newdata = dane.protein[data.split$Fold5,-ncol(dane.protein)])-dane.protein$Y[data.split$Fold5])^2)/length(data.split$Fold5)
    print(c(mse1, mse2, mse3, mse4, mse5))
    mse_new = mean(c(mse1,mse2,mse3,mse4,mse5))
    #mse_new = mse1
    l[[i]]=mse_new
    i = i+1
    print(mse_new)
    plot(1:i, l[1:i], type="l")
    #stop()
  }
  return(l, c(forward.model1, forward.model2, forward.model3, forward.model4, forward.model5))
}

aic.cv = my.cv()

full.formula=terms(Y~.,data=dane.protein)
full.formula=formula(full.formula)
null.model=lm(Y~1, dane.protein)
forward.model = stepAIC(null.model, full.formula, direction = "forward", steps=50)
dane.test=as.data.frame(data.test)
pred.protein=predict(forward.model, newdata = dane.test)


#CANCER
load("/home/jasiek/Dokumenty/Grupa1/cancer.RData")
X = as.matrix(data.train[,-ncol(data.train)])
czy.numeryczne=sapply(X[1,],is.numeric)
length(X[1,czy.numeryczne])==(ncol(data.train)-1)
Y = as.matrix(data.train[,ncol(data.train)])
colnames(data.train)
sum(is.na(data.train))
dane.min = apply(data.train, 2, min)
dane.max = apply(data.train, 2, max)
plot(rep(1:ncol(data.train),2),c(dane.min,dane.max))
ncol(X)
#dane.max[dane.max<2]
dane.cancer=as.data.frame(data.train)
cv.lasy = rfcv(dane.cancer[,-ncol(dane.cancer)], dane.cancer$Y, cv.fold =5 , type="regression")
las2 = randomForest(dane.cancer[,-ncol(dane.cancer)], dane.cancer$Y, mtry=17737, importance=TRUE, ntree = 500)


X.train = as.matrix(dane.cancer[,-ncol(dane.cancer)])
Y.train = as.matrix(dane.cancer$Y)

alfy=seq(0.0, 0.02, length.out = 5)



mses = lapply(alfy, function(x) cv.glmnet(X.train, Y.train, nfolds = 5, standardize=T, alpha=x))

result.cv=cv.glmnet(X.train, Y.train, nfolds = 5, standardize=T, alpha=1)

optimal.result=coef(mses[[2]])
absy = abs(optimal.result)

absy2 = order(absy, decreasing=TRUE)[1:101]
absy3 = (absy2-1)
absy4 = absy3[2:length(absy3)]

dane.final = dane.cancer[,c(ind,ncol(dane.cancer))]

final.model = lm(Y~., data=as.data.frame(dane.final))

pred.cancer = predict(final.model, as.data.frame(data.test[,ind]))


uporzadk = order(absy)

ind = which(optimal.result!=0)-1
ind = ind[2:length(ind)]

vif.cancer = car::vif(final.model)
hist(vif.cancer)


model=gbm.fit(x= X.train, y=Y.train, n.trees=100, distribution = "gaussian")


data.split = createFolds(dane.cancer$Y, k=5, returnTrain = FALSE)
train.set = lapply(data.split, function(x) as.data.frame(dane.cancer[-x,]) )

lambdas=seq(0.001, 0.9, length.out = 20)

trenuj_gbm_i_daj_mse=function(lambda){
  models=lapply(data.split, function(K) gbm.fit(x=X.train[-K,], y=Y.train[-K,], n.trees=100, distribution = "gaussian",shrinkage=lambda))
  mses = lapply(1:5, function(k) sum((predict(models[[k]], newdata = dane.cancer[data.split[[k]],-ncol(dane.cancer)],n.trees=100)-dane.cancer$Y[data.split[[k]]])^2)/length(data.split[[k]]))
  mses = unlist(mses)
  mse = mean(mses)
  print(c(mses, mse, lambda))
  return(mse)
}

mse_finals=sapply(lambdas, trenuj_gbm_i_daj_mse)


vif.protein = car::vif(forward.model)
hist(vif.protein)

predictors.protein = names(forward.model$coefficients[2:6])
predictors.protein
save(pred.protein, predictors.protein, pred.cancer, file = "malinowski.Rdata")
