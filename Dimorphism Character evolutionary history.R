#requeried liberaries
library(ape)
library(maps)
library(phytools)
library(geiger)
library(coda)
install.packages("ggtree")
library(ggtree)
install.packages("plotrix")
library(plotrix)
BiocManager::install("ggtree")a
#uplode newick tree formate and dataset
tree=read.tree("rooted.nwk")
data=read.table("Character_data.csv", sep = ",", header = T, row.names = 1)

#check the consistencey of name in your provieded data
check <- name.check(tree, data)

#Check the topology of ur tree, rooted and ultrametric, also for data numerical values are necessary to express your character states
plot(tree,show.tip.label = T, no.margin = T)
is.rooted(tree)
is.ultrametric(tree)
head(data)
Dimorphism=data$Dimorphism
names(Dimorphism)=row.names(data)
Dimorphism=factor(Dimorphism)
levels(Dimorphism)
levels(Dimorphism)=c("Monomorphic", "Dimorphic")
as.character(Dimorphism)
f<-function(n) c("#4B006E", "#CCCCFF")
plot(tree,show.tip.label = F, no.margin = T)
tiplabels(pie=to.matrix(Dimorphism[tree$tip.label],levels(Dimorphism)),piecol = cols, cex=0.25)

#plot(ladderize(tree),show.tip.label = T, x.lim=2,no.margin = T)
#phydataplot(as.character(Dimorphism),tree,"m",border=NA,offset = 1,width =0.25,funcol = f,legend = "side")

#u_tree <- chronos(tree)
#write.tree(u_tree, file = "rooted.nwk")

dev.off()
cols<-setNames(c("#CCCCFF","#4B006E"),levels(Dimorphism))

#plotTree(tree,type="arc",fsize=0.8,lwd=2,ftype="reg")
#tiplabels(pie=to.matrix(Dimorphism[tree$tip.label],levels(Dimorphism)),piecol = cols,cex=0.1)
#add.simmap.legend(colors = cols,prompt = T)

#Test the most fit model of states changes to our data
fit.ER<-fitMk(tree,Dimorphism,model = "ER")
plot(fit.ER)
print(fit.ER,digits = 6)
fit.SYM<-fitMk(tree,Dimorphism,model = "SYM")
plot(fit.SYM)
print(fit.SYM,digits = 6)

fit.ARD<-fitMk(tree,Dimorphism,model = "ARD")
plot(fit.ARD)
print(fit.ARD,digits = 6)
AIC(fit.ER,fit.SYM,fit.ARD)



#fit.ER.mcmc=mcmcMk(tree,Dimorphism,model = "ER",ngen =1000000)
#plot(fit.ER.mcmc)
#summary(fit.ER.mcmc,burnin=0.2)
#density(fit.ER.mcmc,burnin=0.2)
#plot(density(fit.ER.mcmc,burnin=0.2))

#fit.SYM.mcmc=mcmcMk(tree,Dimorphism,model = "SYM",ngen = 50000)
#plot(fit.SYM.mcmc)
#summary(fit.SYM.mcmc,burnin=0.2)
#density(fit.SYM.mcmc,burnin=0.2)
#plot(density(fit.SYM.mcmc,burnin=0.2))

#fit.ARD.mcmc=mcmcMk(tree,Dimorphism,model = "ARD",ngen = 1000000)
#plot(fit.ARD.mcmc)
#summary(fit.ARD.mcmc,burnin=0.2)
#density(fit.ARD.mcmc,burnin=0.2)
#plot(density(fit.ARD.mcmc,burnin=0.2))


#tree$edge.length <- exp(tree$edge.length)

#Fitting our best model and check the likilhood marginal propabilities for each node including ancestral one
fitARDace<-ace(Dimorphism,tree,model="ER",type="discrete")
fitARDace$lik.anc
dev.off()

#Plot the tree and visulaize the marginal propapilities, i.e. from liklihood approch of the ancestral stateeeeeeee
#plot(tree,show.tip.label = F,no.margin = T, use.edge.length=F)

plot(tree,show.tip.label = F,type = "fan", no.margin=T)
nodelabels(pie=fitARDace$lik.anc,piecol=cols,
           cex=0.25)
tiplabels(pie=to.matrix(Dimorphism[tree$tip.label],levels(Dimorphism)),piecol = cols,
          cex=0.25,offset = 0)
add.simmap.legend(colors = cols,prompt = F,x=0.9*par()$usr[4],y=0.6*par()$usr
                  [3],fsize=1.2)

#Another more robust approach using the posterior probabilities to reconstruct the ancestral state known as stochastic mapping
stochasticTree=make.simmap(tree, Dimorphism, model="ARD",nsim=100000)


#Summarize the stochastic mappings
stochasticTreeSummary <- summary(stochasticTree)

#Check the posterior probabilities for each node in ur tree
stochasticTreeSummary$ace

#plot the tree and visualize the ancestral state from the posterior probabilities approch 
cols2<-setNames(c("#4B006E","#CCCCFF"),levels(Dimorphism))
plot(tree,type="fan",show.tip.label = F,cex = 0.25, no.margin = T)
nodelabels(node=1:tree$Nnode+Ntip(tree),pie = stochasticTreeSummary$ace, piecol = cols2, cex = 0.25)
tiplabels(pie=to.matrix(Dimorphism[tree$tip.label],levels(Dimorphism)),piecol = cols, cex=0.25)
add.simmap.legend(colors = cols,prompt = F,x=0.9*par()$usr[1],y=0.9*par()$usr
                  [4],fsize=1.5)

#To visulize the summery of simulated trees (100000 stochastic trees), use denisty map function
TransitionsOverBranch<-densityMap(stochasticTree, states=levels(Dimorphism)[2:1], plot = F)
dMapDimo<-setMap(TransitionsOverBranch,cols2)
plot(dMapDimo, type="fan",cex = 0.25, no.margin = T,ftype = "off")
nodelabels(node=1:tree$Nnode+Ntip(tree),pie = stochasticTreeSummary$ace, piecol = cols2, cex = 0.25)
tiplabels(pie=to.matrix(Dimorphism[tree$tip.label],levels(Dimorphism)),piecol = cols, cex=0.25)
add.simmap.legend(colors = cols,prompt = F,x=0.9*par()$usr[1],y=0.9*par()$usr
                  [4],fsize=1.5)

#To count the numbers of changes between the character states
dd=density.multiSimmap(stochasticTree,method="changes")
plot(dd)
change.counts <- function(x, name) {
  min <- min(which(x != 0)) - 1
  max <- max(which(x != 0)) + 1
  
  # Calculate maximum value for ylim
  max_val <- max(x[min:max])
  
  # Set graphical parameters for bold labels and title
  par(cex.axis = 1.2, # Size of axis tick labels
      cex.lab = 1.5,  # Size of axis labels
      cex.main = 1.8, # Size of the main title
      font.axis = 2,  # Bold axis tick labels
      font.lab = 2,   # Bold axis labels
      font.main = 2)  # Bold main title
  
  # Create barplot with customized appearance
  barplot(x[min:max], names.arg = c(min:max), 
          xlab = "Number of Changes", 
          ylab = "Total count", 
          main = name, 
          ylim = c(0, max_val * 1.1))  # Extend ylim to be slightly larger than max_val
}

par(mfrow = c(1, 2))
counts <- dd$p$`Monomorphic->Dimorphic`$counts
change.counts(counts, "Monomorphic to Dimorphic")
counts <- dd$p$`Dimorphic->Monomorphic`$counts
change.counts(counts, "Dimorphic to Monomorphic")



#Test the phylogentic signal of dimorphism character
pagelsLambdaTest=fitDiscrete(tree, Dimorphism, model="ARD",transform = "lambda")
pagelsLambdaTest
tree0<-rescale(tree,model = "lambda",0)
tree1<-rescale(tree,model = "lambda",1)
lambda_L0=fitDiscrete(phy = tree0,dat = Dimorphism,model="ARD")
lambda_L1=fitDiscrete(phy = tree1,dat = Dimorphism,model="ARD")
LLR0<- -2*(lambda_L0$opt$lnL-pagelsLambdaTest$opt$lnL)
LLR1<- -2*(lambda_L1$opt$lnL-pagelsLambdaTest$opt$lnL)
pchisq(LLR0,df=1,lower.tail = F)
pchisq(LLR1,df=1,lower.tail = F)

as.Qmatrix.ace<-function(x,...){
  if("index.matrix"%in%names(x)){
    k<-nrow(x$index.matrix)
    Q<-matrix(NA,k,k)
    Q[]<-c(0,x$rates)[x$index.matrix+1]
    rownames(Q)<-colnames(Q)<-colnames(x$lik.anc)
    diag(Q)<--rowSums(Q,na.rm = T)
    class(Q)<-"Qmatrix"
    return(Q)
    } else cat ("\"ace\" object does not appear to contain a Q matrix.\n")
  }
Q=as.Qmatrix.ace(fitARDace)
dev.off()
plot(Q)

#Correlation of dimorphism and lifestyles of our fungi
head(data)

Lifestyle=setNames(data$Lifestyle, row.names(data))
Lifestyle=factor(Lifestyle)
levels(Lifestyle)
levels(Lifestyle)=c("Non-pathogen", "Pathogen", "Ambrosia", "Non-ambrosia")
cols4<-setNames(c("#94d2bd","#001219","#ee9b00", "#9b2226"),levels(Lifestyle))


dimorphismMap<- densityMap(stochasticTree, states=levels(Dimorphism)[2:1], plot = F)
dimorphism_tree<-setMap(dimorphismMap,cols2)

stochasticTree1=make.simmap(tree, Lifestyle, model="ARD",nsim=1000)
lifestyleMap<- densityMap(stochasticTree1, states=levels(Lifestyle)[2:1], plot = F)
lifestyle_tree<-setMap(lifestyleMap,cols4)


layout(matrix(1:3,1,3),widths=c(0.42,0.12,0.42))


plot(dimorphism_tree,ftype="off",legend=F, lwd=4,outline=FALSE)
add.simmap.legend(colors=cols,prompt=F,x=0,y=10,vertical=T,fsize=1.5)
plot.new()
plot(lifestyle_tree,ftype="off",legend=F, lwd=4,outline=FALSE,direction = "leftwards")
add.simmap.legend(colors=cols4,prompt=F,x=0.75,y=25,vertical=T,fsize=1.2)


dev.off()

ambrosialevels=factor(Lifestyle)
levels(ambrosialevels)
levels(ambrosialevels)=c("Non-ambrosia","Non-ambrosia", "Ambrosia", "Non-ambrosia")

fit.xy=fitPagel(tree,Dimorphism,ambrosialevels)

fit.x=fitPagel(tree,Dimorphism,ambrosialevels,dep.var = "x")
fit.y=fitPagel(tree,Dimorphism,ambrosialevels,dep.var = "y")
aic=setNames(c(fit.xy$independent.AIC,fit.x$dependent.AIC,fit.y$dependent.AIC,fit.xy$dependent.AIC),c("Independent", "dependent x", "dependent y", "dependent x&y"))
aic
aic.w(aic)
plot(fit.xy,lwd.by.rate=T)

pathogenlevels=factor(Lifestyle)
levels(pathogenlevels)
levels(pathogenlevels)=c("Non-pathogen","Pathogen", "Non-pathogen", "Non-pathogen")

fit.xy1=fitPagel(tree,Dimorphism,pathogenlevels)

fit.x1=fitPagel(tree,Dimorphism,pathogenlevels,dep.var = "x")
fit.y1=fitPagel(tree,Dimorphism,pathogenlevels,dep.var = "y")
aic1=setNames(c(fit.xy1$independent.AIC,fit.x1$dependent.AIC,fit.y1$dependent.AIC,fit.xy1$dependent.AIC),c("Independent", "dependent x", "dependent y", "dependent x&y"))
aic1
aic.w(aic1)
plot(fit.x1,lwd.by.rate=T)


write.nexus(tree, file = "tree.nexus")

