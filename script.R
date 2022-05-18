library(olsrr)
library(vegan)
library(MKmisc)
library(ResourceSelection)
library(randomForest)

otu = read.table("RS_taxon_table.txt", header=T, sep="\t")
dat = read.table("MappingFile_complete_with_exclusions.txt", header=T, sep="\t")
map = read.table("RS_depression.txt", header=T, sep="\t")
otu = otu[,-c(ncol(otu))]

phe.names = colnames(map[,c(3:ncol(map))])
fin = merge(dat, map, by="RSID")

fin = merge(fin, otu, by="RSID")

tax.names = colnames(otu[,c(3:ncol(otu))])

#############  Taxa #################
results <- data.frame(
            Pheno=as.character(),
            Taxon=as.character(),
            Beta=as.numeric(),
            Se=as.numeric(),
            p=as.numeric(),
            n=as.numeric(),
            resi.corr = as.numeric(),
            stringsAsFactors=FALSE)
for (j in 1:length(tax.names)) {
if (sum(otu[,tax.names[j]] != 0) >= 0.03*dim(otu)[1]){
fin$taxon <- log(fin[,tax.names[j]]+1)
fit <- lm(e5_eiscorewgt ~ taxon + Sex + Age + BMI + Smoking2 + Alcohol_gram_imputed + TimeInMail + Batch + lipidlowering + PPI + metformine + antibiotics , data=fin)
fit.sum <- summary(fit)
tablerow <- data.frame(
Pheno = "e5_eiscorewgt",
Taxon = tax.names[j],
Beta = fit.sum$coefficients[2,1],
Se = fit.sum$coefficients[2,2],
p = fit.sum$coefficients[2,4],
n = dim(fin)[1],
resi.corr = ols_test_correlation(fit),
stringsAsFactors=FALSE)
print(tablerow)
results <- rbind(results, tablerow)
}
}
results$fdr <- p.adjust(results$p, method="fdr")
write.table(results, "RS.LINEAR.MODEL3.TAXA.FINAL.txt", col.names=T, row.names=F, sep="\t", quote=F)

############## Alpha ############# 
genus = fin[, grep("genus.",colnames(fin))]
shan = diversity(genus, "shannon")
invsimp = diversity(genus, "invsimpson")
rich = specnumber(genus)
alphas = cbind(shan, rich, invsimp)
alpha.names = colnames(alphas)
fin = cbind(fin, alphas)
results <- data.frame(
            Pheno=as.character(),
            Alpha=as.character(),
            Beta=as.numeric(),
            Se=as.numeric(),
            p=as.numeric(),
            n=as.numeric(),
            resi.corr = as.numeric(),
            stringsAsFactors=FALSE)
for (j in 1:length(alpha.names)) {
fin$alpha <- fin[,alpha.names[j]]
fit <- lm(e5_eiscorewgt ~ alpha+ Sex + Age + BMI + Smoking2 + Alcohol_gram_imputed + TimeInMail + Batch + lipidlowering + PPI + metformine + antibiotics , data=fin)
fit.sum <- summary(fit)
tablerow <- data.frame(
Pheno = "e5_eiscorewgt",
Alpha = alpha.names[j],
Beta = fit.sum$coefficients[2,1],
Se = fit.sum$coefficients[2,2],
p = fit.sum$coefficients[2,4],
n = dim(fin)[1],
resi.corr = ols_test_correlation(fit),
stringsAsFactors=FALSE)
print(tablerow)
results <- rbind(results, tablerow)
}
results$fdr <- p.adjust(results$p, method="fdr")
write.table(results, "RS.LINEAR.MODEL3.ALPHAS.FINAL.txt", col.names=T, row.names=F, sep="\t", quote=F)

############# Beta ################# 
bc = vegdist(genus, "bray")
ad.rs = adonis(bc ~ Age+Sex+Smoking2+BMI+antibiotics+Alcohol_gram_imputed+TimeInMail + Batch + lipidlowering + PPI + metformine +e5_eiscorewgt, data=fin)


############# RandomForest ############
newRS <- cbind(fin$RSID, fin$e5_eiscorewgt, genus)
newHL = read.table("../RandomForest/HL.DATA.FOR.RF.txt", header=T, check.names=F, row.names=1, sep="\t")

names(newRS)[1] <- "RSID"
names(newRS)[2] <- "eiscorewgt"

common_cols <- intersect(colnames(newRS), colnames(newHL))

hl <- hl.data[,names(newHL) %in% common_cols]
hl$DepScore <- newHL$PHQ9_SS
hl <- hl[,c(357,1:356)]

rs <- newRS[,names(newRS) %in% common_cols]
rs$DepScore <- newRS$eiscorewgt
rs <- rs[,c(357,1:356)]
rownames(rs) <- newRS$RSID

write.table(hl, "HL.DATA.FOR.RF2.txt", col.names=T, row.names=T, quote=F, sep="\t")
write.table(rs, "RS.DATA.FOR.RF2.txt", col.names=T, row.names=T, quote=F, sep="\t")
## include SampleID into first line of each file

rs.data = read.table("RS.DATA.FOR.RF2.txt", header=T, check.names=F, row.names=1, sep="\t")
hl.data = read.table("HL.DATA.FOR.RF2.txt", header=T, check.names=F, row.names=1, sep="\t")

rf.rs.hl <- randomForest(rs.data[-1], rs.data[,1], xtest=hl.data[,-1], ytest=hl.data[,1], ntree=500, importance=TRUE, nPerm=100, do.trace=10, proximity=TRUE, oob.prox=TRUE, keep.forest=TRUE)

pdf("varimpPlot_rs_hl_final.pdf", height=20, width=20)
varImpPlot(rf.rs.hl, n.var=50)
dev.off()

imp <- rf.rs.hl$importance
imp.sd <- rf.rs.hl$importanceSD
output.table <- cbind(imp, imp.sd)
write.table(output.table, "feature_importance_scores_helius_rs_continious_final.txt", sep="\t", quote=F)

tabel(fin$Age)
table(fin$Age)
summary(fin$Age)
sd(fin$Age)
summary(fin$BMI)
sd(fin$BMI)
table(fin$Smoking)
table(fin$Smoking1)
summary(fin$e5_eiscorewgt)
sd(fin$e5_eiscorewgt)
mdd <- fin[which(fin$e5_eiscorewgt >= 21),]
dim(mdd)
savehistory("nov_2020_history.R")
