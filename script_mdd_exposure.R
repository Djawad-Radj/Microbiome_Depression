library(TwoSampleMR)

a<-read_exposure_data(file="PGC_UKB_depression_genome-wide.txt", snp_col = "MarkerName", beta_col = "LogOR", se_col="StdErrLogOR", effect_allele_col = "A1", other_allele_col = "A2",eaf_col="Freq", pval_col="P",ncase_col="ncase", ncontrol_col="ncont", sep="\t")

b<-a[which(a$pval.exposure < 0.0000001),]

mdd_exp<-clump_data(b, clump_r2=0.3)

prob<-c("family.Ruminococcaceae.id.2050","genus.Eggerthella.id.819","genus.Hungatella.id.11306","genus.Lachnoclostridium.id.11308","genus..Eubacteriumventriosumgroup.id.11341","genus.RuminococcaceaeUCG002.id.11360","genus.RuminococcaceaeUCG003.id.11361","genus.Sellimonas.id.14369","genus.Coprococcus3.id.11303","genus.LachnospiraceaeUCG001.id.11321","genus.Subdoligranulum.id.2070","genus..Ruminococcusgauvreauiigroup.id.11342","genus.RuminococcaceaeUCG005.id.11363")

for (i in 1:13) {

outcfile<-paste(prob[i],".summary.txt",sep="")

mb<-read_outcome_data(outcfile, sep="\t", snp_col = "rsID", beta_col = "beta", se_col="SE", effect_allele_col = "eff.allele", other_allele_col = "ref.allele",pval_col="P.weightedSumZ", phenotype_col="bac", eaf_col="freq_eff",samplesize_col = "N")

d<- mb$SNP %in% mdd_exp$SNP

mb_out<-mb[which(d==TRUE),]

dat<-harmonise_data(exposure_dat= mdd_exp, outcome_dat = mb_out)

outfile<-paste("MDD_",prob[i],"_dat.txt",sep="")

df<-add_rsq(dat)

write.table(df,file=outfile, sep="\t",row.names=F,quote=F)

mr_report(dat,output_type="md")}