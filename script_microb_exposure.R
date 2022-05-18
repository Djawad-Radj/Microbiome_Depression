library(TwoSampleMR)

for (i in 1:13) {

list<-c("family.Ruminococcaceae.id.2050","genus.Eggerthella.id.819","genus.Hungatella.id.11306","genus.Lachnoclostridium.id.11308","genus..Eubacteriumventriosumgroup.id.11341","genus.RuminococcaceaeUCG002.id.11360","genus.RuminococcaceaeUCG003.id.11361","genus.Sellimonas.id.14369","genus.Coprococcus3.id.11303","genus.LachnospiraceaeUCG001.id.11321","genus.Subdoligranulum.id.2070","genus.RuminococcaceaeUCG005.id.11363","genus..Ruminococcusgauvreauiigroup.id.11342")

prob<-paste(list[i],".summary.txt",sep="")

a<-read_exposure_data(prob, sep="\t", snp_col = "rsID", beta_col = "beta", se_col="SE", effect_allele_col = "eff.allele", other_allele_col = "ref.allele",pval_col="P.weightedSumZ", eaf_col="freq_eff", phenotype_col= "bac", samplesize_col = "N")

b<-a[which(a$pval.exposure < 0.00001),]

mic_exp<-clump_data(b, clump_r2=0.3)

outc<-read_outcome_data(file="PGC_UKB_depression_genome-wide.txt", snp_col = "MarkerName", beta_col = "LogOR", se_col="StdErrLogOR", effect_allele_col = "A1", other_allele_col = "A2",eaf_col="Freq", pval_col="P", ncase_col="ncase", ncontrol_col="ncont", sep="\t")

d<- outc$SNP %in% mic_exp$SNP

mdd_out<-outc[which(d==TRUE),]

dat<-harmonise_data(exposure_dat= mic_exp, outcome_dat = mdd_out)

outfile<-paste(list[i],"_MDD_dat.txt",sep="")

df<-add_rsq(dat)

write.table(df,file=outfile, sep="\t",row.names=F,quote=F)

mr_report(dat,output_type="md")

}
