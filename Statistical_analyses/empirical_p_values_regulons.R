#Empirical p-values to cell cycle gene sets (S-phase used as an example)

#read MRD data
mrd_combo=read.table("MRD.txt", header=T)
mrd_combo

#read cell cycle scores
score = read.table("cell_cycle_scores.txt", header=T)

#EOI MRD (day29)
#create matrix
cor29 = matrix(nrow = 10000, ncol = 2)
colnames(cor29)= c("S", "G2M")
cor29

#make randomization and calculate the correlation with 10,000 repetitions
for (y in 1:10000){
mrd_ran <- sample(mrd_combo$day29_flow)

for (i in 1:nrow(score)) {
  cor_29=cor.test(mrd_ran, score[i,], method = "spearman")
  cor29[y,i]=cor_29$p.value }
}

write.table(cor29, "cor29_random.txt", sep="\t", row.names = F, quote = F)

#test how many from the randomized results are smaller than the original p-value
sum(cor29[,1] < 0.045303763)
# e.g.--> 453 smaller than the S-phase day 29 p-value 

#calculate the empirical p-value
emp_p <- (453+1)/(10000+1)



#Mid-induction MRD (day15)

#create matrix
cor15 = matrix(nrow = 10000, ncol = 2)
colnames(cor15)= c("S", "G2M")
cor15

#make randomization and calculate the correlation with 10,000 repetitions
for (y in 1:10000){
  mrd_ran15 <- sample(mrd_combo$day15_flow)
  
  for (i in 1:nrow(score)) {
    cor_15=cor.test(mrd_ran15, score[i,], method = "spearman")
    cor15[y,i]=cor_15$p.value }
}

write.table(cor15, "cor15_random.txt.txt", sep="\t", row.names = F, quote = F)

#test how many from the randomized results are smaller than the original p-value
sum(cor15[,1] < 0.0459237634580901)
# e.g.--> 493 smaller than the S-phase day 15 p-value 

#calculate the empirical p-value
emp_p <- (493+1)/(10000+1)
