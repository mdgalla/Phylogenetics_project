setwd("C:/Users/galla/Documents/Masters/Semester_2/Phylogenetics/project")
set.seed(124)

#Packages
library(TreeSim)
library(phyclust)
library(ape)


# Part 1 simulate trees
#test birth and death rates? clock rate? Divergent time?
birth <- 0.1
death <- 0.05
n <- 20
numsim <- 3

trees <- TreeSim::sim.bd.taxa(n, numsim, birth, death, complete=FALSE)

time_tree <- trees[[1]]
plot(time_tree)

# root age
max(TreeSim::getx(time_tree))
write.tree(time_tree, "data/sim_tree")

# Part 2 - simulate clock rates

bl <- length(time_tree$edge.length)

mean_rate1 = 0.2
mean_rate2 = 0.01

rates1 = rlnorm(bl, log(mean_rate1), sdlog = 0.1) # rates1 ~ mean_rate1
rates2 = rlnorm(bl, log(mean_rate2), sdlog = 0.01)

clock_tree1 <- time_tree
clock_tree1$edge.length <- time_tree$edge.length * rates1
plot(clock_tree1)
ape::write.tree(clock_tree1)

clock_tree2 <- time_tree
clock_tree2$edge.length <- time_tree$edge.length * rates2
plot(clock_tree2)

# Part 3 - simulate character data

#make a sequence (nexus file) for clock tree 1
seq1<-phyclust::seqgen(opts="-l500 -mGTR", rooted.tree=clock_tree1)
#make a sequence (nexus file) for clock tree 2 
seq2<-phyclust::seqgen(opts="-l500 -mGTR", rooted.tree=clock_tree2)

final_seq=""
names = c()
for (i in 2:length(seq1)){
  names <- c(names, strsplit(seq1[i], " ")[[1]][1])
  a<-gsub(strsplit(seq1[i], " ")[[1]][1],"", seq1[i])
  b<-gsub(substring(seq1[i],1, 3),"",seq2[i])
  c<-gsub(" ","", b)
  final_seq[i-1]<- gsub(" ", "", paste(a,c))
}

names(final_seq)<-names

write.nexus.data(final_seq, file="data/seq_trial1.nex")


####Repeat####

##make a function to replicate sequence data and save#

mysim.seq<-
  function(tree1=clock_tree1, tree2=clock_tree2, seqfile){
  
#make a sequence (nexus file) for clock tree 1
seq1<-phyclust::seqgen(opts="-l500 -mGTR", rooted.tree=tree1)
#make a sequence (nexus file) for clock tree 2 
seq2<-phyclust::seqgen(opts="-l500 -mGTR", rooted.tree=tree2)

#combine the nexus files 
final_seq=""
names = c()
for (i in 2:length(seq1)){
  names <- c(names, strsplit(seq1[i], " ")[[1]][1])
  a<-gsub(strsplit(seq1[i], " ")[[1]][1],"", seq1[i])
  b<-gsub(substring(seq1[i],1, 3),"",seq2[i])
  c<-gsub(" ","", b)
  final_seq[i-1]<- gsub(" ", "", paste(a,c))
}

names(final_seq)<-names

write.nexus.data(final_seq, file=seqfile)}


###generate sequence data###
set.seed(5689)
mysim.seq(clock_tree1, clock_tree2, "seq_trial2")
set.seed(567)
mysim.seq(clock_tree1, clock_tree2, "seq_trial3")


