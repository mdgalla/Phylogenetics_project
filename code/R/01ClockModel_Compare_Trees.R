setwd("C:/Users/galla/Documents/Masters/Semester_2/Phylogenetics/project")

library(ape)
library(TreeSim)
library(phyclust)
########################## Part1:read in tree ##########################
#simulated tree
time_tree<-read.tree("data/sim_tree")

#generated trees#

    #using relaxed clock
gen_tree1<-read.tree("output/Relaxed_trial1/fixed_tree1_long.tre")
tree1<-gen_tree1[["treeTREE1="]]

gen_tree2<-read.tree("output/relaxed_trial2/fixed_tree2.tre")
tree2<-gen_tree2[["treeTREE1="]]

gen_tree3<-read.tree("output/relaxed_trial3/fixed_tree3.tre")
tree3<-gen_tree3[["treeTREE1="]]

    #using strict clock
gen_tree_strict1<-read.tree("output/fixed_tree1_strict_short.tre")
strict.tree1<-gen_tree_strict1[["treeTREE1="]]

gen_tree_strict3<-read.tree("output/fixed_tree3_strict.tre")
strict.tree3<-gen_tree_strict3[["treeTREE1="]]

gen_tree_strict2<-read.tree("output/fixed_tree2_strict_short.tre")
strict.tree2<-gen_tree_strict2[["treeTREE1="]]

############Part 2: Tree Data#################

#get root age 

  #relaxed clock
real.node.dist<-max(TreeSim::getx(time_tree))
node.dist1<-max(TreeSim::getx(tree1))
node.dist2<-max(TreeSim::getx(tree2))
node.dist3<-max(TreeSim::getx(tree3))

  #strict clock
strict.node.dist1<-max(TreeSim::getx(strict.tree1))
strict.node.dist2<-max(TreeSim::getx(strict.tree2))
strict.node.dist3<-max(TreeSim::getx(strict.tree3))

#relaxed clock
#average branching times 
real.mean<-mean(TreeSim::getx(time_tree))
tree1.mean<-mean(TreeSim::getx(tree1))
tree2.mean<-mean(TreeSim::getx(tree2))
tree3.mean<-mean(TreeSim::getx(tree3))

#range of branching times
range(TreeSim::getx(time_tree))
range(TreeSim::getx(tree1))
range(TreeSim::getx(tree2))
range(TreeSim::getx(tree3))

#branching times
branch1<-(TreeSim::getx(tree1))
branch1<-sort(branch1)

branch2<-(TreeSim::getx(tree2))
branch2<-sort(branch2)

branch3<-(TreeSim::getx(tree3))
branch3<-sort(branch3)

#strict clock branching times
strict.branch1<-TreeSim::getx(strict.tree1)
strict.branch1<-sort(strict.branch1)

strict.branch2<-TreeSim::getx(strict.tree2)
strict.branch2<-sort(strict.branch2)

strict.branch3<-TreeSim::getx(strict.tree3)
strict.branch3<-sort(strict.branch3)

#simulated tree branching times
real_branch<-TreeSim::getx(time_tree)
real_branch<-sort(real_branch)

#dataframe with node ages/branching times
branch_lengths<-data.frame(real_branch, branch1, branch2, branch3, strict.branch1,
                           strict.branch2, strict.branch3, row.names=names(branch1))

for (i in 1:length(names(real_branch))){
  z<-names(real_branch)[i]
  branch<-branch_lengths[z,2:4]
  branch.strict<-branch_lengths[z,5:7]
  branch.mean<-mean(as.numeric(branch))
  branch.mean.strict<-mean(as.numeric(branch.strict))
  branch_lengths$relaxed_branch_means[i]<-branch.mean
  branch_lengths$strict_branch_means[i]<-branch.mean.strict
}


#add percent difference
#relaxed
branch_lengths$relaxed_percent<-branch_lengths$relaxed_branch_means/
  branch_lengths$real_branch

#strict
branch_lengths$strict_percent<-branch_lengths$strict_branch_means/
  branch_lengths$real_branch

#####Dataframe for other tree data####

#make a container#
tree_dat<-data.frame(matrix(nrow=6, ncol=7))
colnames(tree_dat)<-c("trial", "root_age", "root_error",
                      "branch_length_mean", "range_min", "range_max", "branch_std")
rownames<-c("relaxed1", "relaxed2", "relaxed3", "strict1", "strict2", "strict3")

#root ages 
root_ages<-c(node.dist1, node.dist2, node.dist3, 
             strict.node.dist1, strict.node.dist2, strict.node.dist3)
tree_dat$root_age<-root_ages

#stderror
root_errors<-c(0.2907, 0.2051, 0.1786, 0.2004, 0.1568, 0.0568)
tree_dat$root_error<-root_errors

tree_dat$trial<-rownames

#add branch length means 
tree_dat$branch_length_mean<-c(tree1.mean, tree2.mean, tree3.mean, mean(strict.branch1), mean(strict.branch2), mean(strict.branch3))
tree_dat$branch_std<-c(sd(branch1), sd(branch2), sd(branch3), NA, NA, NA)
tree_dat$range_min<-c(range(branch1)[1], range(branch2)[1], range(branch3)[1], NA, NA, NA)
tree_dat$range_max<-c(range(branch1)[2], range(branch2)[2], range(branch3)[2], NA, NA, NA)

##################Part 3: Comparisons##########################

##scatter plot##

  ##relaxed clock
trial_means <- c(root_ages[1:3], real.node.dist)      # Means of different trials
trial_errors <- c(root_errors[1:3], 0)    # Corresponding standard errors
trial_labels <- c("Trial 1", "Trial 2", "Trial 3", "Actual Root Age")

# Create the scatter plot for strict clock
plot(1:length(trial_means), trial_means, ylim=c(0, max(trial_means + trial_errors) + 1), 
     xaxt="n", xlab="Trials", ylab="Root Age", 
     main="Relaxed Clock Root Ages", pch=19)

# Add the x-axis labels
axis(1, at=1:length(trial_means), labels=trial_labels)

# Add error bars
arrows(x0 = 1:length(trial_means), y0 = trial_means - trial_errors, 
       x1 = 1:length(trial_means), y1 = trial_means + trial_errors, 
       angle = 90, code = 3, length = 0.1, col = "black")



  ##strict clock
trial_means <- c(root_ages[4:6], real.node.dist)      # Means of different trials
trial_errors <- c(root_errors[4:6], 0)    # Corresponding standard errors
trial_labels <- c("Trial 1", "Trial 2", "Trial 3", "Actual Root Age")

# Create the scatter plot for strict clock
plot(1:length(trial_means), trial_means, ylim=c(0, max(trial_means + trial_errors) + 1), 
     xaxt="n", xlab="Trials", ylab="Root Age", 
     main="Strict Clock Root Ages", pch=19)

# Add the x-axis labels
axis(1, at=1:length(trial_means), labels=trial_labels)

# Add error bars
arrows(x0 = 1:length(trial_means), y0 = trial_means - trial_errors, 
       x1 = 1:length(trial_means), y1 = trial_means + trial_errors, 
       angle = 90, code = 3, length = 0.1, col = "black")


#root ages mean#
# Create the scatter plot
plot(1:length(trial_means), trial_means, ylim=c(25, max(trial_means + trial_errors) + 1), 
     xaxt="n", xlab="Trials", ylab="Value", 
     main="Comparison of Trials", pch=19)

# Add the x-axis labels
axis(1, at=1:length(trial_means), labels=trial_labels)

# Add error bars
arrows(x0 = 1:length(trial_means), y0 = trial_means - trial_errors, 
       x1 = 1:length(trial_means), y1 = trial_means + trial_errors, 
       angle = 90, code = 3, length = 0.1, col = "black")


#####Checking branch lengths######

  ###box plot###

    #using relaxed clock
boxplot(branch1, branch2, branch3, real_branch, 
        names=c("Relaxed Trial1", "Relaxed Trial 2",
                "Relaxed Trial 3", "Simulated Tree"),
        ylab="Branch Lengths")
title("Relaxed Clock")

    #using strict clock
boxplot(strict.branch1, strict.branch2, strict.branch3, real_branch, 
        names=c("Strict Trial 1", "Strict Trial 2",
                "Strict Trial 3", "Simulated Tree"),
        ylab="Branch Lengths")
title("Strict Clock")



###exporting data as a table####
fin_dat<-data.frame(matrix(nrow=7, ncol=2))
rownames(fin_dat)<-c("Simulated Tree",
           "Relaxed Trial 1","Relaxed Trial 2","Relaxed Trial 3",
           "Strict Trial 1","Strict Trial 2","Strict Trial 3")
colnames(fin_dat)<-c("Root Age", "Branch Lengths Mean")
fin_dat[,1]<-c(round(real.node.dist, digits=2), paste(round(root_ages, digits=2), "+/-", round(root_errors, digits=2)))
fin_dat[,2]<-c(round(real.mean, digits=2) ,round(tree_dat$branch_length_mean, digits=2))

write.csv(fin_dat, file="output/fin.dat.csv")
write.csv(branch_lengths, file="output/branch_lengths.csv")  





