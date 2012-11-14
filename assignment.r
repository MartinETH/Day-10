# Assignment 10 #
# Martin Gubler #

# 1. Consider the state sequence object biofam.seq and the matrix dOM of pairwise OM
# dissimilarities based on the properties matrix considered in the previous assignments.

### CODE COPY-PASTE FROM PREVIOS ASSIGNMENT / SAMPLE SOLUTION

library(TraMineR)
data(biofam)
biofam$cohort <- cut(biofam$birthyr, c(1900,1930,1940,1950,1960),
                     labels=c("1900-1929", "1930-1939", "1940-1949", "1950-1959"),
                     right=FALSE)
bf.states <- c("Parent", "Left", "Married", "Left/Married", "Child",
               "Left/Child", "Left/Married/Child", "Divorced")
bf.shortlab <- c("P","L","M","LM","C","LC", "LMC", "D")
biofam.seq <- seqdef(biofam[,10:25], states=bf.shortlab,
                     labels=bf.states, weights=biofam$wp00tbgs)

properties <- matrix(c(# left, married, child, divorced
  0, 0, 0, 0,  # parent
  1, 0, 0, 0,  # left
  0, 1, .5, 0, # marr
  1, 1, 0, 0,  # left+marr
  0, 0, 1, 0,  # child
  1, 0, 1, 0,  # left+child
  1, 1, 1, 0,  # left+marr+child
  .5, 1, .5, 1 # divorced
), 8, 4, byrow=TRUE)
sm <- as.matrix(dist(properties))
indel <- .5*max(sm)
dOM <- seqdist(biofam.seq, method="OM", indel=indel, sm=sm,
               full.matrix = FALSE)
weight <- attr(biofam.seq, "weight")

# 2. Create a cohort covariate, separating the individual born after the second World War
# from those born before.

biofam$WW2 <- cut(biofam$birthyr, c(1900, 1946,1960),
                     labels=c("Born before end of WW2","Born after WW2"),
                     right=FALSE)

# 3. We are interested in studying the rise of new kind of family trajectories. Plot the
# sequences according to the cohort factor.

seqdplot(biofam.seq, group=biofam$WW2, border = NA)


# 4. To test if the differences are significative, compute the association between the
# sequences and the cohort covariate using dissassoc.

set.seed(1)
WW2.assoc <- dissassoc(dOM, group = biofam$WW2, R = 5000, weights = weight, weight.permutation = "diss")
print(WW2.assoc)

# 5. Some have argued that one of the main change is related with the diversification
# of the family trajectories (this assumption is called destandardisation of the life
# course). To test this hypothesis, we can look at the difference of discrepancies over
# cohorts. What do you think? Are the differences significant?

  # Levene analysis is not significant


# 6. Explore the evolution of the association over time using seqdiff.

biofam.diff <- seqdiff(biofam.seq, group = biofam$WW2, seqdist_arg = list(method="OM", indel = 1, sm = "CONSTANT")
plot(biofam.diff, stat="Levene")
                   
                       
# 7. Build a sequences regression tree using the sex, cohort and plingu02
# (langue of interview) and plot it using seqtreedisplay.

biofam.seqtree <- seqtree(biofam.seq ~ sex + cohort + plingu02, data = biofam, R = 5000, diss = dOM, weight.permutation = "diss")
seqtreedisplay(biofam.seqtree, type = "d", border = NA, showtree = TRUE, showdepth = TRUE, filename = "biofam.seqtree.png")
                       
#### QUESTION:
#### This command (line 74) creates various error messages - but so does the command in the sample solution.
#### For example, it tells me to install www.graphviz.org (which I have now done) but it does not change anything.
#### What is wrong (I think I "copied" the code from slide 52)? 
#### Error message:        Der Befehl "dot" ist entweder falsch geschrieben oder konnte nicht gefunden werden.
####                       Error in disstreedisplayInternal(tree = tree, filename = filename, tmpdisstree = tmpdisstree,  : 
####                       You should install GraphViz to use this function: see http://www.graphviz.org
####                       In addition: Warning message: running command 'C:\Windows\system32\cmd.exe /c dot -Tpng -otmpseqtreebf47540566.png tmpseqtreebf47540566.dot' had status 1            
 