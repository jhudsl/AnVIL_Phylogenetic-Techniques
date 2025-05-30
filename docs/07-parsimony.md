# (PART\*) PARSIMONY {-}



# The parsimony method of tree estimation

This week we begin examining the tree-searching methods of tree estimation. Specifically, we will focus on the parsimony method.

## Parsimony

Parsimony operates on the assumption that the tree containing the fewest sequence changes is the most accurate estimation of what really occurred. Essentially, parsimony assumes that characteristics shared by taxa are shared as the result of common ancestry. Sometimes this assumption is not valid, but often it is, particularly when there is not too much homoplasy (when two or more sequences or structures resemble each other but don't have a common ancestral origin) or variation in rates between sites. Often there will be multiple equally parsimonious trees, in which case you can either make a consensus tree (a combination of all the equally parsimoniously trees), or pick one tree that is illustrative of the point being discussed. If you choose the second, make sure you mention the existence of the other, equally parsimonious trees!

Parsimony is technically not a phylogenetic method because it does not apply a substitution model to the sequence data. A proper parsimony tree does not scale the branches to the number of changes between taxa. Instead, all the taxa are lined up and the branches are however long they need to be in order for the taxa to line up.

<img src="resources/images/07-cladograms.png" title="Major point!! example image" alt="Major point!! example image" style="display: block; margin: auto;" />


We will be scaling the branches in our parsimony trees to reflect the number of changes along each branch, so we will not be working with the parsimony method in its purest form.

The parsimony method is the backbone of the field known as "cladistics." It's largely fallen out of favor as phylogenetics is much more suited to analyzing genetic data. Parsimony still survives as a method of choice for researchers at the American Museum of Natural History and among paleontologists (or anyone who is inferring trees based on non-genetic data). The battle between cladistics proponents and phylogenetics proponents was so intense during the 1980s and 1990s my undergraduate advisor jokingly referred to it as the Clade Wars.


## Tree searching methods

At the heart of the parsimony approach is the idea of tree searching. The parsimony algorithm starts with a random tree (sometimes even a neighbor joining tree, or a tree the researcher supplies). The starting tree is scored, then the tree topology is rearranged in some way to create a new tree. The new tree is scored; if its score is better than the previous tree, the new tree is kept. The process then continues until an optimum tree score is found or it reaches a set end point.

When a tree has only a small number of taxa, an exhaustive search of all possible topologies is possible. However, the number of possible topologies grows exponentially as the number of taxa increases (in his _Phylogenetic Trees Made Easy_ manual, Barry G. Hall estimates that for a tree with just 10 taxa, there are more than 34 million possible rooted topologies). It rapidly becomes impractical to perform an exhaustive topology search to find the best tree. Instead, researchers have developed tree-searching algorithms to more quickly sort through all the possibilities.

### Branch and bound

The branch and bound method, while not exhaustive, can still be time-consuming when dealing with more than ~10 taxa. In this approach, an initial random tree of all the taxa is created and scored.Then, the algorithm starts with a tree containing three taxa. Next, a branch containing a fourth taxon is inserted in all possible locations and each possible tree is scored.  If none of the possible trees scores better than the initial random tree, this particular search is abandoned and a new starting tree of 3 taxa is created. If, however, any of these possible trees has a score smaller than the score of the initial random tree, the search continues and a fifth taxon is added. The process is repeated, with one change - if none of the trees with five taxa score low enough, then the search goes back to the four taxa tree and tries a different fifth taxon. This way, the initial work isn't wasted.

### Nearest Neighbor Interchange (NNI)

Nearest neighbor interchange algorithms take advantage of the fact that there are only three possible ways to connect four groups. A random starting tree is divided into four parts (or subtrees). These four subtrees are connected in all three possible ways. The most parsimonious arrangement is kept, and the new tree is then randomly subdivided into four new subtrees so the process can start again.

### Subtree Pruning and Regrafting (SPR)

Subtree pruning and regrafting has two major steps. First, a random subtree is removed from the original starting tree; then, the subtree is inserted elsewhere on the main tree to create a new node. This rearranged tree is scored. If if is more parsimonious than the starting tree, the new tree is kept, and another subtree is pruned and regrafted.

### Star decomposition

We actually saw the star decomposition method in the neighbor joining section. All taxa are connected together at a single node to form a star tree. Next, a pair of taxa are joined to form a new node and the tree is scored. This is repeated until all possible pairs have been evaluated. The best is kept and then another pair of taxa are joined, and the process continues.

# Building a parsimony tree




We will continue using the `phangorn` package for our parsimony analyses. We are also going to briefly use the `phytools` package to modify the fasta file. The manual for the `phytools` package is [here](https://cran.r-project.org/web/packages/phytools/phytools.pdf).



``` r
install.packages('phylotools')

library(phylotools)
library(phangorn)
```

## Editing the fasta file

Last week we edited the tip labels of our neighbor joining tree so that the tree displayed the common names of our grass samples instead of the GenBank accession numbers. However, we're going to be building multiple trees from the grass fasta file, so it makes sense to edit the fasta file directly so we don't have to manually change the tip labels every time. `Phytools` allows us to edit fasta files using the `rename.fasta` command. Before we can use it, though, we need to create a dataframe with two columns: the first column contains the old sample names (the GenBank accession numbers), while the second column contains the new sample names (the common names for each sample).

A dataframe in R is a table or a two-dimensional array-like structure in which each column contains values of one variable and each row contains one set of values from each column. We can create one by binding (or pasting) two vectors together. The vectors we are binding MUST have the information in the same order (that is, the order of GenBank accession numbers must be the same order as the common names).



``` r
old.labels <- c("JX915632", "EF105403.1", "DQ073553.1",
                "FJ481575.1", "EF204545.1","AJ314771.1",
                "FJ481569.1", "DQ073533.1", "AY804128.1",
                "AY303125.2", "KF887414.1", "D82941.1",
                "JX276655.1")


new.labels <- c('wheat', 'intermediate wheatgrass',
                'mammoth wild rye', 'wheatgrass', 
                'tall wheatgrass', 'rye', 'Asiatic grass', 
                'crested wheatgrass', 'Tauschs goatgrass', 
                'medusahead rye', 'mosquito grass',
                'barley_D-hordein', 'Siberian wild rye_D-hordein' )

sample.df <- cbind(old.labels, new.labels)
sample.df
```

```
##       old.labels   new.labels                   
##  [1,] "JX915632"   "wheat"                      
##  [2,] "EF105403.1" "intermediate wheatgrass"    
##  [3,] "DQ073553.1" "mammoth wild rye"           
##  [4,] "FJ481575.1" "wheatgrass"                 
##  [5,] "EF204545.1" "tall wheatgrass"            
##  [6,] "AJ314771.1" "rye"                        
##  [7,] "FJ481569.1" "Asiatic grass"              
##  [8,] "DQ073533.1" "crested wheatgrass"         
##  [9,] "AY804128.1" "Tauschs goatgrass"          
## [10,] "AY303125.2" "medusahead rye"             
## [11,] "KF887414.1" "mosquito grass"             
## [12,] "D82941.1"   "barley_D-hordein"           
## [13,] "JX276655.1" "Siberian wild rye_D-hordein"
```

In the dataframe `sample.df`, you can see that each row is a sample, with the accession number in column 1 and the common name in column 2.

We can now use the `rename.fasta` command. This command requires three arguments. The first is the fasta file that is being edited, the second is the dataframe that contains both the old and new names, and the third is the name for the newly edited fasta file. The new fasta file will be saved in your working directory (on AnVIL, this is the persistent disk).


``` r
rename.fasta('grass_aligned.fasta', sample.df, 'grass_aligned-renamed.fasta')
```

```
## grass_aligned-renamed.fasta has been saved to  /__w/AnVIL_Phylogenetic-Techniques/AnVIL_Phylogenetic-Techniques
```

## Inferring a parsimony tree

There are several options for inferring parsimony trees using `phangorn`. We can use the branch and bound tree rearrangement method with the `bab` command, or we can apply the nearest neighbor interchange (NNI) or subtree prunning and regrafting (SPR) approaches with either the `optim.parsimony` or `pratchet` (parsimony ratchet) commands. The branch and bound method can be very, very slow if you have more ~10 samples. Let's take a look at the parsimony tree that results from a branch-and-bound search of the grass fasta file.


``` r
g.align <- read.phyDat("grass_aligned-renamed.fasta", format='fasta')
bab.tree <- bab(g.align)

bab.root <- root(bab.tree, 
                 outgroup = c('barley_D-hordein','Siberian wild rye_D-hordein'))
plot(bab.tree, main = 'Parsimony, BAB')
```

<img src="resources/images/07-parsimony_files/figure-html/unnamed-chunk-6-1.png" width="672" />

To test both NNI and SPR approaches for parsimony trees, we will use the `pratchet` command. This particular command has lots of possible parameters we can change.

  - maxit: maximum number of iterations
  
  - minit: minimum number of iterations
  
  - _k_: number of rounds of no improvement after which the ratchet is stopped
  
  - rearrangements: SPR or NNI rearrangements
  
  - trace: how much information to print after each iteration
  
  - all: whether to return all equally parsimonius trees, or just one

The parsimony ratchet is an approach that was developed by KC Nixon in 1999. This approach is a more efficient way to find better trees than by just NNI or SPR rearrangements alone. `phangorn` implements it this way:

1. Create a bootstrapped dataset from the original dataset (more on bootstrapping later).

2. Take the current best tree (or starting tree, for the first go-round) and perform tree rearrangements using the bootstrapped dataset, saving the best bootstrap tree.

3. Use the best bootstrap tree and perform tree rearrangements using the original dataset. Compare the parsimony score (the smallest number of changes necessary to describe the data for a given tree) of the bootstrapped tree and the original best tree. Whichever tree is best (ie, has the lowest parsimony score) is then saved as the new "best" tree.

This process is repeated until either the algorithm reaches the max number of iterations or the _k_ number is reached.


``` r
tree.SPR <- pratchet(g.align, maxit = 10000, minit = 100, k = 10,
                           all = T, rearrangements = 'SPR', trace = 0)

tree.NNI <- pratchet(g.align, maxit = 10000, minit = 100, k = 10,
                           all = T, rearrangements = 'NNI', trace = 0)

#We assign branch lengths to the tree equal to the number of changes

tree.SPR <- acctran(tree.SPR, g.align)
tree.NNI <- acctran(tree.NNI, g.align)

SPR.root <- root(tree.SPR, 
                      outgroup = c('barley_D-hordein','Siberian wild rye_D-hordein'))
NNI.root <- root(tree.NNI, 
                      outgroup = c('barley_D-hordein','Siberian wild rye_D-hordein'))

plot(SPR.root, main = 'Parsimony, SPR')
```


``` r
plot(NNI.root, main = 'Parsimony, NNI')
```

Take a look at the topology of the trees created using the branch-and-bound (BAB), SPR, and NNI approaches. Do you see any differences?


<img src="resources/images/07-parsimony_trees.png" title="Major point!! example image" alt="Major point!! example image" style="display: block; margin: auto;" />


## Parsimony scores

We can compare parsimony trees directly by comparing their parsimony scores. As mentioned above, the parsimony score is just the minimum number of changes necessary to map a dataset onto a particular tree topology. A smaller number is better.


``` r
tree.SPR <- pratchet(g.align, maxit = 10000, minit = 100, k = 10,
                           all = T, rearrangements = 'SPR', trace = 0)

tree.NNI <- pratchet(g.align, maxit = 10000, minit = 100, k = 10,
                           all = T, rearrangements = 'NNI', trace = 0)

#We assign branch lengths to the tree equal to the number of changes

tree.SPR <- acctran(tree.SPR, g.align)
tree.NNI <- acctran(tree.NNI, g.align)

SPR.root <- root(tree.SPR, 
                      outgroup = c('barley_D-hordein','Siberian wild rye_D-hordein'))
NNI.root <- root(tree.NNI, 
                      outgroup = c('barley_D-hordein','Siberian wild rye_D-hordein'))
```


``` r
#for each parsimony command, we need to provide a tree and a phyDat object
parsimony(bab.root, g.align)
```

```
## [1] 1624
```

``` r
parsimony(SPR.root, g.align)
```

```
## [1] 1624
```

``` r
parsimony(NNI.root, g.align)
```

```
## [1] 1624
```

In the case of our grass dataset, all the parsimony scores are the same, so all of the methods performed equally well. We got lucky - this is not always true!

## Bootstrapping

It would be nice to have a way to measure the robustness of each individual branch and clade, not just the overall tree. In phylogenetics, the measurement of choice for determining the "strength" of a branch is the bootstrap value. 

Bootstrapping is a resampling method that attempts to use the original data as a way to judge the strength of phylogenetic inference. For each bootstrap replicate, the program will randomly select a number of sites to create a pseudoalignment. (For example, if the original alignment has 100 bases, the bootstrap algorithm will randomly choose 100 bases to create a new alignment. Some bases might be chosen multiple times, while other bases don't get sampled.) After generating the pseudoalignment, a new tree is built and the relationships are stored. Here is a nice illustration of the process by [www.analyticsvidhya.com](https://www.analyticsvidhya.com/blog/2020/02/what-is-bootstrap-sampling-in-statistics-and-machine-learning/)

<img src="resources/images/07-bootstrap.png" title="Major point!! example image" alt="Major point!! example image" style="display: block; margin: auto;" />



After all the bootstrapping replicates have been analyzed, each branch of the actual tree will be labeled with a value that reflects how frequently that branch was seen among the replicate trees. Higher values mean the branch has more support. In general, researchers consider a branch with a bootstrap value > 0.5 as well-supported.

Luckily, the `pratchet` command automatically does bootstrapping, so we already have the bootstrap information saved. We can see the bootstrap values on the trees using the `plotBS` command.


``` r
plotBS(SPR.root, type = "p", main = 'Parsimony, SPR')
```


``` r
plotBS(NNI.root, type = "p", main = 'Parsimony, NNI')
```


<img src="resources/images/07-bootstrap_trees.png" title="Major point!! example image" alt="Major point!! example image" style="display: block; margin: auto;" />


You may need to resize the Plots window in order to more easily read the bootstrap values for each branch.


::: {.reflection}
can we bootstrap our neighbor joining tree?

In short, yes, although we have to do a bit of R trickery, because neither `ape` nor `phangorn` have the option to bootstrap neighbor joining trees built in. We can get around this by writing our own function.

Basically, we will write a function that creates a distance matrix using the `dist.ml` command. (This command is a more general form of the `dist.dna` command. As of this writing, using `dist.dna` in the function throws an error.) Then we can include both a phyDat object and our function in the `bootstrap.phyDat` command. R will take the phyDat object and pass it to our function, then use the distance matrix in a bootstrap analysis.


``` r
fun <- function(x) nj(dist.ml(x))

bs_nj <- bootstrap.phyDat(g.align,  fun)

dna_dist <- dist.ml(g.align, model = 'JC')
nj.tre <- nj(dna_dist)
nj.root <- root(nj.tre, 
                outgroup = c('barley_D-hordein','Siberian wild rye_D-hordein'))

plotBS(nj.root, bs_nj, 'p', main = 'Neighbor Joining')
```

<img src="resources/images/07-parsimony_files/figure-html/unnamed-chunk-12-1.png" width="672" />

This isn't a perfect approach - for one thing, we can't use the K80 substitution model for the `dist.ml` command. Instead, we simply used the JC model. Also, the bootstrapping values in this approach are from 1 to 100, not 0 to 1 as calculated by the `pratchet` command.

How does the bootstrap support for the neighbor joining tree compare to the bootstrap support for the parsimony trees?

::: 


## Saving files

Make sure to save at least one of the trees you built using the `pratchet` command to the persistent disk.


``` r
write.tree(SPR.root, file = 'spr_grass.tre')
write.tree(NNI.root, file = 'nni_grass.tre')
```





``` r
sessionInfo()
```

```
## R version 4.3.2 (2023-10-31)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Ubuntu 22.04.4 LTS
## 
## Matrix products: default
## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.20.so;  LAPACK version 3.10.0
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## time zone: Etc/UTC
## tzcode source: system (glibc)
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] phangorn_2.11.1  phylotools_0.2.2 ape_5.7-1       
## 
## loaded via a namespace (and not attached):
##  [1] sass_0.4.8       utf8_1.2.4       generics_0.1.3   xml2_1.3.6      
##  [5] lattice_0.21-9   stringi_1.8.3    hms_1.1.3        digest_0.6.34   
##  [9] magrittr_2.0.3   grid_4.3.2       evaluate_0.23    timechange_0.3.0
## [13] bookdown_0.41    fastmap_1.1.1    Matrix_1.6-1.1   rprojroot_2.0.4 
## [17] jsonlite_1.8.8   processx_3.8.3   chromote_0.3.1   ps_1.7.6        
## [21] promises_1.2.1   httr_1.4.7       fansi_1.0.6      ottrpal_1.3.0   
## [25] codetools_0.2-19 jquerylib_0.1.4  cli_3.6.2        rlang_1.1.4     
## [29] cachem_1.0.8     yaml_2.3.8       parallel_4.3.2   tools_4.3.2     
## [33] tzdb_0.4.0       dplyr_1.1.4      fastmatch_1.1-4  vctrs_0.6.5     
## [37] R6_2.5.1         lifecycle_1.0.4  lubridate_1.9.3  snakecase_0.11.1
## [41] stringr_1.5.1    janitor_2.2.0    pkgconfig_2.0.3  pillar_1.9.0    
## [45] bslib_0.6.1      later_1.3.2      glue_1.7.0       Rcpp_1.0.12     
## [49] highr_0.11       xfun_0.48        tibble_3.2.1     tidyselect_1.2.0
## [53] knitr_1.48       igraph_2.0.2     nlme_3.1-164     htmltools_0.5.7 
## [57] websocket_1.4.2  rmarkdown_2.25   webshot2_0.1.1   readr_2.1.5     
## [61] compiler_4.3.2   quadprog_1.5-8   askpass_1.2.0    openssl_2.1.1
```

