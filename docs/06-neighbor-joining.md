# (PART\*) NEIGHBOR JOINING {-}



# The neighbor joining method of tree estimation

There are two basic types of tree estimation methods: distance methods (which includes neighbor joining, discussed here) and tree-searching methods (which we will cover in later sections). The first phylogenetic trees were built using distance methods. Strictly speaking, these methods can be considered phenetic methods, not phylogenetic methods, since these methods group taxa based on similarity instead of attempting to find the most likely tree. For all practical purposes, though, we can include neighbor joining in under the umbrella of phylogenetics.

In distance-based methods like neighbor joining, the genetic data is converted into a distance matrix, which is used to group taxa based on the genetic distance between them. Neighbor joining specifically works via star decomposition. Let's take a look at an example (adapted from (Fred Opperdoes)[http://www.deduveinstitute.be/~opperd/private/neighbor.html]).

We have 5 taxa (A,B,C,D, and E) for which we would like to infer a phylogeny. We've calculated the following distances among all the taxa:

<img src="resources/images/06-distance_1.png" title="Major point!! example image" alt="Major point!! example image" style="display: block; margin: auto;" />

The shortest distance between any taxa is 13, which means that A and B are most closely related to each other and will form our first group. We begin with all the taxa forming a star shape, and then we add a branch that connects A and B.

<img src="resources/images/06-nj_1.png" title="Major point!! example image" alt="Major point!! example image" style="display: block; margin: auto;" />

We keep grouping taxa until the entire star has been resolved into bifurcating branches (nodes that connect two, and only two, branches). For our example tree, we eventually end up with this tree:


<img src="resources/images/06-nj_2.png" title="Major point!! example image" alt="Major point!! example image" style="display: block; margin: auto;" />


The neighbor joining method is very fast, requires very few computational resources, and is statistically consistent. It will also always result in only a single tree topology. Unfortunately, many people will become overly confident that the single neighbor joining tree is the true tree, which is not necessarily true. Nowadays, most researchers might start their analysis with a neighbor joining tree, but they generally go on to infer trees using one of the other methods. 

::: {.dictionary}

There is another type of distance method you might hear about called UPGMA (unweighted pair group method with arithmetic mean). This is a hierarchical clustering method. However, the tree built using UPGMA depends greatly on the order in which the taxa are added, unlike the neighbor joining method (which will result in the same tree no matter what order the samples are added). As a result, the UPGMA method is almost never used in phylogenetics today.

:::


## Models of molecular evolution

In the example above, we kind of glossed over the process of converting DNA sequence into genetic distances. In its simplest form, every mutation is given a value, and the total value of all mutations is then used to calculate a genetic distance. 

But wait, you might be thinking. How do researchers decide what value to assign each mutation? Are all mutations given the same value, or do some changes "cost" more than others? Aren't some mutations considered more unlikely than others? These are all excellent questions that should be considered whenever dealing with phylogenetic analysis. Luckily, we don't have to answer all these questions ourselves. Instead, we can look at published _models of molecular evolution_ (also called substitution models) and decide which model fits our data the best.

### JC69

A model of molecular evolution is a set of rules that determines how much each mutation costs relative to all other mutations. One of the earliest published models is the JC69 model (sometimes written as JC), proposed by Jukes and Cantor in 1969. The JC69 model assumes equal base frequencies (that is, each nucleotide makes up 25% of the bases in the DNA sequence), as well as equal mutation rates. A mutation of an A to a T is the same as the mutation as a T to an A or an A to a C. As a result, the JC69 model has only 1 parameter - the mutation rate. 

### K2P

Many other models of molecular evolution have been published since JC69. One of the more popular ones is the Kimura two-parameter model (also called K80, or K2P). This model still assumes equal base frequencies, but allows for a distinct mutation rate for transitions (mutations between purines or pyrimidines) and a different mutation rate for tranversions (mutations from purine to pyrimidine, or vice versa). As you might have guessed from the name, this particular model has two parameters.

### GTR

One of the most complex models you might come across is the general time reversible model, or GTR. This model is extremely flexible but also parameter-heavy. There are no assumptions about base frequencies, and each possible mutation is given its own mutation rate. (A quick note: the mutation rate from A to G is the same as the mutation rate from G to A in this model.) As a result, this model has 9 parameters: 6 mutation rates + 3 base frequencies. (Because the sum of the 4 base frequencies must equal 1, we only have to describe 3 of them in order to calculate all 4. That's why we only need 3 parameters to model the base frequencies for GTR.)

### Other parameters

Sometimes you will run across a model that has "+ I" added to it (like GTR + I). The "I" refers to a invariant sites parameter. In these models, each base has a certain probability of being invariate, or unlikely to mutate. 

Another option you might see is using a gamma distribution to estimate how many times any given base in a sequence might mutate. We often think of mutation as a one-and-done process, where a nucleotide might mutate only once. This isn't a bad assumption; in fact, given that mutation is so rare, it is extremely unlikely for a nucleotide to mutate multiple times, though it does happen (particularly as divergence times between taxa increase). A solution to this dilemma is to include a parameter that describes how many times a particular nucleotide has undergone mutation. A gamma distribution is used because the highest probability is found in the smallest values of x, with a low probability of larger values of x (which would be multiple mutation hits at the same nucleotide). A model using a gamma distribution with have a "+ G".

Finally, some models with have a "+ SS" added to it. The SS indicates a site-specific variation parameter. These models will estimate different mutation rates for each codon position. Site-specific models are much less common than models that include invariate sites parameters or a \Gamma$ parameter.

### How do I choose a model?

There are enough substitution models out there that it can be a bit overwhelming to pick one for your phylogenetic analyses. (The models discussed above are in no way an exhaustive list.) Many researchers will have a favorite model they use, while others will rely model test programs to tell them which model of molecular evolution best fits their data. Although GTR + I + \Gamma$ might be the most biologically realistic model, it is also parameter-heavy and may not be the best choice when dealing with a sequence alignment that includes only a small number of phylogenetically-informative sites. 


```
## Determining distance matrix based on shared 10-mers:
## ================================================================================
## 
## Time difference of 0.01 secs
## 
## Clustering into groups by similarity:
## ================================================================================
## 
## Time difference of 0 secs
## 
## Aligning Sequences:
## ================================================================================
## 
## Time difference of 0.32 secs
## 
## Iteration 1 of 2:
## 
## Determining distance matrix based on alignment:
## ================================================================================
## 
## Time difference of 0 secs
## 
## Reclustering into groups by similarity:
## ================================================================================
## 
## Time difference of 0 secs
## 
## Realigning Sequences:
## ================================================================================
## 
## Time difference of 0.55 secs
## 
## Iteration 2 of 2:
## 
## Determining distance matrix based on alignment:
## ================================================================================
## 
## Time difference of 0 secs
## 
## Reclustering into groups by similarity:
## ================================================================================
## 
## Time difference of 0 secs
## 
## Realigning Sequences:
## ================================================================================
## 
## Time difference of 0.74 secs
## 
## Refining the alignment:
## ================================================================================
## 
## Time difference of 0.79 secs
```


# Building a neighbor joining tree

We are finally ready to start building trees from our data. 

For these analyses, we will use an R package called `phangorn`. If you'd like to learn more about it, you can find the manual [here](https://cran.r-project.org/web/packages/phangorn/phangorn.pdf).

`phangorn` (Phylogenetic Reconstruction and Analysis) is stored in the CRAN repository, so we will use `install.packages` for the installation.


``` r
install.packages('phangorn')

library(phangorn)
```

## The phyDat object

`phangorn` uses a data structure called phyDat to store information. You can either load the fasta file of your alignment into a phyDat object, or you can directly convert the DNAbin object of your alignment into a phyDat object. You will also want to make sure you have access to your DNAbin object, since we will need that for the model testing.


``` r
#Loading a fasta file

grass.phy <- read.phyDat('grass_aligned.fasta', format = 'fasta', type = 'DNA')

#Converting your DNAbin object

grass.align <- read.dna('grass_aligned.fasta', format = 'fasta')

grass.phy <- phyDat(grass.align)

grass.phy
```


```
## 13 sequences with 3585 character and 1119 different site patterns.
## The states are a c g t
```

The phyDat object is remarkably similar to the DNAbin object - this isn't too surprising, since the same team wrote both `ape` and `phangorn`.

::: {.notice}
**R BASICS**

If you are pausing your R environment on AnVIL (or have not closed an R session), you may still have all the objects you created previously still available. To check this, you can use the `list` command.

ls()

R will print all the objects that are still stored in the workspace. If you see your aligned DNAbin object (which can also bee seen in the Environment tab in the upper right box of RStudio), you can directly convert it to a `phyDat` object instead of loading the fasta file first.
:::

## Choosing a substitution model

`phangorn` can test our data against 24 different substitution models to determine which is the best fit. We first have to create a "guide" tree. (In future releases of `phangorn` creating this tree won't be necessary, but for now it doesn't take much extra time to create a tree.)


``` r
dist.matrix <- dist.dna(grass.align)

dist.matrix
```

```
##              JX915632 EF105403.1 DQ073553.1 FJ481575.1 EF204545.1 AJ314771.1
## EF105403.1 0.01221062                                                       
## DQ073553.1 0.01219682 0.00973376                                            
## FJ481575.1 0.04758637 0.05285651 0.05546604                                 
## EF204545.1 0.05963432 0.06233338 0.06498304 0.05007086                      
## AJ314771.1 0.08255861 0.08255861 0.08248480 0.09450600 0.09638394           
## FJ481569.1 0.06970732 0.07241003 0.07372738 0.08309250 0.09702098 0.11833578
## DQ073533.1 0.09773107 0.09630046 0.09764817 0.11151055 0.11461596 0.14259807
## AY804128.1 0.07268933 0.08095398 0.07950470 0.08464979 0.10731656 0.12910195
## AY303125.2 0.13551646 0.13398741 0.13386896 0.15627157 0.16453062 0.17201317
## KF887414.1 0.10204765 0.10349483 0.10340362 0.11297040 0.13258196 0.14259807
## D82941.1   0.21452588 0.22149463 0.21782031 0.22708927 0.24902239 0.25144294
## JX276655.1 0.20798301 0.21661880 0.21297096 0.22392751 0.23671033 0.24274355
##            FJ481569.1 DQ073533.1 AY804128.1 AY303125.2 KF887414.1   D82941.1
## EF105403.1                                                                  
## DQ073553.1                                                                  
## FJ481575.1                                                                  
## EF204545.1                                                                  
## AJ314771.1                                                                  
## FJ481569.1                                                                  
## DQ073533.1 0.11388614                                                       
## AY804128.1 0.08065836 0.13011175                                            
## AY303125.2 0.14786836 0.11333554 0.16513852                                 
## KF887414.1 0.10968116 0.11385297 0.13050488 0.11776689                      
## D82941.1   0.22628374 0.23078837 0.23717769 0.25882315 0.23332142           
## JX276655.1 0.21794455 0.24374554 0.24113161 0.25554132 0.22798329 0.11452388
```

``` r
str(dist.matrix)
```

```
##  'dist' num [1:78] 0.0122 0.0122 0.0476 0.0596 0.0826 ...
##  - attr(*, "Size")= int 13
##  - attr(*, "Labels")= chr [1:13] "JX915632" "EF105403.1" "DQ073553.1" "FJ481575.1" ...
##  - attr(*, "Upper")= logi FALSE
##  - attr(*, "Diag")= logi FALSE
##  - attr(*, "call")= language dist.dna(x = grass.align)
##  - attr(*, "method")= chr "K80"
```

The distance matrix stores all the distances in the lower half of the table to save space. This particular matrix used the K80 model (also known as the K2P model) to calculate the distances. When we use the `str(dist.matrix)` command, we can see the model used stored in the method attribute.

Next, we make a very quick neighbor joining tree to act as the guide tree. We then run the model test command using our phyDat object and the guide tree.


``` r
tree.guide <- NJ(dist.matrix)

mod.test <- modelTest(grass.phy, tree = tree.guide)
```

```
## Model        df  logLik   AIC      BIC
##           JC 23 -12485.95 25017.89 25160.14 
##         JC+I 24 -12450.63 24949.26 25097.69 
##      JC+G(4) 24 -12440.6 24929.19 25077.62 
##    JC+G(4)+I 25 -12440.6 24931.19 25085.8 
##          F81 26 -12282.14 24616.28 24777.08 
##        F81+I 27 -12246.41 24546.83 24713.81 
##     F81+G(4) 27 -12235.56 24525.12 24692.1 
##   F81+G(4)+I 28 -12235.56 24527.12 24700.29 
##          K80 24 -12259.35 24566.7 24715.13 
##        K80+I 25 -12215.84 24481.69 24636.3 
##     K80+G(4) 25 -12203.6 24457.2 24611.81 
##   K80+G(4)+I 26 -12203.6 24459.2 24620 
##          HKY 27 -12048.85 24151.7 24318.68 
##        HKY+I 28 -12005.72 24067.43 24240.6 
##     HKY+G(4) 28 -11992.71 24041.43 24214.59 
##   HKY+G(4)+I 29 -11992.71 24043.43 24222.78 
##         TrNe 25 -12253.31 24556.62 24711.23 
##       TrNe+I 26 -12210.56 24473.12 24633.91 
##    TrNe+G(4) 26 -12198.44 24448.89 24609.69 
##  TrNe+G(4)+I 27 -12198.44 24450.89 24617.87 
##          TrN 28 -12042.97 24141.94 24315.1 
##        TrN+I 29 -11997.85 24053.69 24233.05 
##     TrN+G(4) 29 -11984.41 24026.83 24206.18 
##   TrN+G(4)+I 30 -11984.41 24028.83 24214.36 
##         TPM1 25 -12258.47 24566.93 24721.54 
##       TPM1+I 26 -12215.02 24482.03 24642.83 
##    TPM1+G(4) 26 -12202.72 24457.45 24618.25 
##  TPM1+G(4)+I 27 -12202.73 24459.45 24626.43 
##          K81 25 -12258.47 24566.93 24721.54 
##        K81+I 26 -12215.02 24482.03 24642.83 
##     K81+G(4) 26 -12202.72 24457.45 24618.25 
##   K81+G(4)+I 27 -12202.73 24459.45 24626.43 
##        TPM1u 28 -12048.78 24153.55 24326.72 
##      TPM1u+I 29 -12005.59 24069.17 24248.52 
##   TPM1u+G(4) 29 -11992.57 24043.14 24222.49 
## TPM1u+G(4)+I 30 -11992.57 24045.14 24230.68 
##         TPM2 25 -12259 24567.99 24722.6 
##       TPM2+I 26 -12215.81 24483.62 24644.42 
##    TPM2+G(4) 26 -12203.6 24459.2 24619.99 
##  TPM2+G(4)+I 27 -12203.6 24461.2 24628.18 
##        TPM2u 28 -12046.75 24149.51 24322.67 
##      TPM2u+I 29 -12004.37 24066.74 24246.09 
##   TPM2u+G(4) 29 -11991.57 24041.14 24220.49 
## TPM2u+G(4)+I 30 -11991.57 24043.14 24228.67 
##         TPM3 25 -12258.58 24567.15 24721.76 
##       TPM3+I 26 -12215.04 24482.08 24642.88 
##    TPM3+G(4) 26 -12202.72 24457.45 24618.25 
##  TPM3+G(4)+I 27 -12202.73 24459.45 24626.43 
##        TPM3u 28 -12038.58 24133.16 24306.33 
##      TPM3u+I 29 -11991.92 24041.83 24221.18 
##   TPM3u+G(4) 29 -11977.59 24013.17 24192.52 
## TPM3u+G(4)+I 30 -11977.59 24015.17 24200.71 
##        TIM1e 26 -12252.42 24556.84 24717.64 
##      TIM1e+I 27 -12209.72 24473.44 24640.43 
##   TIM1e+G(4) 27 -12197.56 24449.12 24616.1 
## TIM1e+G(4)+I 28 -12197.56 24451.12 24624.29 
##         TIM1 29 -12042.9 24143.81 24323.16 
##       TIM1+I 30 -11997.73 24055.47 24241 
##    TIM1+G(4) 30 -11984.28 24028.57 24214.1 
##  TIM1+G(4)+I 31 -11984.28 24030.57 24222.29 
##        TIM2e 26 -12252.94 24557.88 24718.68 
##      TIM2e+I 27 -12210.53 24475.05 24642.03 
##   TIM2e+G(4) 27 -12198.44 24450.89 24617.87 
## TIM2e+G(4)+I 28 -12198.44 24452.89 24626.05 
##         TIM2 29 -12040.85 24139.7 24319.05 
##       TIM2+I 30 -11996.5 24053 24238.53 
##    TIM2+G(4) 30 -11983.33 24026.66 24212.2 
##  TIM2+G(4)+I 31 -11983.33 24028.66 24220.38 
##        TIM3e 26 -12252.45 24556.9 24717.69 
##      TIM3e+I 27 -12209.76 24473.53 24640.51 
##   TIM3e+G(4) 27 -12197.63 24449.26 24616.24 
## TIM3e+G(4)+I 28 -12197.63 24451.26 24624.43 
##         TIM3 29 -12031.72 24121.44 24300.79 
##       TIM3+I 30 -11981.48 24022.96 24208.5 
##    TIM3+G(4) 30 -11966.21 23992.42 24177.95 
##  TIM3+G(4)+I 31 -11966.21 23994.42 24186.14 
##         TVMe 27 -12257.01 24568.02 24735.01 
##       TVMe+I 28 -12214.03 24484.05 24657.22 
##    TVMe+G(4) 28 -12201.77 24459.54 24632.71 
##  TVMe+G(4)+I 29 -12201.77 24461.54 24640.89 
##          TVM 30 -12035.58 24131.15 24316.69 
##        TVM+I 31 -11989.78 24041.55 24233.27 
##     TVM+G(4) 31 -11975.69 24013.38 24205.1 
##   TVM+G(4)+I 32 -11975.69 24015.38 24213.29 
##          SYM 28 -12250.86 24557.72 24730.89 
##        SYM+I 29 -12208.75 24475.5 24654.85 
##     SYM+G(4) 29 -12196.66 24451.33 24630.68 
##   SYM+G(4)+I 30 -12196.66 24453.33 24638.86 
##          GTR 31 -12028.68 24119.37 24311.09 
##        GTR+I 32 -11979.34 24022.69 24220.59 
##     GTR+G(4) 32 -11964.44 23992.87 24190.78 
##   GTR+G(4)+I 33 -11964.44 23994.87 24198.96
```

As the model test runs, R prints out most of the models being tested. For some reason, R doesn't print the basic models. The first three models listed are JC + I, JC + G, and JC + G + I. While these models are being tested, the program is _also_ testing the basic JC model.


``` r
mod.test
```

```
##           Model df    logLik      AIC          AICw     AICc         AICcw
## 1            JC 23 -12485.95 25017.89 8.521162e-224 25018.20 9.676258e-224
## 2          JC+I 24 -12450.63 24949.26 6.831105e-209 24949.59 7.652882e-209
## 3       JC+G(4) 24 -12440.60 24929.19 1.556172e-204 24929.53 1.743378e-204
## 4     JC+G(4)+I 25 -12440.60 24931.19 5.723464e-205 24931.56 6.322241e-205
## 5           F81 26 -12282.14 24616.28 1.379092e-136 24616.67 1.501191e-136
## 6         F81+I 27 -12246.41 24546.83 1.665570e-121 24547.25 1.785617e-121
## 7      F81+G(4) 27 -12235.56 24525.12 8.600271e-117 24525.55 9.220141e-117
## 8    F81+G(4)+I 28 -12235.56 24527.12 3.163049e-117 24527.58 3.337844e-117
## 9           K80 24 -12259.35 24566.70 8.059939e-126 24567.03 9.029544e-126
## 10        K80+I 25 -12215.84 24481.69 2.320082e-107 24482.06 2.562804e-107
## 11     K80+G(4) 25 -12203.60 24457.20 4.825924e-102 24457.56 5.330802e-102
## 12   K80+G(4)+I 26 -12203.60 24459.20 1.774970e-102 24459.59 1.932118e-102
## 13          HKY 27 -12048.85 24151.70  1.051483e-35 24152.13  1.127269e-35
## 14        HKY+I 28 -12005.72 24067.43  2.092988e-17 24067.89  2.208650e-17
## 15     HKY+G(4) 28 -11992.71 24041.43  9.280424e-12 24041.88  9.793276e-12
## 16   HKY+G(4)+I 29 -11992.71 24043.43  3.413362e-12 24043.92  3.543472e-12
## 17         TrNe 25 -12253.31 24556.62 1.242675e-123 24556.99 1.372681e-123
## 18       TrNe+I 26 -12210.56 24473.12 1.687072e-105 24473.51 1.836438e-105
## 19    TrNe+G(4) 26 -12198.44 24448.89 3.076287e-100 24449.28 3.348648e-100
## 20  TrNe+G(4)+I 27 -12198.44 24450.89 1.131444e-100 24451.32 1.212994e-100
## 21          TrN 28 -12042.97 24141.94  1.386096e-33 24142.39  1.462694e-33
## 22        TrN+I 29 -11997.85 24053.69  2.011318e-14 24054.18  2.087985e-14
## 23     TrN+G(4) 29 -11984.41 24026.83  1.373994e-08 24027.32  1.426368e-08
## 24   TrN+G(4)+I 30 -11984.41 24028.83  5.053631e-09 24029.35  5.158084e-09
## 25         TPM1 25 -12258.47 24566.93 7.174516e-126 24567.30 7.925100e-126
## 26       TPM1+I 26 -12215.02 24482.03 1.956068e-107 24482.43 2.129250e-107
## 27    TPM1+G(4) 26 -12202.72 24457.45 4.257453e-102 24457.84 4.634389e-102
## 28  TPM1+G(4)+I 27 -12202.73 24459.45 1.565766e-102 24459.88 1.678619e-102
## 29          K81 25 -12258.47 24566.93 7.174516e-126 24567.30 7.925100e-126
## 30        K81+I 26 -12215.02 24482.03 1.956068e-107 24482.43 2.129250e-107
## 31     K81+G(4) 26 -12202.72 24457.45 4.257453e-102 24457.84 4.634389e-102
## 32   K81+G(4)+I 27 -12202.73 24459.45 1.565766e-102 24459.88 1.678619e-102
## 33        TPM1u 28 -12048.78 24153.55  4.166546e-36 24154.01  4.396796e-36
## 34      TPM1u+I 29 -12005.59 24069.17  8.755752e-18 24069.66  9.089504e-18
## 35   TPM1u+G(4) 29 -11992.57 24043.14  3.934829e-12 24043.63  4.084817e-12
## 36 TPM1u+G(4)+I 30 -11992.57 24045.14  1.447229e-12 24045.67  1.477142e-12
## 37         TPM2 25 -12259.00 24567.99 4.219763e-126 24568.36 4.661226e-126
## 38       TPM2+I 26 -12215.81 24483.62 8.822378e-108 24484.02 9.603473e-108
## 39    TPM2+G(4) 26 -12203.60 24459.20 1.776912e-102 24459.59 1.934232e-102
## 40  TPM2+G(4)+I 27 -12203.60 24461.20 6.535119e-103 24461.62 7.006141e-103
## 41        TPM2u 28 -12046.75 24149.51  3.148425e-35 24149.96  3.322412e-35
## 42      TPM2u+I 29 -12004.37 24066.74  2.954974e-17 24067.23  3.067612e-17
## 43   TPM2u+G(4) 29 -11991.57 24041.14  1.072128e-11 24041.63  1.112996e-11
## 44 TPM2u+G(4)+I 30 -11991.57 24043.14  3.942948e-12 24043.66  4.024444e-12
## 45         TPM3 25 -12258.58 24567.15 6.427868e-126 24567.52 7.100339e-126
## 46       TPM3+I 26 -12215.04 24482.08 1.910411e-107 24482.47 2.079550e-107
## 47    TPM3+G(4) 26 -12202.72 24457.45 4.257475e-102 24457.84 4.634413e-102
## 48  TPM3+G(4)+I 27 -12202.73 24459.45 1.565894e-102 24459.88 1.678757e-102
## 49        TPM3u 28 -12038.58 24133.16  1.113837e-31 24133.62  1.175389e-31
## 50      TPM3u+I 29 -11991.92 24041.83  7.573118e-12 24042.32  7.861790e-12
## 51   TPM3u+G(4) 29 -11977.59 24013.17  1.267532e-05 24013.66  1.315848e-05
## 52 TPM3u+G(4)+I 30 -11977.59 24015.17  4.661767e-06 24015.70  4.758120e-06
## 53        TIM1e 26 -12252.42 24556.84 1.112058e-123 24557.24 1.210515e-123
## 54      TIM1e+I 27 -12209.72 24473.44 1.431998e-105 24473.87 1.535210e-105
## 55   TIM1e+G(4) 27 -12197.56 24449.12 2.737179e-100 24449.55 2.934463e-100
## 56 TIM1e+G(4)+I 28 -12197.56 24451.12 1.006731e-100 24451.58 1.062365e-100
## 57         TIM1 29 -12042.90 24143.81  5.445800e-34 24144.30  5.653383e-34
## 58       TIM1+I 30 -11997.73 24055.47  8.282117e-15 24055.99  8.453299e-15
## 59    TIM1+G(4) 30 -11984.28 24028.57  5.750601e-09 24029.09  5.869459e-09
## 60  TIM1+G(4)+I 31 -11984.28 24030.57  2.115101e-09 24031.13  2.121317e-09
## 61        TIM2e 26 -12252.94 24557.88 6.614142e-124 24558.28 7.199730e-124
## 62      TIM2e+I 27 -12210.53 24475.05 6.407431e-106 24475.48 6.869250e-106
## 63   TIM2e+G(4) 27 -12198.44 24450.89 1.133467e-100 24451.31 1.215162e-100
## 64 TIM2e+G(4)+I 28 -12198.44 24452.89 4.168896e-101 24453.34 4.399276e-101
## 65         TIM2 29 -12040.85 24139.70  4.235978e-33 24140.19  4.397445e-33
## 66       TIM2+I 30 -11996.50 24053.00  2.852468e-14 24053.52  2.911425e-14
## 67    TIM2+G(4) 30 -11983.33 24026.66  1.492272e-08 24027.18  1.523116e-08
## 68  TIM2+G(4)+I 31 -11983.33 24028.66  5.488017e-09 24029.22  5.504145e-09
## 69        TIM3e 26 -12252.45 24556.90 1.083356e-123 24557.29 1.179272e-123
## 70      TIM3e+I 27 -12209.76 24473.53 1.372620e-105 24473.95 1.471552e-105
## 71   TIM3e+G(4) 27 -12197.63 24449.26 2.553916e-100 24449.69 2.737991e-100
## 72 TIM3e+G(4)+I 28 -12197.63 24451.26 9.393268e-101 24451.72 9.912355e-101
## 73         TIM3 29 -12031.72 24121.44  3.921543e-29 24121.93  4.071024e-29
## 74       TIM3+I 30 -11981.48 24022.96  9.472061e-08 24023.49  9.667838e-08
## 75    TIM3+G(4) 30 -11966.21 23992.42  4.067760e-01 23992.94  4.151836e-01
## 76  TIM3+G(4)+I 31 -11966.21 23994.42  1.496128e-01 23994.98  1.500525e-01
## 77         TVMe 27 -12257.01 24568.02 4.153847e-126 24568.45 4.453238e-126
## 78       TVMe+I 28 -12214.03 24484.05 7.122920e-108 24484.51 7.516544e-108
## 79    TVMe+G(4) 28 -12201.77 24459.54 1.494892e-102 24460.00 1.577503e-102
## 80  TVMe+G(4)+I 29 -12201.77 24461.54 5.498279e-103 24462.03 5.707863e-103
## 81          TVM 30 -12035.58 24131.15  3.048796e-31 24131.67  3.111811e-31
## 82        TVM+I 31 -11989.78 24041.55  8.712514e-12 24042.11  8.738119e-12
## 83     TVM+G(4) 31 -11975.69 24013.38  1.141540e-05 24013.94  1.144895e-05
## 84   TVM+G(4)+I 32 -11975.69 24015.38  4.198219e-06 24015.98  4.135045e-06
## 85          SYM 28 -12250.86 24557.72 7.173358e-124 24558.18 7.569770e-124
## 86        SYM+I 29 -12208.75 24475.50 5.134112e-106 24475.99 5.329814e-106
## 87     SYM+G(4) 29 -12196.66 24451.33 9.083953e-101 24451.82 9.430215e-101
## 88   SYM+G(4)+I 30 -12196.66 24453.33 3.341115e-101 24453.85 3.410172e-101
## 89          GTR 31 -12028.68 24119.37  1.102738e-28 24119.93  1.105978e-28
## 90        GTR+I 32 -11979.34 24022.69  1.088435e-07 24023.28  1.072056e-07
## 91     GTR+G(4) 32 -11964.44 23992.87  3.243035e-01 23993.47  3.194235e-01
## 92   GTR+G(4)+I 33 -11964.44 23994.87  1.192745e-01 23995.50  1.153067e-01
##         BIC
## 1  25160.14
## 2  25097.69
## 3  25077.62
## 4  25085.80
## 5  24777.08
## 6  24713.81
## 7  24692.10
## 8  24700.29
## 9  24715.13
## 10 24636.30
## 11 24611.81
## 12 24620.00
## 13 24318.68
## 14 24240.60
## 15 24214.59
## 16 24222.78
## 17 24711.23
## 18 24633.91
## 19 24609.69
## 20 24617.87
## 21 24315.10
## 22 24233.05
## 23 24206.18
## 24 24214.36
## 25 24721.54
## 26 24642.83
## 27 24618.25
## 28 24626.43
## 29 24721.54
## 30 24642.83
## 31 24618.25
## 32 24626.43
## 33 24326.72
## 34 24248.52
## 35 24222.49
## 36 24230.68
## 37 24722.60
## 38 24644.42
## 39 24619.99
## 40 24628.18
## 41 24322.67
## 42 24246.09
## 43 24220.49
## 44 24228.67
## 45 24721.76
## 46 24642.88
## 47 24618.25
## 48 24626.43
## 49 24306.33
## 50 24221.18
## 51 24192.52
## 52 24200.71
## 53 24717.64
## 54 24640.43
## 55 24616.10
## 56 24624.29
## 57 24323.16
## 58 24241.00
## 59 24214.10
## 60 24222.29
## 61 24718.68
## 62 24642.03
## 63 24617.87
## 64 24626.05
## 65 24319.05
## 66 24238.53
## 67 24212.20
## 68 24220.38
## 69 24717.69
## 70 24640.51
## 71 24616.24
## 72 24624.43
## 73 24300.79
## 74 24208.50
## 75 24177.95
## 76 24186.14
## 77 24735.01
## 78 24657.22
## 79 24632.71
## 80 24640.89
## 81 24316.69
## 82 24233.27
## 83 24205.10
## 84 24213.29
## 85 24730.89
## 86 24654.85
## 87 24630.68
## 88 24638.86
## 89 24311.09
## 90 24220.59
## 91 24190.78
## 92 24198.96
```

`phangorn` calculates multiple statistics that can be used to judge the fit of all the models. Which one you choose depends on what you prefer, although the field standard is generally to use the Akaike Information Criterion. This is a single number that combines how well the model fits the data (as determined by the log likelihood score) with a penalty for increasing model complexity. More complex models will almost always fit the data better, but in small datasets it can be difficult to get good estimates of every parameter when the model is complex. With AIC, we might find the GTR + G + I model fits the best, but we only have enough samples to estimates 2 parameters with any confidence. In this case, the AIC score might tell us to use the K80 (K2P) model.

Smaller numbers are better when it comes to AIC. This is also the case for AICc (which is a "second generation" or updated calculate of AIC) and BIC (the Bayesian Information Criterion, a version of AIC that includes a stronger penalty for additional parameters). In RStudio, you canclick on "mod.test" in the Environment tab from the pane in the upper right corner. That will open the results of our model test analysis in a separate window. If we click on "AIC" or "AICc", we can order the models by the AIC values. When we do this, we see the lowest AIC and AICc values are for the GTR + G model, which is what we will use going forward.

## Building the neighbor joining tree

We've already calculated a distance matrix above, but this time we will specify the model we want to use. 


``` r
dist.matrix <- dist.dna(grass.align, model = "GTR + G")

tree <- NJ(dist.matrix)
tree
```

But wait, why are we getting an error message? As it turns out, the `modelTest` command might test 24 different models, but the `dist.dna` command can't use all of those models. This is frustrating but also not uncommon when doing phylogenetic analyses. Unfortunately, not every program supports every substitution model.

`dist.dna` appears to support the following models: RAW, JC69, K80, F81, K81, F84, T92, TN93, GG95, LOGDET, BH87, PARALIN, N, TS, TV, INDEL, INDELBLOCK. We can check which of these models has the lowest AIC from our earlier model test. Surprisingly, it's the K80 model, which is what we used initially.


``` r
dist.matrix <- dist.dna(grass.align, model = "K80")

tree <- NJ(dist.matrix)
```

## Visualizing the neighbor joining tree

Let's take a look at the neighbor joining tree.


``` r
plot(tree, type = "unrooted")
```

<img src="resources/images/06-neighbor-joining_files/figure-html/unnamed-chunk-11-1.png" width="672" />

While we have indeed generated a tree, it's not really the easiest to interpret at the moment. It's hard to see what the relationships among taxa are because the tips are labeled with GenBank accession numbers instead of taxa names. Also, we really ought to declare an outgroup.

First let's change the tip labels. This is simply a matter of replacing each accession number with the taxon name. We do this by creating a vector of taxon names (matching the order of the accession numbers), then replacing the `tip.label` variable in our `tree` object. Because most of us are not plant experts, we'll use the common names for each sample, but it's also acceptable to use the scientific names.


``` r
tree$tip.label
```

```
##  [1] "JX915632"   "EF105403.1" "DQ073553.1" "FJ481575.1" "EF204545.1"
##  [6] "AJ314771.1" "FJ481569.1" "DQ073533.1" "AY804128.1" "AY303125.2"
## [11] "KF887414.1" "D82941.1"   "JX276655.1"
```

``` r
new.labels <- c('wheat', 'intermediate wheatgrass', 'mammoth wild rye', 'wheatgrass', 'tall wheatgrass', 'rye', 'Asiatic grass', 'crested wheatgrass', 'Tauschs goatgrass', 'medusahead rye', 'mosquito grass', 'barley_D-hordein', 'Siberian wild rye_D-hordein' )

tree$tip.label <- new.labels
```

Next, we declare our outgroup (in this case, the two D-hordein samples) and root our tree. (Remember, we did this in our first exercise using R.)


``` r
tree.root <- root(tree, outgroup = c('barley_D-hordein','Siberian wild rye_D-hordein'))

plot(tree.root, type = "phylogram", main = 'Neighbor Joining, Glu-1')
```

<img src="resources/images/06-neighbor-joining_files/figure-html/unnamed-chunk-13-1.png" width="672" />

Now we have a tree that we can begin to make sense of. To an non-botanist, it seems interesting the samples with the common name "wheatgrass" don't appear to be clustering together. Common names can be misleading about phylogenetic relationships!

## Saving your trees

You want to be sure to save the rooted tree and model test results to the persistent disk.


``` r
write.tree(tree.root, file = 'nj_grass.tre')
write.table(mod.test, file = 'grass_model_test', quote=F, sep='\t')
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
## [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] phangorn_2.11.1     DECIPHER_2.30.0     RSQLite_2.3.5      
##  [4] Biostrings_2.70.3   GenomeInfoDb_1.38.8 XVector_0.42.0     
##  [7] IRanges_2.36.0      S4Vectors_0.40.2    BiocGenerics_0.48.1
## [10] ape_5.7-1          
## 
## loaded via a namespace (and not attached):
##  [1] fastmatch_1.1-4         xfun_0.48               bslib_0.6.1            
##  [4] websocket_1.4.2         processx_3.8.3          lattice_0.21-9         
##  [7] tzdb_0.4.0              quadprog_1.5-8          vctrs_0.6.5            
## [10] tools_4.3.2             ps_1.7.6                bitops_1.0-9           
## [13] generics_0.1.3          tibble_3.2.1            fansi_1.0.6            
## [16] highr_0.11              blob_1.2.4              pkgconfig_2.0.3        
## [19] Matrix_1.6-1.1          lifecycle_1.0.4         GenomeInfoDbData_1.2.11
## [22] compiler_4.3.2          stringr_1.5.1           chromote_0.3.1         
## [25] janitor_2.2.0           codetools_0.2-19        snakecase_0.11.1       
## [28] htmltools_0.5.7         sass_0.4.8              RCurl_1.98-1.14        
## [31] yaml_2.3.8              crayon_1.5.2            later_1.3.2            
## [34] pillar_1.9.0            jquerylib_0.1.4         openssl_2.1.1          
## [37] cachem_1.0.8            nlme_3.1-164            webshot2_0.1.1         
## [40] tidyselect_1.2.0        ottrpal_1.3.0           digest_0.6.34          
## [43] stringi_1.8.3           dplyr_1.1.4             bookdown_0.41          
## [46] rprojroot_2.0.4         fastmap_1.1.1           grid_4.3.2             
## [49] cli_3.6.2               magrittr_2.0.3          utf8_1.2.4             
## [52] readr_2.1.5             promises_1.2.1          bit64_4.0.5            
## [55] lubridate_1.9.3         timechange_0.3.0        rmarkdown_2.25         
## [58] httr_1.4.7              igraph_2.0.2            bit_4.0.5              
## [61] askpass_1.2.0           hms_1.1.3               memoise_2.0.1          
## [64] evaluate_0.23           knitr_1.48              rlang_1.1.4            
## [67] Rcpp_1.0.12             DBI_1.2.2               glue_1.7.0             
## [70] xml2_1.3.6              jsonlite_1.8.8          R6_2.5.1               
## [73] zlibbioc_1.48.2
```

