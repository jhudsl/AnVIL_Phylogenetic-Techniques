# (PART\*) ALIGNMENT {-}




# Why do we align sequences?

Sequence alignment is the art of lining up sequences from different samples in such a way that that reflects shared quality. When we perform an alignment in preparation for phylogenetic analyses, we aim to line up our sequences so that the complete alignment reflects the evolutionary relationships among all the samples. When we look at it, this mean that long stretches of the sequences should be fairly similar, with smaller regions of dissimilarity scattered throughout.

As the samples become more distantly related from each other, the regions of similarity will become smaller and the regions of dissimilarity will become larger. How the regions of dissimilarity are arranged can change, depending on our choices of assumptions.

All alignment programs will assign a "price" to each potential alignment the create, then return the least costly alignment as the final result. Most programs create potential alignments using an algorithm that assigns similarity scores to each pairwise comparison. The program then uses these scores to determine the final potential alignment. The algorithm also assigns penalties for alignments that include undesireable features. In general, an alignment algorithm can apply two major costs:

  - gap opening cost: we can apply a penalty for opening (or starting) any gap (indicating an insertion or deletion event)
  - gap extension cost: we can apply a penalty for making a gap longer


Alignments can be created using either the nucleotide sequence or the amino acid sequence. Amino acid sequences can be useful when dealing with more diverse samples where the nucleotide sequence includes lots of regions of dissimilarity and few regions of similarity. Because there are only 4 nucleotides, compared to 20 amino acids, amino acid sequence alignments tend to be less noisy than nucleotide sequence alignments. Amino acid sequences are also slower to change than nucleotide sequences due silent (or synonymous) nucleotide mutations that don't affect the amino acid sequence.
  
Keep in mind that any alignment we use is still just a hypothesis - it may be a well-supported hypothesis that represents our best knowledge, but it may still not be correct. We may never know what the "true" alignment is.





# Creating a multiple sequence alignment

There are several packages available now for performing alignments in R, include `ape` (which we have used before), `Biostrings`, and `mas`. Each of these have their own strengths, and all of them will produce a good alignment. There are also many non-R-based alignment programs. For this class, we will focus on how to use the `DECIPHER` package, which will allow us to create an alignment, as well as look at it.

If you would like to learn more about `DECIPHER`, you can find the manuals and specific examples (vignettes) [here](http://bioconductor.org/packages/release/bioc/html/DECIPHER.html), the the Documentation section.

## Installing the `DECIPHER` package

`DECIPHER` is available via Bioconductor, which means we can install it using the AnVIL::install command. After we install the package, we can then open it via the `library()` command. (If you are prompted to update dependent packages, choose "all".)


```r
#AnVIL::install('DECIPHER') #use this command to install the DECIPHER package

library(DECIPHER)
```

## Loading the fasta file

DECIPHER's commands use a fasta file as their input. We could load the fasta file directly into one of the commands, or we could load the fasta file into an object that we then pass to the commands. The second option gives us more flexibility and ultimately ends up saving time - we can load the file once and use it for multiple calculations. 

We are loading the fasta object we created and saved in our working directory last time. If the fasta file isn't in your working directory, that's okay - you can replace "grass.fasta" with the path to the directory in which it's saved, followed by the file name.


```r
fas <- "grass.fasta"
```


## Creating an alignment

DECIPHER lets us create alignments in a couple of ways, since we downloaded sequence from a protein-coding gene (_Glu-1_). We can align the sequences directly, or we can have the program translate the DNA sequences into amino acid sequences, then align the amino acid sequences. If the sequences you want to align are not protein-coding, then you can only align the sequences directly.

This part has a couple of steps. First, pass the fasta file into an object formatted as a DNAStringSet object using the command `readDNAStringSet`. 


```r
dna <- readDNAStringSet(fas)
dna
```

```
## DNAStringSet object of length 13:
##      width seq                                              names               
##  [1]  1521 ATGGCTAAGCGGCTGGTCCTCTT...ATTGTCGGCCAGCCAGTGATAG JX915632
##  [2]  1506 ATGGCTAAGCGGTTGGTCCTCTT...ATTGTCGGCCAGCCAGTGATAG EF105403.1
##  [3]  1506 ATGGCTAAGCGGCTGGTCCTCTT...ATTGTCGGCTAGCCAGTGATAG DQ073553.1
##  [4]  1749 ATGGCTAAGCGGTTGGTCCTCTT...ATTGTCGGCCAGCCAGTGATAG FJ481575.1
##  [5]  2100 ATGGCTAAGCGGTTGGTCCTCTT...ATTGTCGGCCAGCCAGTGATAG EF204545.1
##  ...   ... ...
##  [9]  2052 ATGGCTAAGCGACTGGTCCTCTT...ATTGTCGGCCAGCCAGTGATAG AY804128.1
## [10]  1650 ATGGCTAAGCGGTTGGTCCTCTT...ATTGTCGGCCAGCCAGTGATAG AY303125.2
## [11]  1665 TCATCACCCACAACACCGAGCAC...CGCATTGTCGGCCAGCCAGTGA KF887414.1
## [12]  2296 CAAAACTAGAGATCAATTCATTG...AAAACGGAAAGCTTCTCCATCC D82941.1
## [13]  2046 ATGGCTAAGCAGCTGGTCCTCTT...CTCCCTGTCGGCCAGCCAGTAG JX276655.1
```

(You might have noticed a warning message about an invalid character type. Sometimes this happens when you download data from GenBank, because people may have made a mistake in the nucleotide sequence when they submitted their samples. Since it is only one base (and the sequences themselves are at least 1500 bp long), so excluding it is unlikely to bias our results.)

Although it might look like the sequences are aligned in the DNAStringSet object, they really aren't. To do so requires another command.

First, let's align our sequences using just the nucleotides.


```r
DNA.no_trans <- AlignSeqs(dna)
```

```
## Determining distance matrix based on shared 10-mers:
## ================================================================================
## 
## Time difference of 0.01 secs
## 
## Clustering into groups by similarity:
## ================================================================================
## 
## Time difference of 0.01 secs
## 
## Aligning Sequences:
## ================================================================================
## 
## Time difference of 0.43 secs
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
## Time difference of 0.63 secs
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
## Time difference of 0.9 secs
## 
## Refining the alignment:
## ================================================================================
## 
## Time difference of 0.72 secs
```

We can also align our sequences after they are first translated. The translated amino acids are aligned, and then the sequences is reverse-translated back to nucleotides.


```r
DNA.trans <- AlignTranslation(dna)
```

```
## Determining distance matrix based on shared 5-mers:
## ================================================================================
## 
## Time difference of 0 secs
## 
## Clustering into groups by similarity:
## ================================================================================
## 
## Time difference of 0 secs
## 
## Aligning Sequences:
## ================================================================================
## 
## Time difference of 0.25 secs
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
## Time difference of 0.24 secs
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
## Time difference of 0.23 secs
```

Translating the nucleotide sequences sped up the alignment process, although both were fast enough it isn't a big deal to directly align the nucleotides. Your choice of which alignment procedure to use will largely come down to whether you are using coding sequence and how divergent your samples are. If you aren't using coding sequence, you will need to align using the nucleotides. If you have samples from deeply divergent species (especially if they come from different phyla), you will generally get a better alignment if you let the program translate your nucleotide sequence to amino acids first.

Now that we've created alignments, it would be helpful to visually check them. (This is possible if your sequences aren't too long, but can become really hard once you start dealing with very long stretches of DNA!) The `BrowseSeqs` command opens a browser window with the aligned sequences. Just remember which window belongs to each alignment!


```r
BrowseSeqs(DNA.no_trans)
BrowseSeqs(DNA.trans)
```

<img src="resources/images/05-alignment_1.png" title="Major point!! example image" alt="Major point!! example image" style="display: block; margin: auto;" />

You can scroll through the browser windows to see the full alignment. DECIPHER automatically color-codes the nucleotides, which makes it easier to pick out when a sequence doesn't match the others. At the very bottom, DECIPHER displays a consensus sequence, which we can also look at to identify which bases have mutations (or gaps).

The first thing that jumps out it how the two alignment methods resulted in alignments of different lengths. The alignment from the `AlignTranslation` method is longer than the alignment from the `AlignSeqs` method. How the gaps are inserted also differ quite a bit.

At this point, the alignment we choose to use is basically a judgment call. You have to decide which alignment seems "better" to you. With the grass sequences, sample AJ314771.1 doesn't quite seem to fit well with the `AlignTranslation` alignment. There are long stretches of the sequence that aren't aligned with anything else. This doesn't seem to be the case for the `AlignSeqs` alignment.

If we chose to go ahead with the `AlignTranslation` alignment, it would make sense to remove sample AJ314771.1 from our fasta file. You may come across a situation where one (or more) sequences don't seem to fit very well with the others. The best option going forward is to remove those sequences, as keeping them will cause issues when it comes time to infer the phylogeny.

The `alignSeqs` command allows us to change the gap opening and gap extension penalties. The default setting (which we used above, since we did not specify anything in the `alignSeqs` or `alignTranslation` commands) is to assign a minimum and maximum cost for each. For gap opening, the cost is -18 for nearly identical sequences and -14 for sequences that are more distantly related. For gap extension, those costs are -3 and -2. Let's see what happens when we change the costs.



```r
DNA.no_trans.1 <- AlignSeqs(dna, gapOpening = c(-20, -10), gapExtension = c(-5, -1))
```

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
## Time difference of 0.44 secs
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
## Time difference of 0.58 secs
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
## Time difference of 0.93 secs
## 
## Refining the alignment:
## ================================================================================
## 
## Time difference of 0.77 secs
```

```r
BrowseSeqs(DNA.no_trans.1)
```

Changing the gap opening and gap extension costs may or may not make a difference. In many cases, the default parameters will work just fine. However, it is always wise to check additional values!

For the grass samples, it would be reasonable to use the alignment generated by the `AlignSeqs` method. This is a judgment call more than anything else - there are large sections of the alignment that include all the sequences with minimal changes. We'll want to trim it to remove the ends with very little sequence coverage. We will keep everything between base 140 and 3373. To trim the sequences, we will write a new fasta file, then use it to create a new `DNAbin` object. If you haven't already loaded the `ape` library in your R session, you should do so now using the `library(ape)` command.

<img src="resources/images/05-alignment_2.png" title="Major point!! example image" alt="Major point!! example image" style="display: block; margin: auto;" />

<img src="resources/images/05-alignment_3.png" title="Major point!! example image" alt="Major point!! example image" style="display: block; margin: auto;" />



```r
writeXStringSet(DNA.no_trans, file="grass_aligned.fasta")

grass.align <- read.dna("grass_aligned.fasta", format="fasta", as.matrix=TRUE)
```

The trick for this step is that we specified the `grass.align` object be loaded as a matrix, not as a list. We can subset our sequences using brackets since `grass.align` is in matrix format.


```r
grass.trimmed <- grass.align[,140:3373]

grass.trimmed
```

```
## 13 DNA sequences in binary format stored in a matrix.
## 
## All sequences of same length: 3234 
## 
## Labels:
## JX915632
## EF105403.1
## DQ073553.1
## FJ481575.1
## EF204545.1
## AJ314771.1
## ...
## 
## Base composition:
##     a     c     g     t 
## 0.310 0.303 0.273 0.115 
## (Total: 42.04 kb)
```

Congratulations! You now have a trimmed alignment, ready for inferring trees. To be on the safe side, we will save it as a fasta file.


```r
write.dna( grass.align, file = 'grass_aligned.fasta', format = 'fasta' )
```

::: {.notice}
**R BASICS**

R has both data types and data structures. Data types (like character (for letters) or numeric (for real or decimal numbers)) can be combined to form data structures. Some of the data structures include:

  * vector
  * list
  * matrix
  * data frame

The data structures differ based on their size and what data types they can accept. A vector is a single dimension collection of one data type. You can have vectors made up of many character data types or many numerical data types, but you can't have a vector of both character and numerical data types together.

For example, a vector of numeric data might look like this:

1   1   2   3   5   

`DNAbin` objects are normally formatted as lists, which can be thought of as a multidimensional collection of vectors where each vector can be of a different type. (Within a vector everything is still a single character type.) In the code above, we specified the `DNAbin` object be formatted as a matrix. The matrix data structure is basically vectors in more than one dimensional, which means all the elements within a matrix must be the same data type (in this case, everything was a character).

A matrix of numeric data like the vector example above might look like this:

      [,1]    [,2]    [,3]    [,4]    [,5]
[1,]    1       1       2       3       5
[2,]    2       3       5       7       11

In a matrix (as well as in a data frame, which we haven't discussed), the position of an entry can be specificed using brackets. The first coordinate in the bracket indicates the row, while the second coordinate indicates the column. For the matrix above, [2,5] points to the number 11, which is the entry in the second row and fifth column. If we want to point to an entire row, we can leave off the second coordinate. [1,] points to the entire first row. If we want to point to an entire column, we have two choices. We could use either [,1] or just [1]. R will interpret both of those commands as pointing to column 1.

By making the `DNAbin` object a matrix, we could subset it like a matrix. That's what we did with the brackets. 

grass.align[,140:3373]

The comma at the beginning of the bracket told R that we wanted to keep all the rows in the matrix `grass.align`. The 140:3373 told R to keep all the columns between column 140 and column 3372 (as well as those two columns).

If we had chosen to remove one of our samples (let's say, the third sample), we could use the brackets to do so as well. The command for that would look like:

kept.sequences <- grass.align[c(1:2,4:13),]

The `c(1:2,4:13)` is telling R that we want to keep all the sequences between the 1st and 2nd rows, as well as all the sequences between the 4th and 13th row. We have to use the c() syntax because we have two sets of rows to keep.

:::



```r
sessionInfo()
```

```
## R version 4.0.2 (2020-06-22)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Ubuntu 20.04.5 LTS
## 
## Matrix products: default
## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/liblapack.so.3
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] stats4    parallel  stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
## [1] DECIPHER_2.18.1     RSQLite_2.2.1       Biostrings_2.58.0  
## [4] XVector_0.30.0      IRanges_2.24.1      S4Vectors_0.28.1   
## [7] BiocGenerics_0.36.1 BiocManager_1.30.10 ape_5.4-1          
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_1.0.10     bslib_0.4.2     compiler_4.0.2  pillar_1.4.6   
##  [5] jquerylib_0.1.4 zlibbioc_1.36.0 tools_4.0.2     bit_4.0.4      
##  [9] digest_0.6.25   memoise_2.0.1   jsonlite_1.7.1  evaluate_0.20  
## [13] lifecycle_1.0.3 tibble_3.0.3    nlme_3.1-149    lattice_0.20-41
## [17] pkgconfig_2.0.3 rlang_1.0.6     DBI_1.1.0       cli_3.6.0      
## [21] yaml_2.2.1      xfun_0.26       fastmap_1.1.1   stringr_1.4.0  
## [25] knitr_1.33      fs_1.5.0        vctrs_0.5.2     sass_0.4.5     
## [29] hms_0.5.3       bit64_4.0.5     grid_4.0.2      R6_2.4.1       
## [33] ottrpal_1.0.1   rmarkdown_2.10  bookdown_0.24   blob_1.2.1     
## [37] readr_1.4.0     magrittr_2.0.3  htmltools_0.5.4 ellipsis_0.3.1 
## [41] stringi_1.5.3   cachem_1.0.7    crayon_1.3.4
```

