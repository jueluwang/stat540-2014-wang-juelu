Semimar01
========================================================


```r
prDat <- read.table("GSE4051_MINI.txt", header = TRUE, row.names = 1)
str(prDat)
```

```
## 'data.frame':	39 obs. of  6 variables:
##  $ sample    : int  20 21 22 23 16 17 6 24 25 26 ...
##  $ devStage  : Factor w/ 5 levels "4_weeks","E16",..: 2 2 2 2 2 2 2 4 4 4 ...
##  $ gType     : Factor w/ 2 levels "NrlKO","wt": 2 2 2 2 1 1 1 2 2 2 ...
##  $ crabHammer: num  10.22 10.02 9.64 9.65 8.58 ...
##  $ eggBomb   : num  7.46 6.89 6.72 6.53 6.47 ...
##  $ poisonFang: num  7.37 7.18 7.35 7.04 7.49 ...
```

```r
prDat <- read.table("GSE4051_MINI.txt")
str(prDat)
```

```
## 'data.frame':	39 obs. of  6 variables:
##  $ sample    : int  20 21 22 23 16 17 6 24 25 26 ...
##  $ devStage  : Factor w/ 5 levels "4_weeks","E16",..: 2 2 2 2 2 2 2 4 4 4 ...
##  $ gType     : Factor w/ 2 levels "NrlKO","wt": 2 2 2 2 1 1 1 2 2 2 ...
##  $ crabHammer: num  10.22 10.02 9.64 9.65 8.58 ...
##  $ eggBomb   : num  7.46 6.89 6.72 6.53 6.47 ...
##  $ poisonFang: num  7.37 7.18 7.35 7.04 7.49 ...
```

```r
`?`(read.table)
```

```
## starting httpd help server ... done
```

```r
nrow(prDat)
```

```
## [1] 39
```

```r
dim(prDat)
```

```
## [1] 39  6
```

```r
dimname(prDat)
```

```
## Error: could not find function "dimname"
```

```r
dimnames(prDat)
```

```
## [[1]]
##  [1] "Sample_20" "Sample_21" "Sample_22" "Sample_23" "Sample_16"
##  [6] "Sample_17" "Sample_6"  "Sample_24" "Sample_25" "Sample_26"
## [11] "Sample_27" "Sample_14" "Sample_3"  "Sample_5"  "Sample_8" 
## [16] "Sample_28" "Sample_29" "Sample_30" "Sample_31" "Sample_1" 
## [21] "Sample_10" "Sample_4"  "Sample_7"  "Sample_32" "Sample_33"
## [26] "Sample_34" "Sample_35" "Sample_13" "Sample_15" "Sample_18"
## [31] "Sample_19" "Sample_36" "Sample_37" "Sample_38" "Sample_39"
## [36] "Sample_11" "Sample_12" "Sample_2"  "Sample_9" 
## 
## [[2]]
## [1] "sample"     "devStage"   "gType"      "crabHammer" "eggBomb"   
## [6] "poisonFang"
```

```r
names(prDat)
```

```
## [1] "sample"     "devStage"   "gType"      "crabHammer" "eggBomb"   
## [6] "poisonFang"
```

```r
levels(prDat$devStage)
```

```
## [1] "4_weeks" "E16"     "P10"     "P2"      "P6"
```

```r
str(prDat$devStage)
```

```
##  Factor w/ 5 levels "4_weeks","E16",..: 2 2 2 2 2 2 2 4 4 4 ...
```

```r
summary(prDat$devStage)
```

```
## 4_weeks     E16     P10      P2      P6 
##       8       7       8       8       8
```

```r
prDat[sample(nrow(prDat), size = 6), ]
```

```
##           sample devStage gType crabHammer eggBomb poisonFang
## Sample_10     10       P6 NrlKO      9.544   6.347      7.252
## Sample_13     13      P10 NrlKO      9.838   7.228      7.459
## Sample_30     30       P6    wt      8.951   6.269      7.274
## Sample_39     39  4_weeks    wt     10.200   7.003      7.320
## Sample_18     18      P10 NrlKO     10.140   7.438      7.363
## Sample_16     16      E16 NrlKO      8.583   6.470      7.494
```

```r
all(sort(prDat$sample))
```

```
## [1] TRUE
```

```r
sort(prDat$sample)
```

```
##  [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
## [24] 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39
```

```r
seq_len(nrow(prDat))
```

```
##  [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
## [24] 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39
```

```r
sort(prDat$sample) == seq_len(nrow(prDat))
```

```
##  [1] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
## [15] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
## [29] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
```

```r
table(prDat$sample, prDat$devStage)
```

```
##     
##      4_weeks E16 P10 P2 P6
##   1        0   0   0  0  1
##   2        1   0   0  0  0
##   3        0   0   0  1  0
##   4        0   0   0  0  1
##   5        0   0   0  1  0
##   6        0   1   0  0  0
##   7        0   0   0  0  1
##   8        0   0   0  1  0
##   9        1   0   0  0  0
##   10       0   0   0  0  1
##   11       1   0   0  0  0
##   12       1   0   0  0  0
##   13       0   0   1  0  0
##   14       0   0   0  1  0
##   15       0   0   1  0  0
##   16       0   1   0  0  0
##   17       0   1   0  0  0
##   18       0   0   1  0  0
##   19       0   0   1  0  0
##   20       0   1   0  0  0
##   21       0   1   0  0  0
##   22       0   1   0  0  0
##   23       0   1   0  0  0
##   24       0   0   0  1  0
##   25       0   0   0  1  0
##   26       0   0   0  1  0
##   27       0   0   0  1  0
##   28       0   0   0  0  1
##   29       0   0   0  0  1
##   30       0   0   0  0  1
##   31       0   0   0  0  1
##   32       0   0   1  0  0
##   33       0   0   1  0  0
##   34       0   0   1  0  0
##   35       0   0   1  0  0
##   36       1   0   0  0  0
##   37       1   0   0  0  0
##   38       1   0   0  0  0
##   39       1   0   0  0  0
```

```r
table(prDat$devStage, prDat$gType)
```

```
##          
##           NrlKO wt
##   4_weeks     4  4
##   E16         3  4
##   P10         4  4
##   P2          4  4
##   P6          4  4
```

```r
with(prDat, table(devStage, gType))
```

```
##          gType
## devStage  NrlKO wt
##   4_weeks     4  4
##   E16         3  4
##   P10         4  4
##   P2          4  4
##   P6          4  4
```

```r
weeDat <- subset(prDat, poisonFang > 7.5)
nrow(weeDat)
```

```
## [1] 9
```

```r
table(weeDat$gType)
```

```
## 
## NrlKO    wt 
##     4     5
```

```r
addmargins(with(weeDat, table(devStage, gType)))
```

```
##          gType
## devStage  NrlKO wt Sum
##   4_weeks     0  0   0
##   E16         0  0   0
##   P10         2  2   4
##   P2          1  3   4
##   P6          1  0   1
##   Sum         4  5   9
```

```r
rownames(prDat)
```

```
##  [1] "Sample_20" "Sample_21" "Sample_22" "Sample_23" "Sample_16"
##  [6] "Sample_17" "Sample_6"  "Sample_24" "Sample_25" "Sample_26"
## [11] "Sample_27" "Sample_14" "Sample_3"  "Sample_5"  "Sample_8" 
## [16] "Sample_28" "Sample_29" "Sample_30" "Sample_31" "Sample_1" 
## [21] "Sample_10" "Sample_4"  "Sample_7"  "Sample_32" "Sample_33"
## [26] "Sample_34" "Sample_35" "Sample_13" "Sample_15" "Sample_18"
## [31] "Sample_19" "Sample_36" "Sample_37" "Sample_38" "Sample_39"
## [36] "Sample_11" "Sample_12" "Sample_2"  "Sample_9"
```

```r
rownames(prDat[prDat$eggBomb < quantile(prDat$eggBomb, 0.1), ])
```

```
## [1] "Sample_25" "Sample_14" "Sample_3"  "Sample_35"
```

```r

```


