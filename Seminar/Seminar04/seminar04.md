Seminar04
========================================================


```r
library(lattice)
prDat <- read.table("GSE4051_data.tsv")
getwd()
```

```
## [1] "F:/STUDY/STAT540/Seminar04"
```

```r
str(prDat, max.level = 0)
```

```
## 'data.frame':	29949 obs. of  39 variables:
```

```r
prDes <- readRDS("GSE4051_design.rds")
str(prDes)
```

```
## 'data.frame':	39 obs. of  4 variables:
##  $ sidChar : chr  "Sample_20" "Sample_21" "Sample_22" "Sample_23" ...
##  $ sidNum  : num  20 21 22 23 16 17 6 24 25 26 ...
##  $ devStage: Factor w/ 5 levels "E16","P2","P6",..: 1 1 1 1 1 1 1 2 2 2 ...
##  $ gType   : Factor w/ 2 levels "wt","NrlKO": 1 1 1 1 2 2 2 1 1 1 ...
```

```r
set.seed(987)
theGene <- sample(1:nrow(prDat), size = 1)
theGene
```

```
## [1] 14294
```

```r
pDat <- data.frame(prDes, gExp = unlist(prDat[theGene, ]))
str(pDat)
```

```
## 'data.frame':	39 obs. of  5 variables:
##  $ sidChar : chr  "Sample_20" "Sample_21" "Sample_22" "Sample_23" ...
##  $ sidNum  : num  20 21 22 23 16 17 6 24 25 26 ...
##  $ devStage: Factor w/ 5 levels "E16","P2","P6",..: 1 1 1 1 1 1 1 2 2 2 ...
##  $ gType   : Factor w/ 2 levels "wt","NrlKO": 1 1 1 1 2 2 2 1 1 1 ...
##  $ gExp    : num  9.88 10.59 10.28 10.22 8.75 ...
```

```r
aggregate(gExp ~ gType, pDat, FUN = mean)
```

```
##   gType  gExp
## 1    wt 9.758
## 2 NrlKO 9.553
```

```r
library(lattice)
stripplot(gType ~ gExp, pDat)
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-11.png) 

```r
t.test(gExp ~ gType, pDat)
```

```
## 
## 	Welch Two Sample t-test
## 
## data:  gExp by gType
## t = 1.479, df = 36.78, p-value = 0.1475
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  -0.07586  0.48605
## sample estimates:
##    mean in group wt mean in group NrlKO 
##               9.758               9.553
```

```r
ttRes <- t.test(gExp ~ gType, pDat)
str(ttRes)
```

```
## List of 9
##  $ statistic  : Named num 1.48
##   ..- attr(*, "names")= chr "t"
##  $ parameter  : Named num 36.8
##   ..- attr(*, "names")= chr "df"
##  $ p.value    : num 0.148
##  $ conf.int   : atomic [1:2] -0.0759 0.4861
##   ..- attr(*, "conf.level")= num 0.95
##  $ estimate   : Named num [1:2] 9.76 9.55
##   ..- attr(*, "names")= chr [1:2] "mean in group wt" "mean in group NrlKO"
##  $ null.value : Named num 0
##   ..- attr(*, "names")= chr "difference in means"
##  $ alternative: chr "two.sided"
##  $ method     : chr "Welch Two Sample t-test"
##  $ data.name  : chr "gExp by gType"
##  - attr(*, "class")= chr "htest"
```

```r
ttRes$statistic
```

```
##     t 
## 1.479
```

```r
ttRes$p.value
```

```
## [1] 0.1475
```

```r
wilcox.test(gType, pDat)
```

```
## Error: object 'gType' not found
```

```r
ks.test(gExp, gType, pDat)
```

```
## Error: object 'gExp' not found
```

```r
str(pDat)
```

```
## 'data.frame':	39 obs. of  5 variables:
##  $ sidChar : chr  "Sample_20" "Sample_21" "Sample_22" "Sample_23" ...
##  $ sidNum  : num  20 21 22 23 16 17 6 24 25 26 ...
##  $ devStage: Factor w/ 5 levels "E16","P2","P6",..: 1 1 1 1 1 1 1 2 2 2 ...
##  $ gType   : Factor w/ 2 levels "wt","NrlKO": 1 1 1 1 2 2 2 1 1 1 ...
##  $ gExp    : num  9.88 10.59 10.28 10.22 8.75 ...
```

```r
x <- subset(pDat, gType = wt)
str(x)
```

```
## 'data.frame':	39 obs. of  5 variables:
##  $ sidChar : chr  "Sample_20" "Sample_21" "Sample_22" "Sample_23" ...
##  $ sidNum  : num  20 21 22 23 16 17 6 24 25 26 ...
##  $ devStage: Factor w/ 5 levels "E16","P2","P6",..: 1 1 1 1 1 1 1 2 2 2 ...
##  $ gType   : Factor w/ 2 levels "wt","NrlKO": 1 1 1 1 2 2 2 1 1 1 ...
##  $ gExp    : num  9.88 10.59 10.28 10.22 8.75 ...
```

```r
wilcox.test(wt, Nrlko)
```

```
## Error: object 'wt' not found
```

```r
getwd()
```

```
## [1] "F:/STUDY/STAT540/Seminar04"
```

```r
kDat <- readRDS("gse4051_MINI.rds")
kMat <- as.matrix(kDat[c("crabHammer", "eggBomb", "poisonFang")])
str(kMat)
```

```
##  num [1:39, 1:3] 10.22 10.02 9.64 9.65 8.58 ...
##  - attr(*, "dimnames")=List of 2
##   ..$ : chr [1:39] "12" "13" "14" "15" ...
##   ..$ : chr [1:3] "crabHammer" "eggBomb" "poisonFang"
```

```r
median(kMat[, 1])
```

```
## [1] 9.611
```

```r
median(kMat[, "eggBomb"])
```

```
## [1] 6.757
```

```r
apply(kMat, 2, median)
```

```
## crabHammer    eggBomb poisonFang 
##      9.611      6.757      7.350
```

```r
apply(kMat, 2, quantile, probs = 0.5)
```

```
## crabHammer    eggBomb poisonFang 
##      9.611      6.757      7.350
```

```r
apply(kMat, 1, min)
```

```
##    12    13    14    15     9    10    11    28    29    30    31    24 
## 7.370 6.890 6.720 6.529 6.470 7.005 6.735 6.587 6.170 6.870 6.800 6.138 
##    25    26    27    36    37    38    39    32    33    34    35    20 
## 6.166 6.269 6.264 6.530 7.100 6.269 6.211 6.286 6.347 6.270 6.188 7.005 
##    21    22    23    16    17    18    19     5     6     7     8     1 
## 7.082 6.757 6.155 7.228 7.226 7.363 7.081 6.993 6.992 6.608 7.003 6.981 
##     2     3     4 
## 7.165 7.075 6.558
```

```r
colnames(kMat)[apply(kMat, 1, which.min)]
```

```
##  [1] "poisonFang" "eggBomb"    "eggBomb"    "eggBomb"    "eggBomb"   
##  [6] "poisonFang" "poisonFang" "eggBomb"    "eggBomb"    "eggBomb"   
## [11] "eggBomb"    "eggBomb"    "eggBomb"    "eggBomb"    "eggBomb"   
## [16] "eggBomb"    "poisonFang" "eggBomb"    "eggBomb"    "eggBomb"   
## [21] "eggBomb"    "eggBomb"    "eggBomb"    "poisonFang" "eggBomb"   
## [26] "eggBomb"    "eggBomb"    "eggBomb"    "eggBomb"    "poisonFang"
## [31] "eggBomb"    "poisonFang" "eggBomb"    "eggBomb"    "eggBomb"   
## [36] "poisonFang" "eggBomb"    "poisonFang" "eggBomb"
```

```r
rowSums(kMat)
```

```
##    12    13    14    15     9    10    11    28    29    30    31    24 
## 25.05 24.09 23.71 23.22 22.55 24.21 24.09 22.96 22.66 23.99 23.26 22.96 
##    25    26    27    36    37    38    39    32    33    34    35    20 
## 22.78 22.60 23.40 22.17 24.80 22.49 22.31 22.58 23.14 22.60 22.75 25.43 
##    21    22    23    16    17    18    19     5     6     7     8     1 
## 24.17 23.86 21.80 24.52 24.76 24.94 24.44 24.82 23.98 23.70 24.52 23.86 
##     2     3     4 
## 23.64 23.93 23.42
```

```r
jRowSums <- rowSums(prDat)
jRowsums
```

```
## Error: object 'jRowsums' not found
```

```r
jRowsums <- apply(prDat, 1, sum)
keepgenes <- c("1431708_a_at", "1424336_at", "1454696_at", "1416119_at", "1432141_x_at", 
    "1429226_at")
miniDat <- subset(prDat, rownames(prDat) %in% keepgenes)
str(miniDat)
```

```
## 'data.frame':	6 obs. of  39 variables:
##  $ Sample_20: num  10.58 7.26 6.4 9.95 6.88 ...
##  $ Sample_21: num  11 7.96 6.65 10.1 6.42 ...
##  $ Sample_22: num  10.85 7.62 6.58 9.83 6.94 ...
##  $ Sample_23: num  10.92 7.93 6.46 9.98 6.5 ...
##  $ Sample_16: num  9.2 7.11 6.92 7.73 7.76 ...
##  $ Sample_17: num  11.01 6.49 6.51 6.85 6.78 ...
##  $ Sample_6 : num  10.9 6.46 6.28 6.89 6.47 ...
##  $ Sample_24: num  10.38 7.49 6.62 9.8 7.31 ...
##  $ Sample_25: num  10.61 7.59 6.58 8.99 7.2 ...
##  $ Sample_26: num  10.25 7.44 6.62 8.99 7.41 ...
##  $ Sample_27: num  9.74 7.33 6.73 9 7.46 ...
##  $ Sample_14: num  10.66 6.93 6.59 7.54 7.05 ...
##  $ Sample_3 : num  10.62 7.01 6.59 7.4 7.36 ...
##  $ Sample_5 : num  10.2 6.65 6.76 7.5 7.31 ...
##  $ Sample_8 : num  9.44 7.14 6.84 7.64 7.43 ...
##  $ Sample_28: num  8.82 7.62 6.76 8.57 7.62 ...
##  $ Sample_29: num  10.32 7.41 6.28 10.64 6.65 ...
##  $ Sample_30: num  10.5 7.74 6.71 9.28 7.15 ...
##  $ Sample_31: num  9.45 7.83 6.74 8.58 7.35 ...
##  $ Sample_1 : num  10.08 6.88 6.63 7.19 7.14 ...
##  $ Sample_10: num  10.67 6.62 6.36 7.27 7.01 ...
##  $ Sample_4 : num  10.16 6.8 6.4 7.53 7.08 ...
##  $ Sample_7 : num  9.11 7.26 6.59 7.75 7.42 ...
##  $ Sample_32: num  10.21 7.8 6.38 10.71 6.68 ...
##  $ Sample_33: num  9.22 7.45 6.72 9.49 7.45 ...
##  $ Sample_34: num  8.99 7.22 6.5 8.95 7.66 ...
##  $ Sample_35: num  9.75 7.97 6.55 8.53 7.3 ...
##  $ Sample_13: num  9.84 6.68 6.55 7.73 7.07 ...
##  $ Sample_15: num  9.55 6.88 6.64 7.69 7.05 ...
##  $ Sample_18: num  10.09 6.59 6.44 7.15 6.71 ...
##  $ Sample_19: num  9.45 7.14 6.65 7.58 7.1 ...
##  $ Sample_36: num  9.58 8.51 6.24 10.05 6.58 ...
##  $ Sample_37: num  8.44 7.79 6.6 9.7 7.27 ...
##  $ Sample_38: num  9.42 7.75 6.43 10.08 6.79 ...
##  $ Sample_39: num  8.83 7.7 6.5 9.87 7.21 ...
##  $ Sample_11: num  9.38 6.58 6.36 7.11 7.01 ...
##  $ Sample_12: num  9 7.19 6.49 8.31 7.52 ...
##  $ Sample_2 : num  9.07 7.06 6.49 7.54 7.05 ...
##  $ Sample_9 : num  10.35 7 6.28 9.58 6.69 ...
```

```r
miniDat <- data.frame(gExp = as.vector(t(as.matrix(miniDat))), gene = factor(rep(rownames(miniDat), 
    each = ncol(miniDat)), levels = keepgenes))
str(miniDat)
```

```
## 'data.frame':	234 obs. of  2 variables:
##  $ gExp: num  10.6 11 10.8 10.9 9.2 ...
##  $ gene: Factor w/ 6 levels "1431708_a_at",..: 4 4 4 4 4 4 4 4 4 4 ...
```

```r
miniDat <- suppressWarnings(data.frame(prDes, miniDat))
str(miniDat)
```

```
## 'data.frame':	234 obs. of  6 variables:
##  $ sidChar : chr  "Sample_20" "Sample_21" "Sample_22" "Sample_23" ...
##  $ sidNum  : num  20 21 22 23 16 17 6 24 25 26 ...
##  $ devStage: Factor w/ 5 levels "E16","P2","P6",..: 1 1 1 1 1 1 1 2 2 2 ...
##  $ gType   : Factor w/ 2 levels "wt","NrlKO": 1 1 1 1 2 2 2 1 1 1 ...
##  $ gExp    : num  10.6 11 10.8 10.9 9.2 ...
##  $ gene    : Factor w/ 6 levels "1431708_a_at",..: 4 4 4 4 4 4 4 4 4 4 ...
```

```r
library(lattice)
stripplot(gType ~ gExp | gene, miniDat, scales = list(x = list(relation = "free")), 
    group = gType, auto.key = TRUE)
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-12.png) 

```r
someDat <- subset(miniDat, gene == keepgenes[1])
t.test(gExp ~ gType, someDat)
```

```
## 
## 	Welch Two Sample t-test
## 
## data:  gExp by gType
## t = 9.838, df = 36.89, p-value = 7.381e-12
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  1.570 2.384
## sample estimates:
##    mean in group wt mean in group NrlKO 
##               9.554               7.578
```

```r
install.packages(pkgs = "plyr")
```

```
## Installing package into 'C:/Users/Pan & Wang/Documents/R/win-library/3.0'
## (as 'lib' is unspecified)
```

```
## Error: trying to use CRAN without setting a mirror
```

```r
library(plyr)
d_ply(miniDat, ~gene, function(x) t.test(gExp ~ gType, x), .print = TRUE)
```

```
## 
## 	Welch Two Sample t-test
## 
## data:  gExp by gType
## t = 9.838, df = 36.89, p-value = 7.381e-12
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  1.570 2.384
## sample estimates:
##    mean in group wt mean in group NrlKO 
##               9.554               7.578 
## 
## 
## 	Welch Two Sample t-test
## 
## data:  gExp by gType
## t = 9.061, df = 36.49, p-value = 7.146e-11
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  0.6239 0.9835
## sample estimates:
##    mean in group wt mean in group NrlKO 
##               7.670               6.867 
## 
## 
## 	Welch Two Sample t-test
## 
## data:  gExp by gType
## t = 8.077, df = 33.49, p-value = 2.278e-09
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  1.402 2.346
## sample estimates:
##    mean in group wt mean in group NrlKO 
##               12.85               10.97 
## 
## 
## 	Welch Two Sample t-test
## 
## data:  gExp by gType
## t = -0.184, df = 36.53, p-value = 0.8551
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  -0.5079  0.4234
## sample estimates:
##    mean in group wt mean in group NrlKO 
##               9.893               9.935 
## 
## 
## 	Welch Two Sample t-test
## 
## data:  gExp by gType
## t = -0.1324, df = 36.31, p-value = 0.8954
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  -0.2413  0.2118
## sample estimates:
##    mean in group wt mean in group NrlKO 
##               7.092               7.107 
## 
## 
## 	Welch Two Sample t-test
## 
## data:  gExp by gType
## t = 0.0983, df = 35.58, p-value = 0.9223
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  -0.1019  0.1123
## sample estimates:
##    mean in group wt mean in group NrlKO 
##               6.551               6.546
```

```r
ttRes <- dlply(miniDat, ~gene, function(x) t.test(gExp ~ gType, x))
ttRes
```

```
## $`1431708_a_at`
## 
## 	Welch Two Sample t-test
## 
## data:  gExp by gType
## t = 9.838, df = 36.89, p-value = 7.381e-12
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  1.570 2.384
## sample estimates:
##    mean in group wt mean in group NrlKO 
##               9.554               7.578 
## 
## 
## $`1424336_at`
## 
## 	Welch Two Sample t-test
## 
## data:  gExp by gType
## t = 9.061, df = 36.49, p-value = 7.146e-11
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  0.6239 0.9835
## sample estimates:
##    mean in group wt mean in group NrlKO 
##               7.670               6.867 
## 
## 
## $`1454696_at`
## 
## 	Welch Two Sample t-test
## 
## data:  gExp by gType
## t = 8.077, df = 33.49, p-value = 2.278e-09
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  1.402 2.346
## sample estimates:
##    mean in group wt mean in group NrlKO 
##               12.85               10.97 
## 
## 
## $`1416119_at`
## 
## 	Welch Two Sample t-test
## 
## data:  gExp by gType
## t = -0.184, df = 36.53, p-value = 0.8551
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  -0.5079  0.4234
## sample estimates:
##    mean in group wt mean in group NrlKO 
##               9.893               9.935 
## 
## 
## $`1432141_x_at`
## 
## 	Welch Two Sample t-test
## 
## data:  gExp by gType
## t = -0.1324, df = 36.31, p-value = 0.8954
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  -0.2413  0.2118
## sample estimates:
##    mean in group wt mean in group NrlKO 
##               7.092               7.107 
## 
## 
## $`1429226_at`
## 
## 	Welch Two Sample t-test
## 
## data:  gExp by gType
## t = 0.0983, df = 35.58, p-value = 0.9223
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  -0.1019  0.1123
## sample estimates:
##    mean in group wt mean in group NrlKO 
##               6.551               6.546 
## 
## 
## attr(,"split_type")
## [1] "data.frame"
## attr(,"split_labels")
##           gene
## 1 1431708_a_at
## 2   1424336_at
## 3   1454696_at
## 4   1416119_at
## 5 1432141_x_at
## 6   1429226_at
```

```r
names(ttRes)
```

```
## [1] "1431708_a_at" "1424336_at"   "1454696_at"   "1416119_at"  
## [5] "1432141_x_at" "1429226_at"
```

```r
ttRes["1454696_at"]
```

```
## $`1454696_at`
## 
## 	Welch Two Sample t-test
## 
## data:  gExp by gType
## t = 8.077, df = 33.49, p-value = 2.278e-09
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  1.402 2.346
## sample estimates:
##    mean in group wt mean in group NrlKO 
##               12.85               10.97
```

```r

```


