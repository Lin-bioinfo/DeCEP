# DeCEP

## A novel method for deciphering context-specific gene programs using single-cell and spatial transcriptomics data

`DeCEP` is a computational framework designed to characterize context-specific gene programs using scRNA-seq and ST data. `DeCEP` leverages functional gene lists and directed graphs to construct functional networks underlying distinct cellular or spatial contexts. It then identifies context-dependent hub genes based on the network topology and assigns gene program activity to individual cells or spatial locations.

## Installation
You can install the released version of **DeCEP** with:

```r
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("Lin-bioinfo/DeCEP")
```

## Getting started
### 1. Static cellular contexts (the test data can be downloaded from this [drive_link](https://drive.google.com/drive/folders/1UzkEEtqDauwOAJrU7JQe-rYR3V20DyDp?usp=drive_link).

```r
library(DeCEP)
library(parallel)
cores = detectCores(logical = FALSE) - 1

load("DeCEP_data_sc_discrete.RData")
dat = pre_sc(x = x, label = label, modes = c("discrete"), genelist = genelist, transition_point = NULL)

cl = makeCluster(cores)
net_A = net_sc(dat$A, modes = c("discrete"), genelist = genelist, select = 1, cluster = cl, seed = 1)
cl = makeCluster(cores)
net_B = net_sc(dat$B, modes = c("discrete"), genelist = genelist, select = 1, cluster = cl, seed = 1)

DeCEP_score_A = ident_sc_score(dat$A, net_A, seed = 1)
DeCEP_score_B = ident_sc_score(dat$B, net_B, seed = 1)
DeCEP_state_A = ident_sc_state(DeCEP_score_A)
DeCEP_state_B = ident_sc_state(DeCEP_score_B)
```

### 2. Dynamic cellular contexts
```r
library(DeCEP)
library(parallel)
cores = detectCores(logical = FALSE) - 1

load("DeCEP_data_sc_continuous.RData")
dat = pre_sc(x = x, label = label, modes = c("continuous"), genelist = genelist, transition_point = 20)

cl = makeCluster(cores)
net = net_sc(dat, modes = c("continuous"), genelist = genelist, select = 1, cluster = cl, seed = 1)

DeCEP_score = ident_sc_score(x, net, seed = 1)
DeCEP_state = ident_sc_state(DeCEP_score)
```

### 3. Spatial contexts
```r
load("DeCEP_data_st.RData")
st_anno = pre_st(x_counts, coords, genelist, ref, cell_state, core_number = 2)
roiID = st_anno[st_anno$DeCEP_state == "High", "spot_ID"]
x_st = net_st_pre(x_data, coords, genelist, roiID, r = 60)

cl = makeCluster(cores)
net_st = net_st(x_st$dat, genelist, select = 1, cluster = cl, seed = 1)
involvedID = c(roiID, x_st$neighborID)
DeCEP_score = ident_st_score(x_impute, involvedID, net_st)   # use data after imputation
```


## License
DeCEP is licensed under the GNU General Public License v3.0.

## Citation
Please cite the following article if you use DeCEP in your research: Genome Research (doi: 10.1101/gr.279689.124
        
        
        
        

## Contact
Should you have any questions, please feel free to write issues for this repository or contact the author directly at [lin.li.bioinfo@gmail.com](lin.li.bioinfo@gmail.com).
