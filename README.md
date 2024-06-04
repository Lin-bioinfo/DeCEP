# CoGPS

## A novel method for deciphering context-specific gene programs using single-cell and spatial transcriptomics data

`CoGPS` is a computational framework for characterizing context-specific gene programs within scRNA-seq and ST datasets. `CoGPS` leverages directed graphs and functional gene lists to construct functional networks underlying distinct cellular and spatial contexts, which further identifies context-dependent hub genes based on the topology of networks and assigns activity levels of gene programs to individual cells or spatial locations.

![CoGPS](https://upload-images.jianshu.io/upload_images/14476738-78ef249acf7cf569.jpg?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

## Installation
You can install the released version of **CoGPS** with:

```r
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("Lin-bioinfo/CoGPS")
```

## Getting started
##### 1. Static cellular contexts (the test data can be downloaded from this [drive_link](https://drive.google.com/drive/folders/1UzkEEtqDauwOAJrU7JQe-rYR3V20DyDp?usp=drive_link).

```r
library(CoGPS)
library(parallel)
cores = detectCores(logical = FALSE) - 1

load("CoGPS_data_sc_discrete.RData")
dat = pre_sc(x = x, label = label, modes = c("discrete"), genelist = genelist, transition_point = NULL)

cl = makeCluster(cores)
net_A = net_sc(dat$A, modes = c("discrete"), genelist = genelist, select = 1, cluster = cl, seed = 1)
cl = makeCluster(cores)
net_B = net_sc(dat$B, modes = c("discrete"), genelist = genelist, select = 1, cluster = cl, seed = 1)

CoGPS_score_A = ident_sc_score(dat$A, net_A, seed = 1)
CoGPS_score_B = ident_sc_score(dat$B, net_B, seed = 1)
CoGPS_state_A = ident_sc_state(CoGPS_score_A)
CoGPS_state_B = ident_sc_state(CoGPS_score_B)
```

##### 2. Dynamic cellular contexts
```r
library(CoGPS)
library(parallel)
cores = detectCores(logical = FALSE) - 1

load("CoGPS_data_sc_continuous.RData")
dat = pre_sc(x = x, label = label, modes = c("continuous"), genelist = genelist, transition_point = 20)

cl = makeCluster(cores)
net = net_sc(dat, modes = c("continuous"), genelist = genelist, select = 1, cluster = cl, seed = 1)

CoGPS_score = ident_sc_score(x, net, seed = 1)
CoGPS_state = ident_sc_state(CoGPS_score)
```

##### 3. Spatial contexts
```r
load("CoGPS_data_st.RData")
st_anno = pre_st(x_counts, coords, genelist, ref, cell_state, core_number = 2)
roiID = st_anno[st_anno$CoGPS_state == "High", "spot_ID"]
x_st = net_st_pre(x_data, coords, genelist, roiID, r = 60)

cl = makeCluster(cores)
net_st = net_st(x_st$dat, genelist, select = 1, cluster = cl, seed = 1)
involvedID = c(roiID, x_st$neighborID)
CoGPS_score = ident_st_score(x_impute, involvedID, net_st)   # use data after imputation
```


## License
CoGPS is licensed under the GNU General Public License v3.0.

## Citation
Please cite the following article if you use CoGPS in your research: 

## Contact
Should you have any questions, please feel free to write issues for this repository or contact the author directly at [lin.li.bioinfo@gmail.com](lin.li.bioinfo@gmail.com).
