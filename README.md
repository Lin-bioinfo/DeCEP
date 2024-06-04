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

## License

CoGPS is licensed under the GNU General Public License v3.0.

## Citation

Please cite the following article if you use CoGPS in your research: [https://](https://)

## Contact

Should you have any questions, please feel free to write issues for this repository or contact the author directly at [lin.li.bioinfo@gmail.com](lin.li.bioinfo@gmail.com).
