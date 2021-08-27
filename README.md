# GOVS
<a href="https://www.r-project.org/" target="_blank"><img src="https://img.shields.io/badge/language-R-orange?style=plastic"></a>
<a href="https://cran.r-project.org/bin/windows/base/old/" target="_blank"><img src="https://img.shields.io/badge/R%20version-%3E%3D%203.6.0-orange?style=plastic"></a>
<a href="https://govs-pack.github.io/" target="_blank"><img src="https://img.shields.io/badge/webpage-ready-green?style=plastic"></a>
![](https://img.shields.io/badge/platform-Win%20%7C%20Linux%20%7C%20MacOS-lightgrey?style=plastic)<br/>

## Overview
__GOVS__ (**G**enome **O**ptimization via **V**irtual **S**imulation) is an integrative R package for maize breeding that streamlines genome optimization via virtual simulation achieve guidance of lines selection and population development. GOVS describes a promising strategy that can help breeders to select materials in a purposeful and directional manner, with the purpose of enabling breeders to combine others technologies with rapidly genetic gain benefit by genome optimization.
<div align="center">
<img src="https://govs-pack.github.io/img/overall.png" alt="drawing" width="900"/>
</div>
<div align="center">
  Overview of GOVS
</div>

## Installation
1.  Github install
```R
## install dependencies and GOVS
install.packages(c("ggplot2","rrBLUP","lsmeans","readr","pbapply","pheatmap","emmeas"))
require("devtools")
install_github("GOVS-pack/GOVS") 
## if you want build vignette in GOVS 
install_github("GOVS-pack/GOVS",build_vignettes = TRUE)
```
2.  Download [.tar.gz package](https://github.com/GOVS-pack/GOVS/raw/master/GOVS_1.0.tar.gz) and install <br/>
```R
## install dependencies and GOVS with bult-in vignette
install.packages(c("ggplot2","rrBLUP","lsmeans","readr","pbapply","pheatmap","emmeas"))
install.packages("DownloadPath/GOVS_1.0.tar.gz")
```
## Links
* GOVS Homepage: https://govs-pack.github.io/
* QuickStart: https://govs-pack.github.io/QuickStart/
* Tutorial: https://govs-pack.github.io/Tutorial/
* Reference Manual: https://github.com/GOVS-pack/GOVS/blob/master/GOVS-User-Manual.pdf
* Download: https://github.com/GOVS-pack/GOVS/raw/master/GOVS_1.0.tar.gz
