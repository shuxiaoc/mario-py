# MARIO: single-cell proteomic data matching and integration pipeline

<img src="https://github.com/shuxiaoc/mario-py/blob/main/media/giphy_mario.gif" width="100" height="100">


## [<img src="https://github.com/shuxiaoc/mario-py/blob/main/media/red.png" width="25" height="25">](https://www.youtube.com/watch?v=2iNKPkTOr5k&ab_channel=13irth) Description

This github repo includes `mario-py` and `mario-R`, which is a Python package for matching and integrating multi-modal single cell data with partially overlapping features. The method is specifically tailored toward proteomic datasets, and for detailed description on the algorithm, including the core methodology, mathmetical ingredients, application on various biological samples, and extensive benchmarking, please refer to the [paper](https://www.biorxiv.org/content/10.1101/2021.12.03.471185v1).

This work has been lead by Shuxiao Chen from [Zongming Lab](http://www-stat.wharton.upenn.edu/~zongming/) @Upenn and Bokai Zhu from [Nolan lab](https://web.stanford.edu/group/nolan/) @Stanford.

<img src="https://github.com/shuxiaoc/mario-py/blob/main/media/mario_flow.png" width="800" height="280">

## <img src="https://github.com/shuxiaoc/mario-py/blob/main/media/green.png" width="25" height="25"> Getting Started

### Dependencies

For easy usage, we suggest builing a ```conda``` virtualenv with ```python = 3.8```.

```{bash}
conda create -n mario python=3.8
```

### Installing

To install ```MARIO```, we can easily install it with ```pip``` function (package name ```pyMARIO```):

```{bash}
python -m pip install pyMARIO
```

## <img src="https://github.com/shuxiaoc/mario-py/blob/main/media/blue.png" width="25" height="25"> How to use

### Quick example:

To use in ```MARIO``` in ```python``` :
```
from mario.match import pipelined_mario
final_matching_lst, embedding_lst = pipelined_mario(data_lst=[df1, df2])
```
Where ```df1``` and ```df2``` are two dataframes for match and integration, with row as cells, columns as features. Remember for shared features, the column names should be identical. Input list can be multiple dataframes, as ```MARIO``` accomodates for multiple dataset match and integration.

The result contains the a matching list (matching), and a embedding list (integration). For detailed usage please refer to the Full tutorial section.  

Similarly, to use in ```MARIO``` in ```R``` (with package ```reticulate```) :
```
library(reticulate)
myenvs=conda_list() # get conda virtualenv list
envname=myenvs$name[12] # specify which virtualenv to use, should use the one for MARIO
use_condaenv(envname, required = TRUE)
mario.match <- import("mario.match") # import main mario-py module

pipelined_res = mario.match$pipelined_mario(data_lst=list(df1, df2))
```
Where the result also contains the matching list and embedding list.

### Full tutorial:
For step by step tutorials on how to use ```MARIO```, with fine-tuned parameters for optimal results and full functionality, please refer to the documents we provided here:

[Python - Jupyter notebook: Match and Integration of Human Bonemarrow datasets](https://github.com/shuxiaoc/mario-py/blob/main/tutorials/mario-py-tutorial-BM.ipynb)

[Python - Jupyter notebook: Match and Integration of multiple Xspecies datasets](https://github.com/shuxiaoc/mario-py/blob/main/tutorials/mario-py-multiple-Xspecies.ipynb)

[R - Rmarkdown: Match and Integration of Human Bonemarrow datasets](https://github.com/shuxiaoc/mario-py/blob/main/tutorials/mario-r-bk.md)


## <img src="https://github.com/shuxiaoc/mario-py/blob/main/media/yellow.png" width="25" height="25"> License and Citation

```MARIO``` is under the [Academic Software License Agreement](https://github.com/shuxiaoc/mario-py/blob/main/LICENSE.txt), please use accordingly.

