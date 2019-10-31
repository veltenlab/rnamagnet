To install from within R, you can use the `install_git` function from the `devtools` package.

```
install.packages("devtools")
devtools::install_github("veltenlab/rnamagnet")
```

In case this fails, download the zip archive from this website, unpack to a directory of your choice (e.g. /path/to) and inside R run
```
devtools::install_local("/path/to/rnamagnet_dir") 
```

**MAGIC is an important requirement of this package, an needs to be installed as a python package.** See https://github.com/KrishnaswamyLab/MAGIC/tree/master/Rmagic#installation

The vignettes requires seurat v3, but the functions work both with seurat v2 and seurat v3 objects.

If you want to build the vignettes yourself or make the data from our manuscript avaialable inside the package, please proceed as follows:
* Download and unpack the zip archive, e.g. to /path/to
* Then download our data bundle (1GB, containing all 10x and LCM data) from https://www.dropbox.com/s/wbnqaebqi74j5ic/RNAMagnetDataBundle.zip
* Move the content of the data bundle to the data/ folder of the package, e.g. /path/to/rnamagnet_dir/data
* Then, run `devtools::install_local("/path/to/rnamagnet_dir")`
* You can find the vignettes in the vignette folder and go through key analyses of the manuscript line by line.
