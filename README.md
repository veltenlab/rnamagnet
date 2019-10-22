To install from within R, you can use the `install_git` function from the `devtools` package.

```
install.packages("devtools")
devtools::install_git("https://git.embl.de/velten/rnamagnet/")
```

In case this fails, download the tar.gz archive, unpack to a directory of your choice (e.g. /path/to) and inside R run
```
devtools::install_local("/path/to/rnamagnet_dir") 
```

**MAGIC is an important requirement of this package, an needs to be installed as a python package.** See https://github.com/KrishnaswamyLab/MAGIC/tree/master/Rmagic#installation

The vignettes now require seurat v3, but the functions work both with seurat v2 and seurat v3 objects.

If you want to build the vignettes yourself or make the demo data from our manuscript avaialable, please proceed as follows:
* Download and unpack the tar.gz archive, e.g. to /path/to
* Then download our data bundle (2GB, containing all 10x and LCM data) from https://www.dropbox.com/s/wbnqaebqi74j5ic/RNAMagnetDataBundle.zip
* Move the content of the data bundle to the data/ folder of the package, e.g. /path/to/rnamagnet_dir/data
* Then, run `devtools::install_local("/path/to/rnamagnet_dir")` or `devtools::install_local("/path/to/rnamagnet_dir", build_vignettes=T)`
