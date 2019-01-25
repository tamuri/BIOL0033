# Software installation

These are instructions for installing software required for the BIOL0033 phylogeny reconstruction practical sessions. We use a number of libraries and packages available for the R programming language.

The practicals require:

1. R programming language
2. RStudio development environment
3. Several specific packages and libraries for phylogenetic analysis

Optional (but recommended) tools to install are:

1. If you are on Windows or Linux, you will need to install [Java Runtime](https://www.java.com/en/download/) for the following programs.
2. [AliView](https://ormbunkar.se/aliview/) - sequence alignment viewer
3. [Figtree](http://tree.bio.ed.ac.uk/software/figtree/) - phylogenetic tree viewer

## Installing the R programming language

1. Download the [R distribution](https://cran.ma.imperial.ac.uk/) for your platform (Windows, Mac and Linux packages available).
2. Run the executable and follow the instructions. The default options are fine.

Warning: you may already have an old version of R installed on your machine. In this case, install the latest R version from (1) above and delete the old version. Some packages we use require the latest version of R.

## Installing RStudio

RStudio is a development environment for the R programming language.

1. Download the free [RStudio Desktop Open Source License](https://www.rstudio.com/products/rstudio/download/) for your platform.
2. Run the executable and follow the instructions. The default options are fine.

## Start and check RStudio

Once you have installed R and RStudio, check that the programs are installed correctly.

1. Start RStudio (should be a program in your list of applications)
2. Find the the 'Console' window. The symbol `>` indicates a prompt waiting for your input.
3. Type the following:

```
version$version.string
```

My output looks like this:

```
> version$version.string
[1] "R version 3.5.2 (2018-12-20)"
```

## Install packages

Researchers and developers have written many packages and libraries for the R platform for phylogenetics. See the "[phylogenetics task view](https://cran.r-project.org/web/views/Phylogenetics.html)" for an overview of many of the available packages. We only install and use some of these packages for our practicals.

Important: during package installation, you may see a message asking "Update all/some/none? [a/s/n]". Type `a` and press return.


1. Find the 'Console' window. The symbol `>` indicates a prompt waiting for your input.
2. Type (or copy and paste) the following into the RStudio console and press return

```
install.packages(c('BiocManager', 'adegenet', 'ape', 
'geiger', 'phangorn', 'phylogram', 'phytools'), dep=TRUE)
```

You will see packages being downloaded and installed. Wait for the installation process to finish (when the output stops and prompt is again displayed). It can take some time, be patient!

3. Install packages for multiple sequence alignment and tree visualisation. Copy & paste the following lines into the RStudio console and press return

```
BiocManager::install("msa", version = "3.8")
BiocManager::install("ggtree", version = "3.8")
```

4. Check that all the packages have been installed and can be loaded. Run the following in the R console

```
packages = c('BiocManager', 'adegenet', 'ape', 
'geiger', 'phangorn', 'phylogram', 'phytools')

package.check <- lapply(packages, FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
        install.packages(x, dependencies = TRUE)
        library(x, character.only = TRUE)
    }
})
```

You should see messages such as "Loading required package: xxxx" for each of the packages. There might be some info messages but no errors.

5. Finally, check the packages are in your search path. Type the command:

```
packages = c('BiocManager', 'adegenet', 'ape', 'geiger', 'phangorn', 'phylogram', 'phytools')
loaded <- search()
found <- loaded[grepl('^package:', loaded)]
found <- unlist(lapply(strsplit(found, split=':'), function(y) y[2]))
packages %in% found
```

The output should be

```
[1] TRUE TRUE TRUE TRUE TRUE TRUE TRUE
```

If it isn't, carefully install the packages again in steps 2 and 3 one at a time and make sure there are no errors

```
install.packages('BiocManager')
```

etc.

