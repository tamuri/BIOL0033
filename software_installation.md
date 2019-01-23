#-------------------------------------------------------------------------------
#  Tue Jan 22 15:55:32 GMT 2019
#-------------------------------------------------------------------------------

# Software installation instructions for molecular evolution practicals

We will use a number of libraries and package available for the R programming language for our phylogeny reconstruction practicals.

You will need to install:

1. The R programming language
2. The RStudio development environment
3. Various packages and libraries for phylogenetic analysis

If you plan on using your UCL@Desktop account, you can skip steps 2 & 3 because R and RStudio is already installed.

Optional (but recommended) tools to install are:

1. AliView - alignment viewer - https://ormbunkar.se/aliview/
2. Figtree - phylogenetic tree viewer - http://tree.bio.ed.ac.uk/software/figtree/

## Installing the R programming language

1. Download the R distribution from https://cran.ma.imperial.ac.uk/ for your computer system (Windows, Mac and Linux packages are available).
2. Run the executable follow the instructions. The default options are fine.

## Installing RStudio

RStudio is a development environment for the R programming language.

1. Download the 'RStudio Desktop Open Source License' from https://www.rstudio.com/products/rstudio/download/ for your platform.
2. Run the executable and follow the instructions. The default options are fine.

## Start RStudio

Once you have installed R and RStudio, check that the programs installed correctly.

1. Start RStudio (should be a program in your list of applications)
2. Find the the 'Console' window. The symbol `>` indicates a prompt waiting for your input.
3. Type the following:

`version$version.string`

4. My final output looks like this:

```
Type 'demo()' for some demos, 'help()' for on-line help, or
'help.start()' for an HTML browser interface to help.
Type 'q()' to quit R.

> version$version.string
[1] "R version 3.5.0 (2018-04-23)"
>
```

## Install packages we need for the practicals

Researchers and developers have written many packages and libraries for the R platform for phylogenetics.

R provides a "task view" for phylogenetics. You can read more about it at https://cran.r-project.org/web/views/Phylogenetics.html

We're only going install some of these packages for our practicals.

During the package installation, if you get a message "Update all/some/none? [a/s/n]" type `a` and continue.

1. Find the 'Console' window. The symbol `>` indicates a prompt waiting for your input.
2. Type the following and press return

`install.packages(c('ape', 'adegenet', 'phangorn'), dep=TRUE)`

You will see packages being downloaded and installed. Wait for the installation process to finish (when the output stops and prompt is again displayed).

3. We now install packages for multiple sequence alignment. Type the following and press return

```
if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")

BiocManager::install("msa", version = "3.8")
```


