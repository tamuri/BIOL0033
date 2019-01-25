# Software installation

These are instructions for installing software required for the BIOL0033 phylogeny reconstruction practical sessions. We use a number of libraries and packages available for the R programming language.

The practicals require:

1. R programming language
2. RStudio development environment
3. Several specific packages and libraries for phylogenetic analysis

If you plan on using [Desktop@UCL Anywhere](https://www.ucl.ac.uk/isd/services/computers/remote-access/desktopucl-anywhere) you can skip steps 2 & 3 because R and RStudio are already installed.

Optional (but recommended) tools to install are:

If you are on Windows or Linux, you will need to install [Java Runtime](https://www.java.com/en/download/) to use the following:

1. [AliView](https://ormbunkar.se/aliview/) - sequence alignment viewer
2. [Figtree](http://tree.bio.ed.ac.uk/software/figtree/) - phylogenetic tree viewer

## Installing the R programming language

1. Download the [R distribution](https://cran.ma.imperial.ac.uk/) for your platform (Windows, Mac and Linux packages available).
2. Run the executable and follow the instructions. The default options are fine.

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

My final output looks like this:

```
> version$version.string
[1] "R version 3.5.2 (2018-12-20)"
>
```

## Install packages

Researchers and developers have written many packages and libraries for the R platform for phylogenetics. See the "[phylogenetics task view](https://cran.r-project.org/web/views/Phylogenetics.html)" for an overview of many of the available packages. We only install and use some of these packages for our practicals.

During the installation of packages, if you see a message asking "Update all/some/none? [a/s/n]" type `a` and press return.

1. Find the 'Console' window. The symbol `>` indicates a prompt waiting for your input.
2. Type (or copy and paste) the following into the RStudio console and press return

```
install.packages(c('BiocManager', 'adegenet', 'ape', 
'geiger', 'phangorn', 'phylogram', 'phytools'), dep=TRUE)
```

You will see packages being downloaded and installed. Wait for the installation process to finish (when the output stops and prompt is again displayed).

3. To install packages for multiple sequence alignment and tree visualation, copy & paste the following into the RStudio console and press return

```
BiocManager::install("msa", version = "3.8")
BiocManager::install("ggtree", version = "3.8")
```



