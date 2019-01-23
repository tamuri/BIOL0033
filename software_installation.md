# Software installation

These are instructions for installing software required for the BIOL0033 practical sessions. Primarily, we use a number of libraries and package available for the R programming language for our phylogeny reconstruction practicals.

The practicals require:

1. The R programming language
2. The RStudio development environment
3. Various packages and libraries for phylogenetic analysis

If you plan on using [Desktop@UCL Anywhere](https://www.ucl.ac.uk/isd/services/computers/remote-access/desktopucl-anywhere) you can skip steps 2 & 3 because R and RStudio are already installed.

Optional (but recommended) tools to install are:

1. [AliView](https://ormbunkar.se/aliview/) - sequence alignment viewer
2. [Figtree](http://tree.bio.ed.ac.uk/software/figtree/) - phylogenetic tree viewer

If you are on Windows or Linux, you will need to install [Java Runtime](https://www.java.com/en/download/) for the above two programs.

## Installing the R programming language

1. Download the [R distribution](https://cran.ma.imperial.ac.uk/) for your computer system (Windows, Mac and Linux packages available).
2. Run the executable and follow the instructions. The default options are fine.

## Installing RStudio

RStudio is a development environment for the R programming language.

1. Download the free [RStudio Desktop Open Source License](https://www.rstudio.com/products/rstudio/download/) for your platform.
2. Run the executable and follow the instructions. The default options are fine.

## Start RStudio

Once you have installed R and RStudio, check that the programs are installed correctly.

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

R provides a [phylogenetics task view](https://cran.r-project.org/web/views/Phylogenetics.html) which gives an overview of many of the available packages.

We're only going install some of these packages for our practicals.

During the package installation, if you get a message "Update all/some/none? [a/s/n]" type `a` and continue.

1. Find the 'Console' window. The symbol `>` indicates a prompt waiting for your input.
2. Type the following and press return

`install.packages(c('ape', 'adegenet', 'phangorn'), dep=TRUE)`

You will see packages being downloaded and installed. Wait for the installation process to finish (when the output stops and prompt is again displayed).

3. To install packages for multiple sequence alignment, copy & paste the following into your RStudio console and press return

```
if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")

BiocManager::install("msa", version = "3.8")
```


