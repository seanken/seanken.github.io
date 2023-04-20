---
title: "Single Cell Level QC"
date: 2023-04-09T21:43:39Z
draft: false
---

## Cell Level QC

One thing that bothered me when I started working with 10X data (roughly 2016, after having worked with Smart-seq2 single nuclei data) is that for many metrics (percent intronic reads, percent intergenic reads, etc), though it was easy to get sample level values, but much less straightforward to get cell level metrics. Since then the situation has gotten better, with much work on using different metrics to help remove low quality cells. In particular, in recent years there has been the realization that % intronic reads in particular give a lot of insight for seperating real cells or nuclei from empty droplets and debris (the first such paper I am aware of was the DropletQC (https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02547-0) paper, with more recent papers from the Kellis lab (https://www.biorxiv.org/content/10.1101/2022.10.21.513315v1) looking at it for single nuclei data, or papers from the Linnarsson lab (https://www.biorxiv.org/content/10.1101/2022.10.24.513487v1) look at it at an atlas scale). 

In particular, in the case of single nuclei data, the more reads that are from the nucleus the higher percent intronic reads one would expect, allowing us to seperate nuclei from empty droplets and other parts of the cell (something a recent paper we were slightly involved with explored to help understand ambient RNA (https://www.biorxiv.org/content/10.1101/2022.11.16.516780v1). This raises two questions: are there other metrics that will give us similiar (or even more) insight but which we arent reporting now? And are there better ways to get the information out after processing 10X data? In particular, for those who use CellRanger, though there is a lot of read level QC information hidden in the output bam, it can be hard to get at. In particular, to get at % intronic levels it is common to use other tools lik STARSolo, Kallisto, or others velocyto. Not to imply they are anything less than great tools (they are all great pieces of software!), but they do not calculate all the metrics we care about, and in many cases we might not want to use them.

Given the above, I decided it might be worth the time to write a script to extract as much cell level QC information at possible from the CellRanger bam. As such, I built a tool I am calling CellLevel_QC (https://github.com/seanken/CellLevel_QC). It is a simple java program, all one has to do is download the jar file and run it on the outs directory from CellRanger (see the github for details). The code does not do anything particularly clever, for the most part just goes line by line in the bam file and extracts information from each bam entry. It returns many cell level QC metrics that might be of interest (% intronic reads, % intergenic, % reads trimmed by CellRanger, etc). 

## Insight From Cell Level QC

To explore this, I used the CellLevel_QC tool to look at one sample (12 week old mouse Hippocampus single nuclei sample taken from a recent prepub of ours, https://www.biorxiv.org/content/10.1101/2022.11.15.516665v1). 

I started by running CellLevel_QC on the sample

```
java -jar SingleCellQC.jar -d outs -o out.txt
```

where `outs` is the output directory from CellRanger (v6.1.2, run with introns included with the mm10 reference). This results in an output file, `out.txt`, into R, and normalizing counts to percent:

```
library(ggplot2)
library(cowplot) //Not strictly needed, just to make look nice
theme_set(theme_cowplot()) //Not strictly needed, just to make look nice

dat=read.table("out.txt",header=T)
for(i in c(2:11,15)){dat[i]=100*dat[,i]/dat[,"total"]}dat=dat[dat$total>10,] //Only look at droplets with at least 10 reads
dat=dat[dat[,1]!="notCell",] //removes non-cell entry
cells=scan("outs/filtered_feature_bc_matrix/barcodes.tsv.gz","") //Thinks cellranger thinks are nuclei
dat["Cell"]=dat$CBC %in% cells
```

The first thing to look at is the nUMI versus the % intronic reads, colored by if a droplet is declared a nuclei by CellRanger:

```
ggplot(dat,aes(x=nUMI,y=intronic,color=Cell))+geom_point()+scale_x_log10()+ylab("Percent Intronic Reads")
```

![Identifying Nuclei](Cell.label.png)

Can see a nice seperation into 3 distinct clusters (with some other droplets here or there) as one might hope/expect, including seeing a cluster that corresponds very well to the CellRanger labelled nuclei (this is not always so clean, sometimes CellRanger definitely sets the cutoff poorly, eleading to loss of nuclei or empty droplets being misidentified as nuclei). We do, however, have many other metrics. We can first look at how the different metric compare using Spearman correlation:

```
library(ComplexHeatmap)

COR=cor(dat[,c(2:14,16)],method="spearman")
Heatmap(COR)
```

![Heatmap](https://github.com/seanken/seanken.github.io/tree/main/static/Images_Single-Cell-Level-Q/Heatmap.png)

For definition of each of the metrics see the CellLevel_QC github.

Can see obvious structure in this plot. The most obvious strucuture comes from total (total number of reads) being highly correlated with nUMI (not suprising), and from the percent antisense being highly correlated to percent intronic (also probably not suprising). One obvious question: can these other metrics be used to improve filtering? For example consider the polyA metric (the percent of reads in a cell that were trimmed to remove a polyA sequence). I particular, we can color the % intronic be nUMI plot by if the cell has >1% of reads being trimmed for polyA:

```
ggplot(dat,aes(x=nUMI,y=intronic,color=polyA>1))+geom_point()+scale_x_log10()+ylab("Percent Intronic Reads")
```

![PolyA Trimming](https://github.com/seanken/seanken.github.io/tree/main/static/Images_Single-Cell-Level-Q/PolyA.label.png)

Can see that most of the droplets with high % polyA trimming are the lower nUMI, lower intronic droplets, which are likely empty droplets. We see something similiar for trimming of the TSO sequence, suggesting these metrics might be useful for better seperating empty droplets and real nuclei (in cases were it is less clear than the above).

## Final Thoughts

Anyways, just wanted to share a tool that we have found useful for some of our analysis--hope others find it useful as well! Wondering what other insight we can get from these metrics. Also: happy to think about adding other metrics as well (assuming I have the time), so let me know if there are any!
