---
title: "Single Cell Level QC"
date: 2023-04-09T21:43:39Z
draft: false
---

## Cell Level QC

One thing that bothered me when I started working with 10X data (roughly 2016, after having worked with Smart-seq2 single nuclei data) is that for many metrics (percent intronic reads, percent intergenic reads, etc), though it was easy to get sample level values, but much less straightforward to get cell level metrics. Since then the situation has gotten better, with much work on using different metrics to help remove low quality cells. In particular, in recent years there has been the realization that % intronic reads in particular give a lot of insight for seperating real cells or nuclei from empty droplets and debris (the first such paper I am aware of was the DropletQC [link] paper, with more recent papers from the Kellis lab [link] looking at it for single nuclei data, or papers from the ... lab [link] look at it at an atlas scale). 

In particular, in the case of single nuclei data, the more reads that are from the nucleus the higher percent intronic reads one would expect, allowing us to seperate nuclei from empty droplets and other parts of the cell (something a recent paper we were slightly involved with explored to help understand ambient RNA [link]). This raises two questions: are there other metrics that will give us similiar (or even more) insight but which we arent reporting now? And are there better ways to get the information out after processing 10X data? In particular, for those who use CellRanger, though there is a lot of read level QC information hidden in the output bam, it can be hard to get at. In particular, to get at % intronic levels it is common to use other tools lik STARSolo, Kallisto, or others velocyto. Not to imply they are anything less than great tools (they are all great pieces of software!), but they do not calculate all the metrics we care about, and in many cases we might not want to use them.

Given the above, I decided it might be worth the time to write a script to extract as much cell level QC information at possible from the CellRanger bam. As such, I built a tool I am calling CellLevel_QC (https://github.com/seanken/CellLevel_QC). It is a simple java program, all one has to do is download the jar file and run it on the outs directory from CellRanger (see the github for details). The code does not do anything particularly clever, for the most part just goes line by line in the bam file and extracts information from each bam entry. It returns many cell level QC metrics that might be of interest (% intronic reads, % intergenic, % reads trimmed by CellRanger, etc). 

## Insight From Cell Level QC

To explore this, I used the CellLevel_QC tool to look at one sample (taken from a recent publication of ours, ...). 

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
for(i in ..){dat[i]=100*dat[,i]/sum(dat[,i])}
dat=dat[dat$total>10,] //Only look at droplets with at least 10 reads
```

The first thing to look at is the nUMI versus the % intronic reads:

```
ggplot(dat,aes(x=nUMI,y=intronic))+geom_point()+scale_x_log10()+ylab("Percent Intronic Reads")
```

[insert plot]

Can see a nice seperation into 3 distinct clusters (with some other droplets here or there) as one might hope/expect. We do, however, have many other metrics. We can first look at how the different metric compare:

```
library(ComplexHeatmap)

COR=cor(dat[,2:dim(dat)[2]])
Heatmap(COR)
```

[insert plot]


Can see obvious structure in this plot, with some 
