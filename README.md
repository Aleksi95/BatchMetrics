# BatchMetrics
Package for evaluating batch correction methods in R

# Installation
Install the package using devtools

```{r}
devtools::install_github("Aleksi95/BatchMetrics")

library(BatchMetrics)
```
Alternatively, you can use the remotes package

# Usage

Load in the dataset of RNA-Seq count data. These commands show description of the loaded datasets


```{r}
str(counts)
```

```{r}
str(data)
```

Show the batch effects

```{r}
print(levels(batch))
```

Show the biological sample groups:

```{r}
print(levels(group))
```

Constructing the Artificial dilution series. This step may take a few minutes:

```{r}
ads = diluteSeries(Data = data, batch = batch, corr_method = "ComBat")
```
This generates a list with datasets with different levels of batch effects removed (0% to 100% with 10% increments, so a list of 11 entries).

Using evaluation metrics to evaluate level of batch effect:

```{r}
#data with full batch effect
evalBatchEffect(ads[[1]], sample_types = group, batch = batch, metric = "Davies-Bouldin")

#data with fully removed batch effect
evalBatchEffect(ads[[11]], sample_types = group, batch = batch, metric = "Davies-Bouldin")
```
The results show three scores: bio, batch and ratio, that refer to the cluster separation level of biological effects, batch effects and their ratio (bio/batch) respectively. In the fully corrected dataset the separation of batch effects should be lower and the ratio of biological effects and batch effects should be higher.

Generating datasets with random variation to evaluate metrics. This may also take some time, depending on the number of iterations:

```{r}
if(!dir.exists("Results")){
  dir.create("Results")
}

setwd("Results")

#select the evaluation metrics
metrics = c("F-score", "Davies-Bouldin", "kNN", "gPCA")

Results = QC_wrapper(counts, batch, group, parallel = TRUE, var_measure = "resampling", iters = 10)
```

Write the results in a text file

```{r}
write_results(res_data = Results, filename = "results.txt")
```

Plot the results in a line plot

```{r}
linePlot2(Results, filename = "resultsLineplot.png")
```

