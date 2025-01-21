# BatchMetrics
Package for evaluating batch correction methods in R

# Installation
Install the package using devtools

```{r}
devtools::install_github("Aleksi95/BatchMetrics")

library(BatchMetrics)
```


# Usage

Load in the dataset of RNA-Seq count data


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

Constructing the Artificial dilution series:

```{r}
ads = diluteSeries(Data = data, batch = batch, corr_method = "ComBat")
```
Using evaluation metrics to evaluate level of batch effect:

```{r}
#data with full batch effect
evalBatchEffect(ads[[1]], sample_types = group, batch = batch, metric = "Davies-Bouldin")

#data with fully removed batch effect
evalBatchEffect(ads[[11]], sample_types = group, batch = batch, metric = "Davies-Bouldin")
```


Generating datasets with random variation to evaluate metrics:

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

