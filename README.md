# BatchMetrics
Package for evaluating batch correction methods in R

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

Generating datasets with random variation to evaluate metrics:

```{r}
if(!dir.exists("Results")){
  dir.create("Results")
}

setwd("InhouseResults")

#select the evaluation metrics
metrics = c("F-score", "Davies-Bouldin", "kNN", "gPCA")

Results = QC_wrapper(counts, batch, group, parallel = TRUE, var_measure = "resampling", iters = 10)
```
