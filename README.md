# "Simulating temporal communities and assessing their detrended synchrony"

This package accompanies the article

Leps, J., Götzenberger, L., Valencia E,. de Bello Francesco (2019) Accounting for long-term directional trends on year-to-year synchrony in species fluctuations. Ecography. 

It provides functions to calculate synchrony indices for temporal ecological communities (i.e. records of co-occurring species abundances through time), and to simulate such communities.

## Loading packages and scripts

Only one additional package, gtools, is required apart from the ones in base R and loaded automatically to run the functions implemented in the script. So if you don't have gtools installed already, do so by
```{r, eval=FALSE}
install.packages("gtools")
```
If it is already installed, you can go ahead and load it, without installing it first.
```{r}
library(gtools)
```

Then, we have to source the scripts from whatever location you saved the script in on your computer. There are two separate scripts, one for the calculation of syncrhony indices, including the t3 version and the decomposition, and one script to simulate various degrees of synchrony among species of a community.

Assuming that your working directory is set to the folder that contains the scripts, we can load them via
```{r}
source("calc_sync.R")
source("syngenr.R")
```
Both these scripts contain a number of functions, but both have one main function, called just like the scripts themselves, i.e. calc_sync and syngenr. For these main functions we also provided two separate help files in pdf format in the appendices.

## Simulating temporal communities with syngenr

We will first have a look at the syngenr function.

As can be seen from the description of the function in the pdf help file, there are a lot of arguments that can be set in syngenr. While this makes the usage of the function slightly more complex, it accomodates all the different kind of fluctuation patterns that we might want our simulated temporal community to carry. To see what all these different arguments do, we can start with the easiest of options, which is to not specify any argument at all. Each of the arguments have default values, so that we can produce a temporal community simply with  

```{r syngenr no arguments}
sim_com <- syngenr()
str(sim_com)
```

The object produced by syngenr is a list with four elements, where the first element contains the actual simulated community. What is in the other three elements then? They contain the set of values that have been used to parametrize the simulation model (fourth element, named param_general), as well as the species specific values for mean abundance, standard deviation of this abundance, the response to the environment, and the trend response, each of which have been calculated based on the aforementioned general parameters. These species specific values are in the third list element, named param_species. There is one more set of generated parameters (named param_years), saved in the second element, which are the environmental cue, and the trend value, for each year. 

We can now modify the arguments one by one, starting with the two most basic ones: Years, and species number. They specify how many years the time series will span, as well as the overall number of species that will be in the simulated community.

```{r}
sim_com <- syngenr(years = 10, n_sp = 5)
```

We just created a temporal community that consists of five species, and contains 10 years of data. Specie are in columns, and years in rows. Let's have a look at that:

```{r}
dim(sim_com[[1]])
```

Now for the more interesting. How do we simulate, for instance, a community with prevailing synchrony between species? For this, we need to set the switch_env argument to "on", to have the
species respond the environmental cue that fluctuates from one year to the other. To generate a synchronous pattern, we also need to set bimodal_env to FALSE, so that the majority of species respond in the same way to the environmental cue. Simulating synchrony this way will not be influenced by any long term trends only if we "turn off" the trend, i.e. by setting switch_trend = "off". Once we have done this, we do not need to worry about the bimodal_trend argument, because it will not be effective anyway under a no trend scenario.

```{r}
set.seed(124)
sim_com_sync <- syngenr(years = 100, n_sp = 5, switch_env = "on", 
                        bimodal_env = FALSE, switch_trend = "off")
```

Looking at the synchrony indices of the temporal community just simulated, we should see (a) relatively high values, indicating synchrony, and not too big of a difference between the standard and t3 versions of the indices, since we did not simulate a long term trend for the species. Of course we might still see some differences, as the simulation model will produce some variation, that might look like a trend to the indices.  

So much for a basic overview of the syngenr function. We will now look at how to use calc_sync to obtain various synchrony measures. 

## Calculating synchrony indices wih calc_sync

To calculate synchrony indices, we use function 'calc_sync'. It produces a suite of indices in their "standard" form (eat, weighted eta, phi, variance ratio, log(variance ration)), and also their t3 versions. (See also the pdf help file for syngenr). So for the index used in the article (log(variance ratio) and its t3 version) we look up log_varrat and log_varrat_t3. Finally, it also produces (when decompose = TRUE) the three indices of the first method in the article, named 'syn_total', 'syn_trend', and 'syn_detrend' (Sdetrend) in the output. These correspond to Stotal, Strend, and Sdetrend in the article, respectively.
```{r}
calc_sync(sim_com_sync[[1]], decompose = TRUE)
```

Playing now with the paremeter settings in the syngenr function, we can turn on the trend, and make it unimodal. The majority of species will have a common long term trend, creating "additional synchrony". In other words, we should detect a pattern via the (standard) syncrhony indices that looks like synchrony according to the values, but is actually, at least in parts, a result of the common trend.
```{r}
set.seed(124)
sim_com_sync_trend <- syngenr(years = 100, n_sp = 5, switch_env = "on", 
                              bimodal_env = FALSE, switch_trend = "on", 
                              bimodal_trend = FALSE)
```
 
We can in fact see that some of the synchrony in our simulated temporal communities is due to a common trend when we compare the t3 versions of the synchrony indices with their standard counter parts, and also when we look at the values of Stotal, Strend, and Sdetrended, the last three entries in the row. 
```{r}
calc_sync(sim_com_sync_trend[[1]], decompose = TRUE)
```

Feel free to play around with all the arguments to create your own temporal communities (see also the pdf help file and description of the simulation model in the article). Let us know if something is not working, or if you want to suggest ideas how to improve the code, or how the mechanisms creating fluctuations could be simulated in better ways. 
