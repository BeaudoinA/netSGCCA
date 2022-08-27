# netSGCCA


<h3> Introduction </h3>
-----------------

netSGCCA is a method based on the integration of multiblock omics data. This method derives from framework of the Generalized Canonical Correlation Analysis, 
and most particulary on the SGCCA (Sparse Generalized Canonical Correlation Analysis), which allows studying the relationship between several groups of variables 
and choosing the most pertinent variables and on the GraphNet penalty, which allows the integration of an information network reflecting the interactions among 
the variables within a given data block. One of the objectives of the method is to remove the “high-frequency” components induced by the GraphNet penalty, meaning 
sharp variations of the loadings between variables that are neighbors in the reference network. 
This method is applied on three data blocks: two genes expression blocks and one clinical data block. 


<h3> Installation </h3>
-----------------
If you want to use this method, you need to install this package with the following commands:
```{r, eval=FALSE }
if (!require(devtools)){install.packages("devtools")}

devtools::install_github("BeaudoinA/netSGCCA")
```



<h3>Parameters Description </h3>
-----------------

We have four parameters to determine. The $\lambda$ which is the regulation parameter between SGCCA and GraphNet penalty. If $\lambda$=0, there is only sgcca 
so the blocks interaction and if $\lambda$=1, we have the graphnet penalty. So we want to find and intermediate $\lambda$. $s_1$ and $s_2$ are the quantity of 
sparsity and $\sigma$ is the eigenvalue filter in order to reduce the high frequencies and give priority to the low frequencies. Because we didn’t have variables 
to help us to optimize these parameters, we defined three criteria in order to choose these parameters. 
