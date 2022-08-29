# netSGCCA package


## Introduction ##

netSGCCA is a method based on the integration of multiblock data. This method derives from framework of the Generalized Canonical Correlation Analysis, and most particulary on the SGCCA (Sparse Generalized Canonical Correlation Analysis) $^2$, which allows studying the relationship between several groups of variables and choosing the most pertinent variables and on the GraphNet $^1$ penalty, which allows the integration of an information network reflecting the interactions among the variables within a given data block. One of the objectives of the method is to remove the “high-frequency” $^3$ components induced by the GraphNet penalty, meaning sharp variations of the loadings between variables that are neighbors in the reference network. 


## Installation ##
If you want to use this method, you need to install this package with the following commands:
```{r, eval=FALSE }
if (!require(devtools)){install.packages("devtools")}

devtools::install_github("BeaudoinA/netSGCCA")
```


## Functions ##

### netSGCCA 
netSGCCA function performs netSGCCA method. It returns a list of weight vectors and the component associated with each block.

### Laplacian_matrix
Laplacian_matrix is a function which calculates the Laplacian matrix of a graph and which apply a filter (sigma) on the eigenvalues of the Laplacien matrix in order to reduce high frequencies induced by the GraphNet penalty.

### explo_param 
explor_param creates a search grid of the parameters in order to find the optimals parameters.
This function creates a dataframe in order to test each combination of the parameters and adds to this
dataframe the different values of the criteria used.
The criteria used to determine the best parameters are : <br/>
    - the sparsity (the number of zero variables),<br/>
    - the adequacy of data (if there are blocks associated with graphs and one block without graph the $R^2$ of a PLS regression is recovered. If there are just one block associated with a graph and a block without a graph, it's the sum of the squared correlation between the block with a graph and the block without a graph. And in the case where all blocks are associated with graphs, we recovered the average of the mean variances explained by each block),<br/>
    - the graph irregularity (The average of the weight vectors variance within communities. We use the Louvain algorithm to find communities).



## Parameters Description ##

We have four parameters to determine. The $\lambda$ which is the regulation parameter between SGCCA and GraphNet penalty. If $\lambda$=0, there is only sgcca so the blocks interaction and if $\lambda$=1, we have the graphnet penalty. So we want to find and intermediate $\lambda$. $s_1$ and $s_2$ are the quantity of sparsity and $\sigma$ is the eigenvalue filter in order to reduce the high frequencies and give priority to the low frequencies. Because we didn’t have variables to help us to optimize these parameters, we defined three criteria in order to choose these parameters. 

Once the optimal parameters are defined, you can apply the function Laplacian_matrix to define the filtered Laplacian matrix, and then the netSGCCA function with the optimal parameters and the filtered Laplacian matrix. Finally, you have weight vectors and components of each block and you can create networks for example.  

Because of the Laplacian matrix, this method doesn't apply on data in large dimension.


## Bibliography
$^1$ Logan Grosenick, Brad Klingenberg, Kiefer Katovich, Brian Knutson, and Jonathan E. Taylor. Interpretable whole-brain prediction analysis with GraphNet. NeuroImage, 72 :304–321, 2013.<br/>
$^2$ Arthur Tenenhaus, Cathy Philippe, Vincent Guillemot, Kim-Anh Le Cao, Jacques Grill, and Vincent Frouin. Variable selection for generalized canonical correlation analysis. Biostatistics (Oxford, England), 15(3) :569–83, jul 2014.<br/>
$^3$ Franck Rapaport, Andrei Zinovyev, Marie Dutreix, Emmanuel Barillot, and Jean- Philippe Vert. Classification of microarray data using gene networks. BMC bioinformatics, 8 :35, feb 2007.<br/>
