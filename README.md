# cdnormbio

### Condition-Decomposition Normalization for Biological Applications

**cdnormbio** is an R package that implements Condition-Decomposition methods 
for data normalization, such as MedianCD normalization and SVCD normalization, 
for the analysis of data from high-throughput biology assays. 

For more details please see:  
Roca C.P., Gomes S.I.L., Amorim M.J.B. & Scott-Fordsmand J.J. 
Variation-preserving normalization unveils blind spots in gene expression 
profiling. *Sci. Rep.* **7**, 42460; 
[doi:10.1038/srep42460](http://dx.doi.org/10.1038/srep42460) 
\(2017\). 


## Installation

To install **cdnormbio** from this GitHub repository, 
use the function `install_github` in the 
[devtools](https://cran.r-project.org/package=devtools) package. 

```R
library( devtools )
install_github( "carlosproca/cdnormbio" )
```


## Help

Use the standard help in R.

```R
? normalize.svcd
? normalize.mediancd
```


## Examples

A minimal example:

```R
gene.n <- 1000
sample.n <- 9
expr.data <- matrix( rnorm( gene.n * sample.n ), nrow = gene.n )
expr.condition <- rep( c( 1, 2, 3 ), each = 3 )

normalize.result <- normalize.svcd( expr.data, expr.condition )

expr.data.normalized <- normalize.result$data
```

Another, a bit more complex, example:

```R
gene.n <- 10000
sample.n <- 9
expr.data <- matrix( rnorm( gene.n * sample.n ), nrow = gene.n )
expr.condition <- rep( c( "treatment.a", "treatment.b", "control" ), each = 3 )
offset.added <- rnorm( sample.n )
expr.data <- sweep( expr.data, 2, offset.added, "+" )

normalize.result <- normalize.svcd( expr.data, expr.condition, 
    stdvec.graph = "stdvec", p.value.graph = "p_value", verbose = TRUE )

expr.data.normalized <- normalize.result$data

boxplot( expr.data )
boxplot( expr.data.normalized )

normalization.offset <- normalize.result$offset
sd( normalize.result$offset - offset.added )

normalization.h0.feature <- normalize.result$h0.feature
length( normalize.result$h0.feature )
```

