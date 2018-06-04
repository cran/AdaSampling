# AdaSampling
## An R implementation of the AdaSampling algorithm for positive unlabeled and label noise learning

### Description
Implements the AdaSampling procedure, a framework for both positive
unlabeled learning and learning with class label noise, which wraps around a
traditional classifying algorithm. See our publication for details, 
documentation and examples.

### References
* **Yang, P.**, Ormerod, J., Liu, W., Ma, C., Zomaya, A., Yang, J.(2018) AdaSampling for positive-unlabeled and label noise learning with bioinformatics applications. _IEEE Transactions on Cybernetics_, doi:10.1109/TCYB.2018.2816984 

* **Yang, P.**, Liu, W., Yang, J. (2017). Positive unlabeled learning via wrapper-based adaptive 
sampling. Proceedings of the 26th International Joint Conference on Artificial 
Intelligence (IJCAI), 3273-3279. [[fulltext](https://doi.org/10.24963/ijcai.2017/457)]

To install this package, use:
```r
devtools::install_github("PengyiYang/AdaSampling", build_vignettes = TRUE)
library(AdaSampling)
```
Current version of this package includes two functions:

- `adaSample()` applies the AdaSampling procedure to reduce noise in the training set, 
and subsequently trains a classifier from the new training set. 
- `adaSvmBenchmark()` which allows the performance of the AdaSampling procedure (with an SVM 
classifier) to be compared against the performance of the SVM classifier on its own. 

In order to see demonstrations of these two functions, see:
```r
browseVignettes("AdaSampling")
```

### Credits
This developmental version of the AdaSampling package is initially created by Kukulege Dinuka Perera.
