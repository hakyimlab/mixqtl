# mixqtl

We provide a unified framework that combines total and allele-specific (at the gene level) read counts and is scalable to large studies with thousands of samples, such as GTEx. 

We develop computationally efficient methods that integrate allele-specific and total read counts to improve cisQTL discovery, fine-mapping, and genetic prediction of the transcriptome. Our log-linear approximation allows scaling-up of computations well-beyond a few hundreds of samples, providing an alternative to current total read count-based methods.

This repository hosts the R implementation of mixQTL, mixFine, and mixPred.

# manuscript 
The manuscript describing the details of the approach can be found [here](https://www.biorxiv.org/content/10.1101/2020.04.22.050666v1)


# Tutorial for applying mixQTL

Find detailed instructions on how to run mixQTL [here](https://github.com/hakyimlab/mixqtl/wiki) via the command line tool or R  functions. 
A computationally efficient GPU-based implementation of mixQTL is embedded in the tensorQTL software [here](https://github.com/broadinstitute/tensorqtl)

