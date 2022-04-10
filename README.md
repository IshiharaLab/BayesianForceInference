# BayesianForceInference
Bayesian force inference in epithelial tissue


## Description

This is a script for performing force inference in epithelial mechanics, proposed in [1] and [2].

Likelihood is accounted by force balance equations among cell membrane tensions and cell pressures at vertices.
Prior is given by the expectation that membrane tensions are positive and uniform. Then the tensions and pressures
are inferred by MAP estimation. Empirical Bayes method is adopted where hyper-parameters are determined to minimize
ABIC.

See [1] and [2] for details.

## Requirement

* PySPQR


## Usage

1. Prepare an input file in the same format as the attached sample (Sample/sample/).
2. Change the variable "filename" in EBayes.py to the input file in step 1.
3. Run Forceinf.py on IDE or IPython.

## Licence

[MIT](https://github.com/tcnksm/tool/blob/master/LICENCE)

## References

1. Shuji Ishihara and Kaoru Sugimura <br>
"Bayesian inference of force dynamics during morphogenesis" <br>
Journal of Theoretical Biology 313, p.201-211 (2012)

2. S Ishihara, K Sugimura, SJ Cox, Isabelle Bonnet, Y Bellaïche, F Graner <br>
Comparative study of non-invasive force and stress inference methods in tissue <br>
The European Physical Journal E 36, p.1-13 (2013)
