# BNP-For-Unmeasured-Confounding: A Bayesian Nonparametric Method to Adjust for Unmeasured Confounding with Negative Controls

### Description

Unmeasured confounding bias is among the largest threats to the validity of observational studies.  Negative control variables have been introduced in the causal inference literature as a promising approach to  account for unmeasured confounding bias by using auxiliary information.   In this repo, we provide tools based on a Bayesian nonparametric method to estimate a causal exposure-response function (CERF) that, under assumptions,  leverages information from negative control variables to address unmeasured confounding. We model the CERF as a  mixture of linear models. A mixture model has the nice feature of capturing the potential nonlinearity of the shape of the CERF while maintaining computational efficiency and leveraging closed-form results that are available under the linear assumption.  

In the *src* folder, we provide functions to fit a gaussian linear mixed model based on Bayesian nonparametric methods with and without using negative control exposures and outcomes.

In the *simulation*, we assess the performance of this method by simulation studies and demonstrate that it can recover the true shape of the CERF in the presence of unmeasured confounding.

In the the *pm25-cardiovasculardisease-exposure* folder We apply this new estimation procedure to adjust for a potential unmeasured confounder when evaluating the relationship between long-term exposure to ambient $PM_{2.5}$ and cardiovascular hospitalization rates among the elderly in the continental U.S. We implement this method in open-source software for reproducibility.

### Getting started
...

### Contact
We welcome contributions and feedback about BNP-For-Unmeasured-Confounding. If you have any suggestions, please open an issue or submit a pull request.

### Documentation
The companion paper is hosted at ...

