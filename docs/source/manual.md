# Manual Multivariate Polynomial Fit

This section describes functions to calculate and use the polynomial coefficients for an arbitrary order polynomial curve fit to a dataset with an arbitrary number of independent variables. Curve fits can be performed with full control over the polynomial terms and custom weighting of datapoints. These functions are related to the manual fitting functions described in the publication

    Ullah, A. H., Fabijanic, C., Estevadeordal, J., Montgomery, Z. S., Hunsaker, D. F., Staiger, J. M., and Joo, J. J., "Experimental and Numerical Evaluation of the Performance of Parabolic Flaps," AIAA Aviation 2019 Forum, June 2019,
    [https://arc.aiaa.org/doi/abs/10.2514/6.2019-2916](https://arc.aiaa.org/doi/abs/10.2514/6.2019-2916)

## Routine Listings

* multivariablePolynomialFit

    function for calculating a curve fit to data with an arbitrary number of independent variables

* multivariablePolynomialFunction

    function for calculating a polynomial with an arbitrary number of independent variables

* multivariableR2

    function calculating the coefficient of determination value, or R^2 value, of a polynomial fit of an arbitrary number of independent variables

* multivariableRMS

    function calculating an RMS (root, mean, squared) error and a custom RMSN (root, mean, squared, normalized) error where normalized means the error is divided by the mean of the absolute value of the dependent variables, for a multidimensional polynomial function.
