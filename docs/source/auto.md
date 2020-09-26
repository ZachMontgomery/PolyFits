# Automatic Multivariate Polynomial Fit

Performs a multivariate polynomial curve fit to a dataset and automatically determines which polynomial terms to use based on a balance between the goodness of the fit and a predictve capabilities measure that attempts to make the model compact.

Based on the method given by: Morelli, E. A., "Global Nonlinear Aerodynamic Modeling using Multivariate Orthogonal Functions," Journal of Aircraft, Vol. 32, Issue 2, 1995, pp. 270-277, [https://arc.aiaa.org/doi/abs/10.2514/3.4](https://arc.aiaa.org/doi/abs/10.2514/3.4)

*module* polyFits.**autoPolyFit**(*X*, *y*, *MaxOrder*=12, *tol*=1.e-12, *sigma*=None, *sigmaMultiplier*=1., *verbose*=True):

#### Parameters

* X : *numpy array*

  Array of shape (N,m). X consists of all the independent variables in the dataset. N is the number of data points in the set and m is the number of independent variables

* y : *list* or *numpy array*

  Array with length N. y is the dependent variable values cooresponding to the independent variables in X

* MaxOrder : *integer, optional*

  Gives the max order of polynomial for any one of the independent varialbes to try. Defaults to 12

* tol : *float, optional*

  Gives the cut-off value for any polynomial coefficient to not be included in the final results. If a coefficient has an absolute value below tol, it won't be included. Defaults to 1e-12

* sigma : *float, optional*

  Value used to determine the trade off between how good of a fit to perform and how many terms to keep. Defaults to None, which causes the function to calculate sigma automatically using the mean squared of the difference of the independent variable values with respect to the mean independent variable value of the dataset

* sigmaMultiplier : *float, optional*

  Term multiplied onto sigma to change it's value. Allows using a multiple of the automatically determined sigma value. Defaults to 1.

* verbose : *boolean, optional*

  Determines the verbosity of the function. Defaults to True.

#### Returns

* *list*

  A list of the polynomial coefficients.

* *list*

  A list of the max polynomial orders for each independent variable. The length of this list is therefore m. This list is comparable to the 'Nvec' object used throughout this module

* *float*

  The coefficient of determination, R^2 value, representing the goodness of the fit


