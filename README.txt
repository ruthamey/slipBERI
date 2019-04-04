slipBERI is a matlab slip inversion code, that uses geodetic data (InSAR/GNSS) to solve for earthquake slip at depth in a Bayesian sense. The primary purpose of this function is to incorporate fractal properties of earthquake slip through the incorporation of von Karman autocorrelation prior. There are also some other inversion techniques included in the function, though these are not thoroughly tested.

For information on the methodology please see:
Amey, R. M. J., Hooper, A., & Walters, R. J. (2018). A Bayesian method for incorporating self‐similarity into earthquake slip inversions. Journal of Geophysical Research: Solid Earth, 123, 6052–6071. https://doi.org/10.1029/2017JB015316
or
Amey, Ruth Mary Joy (2018) The Fractal Nature of Fault Slip and Its Incorporation into Earthquake Slip Inversions. PhD thesis, University of Leeds. http://etheses.whiterose.ac.uk/22137/
or for the trans-dimensional approach
Amey, R. M. J., Hooper, A., & Morishita, Y. ( 2019). Going to any lengths: Solving for fault size and fractal slip for the 2016, Mw 6.2 Central Tottori earthquake, Japan, using a transdimensional inversion scheme. Journal of Geophysical Research: Solid Earth, 124. https://doi.org/10.1029/2018JB016434 
 

And please cite Amey, R. M. J., Hooper, A., & Walters, R. J. (2018) if you use this code in your research, or  Amey, R. M. J., Hooper, A., & Morishita, Y.  (2019) if using the trans-dimensional approach.

This code is expected to work correctly for von Karman regularisation on one fault strand. If using more than one fault strand, some functions may not work correctly.

* * *

Using slipBERI--

Details on all the inputs required for slipBERI can be found in the 'help' section of slipBERI.m or a reader-friendly version can be found here: https://docs.google.com/document/d/1cUXLRxN-oB8Q8kGOueq2c-Zxr3W1vDgWpGUw1MAEx5s/edit#

An example set-up function (make_structure_required_for_slipBERI.m) should be run first, to make structures containing the information required by slipBERI, which can then be edited as appropriate.

slipBERI can then be run in matlab by:
>> slipBERI(fault, data, invert, priors, elastic_params, display, housekeeping )

* * *

Examples--

Two examples can be found in the 'examples' folder.
These datasets contain all the data required to run slipBERI, as well as a function (make_<nameofearthquake>_structure_required_for_slipBERI.m) to be run first, to set up all the inputs required by slipBERI.

1) 2014 Mw 6.0 Napa Valley Earthquake
Details can be found in Amey et al., 2018 (see above)

3) 2016 Mw 6.2 Central Tottori Earthquake
Details can be found in Amey et al., 2019 (see above)



* * *
