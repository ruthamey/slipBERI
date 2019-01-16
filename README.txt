slipBERI is a matlab slip inversion code, that uses geodetic data (InSAR/GNSS) to solve for earthquake slip at depth in a Bayesian sense. The primary purpose of this function is to incorporate fractal properties of earthquake slip through the incorporation of von Karman autocorrelation prior. There are also some other inversion techniques included in the function, though these are not thoroughly tested.

For information on the methodology please see:
Amey, R. M. J., Hooper, A., & Walters, R. J. (2018). A Bayesian method for incorporating self‐similarity into earthquake slip inversions. Journal of Geophysical Research: Solid Earth, 123, 6052–6071. https://doi.org/10.1029/2017JB015316
or
Amey, Ruth Mary Joy (2018) The Fractal Nature of Fault Slip and Its Incorporation into Earthquake Slip Inversions. PhD thesis, University of Leeds. http://etheses.whiterose.ac.uk/22137/

And please cite Amey, R. M. J., Hooper, A., & Walters, R. J. (2018) if you use this code in your research.

This code is expected to work correctly for von Karman regularisation on one fault strand. If using more than one fault strand, some functions may not work correctly.

* * *

Using slipBERI--

Details on all the inputs required for slipBERI can be found in the 'help' section of slipBERI.m

An example set-up function (make_structure_required_for_slipBERI.m) should be run first, to make structures containing the information required by slipBERI.

slipBERI can then be run in matlab by:
>> slipBERI(fault, data, invert, priors, elastic_params, display, housekeeping )

* * *

Examples--

To run this code using examples, please download:
Amey, Ruth M. J. (2018) Napa Valley Earthquake geodetic slip inversion tools and data. University of Leeds. [Dataset] http://archive.researchdata.leeds.ac.uk/446/
or
Amey, Ruth M. J. (2018) Central Tottori Earthquake geodetic slip inversion tools and data. University of Leeds. [Dataset] http://archive.researchdata.leeds.ac.uk/445/

These datasets contain all the data required to run slipBERI, as well as a function (make_<nameofearthquake>_structure_required_for_slipBERI.m) to be run first, to set up all the inputs required by slipBERI.


* * *
