# SimplePlot
A repo collecting simple plotting methods using python

# Structure
The point is to have one simple plotting function, which is found in the main body `simple_plot.py`
Supporting functions can be found in other files, as well as common dictionaries.
This leads to more modular code, with clearer interpretation, and less clutter in the main body.
That way, the main body can be a flexible template for many operations. Additionally,
it becomes easier to extend the code by adding necessary features and dictionaries to other files.

Importantly, the file `XSec.py` should be pulled from the NanoTauAnalysis library as a separate copy
is not maintained here (avoids de-sync errors by forcing a user to pull the up-to-date version from
a second repo they should know about).
