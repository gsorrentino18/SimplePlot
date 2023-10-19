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

To develop in this script, I normally write and define new functions in the main body, and then
move them to a relevant file when they are sufficiently mature and it is evident they could be repurposed.
I try to organize functions by their name, and I avoid using abbreviations in functions. Additionally,
I always write what is imported from where in the most explicit way possible (i.e. no `from x import *` ).
I find doing this leads to more organized code, with clear lines from function call to function implementation.
Finally, since this is meant to be a simple library, I am avoiding using classes. Although they may
have more functionality and modular organization, I find they are overwrought for simply plotting in python.
Maybe this will change as the library evolves.
