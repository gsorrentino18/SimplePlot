# SimplePlot
A repo collecting and combining simple plotting methods using python

# Quickstart
Probably you will want to use this on your local laptop, not lxplus, to plot things faster (< 5 min).
To do that, there is one line to adjust, and an expected file structure to consider. 
Also, XSec.py must be copied from NanoTauAnalysis.

Your final state skims should be structured in the following way

`pwd`
`/Users/ballmond/LocalDesktop/trigger_gain_plotting/Run3FSSplitSamples`
`ls`
`dimuon ditau  etau   mutau`
`ls ditau/`
`DY                ST                TT                VV`
`Data              Signal            TT_AdditionalSkim WJ`

Where each directory has all the subdirectories above (except dimuon which only has Data and DY).
Final state skims can be copied from 
`/eos/user/b/ballmond/NanoTauAnalysis/analysis/HTauTau_2022_fromstep1_FSskimmed/`
Note: TTbar is a huge background and takes quite a bit of time to load on lxplus, so it gets a special
"AdditionalSkim" directory for use on lxplus. For local usage however, use the normal TT directory so that
your sample can be used for more things (the gain in plotting time is not much locally compared to lxplus).

With the above in place, find the variable `home_dir` in `standard_plot.py`, and adjust it to your local
directory. Now... everything should just work. Below is a test command, followed by several useful commands.

`time python3 standard_plot.py --testing`

The default operation (`python3 standard_plot.py`) should result in full 2022F&G plots of mutau data.
The `--testing` flag tells the plotter to use a subset of MCs to test some standard operations, along with
era G Data only. The output plots are not expected to be good quality, and the flag is truly meant to
test developments before making full plots, which could take several minutes. 

To make full ditau plots
`time python3 standard_plot.py --final_state ditau`
For etau, only era G data is used currently, so the luminosity scaling must be adjusted with the following flag
`time python3 standard_plot.py --final_state etau --lumi "2022 G"`

Another useful flag when testing is `--hide_plots` which hides plots. Finally, if you want to make plots without
yields in the legend, simply add the flag `--hide_yields`. This flag hides the yields.

A somewhat untested option is to adjust the "jet\_mode" being plotted. By default the "Inclusive" setting is used,
but other settings are "0j", "1j", "2j", and "GTE2j" (Greater Than or Equal to 2 Jets).

# Tweaking Code
The main reason I wrote this in python was so that it would be easier for me to control and adjust as I
discover mistakes and need to make different kinds of plots with the same data. As such, a lot that is here
already are support functions that aren't meant to be changed, but file lists, cuts, and plot style should
be adjusted freely.

To adjust the cuts to the ditau final state, `grep "make\_ditau\_cut" *.py` and open `cut_and_study_functions.py`.
Here, you should find the function that makes the ditau cuts where you can freely adjust your parameters. 
This is done similarly for other final states.
Adding new variables to plot can be tricky, and I'll try to write more details about this in the future 
(once i understand it better).

Luminosity should be adjusted through the command line (possibly adding new luminosities if necessary), but
files/processed used can be freely adjusted through the `file_maps` dictionaries in `file_functions.py`.
For example, different file maps are used for testing, dimuon, and full sample plots of 2022 F&G.

To adjust plot properties, check out `MC_dictionary.py` which contains the color keycode, label, and other
information for every MC process. There are additional dummy-keys used for grouping subprocesses, or treating
backgrounds estimated from Data as similar objects. Importantly, the August samples that this plotter was built
with are bugged and have to directly use NWEvents calculated from a separate script. This is only done once,
and the value is stored and accessed in `MC_dictionary.py` for each process. This will be removed in the future.

Anything else, just ask OR try to grep a keyword out to find a function to start with.
Additionally, always feel free to make and add new functions, building off the existing structure.

# Structure -- a bit noodly
The point is to have one simple plotting function, which is found in the main body `standard_plot.py`
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
