Blog
---------------------------------------------

**17th international conference on aluminium alloys (ICAA) 2020: AutomAl 6000 poster presentation (TBA)**

The 17th ICAA was converted to an online conference due to the covid situation. I made a brief presentation of AutomAl 6000 (link).
Updating this web-page, commenting source code, writing tutorials and fixing bugs and incomplete features is now a full time job,
as we are slowly starting to promote the software to the community.


**New version (29.09.2020)**

After the master thesis was submitted, I've had the time to do some much needed work on the software. The current
version is available from the source code, and has most features completed. I have not made an executable though, since
this has proved way too time-consuming to get working, but I intend to provide one at a later stage. Main changes:

* Fixed many bugs
* A new scheme for creating statistical models from user data. Check out the youtube tutorials for more details
* Improved column characterization
* Heat map generator (made for a side-project that attempts to identify clustering in 6xxx)
* Can now export svg
* Improved overlay customization
* Automatic bilinear upsampling (the only pre-processing necessary is now the low-pass noise filtering)
* Manually adding or deleting columns now works properly

Compatibility with older save-files has been discontinued with this version, but I intend to uphold backwards compatibility from here on.


**New internal review version in the works (27.02.2020)**

As the work intensifies once again in tandem with the master thesis, there is a new version coming in a few weeks. This
version will carry some changes that are based on some initial feedback on the first internal review version. This
includes, amongst other:

* Fixed some severe bugs and crashes
* A new overlay style customization window
* New keyboard shortcut layouts (See updated webpage)
* New functionality for adding and removing columns manually
* Iterative improvements to the column characterization algorithm

This version will also (somewhat occultly) feature a totally rewritten graph module. This re-write has been done to
better harmonize with the theoretical graph framework that is being developed for the thesis.


**Internal review version available (18.11.2019)**

As promised, I've made an executable available for internal review. Those who were present on the Al-meeting on 15.
Nov 2019, will have received an e-mail. I've termed this version 'alpha 1.0'. The slides from the live demo are included
here in the mail.

Notes on this version:

* I had to disable the plotting and pca module, because matplotlibs's dll's were incurring DLL-hell at build-time. I will try to resolve this for future versions.
* This is not considered a 'public' version! Please see the disclaimers in the slides..


