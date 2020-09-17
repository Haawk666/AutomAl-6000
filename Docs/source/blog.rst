Dev-blog
---------------------------------------------

New internal review version in the works (27.02.2020)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

Haakon

Internal review version available (18.11.2019)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As promised, I've made an executable available for internal review. Those who were present on the Al-meeting on 15.
Nov 2019, will have received an e-mail. I've termed this version 'alpha 1.0'. The slides from the live demo are included
here in the mail.

Notes on this version:

* I had to disable the plotting and pca module, because matplotlibs's dll's were incurring DLL-hell at build-time. I will try to resolve this for future versions.
* This is not considered a 'public' version! Please see the disclaimers in the slides..

Haakon

