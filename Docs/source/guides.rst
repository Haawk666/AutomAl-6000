Guides
---------------------------------------------

Herein we have collected a selection of guides that tries to give a step-by-step introduction to certain tasks that a
user of AutomAl 6000 might want to accomplish.

Installation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To install AutomAl 6000, simply follow the steps below.

    1. Go to the **download** section of this webpage, and download the zip-file provided under **Executable** to a
    desired location on your computer.

    2. Navigate to the downloaded zip-file on your computer, right click and select 'extract all' or something similar,
    depending on your operating system or installed compression software.

    3. (Optional) Navigate into the extracted folder and locate the 'aacc.exe' -file. Right click and select 'Send to
    -> desktop (create shortcut)', if a desktop shortcut is desired.

.. note::

    When starting AutomAl, there will be a significant waiting-time (\~20 sec) before the GUI loads. This is because
    the exe will first build its environment in temporary folders, which takes some time. Unfortunately, as a
    consequence of pyinstaller's --onefile option, during this time there will be no indication that the program is
    running, so be patient before clicking several times! In the future, a full-fledged installer is planned, which will
    eliminate this 'quirk'.

Keyboard shortcuts
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These key-bindigs will be applicable from version Alpha 1.1..

=================   =====================================================   =====================================================
Key                 Function                                                Alternate function
=================   =====================================================   =====================================================
W, A, S, D          Move the central viewport                               Move column by pixel increments (if 'move' enabled)
Space-bar           Toggle permute mode
E                   Align views to current viewport
Q                   Zoom in and centre on currently selected column
Z, X                Move left/right between tabs in the central widget
1                   Set currently selected column to 'Si'
2                   Set currently selected column to 'Cu'
3                   Set currently selected column to 'Al'
4                   Set currently selected column to 'Mg'
5
6
7
8
9
\+                  Toggle currently selected column z-height
P                   Enable move (Enter to accept, P to reset)
R                   Print details about currently selected column
=================   =====================================================   =====================================================

Project workflow with AutomAl's GUI
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There is an implicit logical progression when analyzing images with the AutomAl GUI. When working with images, the
project file will be in certain *stages*, and what state the project is in will affect what you can and/or should do
next. The stages are

    #. A dm3 HAADF-STEM image has been imported as a project file, but no other analysis has taken place yet. the project is in an 'initial' state.

    #. Column detection has been performed, and the project now has information about 2D atomic positions. The project is now in a 'column' state.

    #. Column characterization has been applied, and colums now have information about the probability of its own atomic species, its z-position, its neighbours in the opposite and same crystal plane, etc... The project is in a 'result' stage.

    #. Manual consideration of the data, and manual corrections and control has been performed by the user. This is the final state, and the project in now in a 'control' state.

It is the 'control' state that one would use to analyse data, perform principal component analysis, generate plots
and/or export data. It is important to note though, that these 'states' are only implicit, and is not internally
tracked, and even though the GUI has checks in place to make sure invalid operations are not performed, some of the
software's methods assume a certain state, but can be performed in other states as well, with possibly unpredictable
results. The outline given below, should give a feel for how the GUI is intended to be used.

Initial stage
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

Before importing an .dm3 -file into AutomAl, one will usually find it beneficial to prepare the image in a specific
way. Using a program such as digital micrograph (DMG), one should apply a fft-mask to reduce the noise in the image. In
addition, one could use digital micrograph's scaling algorithm to upscale the image if it is small and/or low
magnification (The scale should typically not be any lower than \~6 pm / pixel). These preparation steps will greatly
increase the effectiveness of the column detection algorithm. In the future, these techniques might be included directly
in the software, but for now, pre-processing in DMG is necessary.

.. Note::

    The filetype of .dm3 must be maintained. It is the only file-type currently supported for import, and contains essential
    metadata. For example, when rescaled in DMG, the 'scale' field of the dm3 metadata is correctly and automatically
    updated.

Now that we have a pre-processed .dm3 file ready, we can open AutomAl, and from the 'file' menu select 'new'. Using the
file-dialog, locate the .dm3 and hit 'open'.

With the program there is a sample image included which can be used to get familiar with the software. This file is
called 'sample.dm3', and is already pre-processed, so can be imported directly. Once 'sample.dm3' is imported, the
project instance that the GUI creates, is now in the 'initial' state, and can now be saved as an AutomAl project file using
'file -> save'. This filetype has no exctension. Typically though, one would proceed directly to column detection from this state.

Column stage
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

We now wish to locate the positions of the columns in the image. To do this, using the built-in centre of mass approach, set the
threshold value under **column detection** to something like 0,3 and hit *start*. This will produce a pop-up; select *Threshold*
from the drop-down menu and press okay. Column detection will now run for some time (1-10 mins depending on the size
and scale of the image). When it's complete, one should evaluate the result from the *atomic positions* tab. If there are too many
columns detected, then the process should be reset with a higher threshold value. If however not all columns where
detected, one should lower the threshold value, and press *start* again. This will continue from the current state, and
will pick up progressively darker columns based on the threshold value. Continue this approach until the number of
undetected columns are approachable by manual intervention. In the current version, you must select columns that are on
the very edge, and press *enable move* and use those to cover columns in the precipitate, if you wish to manually complete the column detection.

.. Note::

    Columns on the very edge of the image will not be considered by the algorithms, so are in effect superfluous. Columns
    to be removed can also be just be moved to the very edge.

.. Note::

    Some manual fiddling is almost always necessary. For a typical image one would expect to have to manually set at most 5-10
    columns depending on properties of the image and precipitate. Additionally one might want to slightly adjust some
    positions, especially columns surrounding Cu or other
    bright columns. All this is due to the crudeness of the column detection. In the future other methods that are
    available, like AtomMap might get integrated as an option for column detection. The column detection algorithm has not been a
    major focus in this work, but it still plays an important part on the end result.

.. Note::

    It is important to get a good result at this stage before proceeding to column characterization, since moving
    columns after column characterization has been performed, is not currently supported.


Result stage
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

To produce an atomic overlay, first set the correct alloy type under *Column characterization* -> *Alloy*. Next, select
a column that is inside the Al-matrix, and manually set its species to Al. This will act as a kind of \'seed\' column.
Then, while said column is still selected,
hit *start* and select *0 - full column characterization*. The algorithm might take anywhere between 1-15 mins, depending on
several factors.

.. Note::

    If no pop-up dialog appears when hitting *start*, it is because no column is selected, or because no project is open.

One can also selectively do the individual steps of the algorithm by selecting the appropriate step in the pop-up menu.
This allows you to review the results at different stages, if for whatever reason. It is not recommended to do this,
unless the user is familiar with the underlying methods.

These and other available sub-steps can also be useful in the manual sub-processing, see next section.

Control stage
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

After the column characterization has run, manual consideration of the result is needed. There are several built-in
tools to aid in this, of which the *atomic graph*, is the central component. See [Master thesis] for details on atomic
graphs and how to interpret them, but here is the TL;DR:

Another tool you you can use to consider the result, is the *info-graph*. This shows...

Generating plots
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

*Coming soon*

Performing built-in principle component analysis (PCA)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

*Coming soon*

Exporting data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

*Coming soon*

Testing the accuracy/effectiveness of the algorithms using the validation data-set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

*Coming soon*

Using core.SuchSoftware as an API without the GUI
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

*Coming soon*

Writing plugins for AutomAl 6000
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

*Coming soon*
