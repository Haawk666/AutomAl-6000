Guides
---------------------------------------------

Herein we have collected a selection of guides that tries to give a step-by-step introduction to certain tasks that a
user of AutomAl 6000 might want to accomplish.


Installation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

AutomAl 6000 can be run from an executable file for windows, or directly from the source code.

Executable
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

To install AutomAl 6000, simply follow the steps below.

    1. Go to the **download** section of this webpage, and download the zip-file provided under **Executable** to a
    desired location on your computer. (Avoid C:Program files)

    2. Navigate to the downloaded zip-file on your computer, right click and select 'extract all' or something similar,
    depending on your operating system or installed compression software.

    3. (Optional) Navigate into the extracted folder and locate the 'AutomAl6000.exe' -file. Right click and select 'Send to
    -> desktop (create shortcut)', if a desktop shortcut is desired.

.. note::

    We hope to make a proper installer in the future.

Source code
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

Using AutomAl 6000 directly from the source code, takes some extra steps to set up. Some of these steps are explained in
more depth in the "quickstart" tutorial available from the **links** section. All though the steps below are explained
for windows users, the software **might** also work on different systems, but this has not been tested as of yet.

    #. Go to `AutomAl 6000 on GitHub <https://github.com/Haawk666/AutomAl-6000>`_, and from the sub-menu of the **clone** button, select **Download ZIP** and select a suitable location on your computer.

    #. Unzip the folder.

    #. Download and install the latest 64-bit version of `python <https://www.python.org/>`_  for windows, and make sure python.exe is added to your system path. (or just use an already existing install if you have python 3.8 or higher).

    #. From the start menu of windows, type ``cmd``, and hit enter. This should bring up a terminal window.

    #. Install **pipenv** by typing ``pip install pipenv``, and hitting enter.

    #. Navigate to the location of aacc.py by first typing the letter of the drive, followed by a colon and hitting enter. Next, type ``cd`` and the path to the unzipped ``AutomAl-6000-master`` folder. For instance, if aacc.py were located in ``C:\programs\automal6000``, one would first type ``C:`` and hit enter, and then ``cd programs\automal6000``.

    #. Type ``pipenv shell`` and hit enter. This will create a virtual environment for the folder and activate it.

    #. Type ``pipenv install numpy`` and hit enter. This will install numpy as well as all other requirements from the pipfile in the folder.

    #. If all went according to plan (see notes below), you can now type ``python aacc.py`` and hit enter to start AutomAl 6000.

.. note::

    On some networks, pip will not be allowed to download from pypi.org. Try a different network.

.. note::

    If pipenv for some reason does not install all the dependencies from the pipfile, manually install ``numpy``, ``scipy``, ``pillow``, ``matplotlib``, ``h5py`` and ``pyqt5``.

.. note::

    A ``setup.py`` is in the works, which will simplify this process.


Keyboard shortcuts
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
\+                  Toggle currently selected column z-height
P                   Enable move (Enter to accept, P to reset)
R                   Print details about currently selected column
=================   =====================================================   =====================================================

GUI familiarization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The GUI has four main components.

* **Control window**. In the left dock widget, there are buttons, readouts and settings used to interact with the software. The controls are grouped by function, and these groups can be expanded or collapsed by double clicking on the group title.
* **Tab view**. The central widget features several tabs. Each tab offers a different viewing function for the image. The tabs can be cycled with the ``z`` and ``x`` keys.
* **Terminal window**. In the right dock widget, AutomAl 6000 modules will output information about what they are doing.
* **System bar**. At the bottom of the GUI, AutomAl 6000 will report current status. When tasks are being performed, the GUI will often be frozen and the system message will typically be ``working...``. If the GUI is ready for inputs, the system message will be ``Ready.``.


Project workflow with AutomAl's GUI
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There is an implicit logical progression when analyzing images with the AutomAl GUI. When working with images, the
project file will be in certain *stages*, and what state the project is in will affect what you can and/or should do
next. The stages are

    #. A dm3 HAADF-STEM image has been imported as a project file, but no other analysis has taken place yet. the project is in an **initial** state.

    #. Column detection has been performed, and the project now has information about 2D atomic positions. The project is now in a **columns** state.

    #. Column characterization has been applied, and columns now have information about the probability of its own atomic species, its z-position, its neighbours in the opposite and same crystal plane, etc... The project is in a **result** stage.

    #. Manual consideration of the data, and manual corrections and control has been performed by the user. This is the final state, and the project in now in a **control** state.

It is the **control** state that one would use to analyse data, generate models
and/or export data. It is important to note though, that these *states* are only implicit, and is not internally
tracked, and even though the GUI has checks in place to make sure invalid operations are not performed, some of the
software's methods assume a certain state, but can be performed in other states as well, with possibly unpredictable
results. The outline given below, should give a feel for how the GUI is intended to be used.

Initial stage
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

Before importing an .dm3 -file into AutomAl 6000, some pre-processing is needed. Below is an excerpt from my master thesis.

.. Note::

    For the column detection to work optimally, images that are to be analyzed using AutomAl 6000
    should be noise filtered using software such as Gatan Microscopy Suite (GMS) by Gatan [13], which
    is a commercial image processor that is widely used in the TEM-field as a constituent of the standard
    software suites on TEM hardware. Applying an appropriate low pass filter on the Fast Fourier Transform (FFT) of the image will eliminate many of the noise frequencies of the image. Filtering out the
    noise in the image is necessary for column detection to work.

    To apply a low pass filter in GMS, start by performing a FFT on the image. Click on the resulting
    FFT with the band pass tool selected, which will produce a donut shaped mask on the FFT. Adjust
    the inner radius of the mask to zero, and the outer radius to approximately 6,7 nm−1, which will
    include the 200 Al reflection, and exclude the 220 Al reflection. This will eliminate features that are
    smaller than 0.15 nm in real space. Finally, perform inverse FFT on the masked FFT to obtain the
    noise filtered image.

    If the scale of the image is greater than 7 pm/pixel, AutomAl 6000 will automatically upsample the
    image so as to double both the width and height of the image. Using bilinear up-scaling has proven
    to have a positive effect on the column detection in images with scales in this high range. This is
    because the circular samples used in the COM calculations becomes over-granulated (non-circular) for
    low scales. AutomAl 6000 uses the resampling method of Scipy’s ndimage module [14].

Now that we have a pre-processed .dm3 file ready, we can open AutomAl 6000, and press **Project->Import**, which will open an import dialogue.

Once the image has been imported, the **species dictionary** dialogue will appear. You can read more about the species dictionary in a later section.

The project can now be saved with the **Project->Save** button.

.. note::

    Save often!

Column stage
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

We now wish to locate the positions of the columns in the image. To do this, using the built-in centre of mass approach, set the
threshold value **Column detection->Detection threshold value, T** to something like 0,3 (between 0 or 1) and hit **Column detection->Start**. This will produce a pop-up; select *Threshold*
from the drop-down menu and press okay. Column detection will now run for some time (5-15 mins depending on the size
and scale of the image). When it's complete, one should evaluate the result with the **atomic positions** tab. If there are too many or too few
columns detected, reconsider the threshold value and hit **Column detection->Start** as before. Column detection will then either roll back detection or continue detection depending on the new threshold value.
If there is only a handful of missing or superfluous columns, this can be corrected manually by using the buttons **Selected column->New** and **Selected column->Delete**.

.. Note::

    Columns on the very edge of the image will not be considered by the algorithms, so are in effect superfluous.

.. Note::

    Some manual fiddling is almost always necessary. For a typical image one would expect to have to manually set at most 5-10
    columns depending on properties of the image and precipitate. Additionally one might want to slightly adjust some
    positions, especially columns surrounding Cu or other
    bright columns. All this is due to the crudeness of the column detection. In the future other methods that are
    available, like AtoMap might get integrated as an option for import. The column detection algorithm has not been a
    major focus in this work, but it still plays an important part on the end result.

.. Note::

    It is important to get a good result at this stage before proceeding to column characterization, since the quality of the column detection might influebnce the quality of the
    column characterization.


Result stage
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

To produce an atomic overlay, first set the correct alloy type by using the species dictionary dialog, accessed from **Project->Species dict**.
Also, if a different model than the default model is to be used, set this under **Project->Associated model**. Next, select
a column that is inside the Al-matrix, and manually set its species to Al by hitting ``3`` on the keyboard, or use **Selected column->Atomic species**. This will act as a kind of \'seed\' column.
Then, while said column is still selected,
hit **Column characterization->Start** and select *0 - full column characterization*. The algorithm might take anywhere between 5-15 mins, depending on
several factors.

.. Note::

    If no pop-up dialog appears when hitting **Column characterization->Start**, it is because no column is selected, or because no project is open.

One can also selectively do the individual steps of the algorithm by selecting the appropriate step in the pop-up menu.
This allows you to review the results at different stages, if for whatever reason. It is not recommended to do this,
unless the user is familiar with the underlying methods.

These and other available sub-steps can also be useful in the manual sub-processing, see next section.

Control stage
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

After the column characterization has run, manual consideration of the result is needed. There are several built-in
tools to aid in this, of which the *atomic graph*, is the central component. See my master thesis, introductory slides, or the youtube tutorials for details on atomic
graphs and how to interpret them and/or manipulate them.


Build new statistical models
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ''Default model'', which is the default statistical model used by AutomAl 6000, is calculated from a wide range of different images. This general approach is not always the most effective though,
and if extended use of AutomAl 6000 is desired, it might be beneficial to build custom statistical models from your own data. Fortunately, this is fairly easy with
AutomAl 6000's model wizard.

Once you have at minimum 4-5 correctly overlayed and graphed images, you can used these to calculate statistical parameters of a multinomial multivariate normal distribution.
To do this, click **Data->Calculate model**, which will bring up the wizard. In principle, one can use
any nominal attribute and any numerical attributes, but the recommended attributes are **Advanced species** for the nominal attribute, and **alpha min**, **alpha max**,
**theta angle mean**, **normalized avg gamma** and **normalized peak gamma** for the numerical attributes. It is also recommended to exclude **edge columns** with the filter settings.
If the files that are used is properly finalized, recalculating graph parameters should **not** be necessary.

Save the model to a convenient location. You can now apply this model on an image by clicking **Project->Associated model**, and selecting the model you saved. This will now be the
model used by AutomAl 6000 the next time column characterization is run. You can also inspect the details of the model by clicking **Data->Model plots**, and then in the dialog which appears,
click **Select model->Load**.


Generating plots
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

*Coming soon*


Exporting data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Exporting data is easy with the export wizard. Click **Project->Export**, and follow the instructions.


Using core.SuchSoftware as an API without the GUI
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

*Coming soon*


