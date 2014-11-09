arduinoFFT
==========

Fast Fourier Transform for Arduino

This is a fork from https://code.google.com/p/makefurt/ which has been abandoned since 2011.

<del>This is a C++ library for Arduino for computing FFT.</del> Now it works both on Arduino and C projects.

Tested on Arduino 1.0.5

Installation on Arduino
--------------------------------------------------------------------------------

To install this library, just place this entire folder as a subfolder in your Arduino installation


When installed, this library should look like:

Arduino\libraries\arduinoFTT              			(this library's folder)
Arduino\libraries\arduinoFTT\arduinoFTT.cpp 			(the library implementation file, uses 32 bits floats vectors)
Arduino\libraries\arduinoFTT\arduinoFTT.h   			(the library description file, uses 32 bits floats vectors)
Arduino\libraries\arduinoFTT\keywords.txt 			(the syntax coloring file)
Arduino\libraries\arduinoFTT\examples     			(the examples in the "open" menu)
Arduino\libraries\arduinoFTT\readme.md   			(this file)

Building on Arduino
--------------------------------------------------------------------------------

After this library is installed, you just have to start the Arduino application.
You may see a few warning messages as it's built.

To use this library in a sketch, go to the Sketch | Import Library menu and
select arduinoFTT.  This will add a corresponding line to the top of your sketch:

`#include <arduinoFTT.h>`

