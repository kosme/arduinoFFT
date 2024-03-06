arduinoFFT
==========

# Fast Fourier Transform for Arduino

This is a fork from https://code.google.com/p/makefurt/ which has been abandoned since 2011.

This is version 2.0 of the library, which has a different [API](#api).

Tested on Arduino 1.8.19 and 2.3.2.

## Installation on Arduino

Use the Arduino Library Manager to install and keep it updated. Just look for arduinoFFT. Only for Arduino 1.5+

## Manual installation on Arduino

To install this library, just place this entire folder as a subfolder in your Arduino installation. When installed, this library should look like:

`Arduino\libraries\arduinoFTT` (this library's folder)
`Arduino\libraries\arduinoFTT\src\arduinoFTT.h` (the library header file. include this in your project)
`Arduino\libraries\arduinoFTT\keywords.txt` (the syntax coloring file)
`Arduino\libraries\arduinoFTT\Examples` (the examples in the "open" menu)
`Arduino\libraries\arduinoFTT\LICENSE` (GPL license file)
`Arduino\libraries\arduinoFTT\README.md` (this file)

## Building on Arduino

After this library is installed, you just have to start the Arduino application.
You may see a few warning messages as it's built.
To use this library in a sketch, go to the Sketch | Import Library menu and
select arduinoFTT.  This will add a corresponding line to the top of your sketch:

`#include <arduinoFTT.h>`

## API

Documentation was moved to the project's [wiki](https://github.com/kosme/arduinoFFT/wiki).
