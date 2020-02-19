arduinoFFT
==========

# Fast Fourier Transform for Arduino

This is a fork from https://code.google.com/p/makefurt/ which has been abandoned since 2011.
~~This is a C++ library for Arduino for computing FFT.~~ Now it works both on Arduino and C projects.  
Tested on Arduino 1.6.11 and 1.8.10.

## Installation on Arduino

Use the Arduino Library Manager to install and keep it updated. Just look for arduinoFFT. Only for Arduino 1.5+

## Manual installation on Arduino

To install this library, just place this entire folder as a subfolder in your Arduino installation. When installed, this library should look like:

`Arduino\libraries\arduinoFTT` (this library's folder)  
`Arduino\libraries\arduinoFTT\arduinoFTT.h` (the library header file, uses 32 bit floats or 64bit doubles)  
`Arduino\libraries\arduinoFTT\keywords.txt` (the syntax coloring file)  
`Arduino\libraries\arduinoFTT\examples` (the examples in the "open" menu)  
`Arduino\libraries\arduinoFTT\LICENSE` (GPL license file)  
`Arduino\libraries\arduinoFTT\README.md` (this file)

## Building on Arduino

After this library is installed, you just have to start the Arduino application.
You may see a few warning messages as it's built.

To use this library in a sketch, go to the Sketch | Import Library menu and
select arduinoFTT.  This will add a corresponding line to the top of your sketch:

`#include <arduinoFTT.h>`

## TODO
* Ratio table for windowing function.
* Document windowing functions advantages and disadvantages.
* Optimize usage and arguments.
* Add new windowing functions.
* ~~Spectrum table?~~

## API

* **ArduinoFFT**(T *vReal, T *vImag, uint_fast16_t samples, T samplingFrequency, T * weighingFactors = nullptr);  
Constructor.
The type `T` can be `float` or `double`. `vReal` and `vImag` are pointers to arrays of real and imaginary data and have to be allocated outside of ArduinoFFT. `samples` is the number of samples in `vReal` and `vImag` and `weighingFactors` (if specified). `samplingFrequency` is the sample frequency of the data. `weighingFactors` can optionally be specified to cache weighing factors for the windowing function. This speeds up repeated calls to **windowing()** significantly.

* **~ArduinoFFT**(void);  
Destructor.
* **complexToMagnitude**();  
Convert complex values to their magnitude and store in vReal.
* **compute**(FFTDirection dir);  
Calcuates the Fast Fourier Transform.
* **dcRemoval**();  
Removes the DC component from the sample data.
* **majorPeak**();  
Looks for and returns the frequency of the biggest spike in the analyzed signal.
* **revision**();  
Returns the library revision.
* **windowing**(FFTWindow windowType, FFTDirection dir, bool withCompensation = false);  
Performs a windowing function on the values array. The possible windowing options are:
  * Rectangle
  * Hamming
  * Hann
  * Triangle
  * Nuttall
  * Blackman
  * Blackman_Nuttall
  * Blackman_Harris
  * Flat_top
  * Welch

  If `withCompensation` == true, the following compensation factors are used:
  * Rectangle: 1.0 * 2.0
  * Hamming: 1.8549343278 * 2.0
  * Hann: 1.8554726898 * 2.0
  * Triangle: 2.0039186079 * 2.0
  * Nuttall: 2.8163172034 * 2.0
  * Blackman: 2.3673474360 * 2.0
  * Blackman Nuttall: 2.7557840395 * 2.0
  * Blackman Harris: 2.7929062517 * 2.0
  * Flat top: 3.5659039231 * 2.0
  * Welch: 1.5029392863 * 2.0

## Special flags

You can define these before including arduinoFFT.h:

* #define FFT_SPEED_OVER_PRECISION  
Define this to use reciprocal multiplication for division and some more speedups that might decrease precision.

* #define FFT_SQRT_APPROXIMATION  
Define this to use a low-precision square root approximation instead of the regular sqrt() call. This might only work for specific use cases, but is significantly faster. Only works if `T == float`.
