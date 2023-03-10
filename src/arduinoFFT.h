/*

	FFT library
	Copyright (C) 2010 Didier Longueville
	Copyright (C) 2014 Enrique Condes

	This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#ifndef ArduinoFFT_h /* Prevent loading library twice */
#define ArduinoFFT_h
#ifdef ARDUINO
#if ARDUINO >= 100
#include "Arduino.h"
#else
#include "WProgram.h" /* This is where the standard Arduino code lies */
#endif
#else
#include <stdio.h>
#include <stdlib.h>

#ifdef __AVR__
#include <avr/io.h>
#include <avr/pgmspace.h>
#endif
#include "defs.h"
#include "types.h"
#include <math.h>

#endif

// Define this to use a low-precision square root approximation instead of the
// regular sqrt() call
// This might only work for specific use cases, but is significantly faster.
// Only works for ArduinoFFT<float>.
// #define FFT_SQRT_APPROXIMATION

#ifdef FFT_SQRT_APPROXIMATION
#include <type_traits>
#else
#define sqrt_internal sqrt
#endif

enum class FFTDirection { Reverse, Forward };

enum class FFTWindow {
  Rectangle,        // rectangle (Box car)
  Hamming,          // hamming
  Hann,             // hann
  Triangle,         // triangle (Bartlett)
  Nuttall,          // nuttall
  Blackman,         // blackman
  Blackman_Nuttall, // blackman nuttall
  Blackman_Harris,  // blackman harris
  Flat_top,         // flat top
  Welch             // welch
};
#define FFT_LIB_REV 0x15
/* Custom constants */
#define FFT_FORWARD FFTDirection::Forward
#define FFT_REVERSE FFTDirection::Reverse

/* Windowing type */
#define FFT_WIN_TYP_RECTANGLE FFTWindow::Rectangle /* rectangle (Box car) */
#define FFT_WIN_TYP_HAMMING FFTWindow::Hamming     /* hamming */
#define FFT_WIN_TYP_HANN FFTWindow::Hann           /* hann */
#define FFT_WIN_TYP_TRIANGLE FFTWindow::Triangle   /* triangle (Bartlett) */
#define FFT_WIN_TYP_NUTTALL FFTWindow::Nuttall     /* nuttall */
#define FFT_WIN_TYP_BLACKMAN FFTWindow::Blackman   /* blackman */
#define FFT_WIN_TYP_BLACKMAN_NUTTALL                                           \
  FFTWindow::Blackman_Nuttall /* blackman nuttall */
#define FFT_WIN_TYP_BLACKMAN_HARRIS                                            \
  FFTWindow::Blackman_Harris                    /* blackman harris*/
#define FFT_WIN_TYP_FLT_TOP FFTWindow::Flat_top /* flat top */
#define FFT_WIN_TYP_WELCH FFTWindow::Welch      /* welch */
/*Mathematial constants*/
#define twoPi 6.28318531
#define fourPi 12.56637061
#define sixPi 18.84955593

#ifdef __AVR__
static const double _c1[] PROGMEM = {
    0.0000000000, 0.7071067812, 0.9238795325, 0.9807852804, 0.9951847267,
    0.9987954562, 0.9996988187, 0.9999247018, 0.9999811753, 0.9999952938,
    0.9999988235, 0.9999997059, 0.9999999265, 0.9999999816, 0.9999999954,
    0.9999999989, 0.9999999997};
static const double _c2[] PROGMEM = {
    1.0000000000, 0.7071067812, 0.3826834324, 0.1950903220, 0.0980171403,
    0.0490676743, 0.0245412285, 0.0122715383, 0.0061358846, 0.0030679568,
    0.0015339802, 0.0007669903, 0.0003834952, 0.0001917476, 0.0000958738,
    0.0000479369, 0.0000239684};
#endif
class arduinoFFT {
public:
  /* Constructor */
  arduinoFFT(void);
  arduinoFFT(double *vReal, double *vImag, uint16_t samples,
             double samplingFrequency);
  /* Destructor */
  ~arduinoFFT(void);
  /* Functions */
  uint8_t Revision(void);
  uint8_t Exponent(uint16_t value);

  void ComplexToMagnitude(double *vReal, double *vImag, uint16_t samples);
  void Compute(double *vReal, double *vImag, uint16_t samples,
               FFTDirection dir);
  void Compute(double *vReal, double *vImag, uint16_t samples, uint8_t power,
               FFTDirection dir);
  void DCRemoval(double *vData, uint16_t samples);
  double MajorPeak(double *vD, uint16_t samples, double samplingFrequency);
  void MajorPeak(double *vD, uint16_t samples, double samplingFrequency,
                 double *f, double *v);
  void Windowing(double *vData, uint16_t samples, FFTWindow windowType,
                 FFTDirection dir);

  void ComplexToMagnitude();
  void Compute(FFTDirection dir);
  void DCRemoval();
  double MajorPeak();
  void MajorPeak(double *f, double *v);
  void Windowing(FFTWindow windowType, FFTDirection dir);

  double MajorPeakParabola();

private:
  /* Variables */
  uint16_t _samples;
  double _samplingFrequency;
  double *_vReal;
  double *_vImag;
  uint8_t _power;
  /* Functions */
  void Swap(double *x, double *y);
  void Parabola(double x1, double y1, double x2, double y2, double x3,
                double y3, double *a, double *b, double *c);
};

#endif
