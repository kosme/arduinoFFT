/*

        FFT library
        Copyright (C) 2010 Didier Longueville
        Copyright (C) 2014 Enrique Condes
        Copyright (C) 2020 Bim Overbohm (template, speed improvements)

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
#include <stdint.h>
#endif
#include "enumsFFT.h"

// This definition uses a low-precision square root approximation instead of the
// regular sqrt() call
// This might only work for specific use cases, but is significantly faster.

#ifndef FFT_SQRT_APPROXIMATION
  #ifndef sqrt_internal
    #define sqrt_internal sqrt
  #endif
#endif

#define FFT_LIB_REV 0x20

template <typename T> class ArduinoFFT {
public:
  ArduinoFFT();
  ArduinoFFT(T *vReal, T *vImag, uint_fast16_t samples, T samplingFrequency,
             bool windowingFactors = false);

  ~ArduinoFFT();

  void complexToMagnitude(void) const;
  void complexToMagnitude(T *vReal, T *vImag, uint_fast16_t samples) const;

  void compute(FFTDirection dir) const;
  void compute(T *vReal, T *vImag, uint_fast16_t samples,
               FFTDirection dir) const;
  void compute(T *vReal, T *vImag, uint_fast16_t samples, uint_fast8_t power,
               FFTDirection dir) const;

  void dcRemoval(void) const;
  void dcRemoval(T *vData, uint_fast16_t samples) const;

  T majorPeak(void) const;
  void majorPeak(T *f, T *v) const;
  T majorPeak(T *vData, uint_fast16_t samples, T samplingFrequency) const;
  void majorPeak(T *vData, uint_fast16_t samples, T samplingFrequency,
                 T *frequency, T *magnitude) const;

  T majorPeakParabola(void) const;
  void majorPeakParabola(T *frequency, T *magnitude) const;
  T majorPeakParabola(T *vData, uint_fast16_t samples,
                      T samplingFrequency) const;
  void majorPeakParabola(T *vData, uint_fast16_t samples, T samplingFrequency,
                         T *frequency, T *magnitude) const;

  uint8_t revision(void);

  void setArrays(T *vReal, T *vImag, uint_fast16_t samples = 0);

  void windowing(FFTWindow windowType, FFTDirection dir,
                 bool withCompensation = false);
  void windowing(T *vData, uint_fast16_t samples, FFTWindow windowType,
                 FFTDirection dir, T *windowingFactors = nullptr,
                 bool withCompensation = false);

private:
  /* Variables */
  static const T _WindowCompensationFactors[11];
#ifdef FFT_SPEED_OVER_PRECISION
  T _oneOverSamples = 0.0;
#endif
  bool _isPrecompiled = false;
  bool _precompiledWithCompensation = false;
  uint_fast8_t _power = 0;
  T *_precompiledWindowingFactors = nullptr;
  uint_fast16_t _samples;
  T _samplingFrequency;
  T *_vImag;
  T *_vReal;
  FFTWindow _windowFunction;
  /* Functions */
  uint_fast8_t exponent(uint_fast16_t value) const;
  void findMaxY(T *vData, uint_fast16_t length, T *maxY,
                uint_fast16_t *index) const;
  void parabola(T x1, T y1, T x2, T y2, T x3, T y3, T *a, T *b, T *c) const;
  void swap(T *a, T *b) const;

#ifdef FFT_SQRT_APPROXIMATION
  float sqrt_internal(float x) const;
  double sqrt_internal(double x) const;
#endif
};

#if defined(__AVR__) && defined(USE_AVR_PROGMEM)
static const float _c1[] PROGMEM = {
    0.0000000000, 0.7071067812, 0.9238795325, 0.9807852804, 0.9951847267,
    0.9987954562, 0.9996988187, 0.9999247018, 0.9999811753, 0.9999952938,
    0.9999988235, 0.9999997059, 0.9999999265, 0.9999999816, 0.9999999954,
    0.9999999989, 0.9999999997};
static const float _c2[] PROGMEM = {
    1.0000000000, 0.7071067812, 0.3826834324, 0.1950903220, 0.0980171403,
    0.0490676743, 0.0245412285, 0.0122715383, 0.0061358846, 0.0030679568,
    0.0015339802, 0.0007669903, 0.0003834952, 0.0001917476, 0.0000958738,
    0.0000479369, 0.0000239684};
#endif

#endif
