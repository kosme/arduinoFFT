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

#include "arduinoFFT.h"

template <typename T> ArduinoFFT<T>::ArduinoFFT() {}

template <typename T>
ArduinoFFT<T>::ArduinoFFT(T *vReal, T *vImag, uint_fast16_t samples,
                          T samplingFrequency, bool windowingFactors)
    : _samples(samples), _samplingFrequency(samplingFrequency), _vImag(vImag),
      _vReal(vReal) {
  if (windowingFactors) {
    _precompiledWindowingFactors = new T[samples / 2];
  }
  _power = exponent(samples);
#ifdef FFT_SPEED_OVER_PRECISION
  _oneOverSamples = 1.0 / samples;
#endif
}

template <typename T> ArduinoFFT<T>::~ArduinoFFT(void) {
  // Destructor
  if (_precompiledWindowingFactors) {
    delete[] _precompiledWindowingFactors;
  }
}

template <typename T> void ArduinoFFT<T>::complexToMagnitude(void) const {
  complexToMagnitude(this->_vReal, this->_vImag, this->_samples);
}

template <typename T>
void ArduinoFFT<T>::complexToMagnitude(T *vReal, T *vImag,
                                       uint_fast16_t samples) const {
  // vM is half the size of vReal and vImag
  for (uint_fast16_t i = 0; i < samples; i++) {
    vReal[i] = sqrt_internal(sq(vReal[i]) + sq(vImag[i]));
  }
}

template <typename T> void ArduinoFFT<T>::compute(FFTDirection dir) const {
  compute(this->_vReal, this->_vImag, this->_samples, exponent(this->_samples),
          dir);
}

template <typename T>
void ArduinoFFT<T>::compute(T *vReal, T *vImag, uint_fast16_t samples,
                            FFTDirection dir) const {
  compute(vReal, vImag, samples, exponent(samples), dir);
}

// Computes in-place complex-to-complex FFT
template <typename T>
void ArduinoFFT<T>::compute(T *vReal, T *vImag, uint_fast16_t samples,
                            uint_fast8_t power, FFTDirection dir) const {
#ifdef FFT_SPEED_OVER_PRECISION
  T oneOverSamples = this->_oneOverSamples;
  if (!this->_oneOverSamples)
    oneOverSamples = 1.0 / samples;
#endif
  // Reverse bits
  uint_fast16_t j = 0;
  for (uint_fast16_t i = 0; i < (samples - 1); i++) {
    if (i < j) {
      swap(&vReal[i], &vReal[j]);
      if (dir == FFTDirection::Reverse)
        swap(&vImag[i], &vImag[j]);
    }
    uint_fast16_t k = (samples >> 1);

    while (k <= j) {
      j -= k;
      k >>= 1;
    }
    j += k;
  }
  // Compute the FFT
  T c1 = -1.0;
  T c2 = 0.0;
  uint_fast16_t l2 = 1;
  for (uint_fast8_t l = 0; (l < power); l++) {
    uint_fast16_t l1 = l2;
    l2 <<= 1;
    T u1 = 1.0;
    T u2 = 0.0;
    for (j = 0; j < l1; j++) {
      for (uint_fast16_t i = j; i < samples; i += l2) {
        uint_fast16_t i1 = i + l1;
        T t1 = u1 * vReal[i1] - u2 * vImag[i1];
        T t2 = u1 * vImag[i1] + u2 * vReal[i1];
        vReal[i1] = vReal[i] - t1;
        vImag[i1] = vImag[i] - t2;
        vReal[i] += t1;
        vImag[i] += t2;
      }
      T z = ((u1 * c1) - (u2 * c2));
      u2 = ((u1 * c2) + (u2 * c1));
      u1 = z;
    }

#if defined(__AVR__) && defined(USE_AVR_PROGMEM)
    c2 = pgm_read_float_near(&(_c2[l]));
    c1 = pgm_read_float_near(&(_c1[l]));
#else
    T cTemp = 0.5 * c1;
    c2 = sqrt_internal(0.5 - cTemp);
    c1 = sqrt_internal(0.5 + cTemp);
#endif

    if (dir == FFTDirection::Forward) {
      c2 = -c2;
    }
  }
  // Scaling for reverse transform
  if (dir == FFTDirection::Reverse) {
    for (uint_fast16_t i = 0; i < samples; i++) {
#ifdef FFT_SPEED_OVER_PRECISION
      vReal[i] *= oneOverSamples;
      vImag[i] *= oneOverSamples;
#else
      vReal[i] /= samples;
      vImag[i] /= samples;
#endif
    }
  }
  // The computation result at position 0 should be as close to 0 as possible.
  // The DC offset on the signal produces a spike on position 0 that should be
  // eliminated to avoid issues.
  vReal[0] = 0;
}

template <typename T> void ArduinoFFT<T>::dcRemoval(void) const {
  dcRemoval(this->_vReal, this->_samples);
}

template <typename T>
void ArduinoFFT<T>::dcRemoval(T *vData, uint_fast16_t samples) const {
  // calculate the mean of vData
  T mean = 0;
  for (uint_fast16_t i = 0; i < samples; i++) {
    mean += vData[i];
  }
  mean /= samples;
  // Subtract the mean from vData
  for (uint_fast16_t i = 0; i < samples; i++) {
    vData[i] -= mean;
  }
}

template <typename T> T ArduinoFFT<T>::majorPeak(void) const {
  return majorPeak(this->_vReal, this->_samples, this->_samplingFrequency);
}

template <typename T> void ArduinoFFT<T>::majorPeak(T *f, T *v) const {
  majorPeak(this->_vReal, this->_samples, this->_samplingFrequency, f, v);
}

template <typename T>
T ArduinoFFT<T>::majorPeak(T *vData, uint_fast16_t samples,
                           T samplingFrequency) const {
  T frequency;
  majorPeak(vData, samples, samplingFrequency, &frequency, nullptr);
  return frequency;
}

template <typename T>
void ArduinoFFT<T>::majorPeak(T *vData, uint_fast16_t samples,
                              T samplingFrequency, T *frequency,
                              T *magnitude) const {
  T maxY = 0;
  uint_fast16_t IndexOfMaxY = 0;
  findMaxY(vData, (samples >> 1) + 1, &maxY, &IndexOfMaxY);

  T delta = 0.5 * ((vData[IndexOfMaxY - 1] - vData[IndexOfMaxY + 1]) /
                   (vData[IndexOfMaxY - 1] - (2.0 * vData[IndexOfMaxY]) +
                    vData[IndexOfMaxY + 1]));
  T interpolatedX = ((IndexOfMaxY + delta) * samplingFrequency) / (samples - 1);
  if (IndexOfMaxY == (samples >> 1)) // To improve calculation on edge values
    interpolatedX = ((IndexOfMaxY + delta) * samplingFrequency) / (samples);
  // returned value: interpolated frequency peak apex
  *frequency = interpolatedX;
  if (magnitude != nullptr) {
#if defined(ESP8266) || defined(ESP32)
    *magnitude = fabs(vData[IndexOfMaxY - 1] - (2.0 * vData[IndexOfMaxY]) +
                      vData[IndexOfMaxY + 1]);
#else
    *magnitude = abs(vData[IndexOfMaxY - 1] - (2.0 * vData[IndexOfMaxY]) +
                     vData[IndexOfMaxY + 1]);
#endif
  }
}

template <typename T> T ArduinoFFT<T>::majorPeakParabola(void) const {
  T freq = 0;
  majorPeakParabola(this->_vReal, this->_samples, this->_samplingFrequency,
                    &freq, nullptr);
  return freq;
}

template <typename T>
void ArduinoFFT<T>::majorPeakParabola(T *frequency, T *magnitude) const {
  majorPeakParabola(this->_vReal, this->_samples, this->_samplingFrequency,
                    frequency, magnitude);
}

template <typename T>
T ArduinoFFT<T>::majorPeakParabola(T *vData, uint_fast16_t samples,
                                   T samplingFrequency) const {
  T freq = 0;
  majorPeakParabola(vData, samples, samplingFrequency, &freq, nullptr);
  return freq;
}

template <typename T>
void ArduinoFFT<T>::majorPeakParabola(T *vData, uint_fast16_t samples,
                                      T samplingFrequency, T *frequency,
                                      T *magnitude) const {
  T maxY = 0;
  uint_fast16_t IndexOfMaxY = 0;
  findMaxY(vData, (samples >> 1) + 1, &maxY, &IndexOfMaxY);

  *frequency = 0;
  if (IndexOfMaxY > 0) {
    // Assume the three points to be on a parabola
    T a, b, c;
    parabola(IndexOfMaxY - 1, vData[IndexOfMaxY - 1], IndexOfMaxY,
             vData[IndexOfMaxY], IndexOfMaxY + 1, vData[IndexOfMaxY + 1], &a,
             &b, &c);

    // Peak is at the middle of the parabola
    T x = -b / (2 * a);

    // And magnitude is at the extrema of the parabola if you want It...
    if (magnitude != nullptr) {
      *magnitude = (a * x * x) + (b * x) + c;
    }

    // Convert to frequency
    *frequency = (x * samplingFrequency) / samples;
  }
}

template <typename T> uint8_t ArduinoFFT<T>::revision(void) {
  return (FFT_LIB_REV);
}

// Replace the data array pointers
template <typename T>
void ArduinoFFT<T>::setArrays(T *vReal, T *vImag, uint_fast16_t samples) {
  _vReal = vReal;
  _vImag = vImag;
  if (samples) {
    _samples = samples;
#ifdef FFT_SPEED_OVER_PRECISION
    _oneOverSamples = 1.0 / samples;
#endif
    if (_precompiledWindowingFactors) {
      delete[] _precompiledWindowingFactors;
    }
    _precompiledWindowingFactors = new T[samples / 2];
  }
}

template <typename T>
void ArduinoFFT<T>::windowing(FFTWindow windowType, FFTDirection dir,
                              bool withCompensation) {
  // The windowing function is the same, precompiled values can be used, and
  // precompiled values exist
  if (this->_precompiledWindowingFactors && this->_isPrecompiled &&
      this->_windowFunction == windowType &&
      this->_precompiledWithCompensation == withCompensation) {
    windowing(this->_vReal, this->_samples, FFTWindow::Precompiled, dir,
              this->_precompiledWindowingFactors, withCompensation);
    // Precompiled values must be generated. Either the function changed or the
    // precompiled values don't exist
  } else if (this->_precompiledWindowingFactors) {
    windowing(this->_vReal, this->_samples, windowType, dir,
              this->_precompiledWindowingFactors, withCompensation);
    this->_isPrecompiled = true;
    this->_precompiledWithCompensation = withCompensation;
    this->_windowFunction = windowType;
    // Don't care about precompiled windowing values
  } else {
    windowing(this->_vReal, this->_samples, windowType, dir, nullptr,
              withCompensation);
  }
}

template <typename T>
void ArduinoFFT<T>::windowing(T *vData, uint_fast16_t samples,
                              FFTWindow windowType, FFTDirection dir,
                              T *windowingFactors, bool withCompensation) {
  // Weighing factors are computed once before multiple use of FFT
  // The weighing function is symmetric; half the weighs are recorded
  if (windowingFactors != nullptr && windowType == FFTWindow::Precompiled) {
    for (uint_fast16_t i = 0; i < (samples >> 1); i++) {
      if (dir == FFTDirection::Forward) {
        vData[i] *= windowingFactors[i];
        vData[samples - (i + 1)] *= windowingFactors[i];
      } else {
#ifdef FFT_SPEED_OVER_PRECISION
        T inverse = 1.0 / windowingFactors[i];
        vData[i] *= inverse;
        vData[samples - (i + 1)] *= inverse;
#else
        vData[i] /= windowingFactors[i];
        vData[samples - (i + 1)] /= windowingFactors[i];
#endif
      }
    }
  } else {
    T samplesMinusOne = (T(samples) - 1.0);
    T compensationFactor;
    if (withCompensation) {
      compensationFactor =
          _WindowCompensationFactors[static_cast<uint_fast8_t>(windowType)];
    }
    for (uint_fast16_t i = 0; i < (samples >> 1); i++) {
      T indexMinusOne = T(i);
      T ratio = (indexMinusOne / samplesMinusOne);
      T weighingFactor = 1.0;
      // Compute and record weighting factor
      switch (windowType) {
      case FFTWindow::Hamming: // hamming
        weighingFactor = 0.54 - (0.46 * cos(twoPi * ratio));
        break;
      case FFTWindow::Hann: // hann
        weighingFactor = 0.54 * (1.0 - cos(twoPi * ratio));
        break;
      case FFTWindow::Triangle: // triangle (Bartlett)
#if defined(ESP8266) || defined(ESP32)
        weighingFactor =
            1.0 - ((2.0 * fabs(indexMinusOne - (samplesMinusOne / 2.0))) /
                   samplesMinusOne);
#else
        weighingFactor =
            1.0 - ((2.0 * abs(indexMinusOne - (samplesMinusOne / 2.0))) /
                   samplesMinusOne);
#endif
        break;
      case FFTWindow::Nuttall: // nuttall
        weighingFactor = 0.355768 - (0.487396 * (cos(twoPi * ratio))) +
                         (0.144232 * (cos(fourPi * ratio))) -
                         (0.012604 * (cos(sixPi * ratio)));
        break;
      case FFTWindow::Blackman: // blackman
        weighingFactor = 0.42323 - (0.49755 * (cos(twoPi * ratio))) +
                         (0.07922 * (cos(fourPi * ratio)));
        break;
      case FFTWindow::Blackman_Nuttall: // blackman nuttall
        weighingFactor = 0.3635819 - (0.4891775 * (cos(twoPi * ratio))) +
                         (0.1365995 * (cos(fourPi * ratio))) -
                         (0.0106411 * (cos(sixPi * ratio)));
        break;
      case FFTWindow::Blackman_Harris: // blackman harris
        weighingFactor = 0.35875 - (0.48829 * (cos(twoPi * ratio))) +
                         (0.14128 * (cos(fourPi * ratio))) -
                         (0.01168 * (cos(sixPi * ratio)));
        break;
      case FFTWindow::Flat_top: // flat top
        weighingFactor = 0.2810639 - (0.5208972 * cos(twoPi * ratio)) +
                         (0.1980399 * cos(fourPi * ratio));
        break;
      case FFTWindow::Welch: // welch
        weighingFactor = 1.0 - sq((indexMinusOne - samplesMinusOne / 2.0) /
                                  (samplesMinusOne / 2.0));
        break;
      default:
        // This is Rectangle windowing which doesn't do anything
        // and Precompiled which shouldn't be selected
        break;
      }
      if (withCompensation) {
        weighingFactor *= compensationFactor;
      }
      if (windowingFactors) {
        windowingFactors[i] = weighingFactor;
      }
      if (dir == FFTDirection::Forward) {
        vData[i] *= weighingFactor;
        vData[samples - (i + 1)] *= weighingFactor;
      } else {
#ifdef FFT_SPEED_OVER_PRECISION
        T inverse = 1.0 / weighingFactor;
        vData[i] *= inverse;
        vData[samples - (i + 1)] *= inverse;
#else
        vData[i] /= weighingFactor;
        vData[samples - (i + 1)] /= weighingFactor;
#endif
      }
    }
  }
}

// Private functions

template <typename T>
uint_fast8_t ArduinoFFT<T>::exponent(uint_fast16_t value) const {
  // Calculates the base 2 logarithm of a value
  uint_fast8_t result = 0;
  while (value >>= 1)
    result++;
  return result;
}

template <typename T>
void ArduinoFFT<T>::findMaxY(T *vData, uint_fast16_t length, T *maxY,
                             uint_fast16_t *index) const {
  *maxY = 0;
  *index = 0;
  // If sampling_frequency = 2 * max_frequency in signal,
  // value would be stored at position samples/2
  for (uint_fast16_t i = 1; i < length; i++) {
    if ((vData[i - 1] < vData[i]) && (vData[i] > vData[i + 1])) {
      if (vData[i] > vData[*index]) {
        *index = i;
      }
    }
  }
  *maxY = vData[*index];
}

template <typename T>
void ArduinoFFT<T>::parabola(T x1, T y1, T x2, T y2, T x3, T y3, T *a, T *b,
                             T *c) const {
  // const T reversed_denom = 1 / ((x1 - x2) * (x1 - x3) * (x2 - x3));
  // This is a special case in which the three X coordinates are three positive,
  // consecutive integers. Therefore the reverse denominator will always be -0.5
  const T reversed_denom = -0.5;

  *a = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2)) * reversed_denom;
  *b = (x3 * x3 * (y1 - y2) + x2 * x2 * (y3 - y1) + x1 * x1 * (y2 - y3)) *
       reversed_denom;
  *c = (x2 * x3 * (x2 - x3) * y1 + x3 * x1 * (x3 - x1) * y2 +
        x1 * x2 * (x1 - x2) * y3) *
       reversed_denom;
}

template <typename T> void ArduinoFFT<T>::swap(T *a, T *b) const {
  T temp = *a;
  *a = *b;
  *b = temp;
}

#ifdef FFT_SQRT_APPROXIMATION
// Fast inverse square root aka "Quake 3 fast inverse square root", multiplied
// by x. Uses one iteration of Halley's method for precision. See:
// https://en.wikipedia.org/wiki/Methods_of_computing_square_roots#Iterative_methods_for_reciprocal_square_roots
// And: https://github.com/HorstBaerbel/approx
template <typename T> float ArduinoFFT<T>::sqrt_internal(float x) const {
  union // get bits for floating point value
  {
    float x;
    int32_t i;
  } u;
  u.x = x;
  u.i = 0x5f375a86 - (u.i >> 1); // gives initial guess y0.
  float xu = x * u.x;
  float xu2 = xu * u.x;
  // Halley's method, repeating increases accuracy
  u.x = (0.125 * 3.0) * xu * (5.0 - xu2 * ((10.0 / 3.0) - xu2));
  return u.x;
}

template <typename T> double ArduinoFFT<T>::sqrt_internal(double x) const {
  // According to HosrtBaerbel, on the ESP32 the approximation is not faster, so
  // we use the standard function
#ifdef ESP32
  return sqrt(x);
#else
  union // get bits for floating point value
  {
    double x;
    int64_t i;
  } u;
  u.x = x;
  u.i = 0x5fe6ec85e7de30da - (u.i >> 1); // gives initial guess y0.
  double xu = x * u.x;
  double xu2 = xu * u.x;
  // Halley's method, repeating increases accuracy
  u.x = (0.125 * 3.0) * xu * (5.0 - xu2 * ((10.0 / 3.0) - xu2));
  return u.x;
#endif
}
#endif

template <typename T>
const T ArduinoFFT<T>::_WindowCompensationFactors[11] = {
    1.0000000000 * 2.0, // rectangle (Box car)
    1.8549343278 * 2.0, // hamming
    1.8554726898 * 2.0, // hann
    2.0039186079 * 2.0, // triangle (Bartlett)
    2.8163172034 * 2.0, // nuttall
    2.3673474360 * 2.0, // blackman
    2.7557840395 * 2.0, // blackman nuttall
    2.7929062517 * 2.0, // blackman harris
    3.5659039231 * 2.0, // flat top
    1.5029392863 * 2.0, // welch
    // This is added as a precaution, since this index should never be
    // accessed under normal conditions
    1.0 // Custom, precompiled value.
};

template class ArduinoFFT<double>;
template class ArduinoFFT<float>;
