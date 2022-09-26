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

template <class T>
ArduinoFFT<T>::ArduinoFFT(T *vReal, T *vImag, uint_fast16_t samples, T samplingFrequency, T *windowWeighingFactors = nullptr)
	: _vReal(vReal), _vImag(vImag), _samples(samples)
#ifdef FFT_SPEED_OVER_PRECISION
	  ,
	  _oneOverSamples(1.0 / samples)
#endif
	  ,
	  _samplingFrequency(samplingFrequency), _windowWeighingFactors(windowWeighingFactors)
{
	// Calculates the base 2 logarithm of sample count
	_power = 0;
	while (((samples >> _power) & 1) != 1)
	{
		_power++;
	}
}

template <class T>
ArduinoFFT<T>::~ArduinoFFT(void)
{
}

template <class T>
uint8_t ArduinoFFT<T>::revision(void)
{
	return 0x19;
}

template <class T>
void ArduinoFFT<T>::setArrays(T *vReal, T *vImag)
{
	_vReal = vReal;
	_vImag = vImag;
}

template <class T>
void ArduinoFFT<T>::compute(FFTDirection dir) const
{
	// Reverse bits /
	uint_fast16_t j = 0;
	for (uint_fast16_t i = 0; i < (this->_samples - 1); i++)
	{
		if (i < j)
		{
			Swap(this->_vReal[i], this->_vReal[j]);
			if (dir == FFTDirection::Reverse)
			{
				Swap(this->_vImag[i], this->_vImag[j]);
			}
		}
		uint_fast16_t k = (this->_samples >> 1);
		while (k <= j)
		{
			j -= k;
			k >>= 1;
		}
		j += k;
	}
	// Compute the FFT
#ifdef __AVR__
	uint_fast8_t index = 0;
#endif
	T c1 = -1.0;
	T c2 = 0.0;
	uint_fast16_t l2 = 1;
	for (uint_fast8_t l = 0; (l < this->_power); l++)
	{
		uint_fast16_t l1 = l2;
		l2 <<= 1;
		T u1 = 1.0;
		T u2 = 0.0;
		for (j = 0; j < l1; j++)
		{
			for (uint_fast16_t i = j; i < this->_samples; i += l2)
			{
				uint_fast16_t i1 = i + l1;
				T t1 = u1 * this->_vReal[i1] - u2 * this->_vImag[i1];
				T t2 = u1 * this->_vImag[i1] + u2 * this->_vReal[i1];
				this->_vReal[i1] = this->_vReal[i] - t1;
				this->_vImag[i1] = this->_vImag[i] - t2;
				this->_vReal[i] += t1;
				this->_vImag[i] += t2;
			}
			T z = ((u1 * c1) - (u2 * c2));
			u2 = ((u1 * c2) + (u2 * c1));
			u1 = z;
		}
#ifdef __AVR__
		c2 = pgm_read_float_near(&(_c2[index]));
		c1 = pgm_read_float_near(&(_c1[index]));
		index++;
#else
		T cTemp = 0.5 * c1;
		c2 = sqrt_internal(0.5 - cTemp);
		c1 = sqrt_internal(0.5 + cTemp);
#endif
		c2 = dir == FFTDirection::Forward ? -c2 : c2;
	}
	// Scaling for reverse transform
	if (dir != FFTDirection::Forward)
	{
		for (uint_fast16_t i = 0; i < this->_samples; i++)
		{
#ifdef FFT_SPEED_OVER_PRECISION
			this->_vReal[i] *= _oneOverSamples;
			this->_vImag[i] *= _oneOverSamples;
#else
			this->_vReal[i] /= this->_samples;
			this->_vImag[i] /= this->_samples;
#endif
		}
	}
}

template <class T>
void ArduinoFFT<T>::complexToMagnitude() const
{
	// vM is half the size of vReal and vImag
	for (uint_fast16_t i = 0; i < this->_samples; i++)
	{
		this->_vReal[i] = sqrt_internal(sq(this->_vReal[i]) + sq(this->_vImag[i]));
	}
}

template <class T>
void ArduinoFFT<T>::dcRemoval() const
{
	// calculate the mean of vData
	T mean = 0;
	for (uint_fast16_t i = 1; i < ((this->_samples >> 1) + 1); i++)
	{
		mean += this->_vReal[i];
	}
	mean /= this->_samples;
	// Subtract the mean from vData
	for (uint_fast16_t i = 1; i < ((this->_samples >> 1) + 1); i++)
	{
		this->_vReal[i] -= mean;
	}
}

template <class T>
void ArduinoFFT<T>::windowing(FFTWindow windowType, FFTDirection dir, bool withCompensation = false)
{
	// check if values are already pre-computed for the correct window type and compensation
	if (_windowWeighingFactors && _weighingFactorsComputed &&
		_weighingFactorsFFTWindow == windowType &&
		_weighingFactorsWithCompensation == withCompensation)
	{
		// yes. values are precomputed
		if (dir == FFTDirection::Forward)
		{
			for (uint_fast16_t i = 0; i < (this->_samples >> 1); i++)
			{
				this->_vReal[i] *= _windowWeighingFactors[i];
				this->_vReal[this->_samples - (i + 1)] *= _windowWeighingFactors[i];
			}
		}
		else
		{
			for (uint_fast16_t i = 0; i < (this->_samples >> 1); i++)
			{
#ifdef FFT_SPEED_OVER_PRECISION
				// on many architectures reciprocals and multiplying are much faster than division
				T oneOverFactor = 1.0 / _windowWeighingFactors[i];
				this->_vReal[i] *= oneOverFactor;
				this->_vReal[this->_samples - (i + 1)] *= oneOverFactor;
#else
				this->_vReal[i] /= _windowWeighingFactors[i];
				this->_vReal[this->_samples - (i + 1)] /= _windowWeighingFactors[i];
#endif
			}
		}
	}
	else
	{
		// no. values need to be pre-computed or applied
		T samplesMinusOne = (T(this->_samples) - 1.0);
		T compensationFactor = _WindowCompensationFactors[static_cast<uint_fast8_t>(windowType)];
		for (uint_fast16_t i = 0; i < (this->_samples >> 1); i++)
		{
			T indexMinusOne = T(i);
			T ratio = (indexMinusOne / samplesMinusOne);
			T weighingFactor = 1.0;
			// Compute and record weighting factor
			switch (windowType)
			{
			case FFTWindow::Rectangle: // rectangle (box car)
				weighingFactor = 1.0;
				break;
			case FFTWindow::Hamming: // hamming
				weighingFactor = 0.54 - (0.46 * cos(TWO_PI * ratio));
				break;
			case FFTWindow::Hann: // hann
				weighingFactor = 0.54 * (1.0 - cos(TWO_PI * ratio));
				break;
			case FFTWindow::Triangle: // triangle (Bartlett)
				weighingFactor = 1.0 - ((2.0 * abs(indexMinusOne - (samplesMinusOne / 2.0))) / samplesMinusOne);
				break;
			case FFTWindow::Nuttall: // nuttall
				weighingFactor = 0.355768 - (0.487396 * (cos(TWO_PI * ratio))) + (0.144232 * (cos(FOUR_PI * ratio))) - (0.012604 * (cos(SIX_PI * ratio)));
				break;
			case FFTWindow::Blackman: // blackman
				weighingFactor = 0.42323 - (0.49755 * (cos(TWO_PI * ratio))) + (0.07922 * (cos(FOUR_PI * ratio)));
				break;
			case FFTWindow::Blackman_Nuttall: // blackman nuttall
				weighingFactor = 0.3635819 - (0.4891775 * (cos(TWO_PI * ratio))) + (0.1365995 * (cos(FOUR_PI * ratio))) - (0.0106411 * (cos(SIX_PI * ratio)));
				break;
			case FFTWindow::Blackman_Harris: // blackman harris
				weighingFactor = 0.35875 - (0.48829 * (cos(TWO_PI * ratio))) + (0.14128 * (cos(FOUR_PI * ratio))) - (0.01168 * (cos(SIX_PI * ratio)));
				break;
			case FFTWindow::Flat_top: // flat top
				weighingFactor = 0.2810639 - (0.5208972 * cos(TWO_PI * ratio)) + (0.1980399 * cos(FOUR_PI * ratio));
				break;
			case FFTWindow::Welch: // welch
				weighingFactor = 1.0 - sq((indexMinusOne - samplesMinusOne / 2.0) / (samplesMinusOne / 2.0));
				break;
			}
			if (withCompensation)
			{
				weighingFactor *= compensationFactor;
			}
			if (_windowWeighingFactors)
			{
				_windowWeighingFactors[i] = weighingFactor;
			}
			if (dir == FFTDirection::Forward)
			{
				this->_vReal[i] *= weighingFactor;
				this->_vReal[this->_samples - (i + 1)] *= weighingFactor;
			}
			else
			{
#ifdef FFT_SPEED_OVER_PRECISION
				// on many architectures reciprocals and multiplying are much faster than division
				T oneOverFactor = 1.0 / weighingFactor;
				this->_vReal[i] *= oneOverFactor;
				this->_vReal[this->_samples - (i + 1)] *= oneOverFactor;
#else
				this->_vReal[i] /= weighingFactor;
				this->_vReal[this->_samples - (i + 1)] /= weighingFactor;
#endif
			}
		}
		// mark cached values as pre-computed
		_weighingFactorsFFTWindow = windowType;
		_weighingFactorsWithCompensation = withCompensation;
		_weighingFactorsComputed = true;
	}
}

template <class T>
T ArduinoFFT<T>::majorPeak() const
{
	T maxY = 0;
	uint_fast16_t IndexOfMaxY = 0;
	//If sampling_frequency = 2 * max_frequency in signal,
	//value would be stored at position samples/2
	for (uint_fast16_t i = 1; i < ((this->_samples >> 1) + 1); i++)
	{
		if ((this->_vReal[i - 1] < this->_vReal[i]) && (this->_vReal[i] > this->_vReal[i + 1]))
		{
			if (this->_vReal[i] > maxY)
			{
				maxY = this->_vReal[i];
				IndexOfMaxY = i;
			}
		}
	}
	T delta = 0.5 * ((this->_vReal[IndexOfMaxY - 1] - this->_vReal[IndexOfMaxY + 1]) / (this->_vReal[IndexOfMaxY - 1] - (2.0 * this->_vReal[IndexOfMaxY]) + this->_vReal[IndexOfMaxY + 1]));
	T interpolatedX = ((IndexOfMaxY + delta) * this->_samplingFrequency) / (this->_samples - 1);
	if (IndexOfMaxY == (this->_samples >> 1))
	{
		//To improve calculation on edge values
		interpolatedX = ((IndexOfMaxY + delta) * this->_samplingFrequency) / (this->_samples);
	}
	// returned value: interpolated frequency peak apex
	return interpolatedX;
}

template <class T>
void ArduinoFFT<T>::majorPeak(T &frequency, T &value) const
{
	T maxY = 0;
	uint_fast16_t IndexOfMaxY = 0;
	//If sampling_frequency = 2 * max_frequency in signal,
	//value would be stored at position samples/2
	for (uint_fast16_t i = 1; i < ((this->_samples >> 1) + 1); i++)
	{
		if ((this->_vReal[i - 1] < this->_vReal[i]) && (this->_vReal[i] > this->_vReal[i + 1]))
		{
			if (this->_vReal[i] > maxY)
			{
				maxY = this->_vReal[i];
				IndexOfMaxY = i;
			}
		}
	}
	T delta = 0.5 * ((this->_vReal[IndexOfMaxY - 1] - this->_vReal[IndexOfMaxY + 1]) / (this->_vReal[IndexOfMaxY - 1] - (2.0 * this->_vReal[IndexOfMaxY]) + this->_vReal[IndexOfMaxY + 1]));
	T interpolatedX = ((IndexOfMaxY + delta) * this->_samplingFrequency) / (this->_samples - 1);
	if (IndexOfMaxY == (this->_samples >> 1))
	{
		//To improve calculation on edge values
		interpolatedX = ((IndexOfMaxY + delta) * this->_samplingFrequency) / (this->_samples);
	}
	// returned value: interpolated frequency peak apex
	frequency = interpolatedX;
	value = abs(this->_vReal[IndexOfMaxY - 1] - (2.0 * this->_vReal[IndexOfMaxY]) + this->_vReal[IndexOfMaxY + 1]);
}