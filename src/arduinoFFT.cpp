/*

	FFT libray
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

#include "arduinoFFT.h"

arduinoFFT::arduinoFFT(void)
{
/* Constructor */
}

arduinoFFT::~arduinoFFT(void)
{
/* Destructor */
}

uint8_t arduinoFFT::Revision(void)
{
	return(FFT_LIB_REV);
}

void arduinoFFT::Compute(double *vReal, double *vImag, uint16_t samples, uint8_t dir)
{
	Compute(vReal, vImag, samples, Exponent(samples), dir);
}

void arduinoFFT::Compute(double *vReal, double *vImag, uint16_t samples, uint8_t power, uint8_t dir)
{
/* Computes in-place complex-to-complex FFT */
	/* Reverse bits */
	uint16_t j = 0;
	for (uint16_t i = 0; i < (samples - 1); i++) {
		if (i < j) {
			Swap(&vReal[i], &vReal[j]);
			Swap(&vImag[i], &vImag[j]);
		}
		uint16_t k = (samples >> 1);
		while (k <= j) {
			j -= k;
			k >>= 1;
		}
		j += k;
	}
	/* Compute the FFT  */
	double c1 = -1.0;
	double c2 = 0.0;
	uint16_t l2 = 1;
	for (uint8_t l = 0; (l < power); l++) {
		uint16_t l1 = l2;
		l2 <<= 1;
		double u1 = 1.0;
		double u2 = 0.0;
		for (j = 0; j < l1; j++) {
			 for (uint16_t i = j; i < samples; i += l2) {
					uint16_t i1 = i + l1;
					double t1 = u1 * vReal[i1] - u2 * vImag[i1];
					double t2 = u1 * vImag[i1] + u2 * vReal[i1];
					vReal[i1] = vReal[i] - t1;
					vImag[i1] = vImag[i] - t2;
					vReal[i] += t1;
					vImag[i] += t2;
			 }
			 double z = ((u1 * c1) - (u2 * c2));
			 u2 = ((u1 * c2) + (u2 * c1));
			 u1 = z;
		}
		c2 = sqrt((1.0 - c1) / 2.0);
		if (dir == FFT_FORWARD) {
			c2 = -c2;
		}
		c1 = sqrt((1.0 + c1) / 2.0);
	}
	/* Scaling for reverse transform */
	if (dir != FFT_FORWARD) {
		for (uint16_t i = 0; i < samples; i++) {
			 vReal[i] /= samples;
			 vImag[i] /= samples;
		}
	}
}

void arduinoFFT::ComplexToMagnitude(double *vReal, double *vImag, uint16_t samples)
{
/* vM is half the size of vReal and vImag */
	for (uint16_t i = 0; i < samples; i++) {
		vReal[i] = sqrt(sq(vReal[i]) + sq(vImag[i]));
	}
}

void arduinoFFT::Windowing(double *vData, uint16_t samples, uint8_t windowType, uint8_t dir)
{
/* Weighing factors are computed once before multiple use of FFT */
/* The weighing function is symetric; half the weighs are recorded */
	double samplesMinusOne = (double(samples) - 1.0);
	for (uint16_t i = 0; i < (samples >> 1); i++) {
		double indexMinusOne = double(i);
		double ratio = (indexMinusOne / samplesMinusOne);
		double weighingFactor = 1.0;
		/* Compute and record weighting factor */
		switch (windowType) {
		case FFT_WIN_TYP_RECTANGLE: /* rectangle (box car) */
			weighingFactor = 1.0;
			break;
		case FFT_WIN_TYP_HAMMING: /* hamming */
			weighingFactor = 0.54 - (0.46 * cos(twoPi * ratio));
			break;
		case FFT_WIN_TYP_HANN: /* hann */
			weighingFactor = 0.54 * (1.0 - cos(twoPi * ratio));
			break;
		case FFT_WIN_TYP_TRIANGLE: /* triangle (Bartlett) */
			weighingFactor = 1.0 - ((2.0 * abs(indexMinusOne - (samplesMinusOne / 2.0))) / samplesMinusOne);
			break;
		case FFT_WIN_TYP_BLACKMAN: /* blackmann */
			weighingFactor = 0.42323 - (0.49755 * (cos(twoPi * ratio))) + (0.07922 * (cos(fourPi * ratio)));
			break;
		case FFT_WIN_TYP_FLT_TOP: /* flat top */
			weighingFactor = 0.2810639 - (0.5208972 * cos(twoPi * ratio)) + (0.1980399 * cos(fourPi * ratio));
			break;
		case FFT_WIN_TYP_WELCH: /* welch */
			weighingFactor = 1.0 - sq((indexMinusOne - samplesMinusOne / 2.0) / (samplesMinusOne / 2.0));
			break;
		}
		if (dir == FFT_FORWARD) {
			vData[i] *= weighingFactor;
			vData[samples - (i + 1)] *= weighingFactor;
		}
		else {
			vData[i] /= weighingFactor;
			vData[samples - (i + 1)] /= weighingFactor;
		}
	}
}

double arduinoFFT::MajorPeak(double *vD, uint16_t samples, double samplingFrequency)
{
	double maxY = 0;
	uint16_t IndexOfMaxY = 0;
	for (uint16_t i = 1; i < ((samples >> 1) - 1); i++) {
		if ((vD[i-1] < vD[i]) && (vD[i] > vD[i+1])) {
			if (vD[i] > maxY) {
				maxY = vD[i];
				IndexOfMaxY = i;
			}
		}
	}
	double delta = 0.5 * ((vD[IndexOfMaxY-1] - vD[IndexOfMaxY+1]) / (vD[IndexOfMaxY-1] - (2.0 * vD[IndexOfMaxY]) + vD[IndexOfMaxY+1]));
	double interpolatedX = ((IndexOfMaxY + delta)  * samplingFrequency) / (samples-1);
	/* retuned value: interpolated frequency peak apex */
	return(interpolatedX);
}

/* Private functions */

void arduinoFFT::Swap(double *x, double *y)
{
	double temp = *x;
	*x = *y;
	*y = temp;
}

uint8_t arduinoFFT::Exponent(uint16_t value)
{
	/* Calculates the base 2 logarithm of a value */
	uint8_t result = 0;
	while (((value >> result) & 1) != 1) result++;
	return(result);
}
