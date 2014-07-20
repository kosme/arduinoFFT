/*

	Example of use of the FFT libray
	Copyright (C) 2011 Didier Longueville

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

#include "PlainFFT.h"

PlainFFT FFT = PlainFFT(); /* Create FFT object */
/* 
These values can be changed in order to evaluate the functions 
*/
const uint16_t samples = 64;
double signalFrequency = 1000;
double samplingFrequency = 5000;
uint8_t signalIntensity = 100;
/* 
These are the input and output vectors 
Input vectors receive computed results from FFT
*/
double vReal[samples]; 
double vImag[samples];

#define SCL_INDEX 0x00
#define SCL_TIME 0x01
#define SCL_FREQUENCY 0x02

void setup(){
	Serial.begin(115200);
	Serial.println("Ready");
}

void loop() 
{
	/* Build raw data */
	double cycles = (((samples-1) * signalFrequency) / samplingFrequency);
	for (uint8_t i = 0; i < samples; i++) {
		vReal[i] = uint8_t((signalIntensity * (sin((i * (6.2831 * cycles)) / samples) + 1.0)) / 2.0);
	}
	PrintVector(vReal, samples, SCL_TIME);
	FFT.Windowing(vReal, samples, FFT_WIN_TYP_HAMMING, FFT_FORWARD);	/* Weigh data */
	PrintVector(vReal, samples, SCL_TIME);
	FFT.Compute(vReal, vImag, samples, FFT_FORWARD); /* Compute FFT */
	PrintVector(vReal, samples, SCL_INDEX);
	PrintVector(vImag, samples, SCL_INDEX);
	FFT.ComplexToMagnitude(vReal, vImag, samples); /* Compute magnitudes */
	PrintVector(vReal, (samples >> 1), SCL_FREQUENCY);	
	double x = FFT.MajorPeak(vReal, samples, samplingFrequency);
	Serial.println(x, 6);
	while(1); /* Run Once */
	// delay(2000); /* Repeat after delay */
}

void PrintVector(double *vData, uint8_t bufferSize, uint8_t scaleType) 
{	
	for (uint16_t i = 0; i < bufferSize; i++) {
		double abscissa;
		/* Print abscissa value */
		switch (scaleType) {
		case SCL_INDEX:
			abscissa = (i * 1.0);
			break;
		case SCL_TIME:
			abscissa = ((i * 1.0) / samplingFrequency);
			break;
		case SCL_FREQUENCY:
			abscissa = ((i * 1.0 * samplingFrequency) / samples);
			break;
		}
		Serial.print(abscissa, 6);
		Serial.print(" ");
		Serial.print(vData[i], 4);
		Serial.println();
	}
	Serial.println();
}