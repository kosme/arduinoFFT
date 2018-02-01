/*

	Example of use of the FFT libray to compute FFT for a signal sampled through the ADC.
        Copyright (C) 2017 Enrique Condes

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

arduinoFFT FFT = arduinoFFT(); /* Create FFT object */
/*
These values can be changed in order to evaluate the functions
*/
#define CHANNEL A0
const uint16_t samples = 64; //This value MUST ALWAYS be a power of 2
const double samplingFrequency = 200;

unsigned int delayTime = 0;

/*
These are the input and output vectors
Input vectors receive computed results from FFT
*/
double vReal[samples];
double vImag[samples];

#define SCL_INDEX 0x00
#define SCL_TIME 0x01
#define SCL_FREQUENCY 0x02

void setup()
{
  if(samplingFrequency<=1000)
    delayTime = 1000/samplingFrequency;
  else
    delayTime = 1000000/samplingFrequency;
  Serial.begin(115200);
  Serial.println("Ready");
}

void loop()
{
  for(uint16_t i =0;i<samples;i++)
  {
    vReal[i] = double(analogRead(CHANNEL));
    vImag[i] = 0.0; //Imaginary part must be zeroed in case of looping to avoid wrong calculations and overflows
    if(samplingFrequency<=1000)
      delay(delayTime);
    else
      delayMicroseconds(delayTime);
  }
  /* Print the results of the sampling according to time */
  Serial.println("Data:");
  FFT.PrintSignal(vReal, samples, samplingFrequency);
  FFT.Windowing(vReal, samples, FFT_WIN_TYP_HAMMING, FFT_FORWARD);	/* Weigh data */
  Serial.println("Weighed data:");
  FFT.PrintSignal(vReal, samples, samplingFrequency);
  FFT.Compute(vReal, vImag, samples, FFT_FORWARD); /* Compute FFT */
  Serial.println("Computed Real values:");
  FFT.PrintVector(vReal, samples, samplingFrequency);
  Serial.println("Computed Imaginary values:");
  FFT.PrintVector(vImag, samples, samplingFrequency);
  FFT.ComplexToMagnitude(vReal, vImag, samples); /* Compute magnitudes */
  Serial.println("Computed magnitudes:");
  FFT.PrintSpectrum(vReal, samples, samplingFrequency);
  double x = FFT.MajorPeak(vReal, samples, samplingFrequency);
  Serial.println(x, 6);
  while(1); /* Run Once */
  // delay(2000); /* Repeat after delay */
}
