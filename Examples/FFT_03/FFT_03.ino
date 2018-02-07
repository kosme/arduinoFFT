/*

	Example of use of the FFT libray to compute FFT for a signal sampled through the ADC.
        Copyright (C) 2018 Enrique Condés and Ragnar Ranøyen Homb

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
const uint16_t SAMPLES = 64; //This value MUST ALWAYS be a power of 2
const double samplingFrequency = 100; //Hz, must be less than 10000 due to ADC

unsigned int sampling_period_us;
unsigned long microseconds;

/*
These are the input and output vectors
Input vectors receive computed results from FFT
*/
double vReal[samples];
double vImag[samples];

void setup()
{
  sampling_period_us = round(1000000*(1.0/SAMPLING_FREQUENCY));
  Serial.begin(115200);
  Serial.println("Ready");
}

void loop()
{
  /*SAMPLING*/
  for(int i=0; i<SAMPLES; i++)
  {
      microseconds = micros();    //Overflows after around 70 minutes!

      vReal[i] = analogRead(CHANNEL);
      vImag[i] = 0;
      while(micros() < (microseconds + sampling_period_us)){
        //empty loop
      }
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
  Serial.println(x, 6); //Print out what frequency is the most dominant.
  while(1); /* Run Once */
  // delay(2000); /* Repeat after delay */
}
