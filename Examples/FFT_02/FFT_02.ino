/*

	Example of use of the FFT libray to compute FFT for several signals over a range of frequencies.
        The exponent is calculated once before the excecution since it is a constant.
        This saves resources during the excecution of the sketch and reduces the compiled size.
        The sketch shows the time that the computing is taking.
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

arduinoFFT FFT = arduinoFFT(); /* Create FFT object */
/*
These values can be changed in order to evaluate the functions
*/

const uint16_t samples = 64;
const double sampling = 40;
const uint8_t amplitude = 4;
uint8_t exponent;
const double startFrequency = 2;
const double stopFrequency = 16.4;
const double step_size = 0.1;

/*
These are the input and output vectors
Input vectors receive computed results from FFT
*/
double vReal[samples];
double vImag[samples];

unsigned long time;

void setup()
{
  Serial.begin(115200);
  Serial.println("Ready");
  exponent = FFT.Exponent(samples);
}

void loop()
{
  Serial.println("Frequency\tDetected\ttakes (ms)");
  Serial.println("=======================================\n");
  for(double frequency = startFrequency; frequency<=stopFrequency; frequency+=step_size)
  {
    /* Build raw data */
    double cycles = (((samples-1) * frequency) / sampling);
    for (uint16_t i = 0; i < samples; i++)
    {
      vReal[i] = int8_t((amplitude * (sin((i * (twoPi * cycles)) / samples))) / 2.0);
      vImag[i] = 0; //Reset the imaginary values vector for each new frequency
    }
    /*Serial.println("Data:");
    FFT.PrintSignal(vReal, samples, samplingFrequency);*/
    time=millis();
    FFT.Windowing(vReal, samples, FFT_WIN_TYP_HAMMING, FFT_FORWARD);	/* Weigh data */
    /*Serial.println("Weighed data:");
    FFT.PrintSignal(vReal, samples, samplingFrequency);*/
    FFT.Compute(vReal, vImag, samples, exponent, FFT_FORWARD); /* Compute FFT */
    /*Serial.println("Computed Real values:");
    FFT.PrintVector(vReal, samples, samplingFrequency);
    Serial.println("Computed Imaginary values:");
    FFT.PrintVector(vImag, samples, samplingFrequency);*/
    FFT.ComplexToMagnitude(vReal, vImag, samples); /* Compute magnitudes */
    /*Serial.println("Computed magnitudes:");
    FFT.PrintSpectrum(vReal, samples, samplingFrequency);*/
    double x = FFT.MajorPeak(vReal, samples, sampling);
    Serial.print(frequency);
    Serial.print(": \t\t");
    Serial.print(x, 4);
    Serial.print("\t\t");
    Serial.print(millis()-time);
    Serial.println(" ms");
    // delay(2000); /* Repeat after delay */
  }
  while(1); /* Run Once */
}
