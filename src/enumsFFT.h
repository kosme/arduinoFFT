#ifndef enumsFFT_h
#define enumsFFT_h
/* Custom constants */
/* These defines keep compatibility with pre 2.0 code */
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
/* End of compatibility defines */

/* Mathematial constants */
#define twoPi 6.28318531
#define fourPi 12.56637061
#define sixPi 18.84955593

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
  Welch,            // welch
  Precompiled       // Placeholder for using custom or precompiled window values
};

enum class FFTDirection { Forward, Reverse };
#endif