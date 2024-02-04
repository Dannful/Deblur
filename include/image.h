#ifndef IMAGE_H
#define IMAGE_H

#include <vector>
#include <complex>
#include <cmath>

constexpr unsigned short DESIRED_CHANNELS = 3;

typedef unsigned char PixelValue;

typedef double Real;
typedef std::complex<Real> ComplexNumber;
typedef std::vector<ComplexNumber> ComplexVector;
typedef std::vector<double> DoubleVector;

class FrequencyDomain {
public:
    FrequencyDomain(ComplexVector data[DESIRED_CHANNELS], int originalImageWidth, int originalImageHeight);

    void convolveCircular(const ComplexVector &x);

    void wienerFilter(const std::vector<DoubleVector> &kernel, double K);

    void deconvolveCircular(const ComplexVector &x);

    void convolveKernel(const std::vector<DoubleVector> &kernel);

    void deconvolveKernel(const std::vector<DoubleVector> &kernel);

    void addBrightness(int value);

    void adjustContrast(float value);

    [[nodiscard]] int getImageWidth() const;

    [[nodiscard]] int getImageHeight() const;

    ComplexVector *getData();

    static void fftshift(ComplexVector &x);

    static void bluestein_fft(ComplexVector &x, bool inverse);

private
:
    ComplexVector data[DESIRED_CHANNELS];
    int originalImageWidth;
    int originalImageHeight;

    static void fft(ComplexVector &x);

    static void ifft(ComplexVector &x);

    static void ifftRec(ComplexVector &x);

    static ComplexVector circularConvolution(ComplexVector x, ComplexVector y);

    [[nodiscard]] ComplexVector fourierKernel(const std::vector<DoubleVector> &kernel) const;

};

class Image {
public:
    explicit Image(const char *filename);

    ~Image();

    [[nodiscard]] PixelValue **getPixels() const;

    [[nodiscard]] int getWidth() const;

    [[nodiscard]] int getHeight() const;

    [[nodiscard]] int getChannels() const;

    [[nodiscard]] FrequencyDomain toFrequencyDomain() const;

    void restoreFFT(FrequencyDomain &domain, bool log = false) const;

    void save(const char *filename, bool log = false) const;

    void equalizeHistogram() const;

    [[nodiscard]] std::vector<int> getHistogram() const;

private:
    int width;
    int height;
    int channels;
    PixelValue **pixels;
    const unsigned char *rawData;

    static ComplexVector convertToComplex(const unsigned char *data, int width, int height,
                                          unsigned short channels);
};

#endif
