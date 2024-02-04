#include <iostream>
#include <algorithm>
#include <bits/stl_algo.h>

#include <image.h>

#define STB_IMAGE_IMPLEMENTATION
#include <stb_image.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>


int Image::getChannels() const {
    return this->channels;
}

int Image::getHeight() const {
    return this->height;
}

int Image::getWidth() const {
    return this->width;
}

PixelValue **Image::getPixels() const {
    return this->pixels;
}


Image::Image(const char *filename): pixels(new PixelValue *[DESIRED_CHANNELS]) {
    width = height = channels = 0;
    this->rawData = stbi_load(filename, &width, &height, &channels, DESIRED_CHANNELS);
    std::cout << "Loading image " << filename << std::endl;
    if (this->rawData == nullptr) {
        std::cerr << "Error: could not load image " << filename << std::endl;
        exit(1);
    }
    std::cout << "Image loaded successfully with " << width << "x" << height << " pixels and " << DESIRED_CHANNELS
            << " channel(s)!" << std::endl;
    for (int channel = 0; channel < DESIRED_CHANNELS; channel++) {
        this->pixels[channel] = new PixelValue[width * height];
        for (int pixel = 0; pixel < width * height; pixel++)
            this->pixels[channel][pixel] = this->rawData[DESIRED_CHANNELS * pixel + channel];
    }
}

Image::~Image() {
    for (int i = 0; i < DESIRED_CHANNELS; i++)
        delete[] pixels[i];
    delete[] pixels;
    stbi_image_free(const_cast<PixelValue *>(this->rawData));
}

ComplexVector
Image::convertToComplex(const PixelValue *data, const int width, const int height,
                        const unsigned short channels = DESIRED_CHANNELS) {
    ComplexVector result(width * height * channels, ComplexNumber(0, 0));
    for (int i = 0; i < width * height * channels; i++)
        result[i] = ComplexNumber(data[i], 0);
    return result;
}

void FrequencyDomain::fft(ComplexVector &x) {
    const auto size = x.size();
    if (size <= 1)
        return;
    ComplexVector even(size / 2), odd(size / 2);
    for (int i = 0; i < size / 2; i++) {
        even[i] = x[2 * i];
        odd[i] = x[2 * i + 1];
    }
    fft(even);
    fft(odd);
    for (int i = 0; i < size / 2; i++) {
        ComplexNumber t = std::exp(ComplexNumber(0, -2 * M_PI * i / static_cast<Real>(size))) * odd[i];
        x[i] = even[i] + t;
        x[i + size / 2] = even[i] - t;
    }
}

void FrequencyDomain::ifftRec(ComplexVector &x) {
    const auto size = x.size();
    if (size <= 1)
        return;
    ComplexVector even(size / 2), odd(size / 2);
    for (int i = 0; i < size / 2; i++) {
        even[i] = x[2 * i];
        odd[i] = x[2 * i + 1];
    }
    ifftRec(even);
    ifftRec(odd);
    for (int i = 0; i < size / 2; i++) {
        ComplexNumber t = std::exp(ComplexNumber(0, 2 * M_PI * i / static_cast<Real>(size))) * odd[i];
        x[i] = even[i] + t;
        x[i + size / 2] = even[i] - t;
    }
}

void FrequencyDomain::ifft(ComplexVector &x) {
    ifftRec(x);
    for (int i = 0; i < x.size(); i++)
        x[i] /= static_cast<Real>(x.size());
}

ComplexVector FrequencyDomain::circularConvolution(ComplexVector x, ComplexVector y) {
    const int n = static_cast<int>(x.size());
    if (n != static_cast<int>(y.size()))
        throw std::invalid_argument("Vectors must have the same size!");
    fft(x);
    fft(y);
    auto result = ComplexVector(n);
    for (int i = 0; i < n; i++)
        result[i] = x[i] * y[i];
    ifft(result);
    return result;
}

Real roundToNDecimals(const Real value, const int decimals) {
    const Real factor = std::pow(10, decimals);
    return std::round(value * factor) / factor;
}

void FrequencyDomain::bluestein_fft(ComplexVector &x, const bool inverse) {
    const size_t n = x.size();
    size_t m = 1;
    while (m / 2 <= n) {
        if (m > SIZE_MAX / 2)
            throw std::length_error("Vector too large for Bluestein FFT!");
        m *= 2;
    }

    ComplexVector expTable(n);
    for (size_t i = 0; i < n; i++) {
        uintmax_t temp = i * i;
        temp %= n * 2;
        double angle = (inverse ? M_PI : -M_PI) * static_cast<Real>(temp) / static_cast<Real>(n);
        expTable[i] = std::polar(1.0, angle);
    }

    ComplexVector avec(m);
    for (size_t i = 0; i < n; i++)
        avec[i] = x[i] * expTable[i];
    ComplexVector bvec(m);
    bvec[0] = expTable[0];
    for (size_t i = 1; i < n; i++)
        bvec[i] = bvec[m - i] = std::conj(expTable[i]);

    const ComplexVector cvec = circularConvolution(std::move(avec), std::move(bvec));

    for (size_t i = 0; i < n; i++) {
        x[i] = cvec[i] * expTable[i] / (inverse ? static_cast<double>(n) : 1.0);
    }
}

ComplexVector *FrequencyDomain::getData() {
    return data;
}

FrequencyDomain::FrequencyDomain(ComplexVector data[DESIRED_CHANNELS], const int originalImageWidth,
                                 const int originalImageHeight) {
    this->originalImageWidth = originalImageWidth;
    this->originalImageHeight = originalImageHeight;
    for (int i = 0; i < DESIRED_CHANNELS; i++)
        this->data[i] = data[i];
}

FrequencyDomain Image::toFrequencyDomain() const {
    std::cout << "Computing FFT..." << std::endl;
    ComplexVector channels[DESIRED_CHANNELS];
    for (int channel = 0; channel < DESIRED_CHANNELS; channel++) {
        auto fourierChannel = ComplexVector(getWidth() * getHeight());
        for (int pixel = 0; pixel < this->width * this->height; pixel++)
            fourierChannel[pixel] = this->rawData[DESIRED_CHANNELS * pixel + channel];
        FrequencyDomain::bluestein_fft(fourierChannel, false);
        channels[channel] = fourierChannel;
    }
    return FrequencyDomain(channels, getWidth(), getHeight());
}

void Image::restoreFFT(FrequencyDomain &domain, const bool log) const {
    if (log)
        std::cout << "Computing IFFT..." << std::endl;
    for (int channel = 0; channel < DESIRED_CHANNELS; channel++) {
        FrequencyDomain::bluestein_fft(domain.getData()[channel], true);
        for (int pixel = 0; pixel < this->width * this->height; pixel++) {
            constexpr Real upperBound = 255.0;
            constexpr Real lowerBound = 0.0;
            pixels[channel][pixel] = static_cast<PixelValue>(std::clamp(
                round(domain.getData()[channel][pixel].real()), lowerBound,
                upperBound));
        }
    }
    if (log)
        std::cout << "IFFT computed successfully!" << std::endl;
}

void FrequencyDomain::convolveCircular(const ComplexVector &x) {
    if (x.size() != getData()->size())
        throw std::invalid_argument("Vector must have the same size as the image!");
    for (int channel = 0; channel < DESIRED_CHANNELS; channel++)
        for (int coef = 0; coef < getData()->size(); coef++)
            getData()[channel][coef] *= x[coef];
}

void FrequencyDomain::wienerFilter(const std::vector<DoubleVector> &kernel, const double K) {
    const ComplexVector H = fourierKernel(kernel);

    ComplexVector H2(H.size());
    for (int i = 0; i < H.size(); i++)
        H2[i] = std::norm(H[i]);

    ComplexVector Wiener(H.size());
    for (int i = 0; i < H.size(); i++)
        Wiener[i] = (1.0 / H[i]) * (H2[i] / (H2[i] + K));

    for (int channel = 0; channel < DESIRED_CHANNELS; channel++)
        for (int i = 0; i < getData()[channel].size(); i++)
            getData()[channel][i] *= Wiener[i];
}

void FrequencyDomain::deconvolveCircular(const ComplexVector &x) {
    if (x.size() != getData()->size())
        throw std::invalid_argument("Vector must have the same size as the image!");
    for (int channel = 0; channel < DESIRED_CHANNELS; channel++)
        for (int coef = 0; coef < getData()->size(); coef++)
            getData()[channel][coef] /= x[coef];
}

int FrequencyDomain::getImageWidth() const {
    return this->originalImageWidth;
}

int FrequencyDomain::getImageHeight() const {
    return this->originalImageHeight;
}


void FrequencyDomain::convolveKernel(const std::vector<DoubleVector> &kernel) {
    const auto newKernel = fourierKernel(kernel);
    convolveCircular(newKernel);
}

void FrequencyDomain::deconvolveKernel(const std::vector<DoubleVector> &kernel) {
    const auto newKernel = fourierKernel(kernel);
    deconvolveCircular(newKernel);
}

void Image::save(const char *filename, const bool log) const {
    if (log)
        std::cout << "Saving image to " << filename << "..." << std::endl;
    const auto data = new PixelValue[DESIRED_CHANNELS * this->width * this->height];
    for (int pixel = 0; pixel < this->width * this->height; pixel++)
        for (int channel = 0; channel < DESIRED_CHANNELS; channel++)
            data[DESIRED_CHANNELS * pixel + channel] = this->pixels[channel][pixel];
    stbi_write_png(filename, this->width, this->height, DESIRED_CHANNELS, data, DESIRED_CHANNELS * this->width);
    delete[] data;
    if (log)
        std::cout << "Image saved successfully!" << std::endl;
}

void Image::equalizeHistogram() const {
    const auto histogram = getHistogram();
    int acHistogram[256];
    const float ratio = 255.0f / static_cast<float>(getWidth() * getHeight());
    acHistogram[0] = static_cast<int>(round(ratio * static_cast<float>(histogram[0])));
    for (int i = 1; i < 256; i++)
        acHistogram[i] = acHistogram[i - 1] + static_cast<int>(round(ratio * static_cast<float>(histogram[i])));

    for (int pixel = 0; pixel < getWidth() * getHeight(); pixel++)
        for (int channel = 0; channel < DESIRED_CHANNELS; channel++)
            pixels[channel][pixel] = static_cast<PixelValue>(std::clamp(
                static_cast<float>(acHistogram[pixels[channel][pixel]]), 0.0f,
                255.0f));
}

void FrequencyDomain::addBrightness(const int value) {
    for (int channel = 0; channel < DESIRED_CHANNELS; channel++)
        getData()[channel][0] += ComplexNumber(getImageWidth() * getImageHeight() * value, 0);
}

void FrequencyDomain::adjustContrast(const float value) {
    for (auto &channel: data)
        for (int pixel = 0; pixel < getImageWidth() * getImageHeight(); pixel++)
            channel[pixel] *= ComplexNumber(value, 0);
}

std::vector<int> Image::getHistogram() const {
    auto histogram = std::vector(256, 0);
    for (int i = 0; i < getWidth() * getHeight(); i++) {
        const double value = DESIRED_CHANNELS == 1
                                 ? pixels[0][i]
                                 : 0.299 * pixels[0][i] + 0.587 * pixels[1][i] + 0.114 * pixels[2][i];
        histogram[static_cast<size_t>(value)]++;
    }
    return histogram;
}

void FrequencyDomain::fftshift(ComplexVector &x) {
    const auto size = x.size();
    for (int i = 0; i < size / 2; i++)
        std::swap(x[i], x[i + size / 2]);
}

ComplexVector FrequencyDomain::fourierKernel(const std::vector<DoubleVector> &kernel) const {
    ComplexVector centeredKernel(getImageWidth() * getImageHeight(), ComplexNumber(0, 0));

    const auto sizeY = kernel.size();
    const auto sizeX = kernel[0].size();
    const int centerX = getImageWidth() / 2;
    const int centerY = getImageHeight() / 2;
    const int kernelCenterX = static_cast<int>(sizeX) / 2;
    const int kernelCenterY = static_cast<int>(sizeY) / 2;

    for (int i = 0; i < sizeY; i++) {
        for (int j = 0; j < sizeX; j++) {
            centeredKernel[(centerY - kernelCenterY + i) * getImageWidth() + centerX - kernelCenterX + j] = kernel[i][
                j];
        }
    }
    fftshift(centeredKernel);
    bluestein_fft(centeredKernel, false);
    return centeredKernel;
}
