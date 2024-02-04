#include <image.h>
#include <iostream>

std::vector<DoubleVector> gaussianKernel() {
    std::vector kernel(3, DoubleVector(3));
    kernel[0][0] = 0.0625;
    kernel[0][1] = 0.125;
    kernel[0][2] = 0.0625;
    kernel[1][0] = 0.125;
    kernel[1][1] = 0.25;
    kernel[1][2] = 0.125;
    kernel[2][0] = 0.0625;
    kernel[2][1] = 0.125;
    kernel[2][2] = 0.0625;
    return kernel;
}

std::vector<DoubleVector> prewittHxKernel() {
    std::vector kernel(3, DoubleVector(3));
    kernel[0][0] = -1;
    kernel[0][1] = 0;
    kernel[0][2] = 1;
    kernel[1][0] = -1;
    kernel[1][1] = 0;
    kernel[1][2] = 1;
    kernel[2][0] = -1;
    kernel[2][1] = 0;
    kernel[2][2] = 1;
    return kernel;
}

std::vector<DoubleVector> sobelHxKernel() {
    std::vector kernel(3, DoubleVector(3));
    kernel[0][0] = -1;
    kernel[0][1] = 0;
    kernel[0][2] = 1;
    kernel[1][0] = -2;
    kernel[1][1] = 0;
    kernel[1][2] = 2;
    kernel[2][0] = -1;
    kernel[2][1] = 0;
    kernel[2][2] = 1;
    return kernel;
}

std::vector<DoubleVector> sobelHyKernel() {
    std::vector kernel(3, DoubleVector(3));
    kernel[0][0] = -1;
    kernel[0][1] = -2;
    kernel[0][2] = -1;
    kernel[1][0] = -0;
    kernel[1][1] = 0;
    kernel[1][2] = 0;
    kernel[2][0] = 1;
    kernel[2][1] = 2;
    kernel[2][2] = 1;
    return kernel;
}

std::vector<DoubleVector> sqrtKernel() {
    constexpr auto count = 3;
    std::vector kernel(count, DoubleVector(count));
    for (int i = 0; i < count; i++)
        for (int j = 0; j < count; j++)
            kernel[i][j] = 1.0 / (i * i + j * j + 1);
    return kernel;
}

std::vector<DoubleVector> horizontalBlurKernel(const int count) {
    std::vector kernel(1, DoubleVector(count));
    for (int i = 0; i < count; i++)
        kernel[0][i] = 1.0 / count;
    return kernel;
}

int main() {
    std::cout << "Enter the path to the image: " << std::endl;
    std::string path;
    std::cin >> path;
    const auto image = Image(path.c_str());
    auto fourier = image.toFrequencyDomain();
    bool equalize = false;

    while (true) {
        std::cout << "Apply Gaussian filter (0)" << std::endl;
        std::cout << "Apply Prewitt filter (1)" << std::endl;
        std::cout << "Apply Sobel filter (2)" << std::endl;
        std::cout << "Apply horizontal deblur filter (3)" << std::endl;
        std::cout << "Apply gaussian deblur filter (4)" << std::endl;
        std::cout << "Remove blur with Wiener (5)" << std::endl;
        std::cout << "Equalize histogram (6)" << std::endl;
        std::cout << "Adjust contrast (7)" << std::endl;
        std::cout << "Adjust brightness (8)" << std::endl;
        std::cout << "Save and exit (9)" << std::endl;

        char answer;
        std::cin >> answer;

        switch (answer) {
            case '0':
                fourier.convolveKernel(gaussianKernel());
                break;
            case '1':
                fourier.convolveKernel(prewittHxKernel());
                break;
            case '2':
                fourier.convolveKernel(sobelHxKernel());
                break;
            case '3':
                std::cout << "Enter the desired delta: " << std::endl;
                int delta;
                std::cin >> delta;
                fourier.deconvolveKernel(horizontalBlurKernel(delta));
                break;
            case 4:
                fourier.deconvolveKernel(gaussianKernel());
                break;
            case '5':
                std::cout << "Enter the value of K: " << std::endl;
                double K;
                std::cin >> K;
                fourier.wienerFilter(sqrtKernel(), K);
                break;
            case '6':
                equalize = true;
                break;
            case '7':
                std::cout << "Enter the value of contrast: " << std::endl;
                float contrast;
                std::cin >> contrast;
                fourier.adjustContrast(contrast);
                break;
            case '8':
                std::cout << "Enter the value of brightness: " << std::endl;
                int brightness;
                std::cin >> brightness;
                fourier.addBrightness(brightness);
                break;
            case '9':
                std::cout << "Enter the path to save the image: " << std::endl;
                std::cin >> path;
                image.restoreFFT(fourier, true);
                if (equalize)
                    image.equalizeHistogram();
                image.save(path.c_str(), true);
                return 0;
            default: break;
        }
    }
}
