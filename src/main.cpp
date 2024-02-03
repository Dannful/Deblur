#include <image.h>

double **gaussianKernel() {
    const auto kernel = new double *[3];
    for (int i = 0; i < 3; i++)
        kernel[i] = new double[3];
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

double **prewittHxKernel() {
    const auto kernel = new double *[3];
    for (int i = 0; i < 3; i++)
        kernel[i] = new double[3];
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

double **sobelHxKernel() {
    const auto kernel = new double *[3];
    for (int i = 0; i < 3; i++)
        kernel[i] = new double[3];
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

double **sobelHyKernel() {
    const auto kernel = new double *[3];
    for (int i = 0; i < 3; i++)
        kernel[i] = new double[3];
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

double **createDeblurKernel() {
    constexpr int size = 5;  // Size of the kernel
    const auto kernel = new double *[size];
    for (int i = 0; i < size; i++)
        kernel[i] = new double[size];

    for (int i = 0; i < size; i++)
        for (int j = 0; j < size; j++)
            kernel[i][j] = 0;

    for (int i = 0; i < size; i++)
        kernel[i][size / 2] = 1;

    return kernel;
}

int main() {
    const auto image = Image("../files/Car.jpg");
    auto fourier = image.toFrequencyDomain();
    fourier.convolveKernel(sobelHxKernel(), 3, 3);
    image.restoreFFT(fourier);
    image.adjustContrast(2);
    image.save("../files/output.png");
}
