import cmath
import numpy as np
import matplotlib.pyplot as plt

period = 2 * np.pi
count = 0
# N = 64


class FourierTransformation(object):
    @staticmethod
    def DFT(ValueVector, dir):
        if dir != 1 and dir != -1:
            return 0

        N = len(ValueVector)
        output = [0] * N

        for k in xrange(N):
            for m in xrange(N):
                output[k] += ValueVector[m] * np.exp(-dir * 2j * np.pi * k * m / N)
                # count += 1
            if dir == 1:
                output[k] /= N

        return output

    @staticmethod
    def FFFT(ValueVector, dir):
        result = FourierTransformation.__FFT(ValueVector, dir)
        return FourierTransformation.__FFFTReOrder(result)

    # private

    @staticmethod
    def __FFT(a, dir):
        N = len(a)
        if N == 1:
            return a

        y = [0] * N
        left = [0] * (N / 2)
        right = [0] * (N / 2)
        wN = np.exp(-dir * 2 * np.pi * 1j / N)
        w = 1

        for j in xrange(N / 2):
            left[j] = a[j] + a[j + N / 2]  # LEFT
            right[j] = (a[j] - a[j + N / 2]) * w  # RIGHT
            w *= wN
            # count +=1

        b_left = FourierTransformation.__FFT(left, dir)
        b_right = FourierTransformation.__FFT(right, dir)

        if dir == 1:
            for j in xrange(N / 2):
                y[j] = b_left[j] / 2
                y[j + N / 2] = b_right[j] / 2
        else:
            for j in xrange(N / 2):
                y[j] = b_left[j]
                y[j + N / 2] = b_right[j]
        return y

    @staticmethod
    def __FFFTReOrder(result):
        length = len(result)
        r = 0

        if length <= 2:
            return result

        for x in range(1, length):
            r = FourierTransformation.__rev_next(r, length)
            if r > x:
                result[x], result[r] = result[r], result[x]
        return result

    @staticmethod
    def __rev_next(r, length):
        while True:
            length >>= 1
            r ^= length
            if (r & length) == 0:
                break
        return r


def my_func(x):
    return (np.sin(2 * x) + np.cos(7 * x))

n = 64
x = np.arange(0, period, period / n)
y = my_func(x)

y_dft = FourierTransformation.DFT(y, 1)
dft_abs = map(abs, y_dft)
dft_phase = map(cmath.phase, y_dft)
y_idft = np.real(FourierTransformation.DFT(y_dft, -1))

y_fft = FourierTransformation.FFFT(y, 1)
fft_abs = map(abs, y_fft)
fft_phase = map(cmath.phase, y_fft)
y_ifft = np.real(np.divide(FourierTransformation.FFFT(y_fft, -1), n))

fig, axes = plt.subplots(nrows=2, ncols=4, figsize=(6, 6))

axes[0, 0].plot(x, y, 'r')
axes[0, 1].plot(x, dft_abs, 'g')
axes[0, 2].plot(x, dft_phase, 'blue')
axes[0, 3].plot(x, y_idft, 'black')

axes[1, 0].plot(x, y, 'r')
axes[1, 1].plot(x, fft_abs, 'g')
axes[1, 2].plot(x, fft_phase, 'blue')
axes[1, 3].plot(x, y_ifft, 'black')

plt.show()
# ==============================================