#include <cmath>
#include <immintrin.h>
const int VECTOR_WIDTH = 8; // avx2寄存器宽度为8个float

void sqrtAVX2(int N, float initialGuess, float values[], float output[]) {
  static const __m256 kThreshold = _mm256_set1_ps(0.00001f);

  for (int i = 0; i < N; i += VECTOR_WIDTH) {
    __m256 x = _mm256_loadu_ps(values + i);      // float x = values[i]
    __m256 guess = _mm256_set1_ps(initialGuess); // float guess = initialGuess;

    // float error = fabs(guess * guess * x - 1.f);
    //这里先对-0.0f
    //按位取反，然后将两个操作数按位与，即将浮点数的符号位设置为0，取其绝对值
    __m256 error = _mm256_andnot_ps(
        _mm256_set1_ps(-0.0f), // 将符号位设为1，其余为0
        _mm256_sub_ps(_mm256_mul_ps(_mm256_mul_ps(guess, guess), x),
                      _mm256_set1_ps(1.f)));

    // while (error > 0.00001f),
    // _mm256_cmp_ps先获取有序大于等于0.00001f的掩码，返回256位掩码值中，如果大于0.00001f，则32位全1，否者32位全0
    // _mm256_movemask_ps提取每个32位的掩码，将其合并成一个整数进行返回
    while (_mm256_movemask_ps(_mm256_cmp_ps(error, kThreshold, _CMP_GT_OS))) {
      // guess = (3.f * guess - x * guess * guess * guess) * 0.5f;
      // 循环中得计算是为了保证向量中得每个元素都达到计算精度，如果向量中某个元素已经提前达到计算精度，继续计算也不影响最终结果正确性
      guess = _mm256_mul_ps(
          _mm256_sub_ps(_mm256_mul_ps(_mm256_set1_ps(3.0f), guess),
                        _mm256_mul_ps(_mm256_mul_ps(x, guess),
                                      _mm256_mul_ps(guess, guess))),
          _mm256_set1_ps(0.5f));
      // error = fabs(guess * guess * x - 1.f);
      error = _mm256_andnot_ps(
          _mm256_set1_ps(-0.0f),
          _mm256_sub_ps(_mm256_mul_ps(_mm256_mul_ps(guess, guess), x),
                        _mm256_set1_ps(1.f)));
    }

    // output[i] = x * guess;
    _mm256_storeu_ps(output + i, _mm256_mul_ps(x, guess));
  }
}