#include "CS149intrin.h"
#include "logger.h"
#include <algorithm>
#include <getopt.h>
#include <math.h>
#include <stdio.h>

#include "../common/CycleTimer.h"

using namespace std;

#define EXP_MAX 10 // 生成随机的指数值

Logger CS149Logger;

void usage(const char *progname);
void initValue(float *values, int *exponents, float *output, float *gold,
               unsigned int N);
void absSerial(float *values, float *output, int N);
void absVector(float *values, float *output, int N);
void clampedExpSerial(float *values, int *exponents, float *output, int N);
void clampedExpVector(float *values, int *exponents, float *output, int N);
float arraySumSerial(float *values, int N);
float arraySumVector(float *values, int N);
bool verifyResult(float *values, int *exponents, float *output, float *gold,
                  int N);

int main(int argc, char *argv[]) {
  int N = 16;
  bool printLog = false;

  // parse commandline options ////////////////////////////////////////////
  int opt;
  static struct option long_options[] = {{"size", 1, 0, 's'},
                                         {"log", 0, 0, 'l'},
                                         {"help", 0, 0, '?'},
                                         {0, 0, 0, 0}};

  while ((opt = getopt_long(argc, argv, "s:l?", long_options, NULL)) != EOF) {

    switch (opt) {
    case 's':
      N = atoi(optarg);
      if (N <= 0) {
        printf("Error: Workload size is set to %d (<0).\n", N);
        return -1;
      }
      break;
    case 'l':
      printLog = true;
      break;
    case '?':
    default:
      usage(argv[0]);
      return 1;
    }
  }

  float *values = new float[N + VECTOR_WIDTH];
  int *exponents = new int[N + VECTOR_WIDTH];
  float *output = new float[N + VECTOR_WIDTH];
  float *gold = new float[N + VECTOR_WIDTH];
  initValue(values, exponents, output, gold, N);

  double startTime;
  double endTime;

  // 基本上只能算一个模拟
  // startTime = CycleTimer::currentSeconds();
  clampedExpSerial(values, exponents, gold, N);
  // endTime = CycleTimer::currentSeconds();
  // printf("serial time:\t\t[%f] ms\n", (endTime - startTime));

  // startTime = CycleTimer::currentSeconds();
  clampedExpVector(values, exponents, output, N);
  // endTime = CycleTimer::currentSeconds();
  // printf("vector time:\t\t[%f] ms\n", (endTime - startTime));

  printf("\e[1;31mCLAMPED EXPONENT\e[0m (required) \n");
  bool clampedCorrect = verifyResult(values, exponents, output, gold, N);
  if (printLog)
    CS149Logger.printLog();
  CS149Logger.printStats();

  printf("************************ Result Verification "
         "*************************\n");
  if (!clampedCorrect) {
    printf("@@@ Failed!!!\n");
  } else {
    printf("Passed!!!\n");
  }

  printf("\n\e[1;31mARRAY SUM\e[0m (bonus) \n");
  if (N % VECTOR_WIDTH == 0) {
    float sumGold = arraySumSerial(values, N);
    float sumOutput = arraySumVector(values, N);
    float epsilon = 0.1;
    bool sumCorrect = abs(sumGold - sumOutput) < epsilon * 2;
    if (!sumCorrect) {
      printf("Expected %f, got %f\n.", sumGold, sumOutput);
      printf("@@@ Failed!!!\n");
    } else {
      printf("Passed!!!\n");
    }
  } else {
    printf("Must have N %% VECTOR_WIDTH == 0 for this problem (VECTOR_WIDTH is "
           "%d)\n",
           VECTOR_WIDTH);
  }

  delete[] values;
  delete[] exponents;
  delete[] output;
  delete[] gold;

  return 0;
}

void usage(const char *progname) {
  printf("Usage: %s [options]\n", progname);
  printf("Program Options:\n");
  printf("  -s  --size <N>     Use workload size N (Default = 16)\n");
  printf("  -l  --log          Print vector unit execution log\n");
  printf("  -?  --help         This message\n");
}

// 生成随机测试数据
void initValue(float *values, int *exponents, float *output, float *gold,
               unsigned int N) {

  for (unsigned int i = 0; i < N + VECTOR_WIDTH; i++) {
    // random input values
    values[i] = -1.f + 4.f * static_cast<float>(rand()) / RAND_MAX;
    exponents[i] = rand() % EXP_MAX; // 指数数组
    output[i] = 0.f;                 // 初始化输出结果缓冲区
    gold[i] = 0.f;                   // 初始化预期结果缓冲区
  }
}

// 验证向量化计算结果与串行结果的一致性
bool verifyResult(float *values, int *exponents, float *output, float *gold,
                  int N) {
  int incorrect = -1;
  float epsilon = 0.00001;
  for (int i = 0; i < N + VECTOR_WIDTH; i++) {
    if (abs(output[i] - gold[i]) > epsilon) {
      incorrect = i;
      break;
    }
  }

  if (incorrect != -1) {
    if (incorrect >= N)
      printf("You have written to out of bound value!\n");

    printf("Wrong calculation at value[%d]!\n", incorrect);

    printf("value  = ");

    for (int i = 0; i < N; i++) {
      printf("% f ", values[i]);
    }
    printf("\n");

    printf("exp    = ");
    for (int i = 0; i < N; i++) {
      printf("% 9d ", exponents[i]);
    }
    printf("\n");

    printf("output = ");
    for (int i = 0; i < N; i++) {
      printf("% f ", output[i]);
    }
    printf("\n");

    printf("gold   = ");
    for (int i = 0; i < N; i++) {
      printf("% f ", gold[i]);
    }
    printf("\n");
    return false;
  } //  if(incorrect != -1)

  printf("Results matched with answer!\n");
  return true;
}

// computes the absolute value of all elements in the input array
// values, stores result in output
// 串行化计算数组的绝对值
void absSerial(float *values, float *output, int N) {
  for (int i = 0; i < N; i++) {
    float x = values[i];
    if (x < 0) {
      output[i] = -x;
    } else {
      output[i] = x;
    }
  }
}

// implementation of absSerial() above, but it is vectorized using CS149
// intrinsics
// 向量化计算数组元素的绝对值，使用CS149向量指令集
void absVector(float *values, float *output, int N) {
  __cs149_vec_float x;
  __cs149_vec_float result;
  __cs149_vec_float zero = _cs149_vset_float(0.f);
  __cs149_mask maskAll, maskIsNegative, maskIsNotNegative;

  //  Note: Take a careful look at this loop indexing.  This example
  //  code is not guaranteed to work when (N % VECTOR_WIDTH) != 0.
  //  Why is that the case?
  for (int i = 0; i < N; i += VECTOR_WIDTH) {

    // All ones
    maskAll = _cs149_init_ones();

    // All zeros
    maskIsNegative = _cs149_init_ones(0);

    // Load vector of values from contiguous memory addresses
    _cs149_vload_float(x, values + i, maskAll); // x = values[i];

    // Set mask according to predicate
    // 判断x的值是否小于0，如果是则将maskIsNegative对应位置设置为1
    _cs149_vlt_float(maskIsNegative, x, zero, maskAll); // if (x < 0) {

    // Execute instruction using mask ("if" clause)
    // 根据掩码的有效位，对两个向量进行剑法操作（根据掩码有效位得到对应的-x值）
    _cs149_vsub_float(result, zero, x, maskIsNegative); //   output[i] = -x;

    // Inverse maskIsNegative to generate "else" mask
    // maskIsNotNegative = ~maskIsNegative;
    maskIsNotNegative = _cs149_mask_not(maskIsNegative); // } else {

    // Execute instruction ("else" clause)
    // 根据mask的掩码有效值将values中的值存入result
    _cs149_vload_float(result, values + i,
                       maskIsNotNegative); //   output[i] = x; }

    // Write results back to memory
    _cs149_vstore_float(output + i, result, maskAll);
  }
}

// accepts an array of values and an array of exponents
//
// For each element, compute values[i]^exponents[i] and clamp value to
// 9.999.  Store result in output.
// 串行计算数组元素的指数幂，并将结果存储在output中，并将结果限制在9.999999以下
void clampedExpSerial(float *values, int *exponents, float *output, int N) {
  for (int i = 0; i < N; i++) {

    float x = values[i];
    int y = exponents[i];

    if (y == 0) {
      output[i] = 1.f;
    } else {
      float result = x;

      int count = y - 1;
      while (count > 0) {
        result *= x;
        count--;
      }

      if (result > 9.999999f) {
        result = 9.999999f;
      }

      output[i] = result;
    }
  }
}

// 向量化计算数组元素的指数幂，并将结果存储在output中
void clampedExpVector(float *values, int *exponents, float *output, int N) {

  //
  // CS149 STUDENTS TODO: Implement your vectorized version of
  // clampedExpSerial() here.
  //
  // Your solution should work for any value of
  // N and VECTOR_WIDTH, not just when VECTOR_WIDTH divides N
  //

  __cs149_vec_int zero_i = _cs149_vset_int(0);
  __cs149_vec_int one_i = _cs149_vset_int(1);
  __cs149_vec_float nine_999 = _cs149_vset_float(9.999999f);
  __cs149_mask mask_one = _cs149_init_ones(VECTOR_WIDTH); // mask_one = 1111
  __cs149_vec_float s_x;                                  // 输入值
  __cs149_vec_int s_e;                                    // 指数值
  __cs149_vec_float results;

  __cs149_mask mask_lt_n;       // 标记索引小于n的元素
  __cs149_mask mask_unfinished; // 索引有效位
  __cs149_mask mask_zero_exp;   // 指数为0的mask
  __cs149_mask gtn_999;         // 标记大于9.99的元素

  for (int i = 0; i < N; i += VECTOR_WIDTH) {

    // 每次循环处理VECTOR_WIDTH个元素，并标记小于N的元素
    mask_lt_n = _cs149_init_ones(min(abs(N - i), VECTOR_WIDTH));

    // 根据掩码有效位， s_x = values[i]
    _cs149_vload_float(s_x, values + i, mask_lt_n);
    // 根据掩码有效位， s_e = exponents[i],加载指数
    _cs149_vload_int(s_e, exponents + i, mask_lt_n);
    // results = 1.f，设置初始值
    _cs149_vset_float(results, 1.f, mask_lt_n);

    // 初始情况下，只要索引有效就表示计算未结束
    mask_unfinished = mask_lt_n;
    // 先将mask清空,初始化为0
    gtn_999 = _cs149_init_ones(0);

    // 如果 exponent == 0, 则标记为 1，否则标记为 0
    _cs149_veq_int(mask_zero_exp, s_e, zero_i, mask_one);
    // 标记 exponent != 0的位置，即指数为0的位置为0，指数不为0的位置为1
    mask_zero_exp = _cs149_mask_not(mask_zero_exp);

    // 开始计算前判断初始指数值是否为0
    mask_unfinished = _cs149_mask_and(mask_unfinished, mask_zero_exp);
    // 和mask_lt_n取交集，保证和索引无效的部分为0
    mask_unfinished = _cs149_mask_and(mask_unfinished, mask_lt_n);

    // 如果当前mask中存在元素为完成指数计算，继续计算
    while (_cs149_cntbits(mask_unfinished) > 0) {
      // 计算幂，将结果存入results中
      _cs149_vmult_float(results, results, s_x, mask_unfinished);

      // 根据掩码值的有效位，将对应的指数减1
      _cs149_vsub_int(s_e, s_e, one_i, mask_unfinished);

      // 标记大于9.999999的元素，设置到gtn_999中
      _cs149_vgt_float(gtn_999, results, nine_999, mask_unfinished);
      // 根据gtn_999中掩码有效位，如果结果大于9.999999，将其设为9.999999
      _cs149_vset_float(results, 9.999999f, gtn_999);

      // 计算过一次乘法，重新标记指数为0的位置，表示这些位置计算结束
      mask_zero_exp = _cs149_init_ones(0);
      _cs149_veq_int(mask_zero_exp, s_e, zero_i, mask_one);

      // 如果指数为 0 或结果大于 9.999, 则表示计算结束，用 1 表示
      __cs149_mask new_finished = _cs149_mask_or(mask_zero_exp, gtn_999);

      // 反转，计算结束的位置用 0 表示，计算未结束的位置用 1 表示
      mask_zero_exp = _cs149_mask_not(new_finished);
      // 如果指数为 0, 则表示计算结束，设置索引有效位
      mask_unfinished = _cs149_mask_and(mask_unfinished, mask_zero_exp);
      // 保险起见，再和 mask_lt_N 取交集，保证索引无效的部分为 0
      mask_unfinished = _cs149_mask_and(mask_unfinished, mask_lt_n);
    }

    // 结果存储到output中
    _cs149_vstore_float(output + i, results, mask_lt_n);
  }
}

// returns the sum of all elements in values
float arraySumSerial(float *values, int N) {
  float sum = 0;
  for (int i = 0; i < N; i++) {
    sum += values[i];
  }

  return sum;
}

// returns the sum of all elements in values
// You can assume N is a multiple of VECTOR_WIDTH
// You can assume VECTOR_WIDTH is a power of 2
float arraySumVector(float *values, int N) {

  //
  // CS149 STUDENTS TODO: Implement your vectorized version of arraySumSerial

  __cs149_vec_float v_sum; // 结果寄存器
  __cs149_mask mask_one =
      _cs149_init_ones(min(N, VECTOR_WIDTH)); // 掩码初始化为1

  _cs149_vset_float(v_sum, 0.f, mask_one); // 初始化为0
  __cs149_vec_float v_temp;

  //

  for (int i = 0; i < N; i += VECTOR_WIDTH) {
    mask_one = _cs149_init_ones(min(abs(N - i), VECTOR_WIDTH)); // 掩码初始化为1
    _cs149_vload_float(v_temp, values + i, mask_one); // 加载数据到寄存器
    _cs149_vadd_float(v_sum, v_sum, v_temp, mask_one); // 计算和
  }

  float sum = 0.f;
  for (int i = 0; i < VECTOR_WIDTH; i++) {
    sum += v_sum.value[i]; // 将结果存储到sum中
  }

  return sum;
}
