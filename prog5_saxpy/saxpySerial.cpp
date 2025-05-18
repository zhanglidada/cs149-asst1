
void saxpySerial(int N, float scale, float X[], float Y[], float result[]) {

  for (int i = 0; i < N; i++) {
    // 对每三个元素执行两个数学运算，显然可并行化计算，具有可预测，规律的数据访问和可预测的执行成本
    result[i] = scale * X[i] + Y[i];
  }
}
