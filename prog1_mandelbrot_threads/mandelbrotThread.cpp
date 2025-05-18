/**
 * @file mandelbrotThread.cpp
 * @author zhangli
 * @brief
 * @version 0.1
 * @date 2025-01-05
 *
 * @copyright Copyright (c) 2025
 *
 * 实现Mandelbrot集图像生成的多线程版本
 *
 */

#include <stdio.h>
#include <thread>

#include "../common/CycleTimer.h"
//#include "CycleTimer.h"

// 工作线程对应参数
typedef struct {
  float x0;
  float x1;
  float y0;
  float y1;
  unsigned int width;
  unsigned int height;
  int maxIterations;
  int *output;
  int threadId;
  int numThreads;
} WorkerArgs;

extern void mandelbrotSerial(float x0, float y0, float x1, float y1, int width,
                             int height, int startRow, int numRows,
                             int maxIterations, int output[]);

// 工作线程入口函数，每个线程使用串行版本计算部分数据
void workerThreadStart(WorkerArgs *const args) {

  // TODO FOR CS149 STUDENTS: Implement the body of the worker
  // thread here. Each thread should make a call to mandelbrotSerial()
  // to compute a part of the output image.  For example, in a
  // program that uses two threads, thread 0 could compute the top
  // half of the image and thread 1 could compute the bottom half.

  // 单个线程的计算行数
  int row_num_per_t = args->height / args->numThreads;

  // 每个线程的计算区域（根据线程的threadid计算当前线程需要计算的图像区域）
  int startRow = args->threadId * row_num_per_t;
  int endRow = startRow + row_num_per_t;

  if (args->threadId == args->numThreads - 1) {
    endRow = args->height; // 最后一个线程计算剩余的行
  }

  double t_startTime = CycleTimer::currentSeconds();

  // 调用串行版本的mandelbrotSerial函数，计算每个线程对应的图像区域
  mandelbrotSerial(args->x0, args->y0, args->x1, args->y1, args->width,
                   args->height, startRow, endRow - startRow,
                   args->maxIterations, args->output);

  double t_endTime = CycleTimer::currentSeconds();

  printf("Thread %d: [%.3f] ms\n", args->threadId,
         (t_endTime - t_startTime) * 1000);
}

// 优化工作线程，对每个线程的计算粒度细分，每个线程从上到下按照一定间隔计算多块区域
void workerThreadStart_ex(WorkerArgs *const args) {

  double t_startTime = CycleTimer::currentSeconds();

  const unsigned int chunk_size = 16; // 每个线程计算的块大小

  // 从当前线程id开始，每隔numThreads*chunk_size计算一块区域
  for (unsigned int cur_row = args->threadId * chunk_size;
       cur_row < args->height; cur_row += args->numThreads * chunk_size) {
    // 防止越界，取较小值
    int numRows = std::min(chunk_size, args->height - cur_row);

    // 调用串行版本的mandelbrotSerial函数，计算每个线程对应的图像区域
    mandelbrotSerial(args->x0, args->y0, args->x1, args->y1, args->width,
                     args->height, cur_row, numRows, args->maxIterations,
                     args->output);
  }

  double t_endTime = CycleTimer::currentSeconds();

  printf("Thread %d: [%.3f] ms\n", args->threadId,
         (t_endTime - t_startTime) * 1000);
}

// 将操作分配到多个线程中去执行
void mandelbrotThread(int numThreads, // number of threads to use
                      float x0, float y0, float x1, float y1, int width,
                      int height, int maxIterations, int output[]) {

  static constexpr int MAX_THREADS =
      12; // 限制最大线程个数（本地cpu只有6核12线程）

  if (numThreads > MAX_THREADS) {
    fprintf(stderr, "Error: Max allowed threads is %d\n", MAX_THREADS);
    exit(1);
  }

  // Creates thread objects that do not yet represent a thread.
  std::thread workers[MAX_THREADS];
  WorkerArgs args[MAX_THREADS];

  // 设置每个线程的初始化arg值
  for (int i = 0; i < numThreads; i++) {
    args[i].x0 = x0;
    args[i].y0 = y0;
    args[i].x1 = x1;
    args[i].y1 = y1;
    args[i].width = width;
    args[i].height = height;
    args[i].maxIterations = maxIterations;
    args[i].numThreads = numThreads;
    args[i].output = output;

    args[i].threadId = i;
  }

  // 创建线程数组，每个线程加进去
  //   for (int i = 1; i < numThreads; i++) {
  //     workers[i] = std::thread(workerThreadStart, &args[i]);
  //   }

  //   workerThreadStart(&args[0]);

  // 使用优化后的算法，每个线程计算一块区域
  for (int i = 1; i < numThreads; i++) {
    workers[i] = std::thread(workerThreadStart_ex, &args[i]);
  }

  workerThreadStart_ex(&args[0]);

  // join worker threads
  for (int i = 1; i < numThreads; i++) {
    workers[i].join();
  }
}
