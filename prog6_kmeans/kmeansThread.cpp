#include <algorithm>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <thread>

#include "CycleTimer.h"

using namespace std;

typedef struct {
  // Control work assignments
  int start, end; // 每个工作任务计算数据范围

  // 下面的属性指向所有arg共享的数据
  double *data;
  double *clusterCentroids;
  int *clusterAssignments;
  double *currCost;
  int M, N, K;
} WorkerArgs;

/**
 * Checks if the algorithm has converged.
 *
 * @param prevCost Pointer to the K dimensional array containing cluster costs
 *    from the previous iteration.
 * @param currCost Pointer to the K dimensional array containing cluster costs
 *    from the current iteration.
 * @param epsilon Predefined hyperparameter which is used to determine when
 *    the algorithm has converged.
 * @param K The number of clusters.
 *
 * NOTE: DO NOT MODIFY THIS FUNCTION!!!
 */
// 此接口不可调整
static bool stoppingConditionMet(double *prevCost, double *currCost,
                                 double epsilon, int K) {
  for (int k = 0; k < K; k++) {
    if (abs(prevCost[k] - currCost[k]) > epsilon)
      return false;
  }
  return true;
}

/**
 * Computes L2 distance between two points of dimension nDim.
 *
 * @param x Pointer to the beginning of the array representing the first
 *     data point.
 * @param y Poitner to the beginning of the array representing the second
 *     data point.
 * @param nDim The dimensionality (number of elements) in each data point
 *     (must be the same for x and y).
 */
// 可优化函数
double dist(double *x, double *y, int nDim) {
  double accum = 0.0;
  for (int i = 0; i < nDim; i++) {
    accum += pow((x[i] - y[i]), 2);
  }
  return sqrt(accum);
}

/**
 * Assigns each data point to its "closest" cluslter centroid.
 */
// 可优化函数
// 将每个数据点添加到最近的聚类中心中
void computeAssignments(WorkerArgs *const args) {
  double *minDist = new double[args->M];

  // Initialize arrays
  for (int m = args->start; m < args->end; m++) {
    minDist[m] = 1e30;
    args->clusterAssignments[m] = -1;
  }

  // Assign datapoints to closest centroids
  for (int m = args->start; m < args->end; m++) {
    for (int k = 0; k < args->K; k++) {
      double d = dist(&args->data[m * args->N],
                      &args->clusterCentroids[k * args->N], args->N);
      if (d < minDist[m]) {
        minDist[m] = d;
        args->clusterAssignments[m] = k;
      }
    }
  }

  free(minDist);
}

/**
 * Given the cluster assignments, computes the new centroid locations for
 * each cluster.
 */
// 可优化函数
void computeCentroids(WorkerArgs *const args) {
  int *counts = new int[args->K];

  // Zero things out
  for (int k = 0; k < args->K; k++) {
    counts[k] = 0;
    for (int n = 0; n < args->N; n++) {
      args->clusterCentroids[k * args->N + n] = 0.0;
    }
  }

  // Sum up contributions from assigned examples
  for (int m = 0; m < args->M; m++) {
    int k = args->clusterAssignments[m];
    for (int n = 0; n < args->N; n++) {
      args->clusterCentroids[k * args->N + n] += args->data[m * args->N + n];
    }
    counts[k]++;
  }

  // Compute means
  for (int k = 0; k < args->K; k++) {
    counts[k] = max(counts[k], 1); // prevent divide by 0
    for (int n = 0; n < args->N; n++) {
      args->clusterCentroids[k * args->N + n] /= counts[k];
    }
  }

  free(counts);
}

/**
 * Computes the per-cluster cost. Used to check if the algorithm has converged.
 */
// 可优化函数
void computeCost(WorkerArgs *const args) {
  double *accum = new double[args->K];

  // Zero things out
  for (int k = 0; k < args->K; k++) {
    accum[k] = 0.0;
  }

  // Sum cost for all data points assigned to centroid
  for (int m = 0; m < args->M; m++) {
    int k = args->clusterAssignments[m];
    accum[k] += dist(&args->data[m * args->N],
                     &args->clusterCentroids[k * args->N], args->N);
  }

  // Update costs
  for (int k = args->start; k < args->end; k++) {
    args->currCost[k] = accum[k];
  }

  free(accum);
}

/**
 * Computes the K-Means algorithm, using std::thread to parallelize the work.
 *
 * @param data Pointer to an array of length M*N representing the M different N
 *     dimensional data points clustered. The data is layed out in a "data point
 *     major" format, so that data[i*N] is the start of the i'th data point in
 *     the array. The N values of the i'th datapoint are the N values in the
 *     range data[i*N] to data[(i+1) * N].
 * @param clusterCentroids Pointer to an array of length K*N representing the K
 *     different N dimensional cluster centroids. The data is laid out in
 *     the same way as explained above for data.
 * @param clusterAssignments Pointer to an array of length M representing the
 *     cluster assignments of each data point, where clusterAssignments[i] = j
 *     indicates that data point i is closest to cluster centroid j.
 * @param M The number of data points to cluster.
 * @param N The dimensionality of the data points.
 * @param K The number of cluster centroids.
 * @param epsilon The algorithm is said to have converged when
 *     |currCost[i] - prevCost[i]| < epsilon for all i where i = 0, 1, ..., K-1
 */
/*
  data: 指向M*N维的数据起始地址
  clusterCentroids: 指向一个长度为 K * N 的数组的指针，该数组表示 K 个不同的 N
  维聚类质心。数据的排列方式与上述数据的排列方式相同。
  clusterAssignments:
  指向一个长度为M的数组的指针，该数组表示每个数据点的聚类分配情况。
  其中，clusterAssignments[i] = j 表明数据点 i 距离聚类质心 j 最近。
  M: 集合中数据个数。
  N: 数据点的维度。
  K: 聚类中心的数量。
  epsilon: 当所有数据点都满足|currCost[i] - prevCost[i]| < epsilon
  的时候，说明数据已经收敛
*/
void kMeansThread(double *data, double *clusterCentroids,
                  int *clusterAssignments, int M, int N, int K,
                  double epsilon) {

  // Used to track convergence
  double *prevCost = new double[K];
  double *currCost = new double[K];

  // The WorkerArgs array is used to pass inputs to and return output from
  // 初始化共享数据
  WorkerArgs args;
  args.data = data; // 每个arg都指向全部数据的起始位置
  args.clusterCentroids = clusterCentroids;
  args.clusterAssignments = clusterAssignments;
  args.currCost = currCost;
  args.M = M;
  args.N = N;
  args.K = K;

  // Initialize arrays to track cost
  for (int i = 0; i < K; i++) {
    prevCost[i] = 1e30;
    currCost[i] = 0.0;
  }

  // 想到最简单提升效率的办法就是多线程
  const int threadnum = 12;
  // 暂时先定义和物理线程数一样的线程个数
  auto workers = new std::thread[threadnum];
  // 定义同样数量work arg
  auto workargs = new WorkerArgs[threadnum];

  /* Main K-Means Algorithm Loop */
  int iter = 0;
  while (!stoppingConditionMet(prevCost, currCost, epsilon, K)) {
    // Update cost arrays (for checking convergence criteria)
    for (int i = 0; i < K; i++) {
      prevCost[i] = currCost[i];
    }

    // Setup args struct
    args.start = 0;
    args.end = K;

    // 分配多线程
    for (int i = 1; i < threadnum; i++) {
      workargs[i] = args;
      workargs[i].start = M / threadnum * i;
      workargs[i].end = min(M, workargs[i].start + M / threadnum);
      workers[i] = std::thread(computeAssignments, &workargs[i]);
    }

    workargs[0] = args;
    workargs[0].start = 0;
    workargs[0].end = M / threadnum;
    computeAssignments(&workargs[0]);

    // 启动各个线程
    for (int i = 1; i < threadnum; i++) {
      workers[i].join();
    }

    computeCentroids(&args);
    computeCost(&args);

    iter++;
  }

  delete[] currCost;
  delete[] prevCost;
  delete[] workers;
  delete[] workargs;
}
