/*

  Note: This code was modified from example code
  originally provided by Intel.  To comply with Intel's open source
  licensing agreement, their copyright is retained below.

  -----------------------------------------------------------------

  Copyright (c) 2010-2011, Intel Corporation
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are
  met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    * Neither the name of Intel Corporation nor the names of its
      contributors may be used to endorse or promote products derived from
      this software without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
   IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
   TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
   PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
   OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

static inline int mandel(float c_re, float c_im, int count) {
  float z_re = c_re, z_im = c_im;
  int i;
  // 迭代count次
  for (i = 0; i < count; ++i) {

    if (z_re * z_re + z_im * z_im > 4.f)
      break;

    float new_re = z_re * z_re - z_im * z_im;
    float new_im = 2.f * z_re * z_im;
    z_re = c_re + new_re;
    z_im = c_im + new_im;
  }

  return i;
}

//
// MandelbrotSerial --
//
// Compute an image visualizing the mandelbrot set.  The resulting
// array contains the number of iterations required before the complex
// number corresponding to a pixel could be rejected from the set.
//
// * x0, y0, x1, y1 describe the complex coordinates mapping
//   into the image viewport.
// * width, height describe the size of the output image
// * startRow, totalRows describe how much of the image to compute
/**
 * @brief
 * 这四个参数定义了要计算的Mandelbrot集合的矩形区域在复数平面上的左下角和右上角的坐标。
   这里的坐标是以实数部分和虚数部分分别表示的，即(x0, y0)是左下角，(x1,
 y1)是右上角。
 * @param x0
 * @param y0
 * @param x1
 * @param y1
 * @param width         输出图像的宽度
 * @param height        输出图像的高度
 * @param startRow      用于支持并行计算，指定从哪一行开始计算
 * @param totalRows     用于支持并行计算，指定要计算的行数
 * @param maxIterations 迭代次数的上限
 * @param output
 输出数组，用于存储每个像素（或小块）是否属于Mandelbrot集合的结果。
                        通常这个数组会被映射到图像的颜色上，以可视化Mandelbrot集合
 */
void mandelbrotSerial(float x0, float y0, float x1, float y1, int width,
                      int height, int startRow, int totalRows,
                      int maxIterations, int output[]) {
  // 计算x轴和y轴的步长
  float dx = (x1 - x0) / width;
  float dy = (y1 - y0) / height;

  // 计算结束行号
  int endRow = startRow + totalRows;

  // 遍历每一行
  for (int j = startRow; j < endRow; j++) {
    // 遍历每一列
    for (int i = 0; i < width; ++i) {
      // 计算当前点的x和y坐标
      float x = x0 + i * dx;
      float y = y0 + j * dy;

      // 计算当前点在输出数组中的索引
      int index = (j * width + i);

      // 计算当前点的Mandelbrot集合值并存储到输出数组中
      output[index] = mandel(x, y, maxIterations);
    }
  }
}
