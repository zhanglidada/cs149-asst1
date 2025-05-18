// Define vector unit width here
#define VECTOR_WIDTH 2

#ifndef CS149INTRIN_H_
#define CS149INTRIN_H_

#include "logger.h"
#include <cmath>
#include <cstdlib>

//*******************
//* Type Definition *
//*******************

extern Logger CS149Logger;

template <typename T> struct __cs149_vec { T value[VECTOR_WIDTH]; };

// Declare a mask with __cs149_mask
struct __cs149_mask : __cs149_vec<bool> {};

// Declare a floating point vector register with __cs149_vec_float
#define __cs149_vec_float __cs149_vec<float>

// Declare an integer vector register with __cs149_vec_int
#define __cs149_vec_int __cs149_vec<int>

//***********************
//* Function Definition *
//***********************

// Return a mask initialized to 1 in the first N lanes and 0 in the others
// 返回一个掩码，前n个初始化为1，后面都为0（参数默认值为VECTOR_WIDTH）
__cs149_mask _cs149_init_ones(int first = VECTOR_WIDTH);

// 将maska中的每个元素取反，返回一个新的掩码
__cs149_mask _cs149_mask_not(__cs149_mask &maska);

// Return (maska | maskb) 返回掩码 或
__cs149_mask _cs149_mask_or(__cs149_mask &maska, __cs149_mask &maskb);

// Return (maska & maskb) 返回掩码 与
__cs149_mask _cs149_mask_and(__cs149_mask &maska, __cs149_mask &maskb);

// Count the number of 1s in maska 返回掩码 a 中1的个数
int _cs149_cntbits(__cs149_mask &maska);

// Set register to value if vector lane is active，otherwise keep the old value
// 如果mask表示vector中值有效，将其设置为value，否则vector寄存器（数组）保持原值
void _cs149_vset_float(__cs149_vec_float &vecResult, float value,
                       __cs149_mask &mask);
void _cs149_vset_int(__cs149_vec_int &vecResult, int value, __cs149_mask &mask);

// For user's convenience, returns a vector register with all lanes initialized
// to value 为了方便用户，返回一个所有元素初始化为value的向量寄存器
__cs149_vec_float _cs149_vset_float(float value);
__cs149_vec_int _cs149_vset_int(int value);

// Copy values from vector register src to vector register dest if vector lane
// active， otherwise keep the old value
// 根据mask表示的有效位，将src中对应有效的value复制到dst中
void _cs149_vmove_float(__cs149_vec_float &dest, __cs149_vec_float &src,
                        __cs149_mask &mask);
void _cs149_vmove_int(__cs149_vec_int &dest, __cs149_vec_int &src,
                      __cs149_mask &mask);

// Load values from array src to vector register dest if vector lane active，
// otherwise keep the old value
// 根据mask表示的有效位，将src中对应有效的value加载到dst中
void _cs149_vload_float(__cs149_vec_float &dest, float *src,
                        __cs149_mask &mask);
void _cs149_vload_int(__cs149_vec_int &dest, int *src, __cs149_mask &mask);

// Store values from vector register src to array dest if vector lane active,
// otherwise keep the old value
// 根据mask表示的有效位，将src中对应有效的value存储到dst中
void _cs149_vstore_float(float *dest, __cs149_vec_float &src,
                         __cs149_mask &mask);
void _cs149_vstore_int(int *dest, __cs149_vec_int &src, __cs149_mask &mask);

// Return calculation of (veca + vecb) if vector lane active , otherwise keep
// the old value 根据mask的有效位，返回veca + vecb的结果，否则保持原值
void _cs149_vadd_float(__cs149_vec_float &vecResult, __cs149_vec_float &veca,
                       __cs149_vec_float &vecb, __cs149_mask &mask);
void _cs149_vadd_int(__cs149_vec_int &vecResult, __cs149_vec_int &veca,
                     __cs149_vec_int &vecb, __cs149_mask &mask);

// Return calculation of (veca - vecb) if vector lane active， otherwise keep
// the old value 根据mask的有效位，返回veca - vecb的结果，否则保持原值
void _cs149_vsub_float(__cs149_vec_float &vecResult, __cs149_vec_float &veca,
                       __cs149_vec_float &vecb, __cs149_mask &mask);
void _cs149_vsub_int(__cs149_vec_int &vecResult, __cs149_vec_int &veca,
                     __cs149_vec_int &vecb, __cs149_mask &mask);

// Return calculation of (veca * vecb) if vector lane active, otherwise keep the
// old value 根据mask的有效位，返回veca * vecb的结果，否则保持原值
void _cs149_vmult_float(__cs149_vec_float &vecResult, __cs149_vec_float &veca,
                        __cs149_vec_float &vecb, __cs149_mask &mask);
void _cs149_vmult_int(__cs149_vec_int &vecResult, __cs149_vec_int &veca,
                      __cs149_vec_int &vecb, __cs149_mask &mask);

// Return calculation of (veca / vecb) if vector lane active, otherwise keep the
// old value 根据mask的有效位，返回veca / vecb的结果，否则保持原值
void _cs149_vdiv_float(__cs149_vec_float &vecResult, __cs149_vec_float &veca,
                       __cs149_vec_float &vecb, __cs149_mask &mask);
void _cs149_vdiv_int(__cs149_vec_int &vecResult, __cs149_vec_int &veca,
                     __cs149_vec_int &vecb, __cs149_mask &mask);

// Return calculation of absolute value abs(veca) if vector lane active,
// otherwise keep the old value 根据mask的有效位，返回 veca
// 的绝对值结果，否则保持原值
void _cs149_vabs_float(__cs149_vec_float &vecResult, __cs149_vec_float &veca,
                       __cs149_mask &mask);
void _cs149_vabs_int(__cs149_vec_int &vecResult, __cs149_vec_int &veca,
                     __cs149_mask &mask);

// Return a mask of (veca > vecb) if vector lane active, otherwise keep the old
// value 根据mask的有效位，返回veca > vecb的结果，否则保持原值
void _cs149_vgt_float(__cs149_mask &vecResult, __cs149_vec_float &veca,
                      __cs149_vec_float &vecb, __cs149_mask &mask);
void _cs149_vgt_int(__cs149_mask &vecResult, __cs149_vec_int &veca,
                    __cs149_vec_int &vecb, __cs149_mask &mask);

// Return a mask of (veca < vecb) if vector lane active，otherwise keep the old
// value 根据mask的有效位，返回 veca < vecb的结果，否则保持原值
void _cs149_vlt_float(__cs149_mask &vecResult, __cs149_vec_float &veca,
                      __cs149_vec_float &vecb, __cs149_mask &mask);
void _cs149_vlt_int(__cs149_mask &vecResult, __cs149_vec_int &veca,
                    __cs149_vec_int &vecb, __cs149_mask &mask);

// Return a mask of (veca == vecb) if vector lane active, otherwise keep the old
// value 根据mask的有效位，返回 veca == vecb的结果，否则保持原值
void _cs149_veq_float(__cs149_mask &vecResult, __cs149_vec_float &veca,
                      __cs149_vec_float &vecb, __cs149_mask &mask);
void _cs149_veq_int(__cs149_mask &vecResult, __cs149_vec_int &veca,
                    __cs149_vec_int &vecb, __cs149_mask &mask);

// 从前往后，将相邻的元素对相加, 表现如下：
// [0 1 2 3 4 5] -> [0+1 0+1 2+3 2+3 4+5 4+5]
void _cs149_hadd_float(__cs149_vec_float &vecResult, __cs149_vec_float &vec);

// Performs an even-odd interleaving where all even-indexed elements move to
// front half， of the array and odd-indexed to the back half, so
//  [0 1 2 3 4 5 6 7] -> [0 2 4 6 1 3 5 7]
void _cs149_interleave_float(__cs149_vec_float &vecResult,
                             __cs149_vec_float &vec);

// Add a customized log to help debugging
void addUserLog(const char *logStr);

#endif
