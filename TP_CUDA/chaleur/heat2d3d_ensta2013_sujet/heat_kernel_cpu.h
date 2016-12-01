#ifndef HEAT_KERNEL_CPU_H_
#define HEAT_KERNEL_CPU_H_

// 2D
void heat2d_ftcs_cpu_order2 (  real_t* data, real_t* dataNext );
void heat2d_ftcs_cpu_order2b(  real_t* data, real_t* dataNext );
void heat2d_ftcs_cpu_order4 (  real_t* data, real_t* dataNext );

void heat2d_ftcs_cpu_order2_with_mask (real_t* data, real_t* dataNext, int* mask );

// 3D
void heat3d_ftcs_cpu_order2 (  real_t* data, real_t* dataNext );
void heat3d_ftcs_cpu_order2b(  real_t* data, real_t* dataNext );
void heat3d_ftcs_cpu_order4 (  real_t* data, real_t* dataNext );

#endif // HEAT_KERNEL_CPU_H_
