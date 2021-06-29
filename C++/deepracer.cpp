/*-----------------------------------------------------------------------
 * File: deepracer.cpp
 * Date: 01.06.2021
 * Athr: M. Khaled
 *-----------------------------------------------------------------------*/

#include <iostream>
#include <math.h>

#ifndef M_PI
#define M_PI 3.1415926f
#endif

#define ssDim 4
#define isDim 2

// utilitity functions
#define M_2PI ((2.0f*(M_PI)))
inline void wrapToPi(float* rad){
	*rad -= M_2PI * floor((*rad + M_PI) * (1.0f/M_2PI));
}

// some maps
inline float map_steering(const float angle_in){
    float p1 = -0.1167;
    float p2 = 0.01949;
    float p3 = 0.3828;
    float p4 = -0.0293;
    float x = angle_in;
    return p1*x*x*x + p2*x*x + p3*x + p4;
}
inline float map_speed(const int speed_in){
	switch(speed_in) {
		case  6: return  0.70f;
		case  5: return  0.65f;
		case  4: return  0.60f;
		case  3: return  0.55f;
		case  2: return  0.50f;
		case  1: return  0.45f;
		case  0: return  0.00f;
		case -1: return -0.45f;
		case -2: return -0.50f;
		case -3: return -0.55f;
		case -4: return -0.60f;
		case -5: return -0.65f;
		case -6: return -0.70f;
	}
	return 0.0f;
}
inline void get_v_params(float u_speed, float* a, float* b){
	float K,T;
    if (u_speed == 0.0f){
        K = 0.0f; 
		T = 0.25f;
	}
	else if (u_speed == 0.45f){
        K = 1.9953f; 
		T = 0.9933f;
    }
	else if (u_speed == 0.50f){
        K = 2.3567f;
		T = 0.8943f;
    }
	else if (u_speed == 0.55f){
        K = 3.0797f; 
		T = 0.88976f;
    }
	else if (u_speed == 0.60f){
        K = 3.2019f; 
		T = 0.87595f;
    }
	else if (u_speed == 0.65f){
        K = 3.3276f; 
		T = 0.89594f;
    }
	else if (u_speed == 0.70f){
        K = 3.7645f; 
		T = 0.92501f;
    }
	else if (u_speed == -0.45f){
        K = 1.8229f; 
		T = 1.8431f;
    }
	else if (u_speed == -0.50f){
        K = 2.3833f; 
		T = 1.2721f;
    }
	else if (u_speed == -0.55f){
        K = 2.512f; 
		T = 1.1403f;
    }
	else if (u_speed == -0.60f){
        K = 3.0956f; 
		T = 1.1278f;
    }
	else if (u_speed == -0.65f){
        K = 3.55f; 
		T = 1.1226f;
    }
	else if (u_speed == -0.70f){
        K = 3.6423f; 
		T = 1.1539f;
	}
    else{
        printf("get_v_params: Invalid input !\n");
		K = 0.0f;
		T = 0.0f;
	}
    *a = -1.0f/T;
    *b = K/T;
}

// this function should get called after compuing the solution to the ODE
// we use it to modiy the angle to be withtin the range -pi:pi
void post_ODE(float* xx, const float* x, const float* r, const float* u) {
	wrapToPi(&xx[2]);
}


// RHS of the ODE
// u[0] steering angle in [-1,1] (must be mapped -> [-2.8,2.8])
// u[1] speed in [-6;6] (must be mapped -> [-0.7,-0.45]U{0.0}U[0.45,0.7])
// x[0] x-pos
// x[1] y-pos
// x[2] theta (orientation)
// x[3] forward velocity
void ode_rhs(float* xx, const float* x, const float* u);
void ode_rhs(float* xx, const float* x, const float* u) {
	
	float u_steer = map_steering(u[0]);
	float u_speed = map_speed((int)u[1]);
	float L = 0.165f;
	float a,b;
	get_v_params(u_speed, &a, &b);

	xx[0] = x[3]*cos(x[2]);
	xx[1] = x[3]*sin(x[2]);
	xx[2] = (x[3]/L)*tan(u_steer);
	xx[3] = a*x[3] + b*u_speed;
}


// Runge-Kutta ODE Solver
#define SAMPLING_PERIOD 0.5
#define RK4_NINT 5
#define RK4_H ((((float)SAMPLING_PERIOD)/((float)RK4_NINT)))
void rk4OdeSolver(float* xx, const float* x, const float* u);
void rk4OdeSolver(float* xx, const float* x, const float* u) {
	float k0[ssDim];
	float k1[ssDim];
	float k2[ssDim];
	float k3[ssDim];
	float tmp[ssDim];

	for (unsigned int i = 0; i < ssDim; i++)
		xx[i] = x[i];

	for (unsigned int t = 0; t < RK4_NINT; t++) {

        ode_rhs(k0, xx, u);

		for (unsigned int i = 0; i < ssDim; i++)
			tmp[i] = xx[i] + RK4_H / 2.0f*k0[i];

        ode_rhs(k1, tmp, u);

		for (unsigned int i = 0; i < ssDim; i++)
			tmp[i] = xx[i] + RK4_H / 2.0f*k1[i];

        ode_rhs(k2, tmp, u);

		for (unsigned int i = 0; i < ssDim; i++)
			tmp[i] = xx[i] + RK4_H * k2[i];

        ode_rhs(k3, tmp, u);

		for (unsigned int i = 0; i < ssDim; i++)
			xx[i] = xx[i] + (RK4_H / 6.0f)*(k0[i] + 2.0f*k1[i] + 2.0f*k2[i] + k3[i]);

	}
}


int main() {

    // some required vars
	float u[isDim] = {0.8, 3.0};
    float x_init[ssDim] = {0.0f, 0.0f, 0.0f, 0.0f};
	float x_final[ssDim];
    

    // solve the ODEs and compute the growth-bound
    rk4OdeSolver(x_final, x_init, u);

	// wrap the orrinetation angle xx[2] in [-pi,pi]
	wrapToPi(&x_final[2]);

	// print post state
    std::cout << "Post state: ";
	for (unsigned int i = 0; i<ssDim; i++) {
		std::cout << x_final[i] << ", ";
	}    
    std::cout << std::endl;
	return 0;
}

