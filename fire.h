#pragma once

#include <stdio.h>
#define N 30
#define INDEX(x,y,z) ((x)+(N+2)*(y) + (N+2)*(N+2)*(z))
#define SIZE ((N+2)*(N+2)*(N+2))

#define SWAP(x0,x) {float *tmp=x0;x0=x;x=tmp;}

class Fire {

public:
	float buffers[10][SIZE];
	float* d, * d0;
	float* T, * T0;	
	float* u, * u0;			
	float* v, * v0;			
	float* w, * w0;

protected:
	void add_source(float* src, float* dst, float dt);
	void add_buoyancy(float dt);
	void set_boundaries(int b, float* x);
	void diffuse(int b, float* x0, float* x, float diff, float dt);
	void advect(int b, float* x0, float* x, float* uu, float* vv, float* ww, float dt);
	void advect_cool(int b, float* x0, float* x, float* y0, float* y, float* uu, float* vv, float* ww, float dt);
	void set_equilibrium(void);
	void vorticity_confinement(float dt);

	void vel_step(float dt);
	void dens_step(float dt);
	void dens_temp_step(float dt);

	void clear_buffer(float* x);
	void clear_sources(void);

public:
	float sd[SIZE], su[SIZE], sv[SIZE], sw[SIZE], sT[SIZE];	
	float diffusion, viscosity, cooling, buoyancy, vc_eps;

	Fire(void);
	~Fire(void);

	void step(float dt);

};