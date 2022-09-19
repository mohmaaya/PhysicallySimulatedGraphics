#include <iostream>
#include "fire.h"
#include <math.h>

Fire::Fire() {

	for (int i = 0; i < 10; i++)
		clear_buffer(buffers[i]);

	int i = 0;
	d = buffers[i++]; d0 = buffers[i++];
	T = buffers[i++]; T0 = buffers[i++];
	u = buffers[i++]; u0 = buffers[i++];
	v = buffers[i++]; v0 = buffers[i++];
	w = buffers[i++]; w0 = buffers[i++];

	clear_sources();

	int size = (N + 2) * (N + 2) * (N + 2);
	for (int i = 0; i < size; i++)
		v[i] = -0.5f;
}


Fire::~Fire() {}

void Fire::set_boundaries(int b, float* x)
{
	int i, j;
		for (i=1; i<=N; i++)
		{
			for (j=1; j<=N; j++) {
				x[INDEX(0,i,j)]    = (b==1) ? -x[INDEX(1,i,j)] : x[INDEX(1,i,j)];
				x[INDEX(N+1,i,j)]  = (b==1) ? -x[INDEX(N,i,j)] : x[INDEX(N,i,j)];
				x[INDEX(i,0,j)]    = (b==2) ? -x[INDEX(i,1,j)] : x[INDEX(i,1,j)];
				x[INDEX(i,N+1,j)]  = (b==2) ? -x[INDEX(i,N,j)] : x[INDEX(i,N,j)];
				x[INDEX(i,j,0)]    = (b==3) ? -x[INDEX(i,j,1)] : x[INDEX(i,j,1)];
				x[INDEX(i,j,N+1)]  = (b==3) ? -x[INDEX(i,j,N)] : x[INDEX(i,j,N)];
			}
		}

		x[INDEX(0,0,0)]       = (x[INDEX(1,0,0)]    +x[INDEX(0,1,0)]    +x[INDEX(0,0,1)])    /3;
		x[INDEX(0,N+1,0)]     = (x[INDEX(1,N+1,0)]  +x[INDEX(0,N,0)]    +x[INDEX(0,N+1,1)])  /3;
		x[INDEX(N+1,0,0)]     = (x[INDEX(N,0,0)]    +x[INDEX(N+1,1,0)]  +x[INDEX(N+1,0,1)])  /3;
		x[INDEX(N+1,N+1,0)]   = (x[INDEX(N,N+1,0)]  +x[INDEX(N+1,N,0)]  +x[INDEX(N+1,N+1,1)])/3;
		x[INDEX(0,0,N+1)]     = (x[INDEX(1,0,N+1)]  +x[INDEX(0,1,N+1)]  +x[INDEX(0,0,N)])    /3;
		x[INDEX(0,N+1,N+1)]   = (x[INDEX(1,N+1,N+1)]+x[INDEX(0,N,N+1)]  +x[INDEX(0,N+1,N)])  /3;
		x[INDEX(N+1,0,N+1)]   = (x[INDEX(N,0,N+1)]  +x[INDEX(N+1,1,N+1)]+x[INDEX(N+1,0,N)])  /3;
		x[INDEX(N+1,N+1,N+1)] = (x[INDEX(N,N+1,N+1)]+x[INDEX(N+1,N,N+1)]+x[INDEX(N+1,N+1,N)])/3;
}

void Fire::add_source(float* src, float* dst, float dt)
{
	int i, size = (N + 2) * (N + 2) * (N + 2);

	for (i = 0; i < size; i++)
		dst[i] += src[i] * dt;
}

void Fire::add_buoyancy(float dt)
{
	int i, size = (N + 2) * (N + 2) * (N + 2);

	for (i = 0; i < size; i++)
		v[i] += -T[i] * buoyancy * dt;
}

inline void Fire::diffuse(int b, float* x0, float* x, float diff, float dt)
{
	int i, j, k, l;
	float a = dt * diff * N * N * N;
	for (l = 0; l < 20; l++)
	{
		for (k = 1; k <= N; k++)
		{
			for (j = 1; j <= N; j++)
			{
				for (i = 1; i <= N; i++)
				{
					x[INDEX(i, j, k)] = (x0[INDEX(i, j, k)] + a * (
						x[INDEX(i - 1, j, k)] + x[INDEX(i + 1, j, k)] +
						x[INDEX(i, j - 1, k)] + x[INDEX(i, j + 1, k)] +
						x[INDEX(i, j, k - 1)] + x[INDEX(i, j, k + 1)])) / (1 + 6 * a);
				}
			}
		}
		set_boundaries(b, x);
	}
}

inline void Fire::advect(int b, float* x0, float* x, float* uu, float* vv, float* ww, float dt)
{
	int i, j, k, i0, j0, k0, i1, j1, k1;
	float sx0, sx1, sy0, sy1, sz0, sz1, v0, v1;
	float xx, yy, zz, dt0;
	dt0 = dt * N;
	for (k = 1; k <= N; k++)
	{
		for (j = 1; j <= N; j++)
		{
			for (i = 1; i <= N; i++)
			{
				xx = i - dt0 * uu[INDEX(i, j, k)];
				yy = j - dt0 * vv[INDEX(i, j, k)];
				zz = k - dt0 * ww[INDEX(i, j, k)];
				if (xx < 0.5) xx = 0.5f; if (xx > N + 0.5) xx = N + 0.5f; i0 = (int)xx; i1 = i0 + 1;
				if (yy < 0.5) yy = 0.5f; if (yy > N + 0.5) yy = N + 0.5f; j0 = (int)yy; j1 = j0 + 1;
				if (zz < 0.5) zz = 0.5f; if (zz > N + 0.5) zz = N + 0.5f; k0 = (int)zz; k1 = k0 + 1;
				sx1 = xx - i0; sx0 = 1 - sx1;
				sy1 = yy - j0; sy0 = 1 - sy1;
				sz1 = zz - k0; sz0 = 1 - sz1;
				v0 = sx0 * (sy0 * x0[INDEX(i0, j0, k0)] + sy1 * x0[INDEX(i0, j1, k0)]) + sx1 * (sy0 * x0[INDEX(i1, j0, k0)] + sy1 * x0[INDEX(i1, j1, k0)]);
				v1 = sx0 * (sy0 * x0[INDEX(i0, j0, k1)] + sy1 * x0[INDEX(i0, j1, k1)]) + sx1 * (sy0 * x0[INDEX(i1, j0, k1)] + sy1 * x0[INDEX(i1, j1, k1)]);
				x[INDEX(i, j, k)] = sz0 * v0 + sz1 * v1;
			}
		}
	}
	set_boundaries(b, d);
}

inline void Fire::advect_cool(int b, float* x0, float* x, float* y0, float* y, float* uu, float* vv, float* ww, float dt)
{
	int i, j, k, i0, j0, k0, i1, j1, k1;
	float sx0, sx1, sy0, sy1, sz0, sz1, v0, v1;
	float xx, yy, zz, dt0, c0;
	dt0 = dt * N;
	c0 = 1.0f - cooling * dt;
	for (k = 1; k <= N; k++)
	{
		for (j = 1; j <= N; j++)
		{
			for (i = 1; i <= N; i++)
			{
				xx = i - dt0 * uu[INDEX(i, j, k)];
				yy = j - dt0 * vv[INDEX(i, j, k)];
				zz = k - dt0 * ww[INDEX(i, j, k)];
				if (xx < 0.5) xx = 0.5f; if (xx > N + 0.5) xx = N + 0.5f; i0 = (int)xx; i1 = i0 + 1;
				if (yy < 0.5) yy = 0.5f; if (yy > N + 0.5) yy = N + 0.5f; j0 = (int)yy; j1 = j0 + 1;
				if (zz < 0.5) zz = 0.5f; if (zz > N + 0.5) zz = N + 0.5f; k0 = (int)zz; k1 = k0 + 1;
				sx1 = xx - i0; sx0 = 1 - sx1;
				sy1 = yy - j0; sy0 = 1 - sy1;
				sz1 = zz - k0; sz0 = 1 - sz1;
				v0 = sx0 * (sy0 * x0[INDEX(i0, j0, k0)] + sy1 * x0[INDEX(i0, j1, k0)]) + sx1 * (sy0 * x0[INDEX(i1, j0, k0)] + sy1 * x0[INDEX(i1, j1, k0)]);
				v1 = sx0 * (sy0 * x0[INDEX(i0, j0, k1)] + sy1 * x0[INDEX(i0, j1, k1)]) + sx1 * (sy0 * x0[INDEX(i1, j0, k1)] + sy1 * x0[INDEX(i1, j1, k1)]);
				x[INDEX(i, j, k)] = sz0 * v0 + sz1 * v1;
				v0 = sx0 * (sy0 * y0[INDEX(i0, j0, k0)] + sy1 * y0[INDEX(i0, j1, k0)]) + sx1 * (sy0 * y0[INDEX(i1, j0, k0)] + sy1 * y0[INDEX(i1, j1, k0)]);
				v1 = sx0 * (sy0 * y0[INDEX(i0, j0, k1)] + sy1 * y0[INDEX(i0, j1, k1)]) + sx1 * (sy0 * y0[INDEX(i1, j0, k1)] + sy1 * y0[INDEX(i1, j1, k1)]);
				y[INDEX(i, j, k)] = (sz0 * v0 + sz1 * v1) * c0;
			}
		}
	}
	set_boundaries(b, d);
}

void Fire::set_equilibrium(void)
{
	float* p = u0;	float* div = v0;	
	int i, j, k, l;
	float h;
	h = 1.0f / N;
	for (k = 1; k <= N; k++) {
		for (j = 1; j <= N; j++) {
			for (i = 1; i <= N; i++) {
				div[INDEX(i, j, k)] = -h * (
					u[INDEX(i + 1, j, k)] - u[INDEX(i - 1, j, k)] +
					v[INDEX(i, j + 1, k)] - v[INDEX(i, j - 1, k)] +
					w[INDEX(i, j, k + 1)] - w[INDEX(i, j, k - 1)]) / 3;
				p[INDEX(i, j, k)] = 0;
			}
		}
	}
	set_boundaries(0, div); set_boundaries(0, p);
	for (l = 0; l < 20; l++)
	{
		for (k = 1; k <= N; k++) {
			for (j = 1; j <= N; j++) {
				for (i = 1; i <= N; i++) {
					p[INDEX(i, j, k)] = (div[INDEX(i, j, k)] +
						p[INDEX(i - 1, j, k)] + p[INDEX(i + 1, j, k)] +
						p[INDEX(i, j - 1, k)] + p[INDEX(i, j + 1, k)] +
						p[INDEX(i, j, k - 1)] + p[INDEX(i, j, k + 1)]) / 6;
				}
			}
		}
		set_boundaries(0, p);
	}
	for (k = 1; k <= N; k++) {
		for (j = 1; j <= N; j++) {
			for (i = 1; i <= N; i++) {
				u[INDEX(i, j, k)] -= (p[INDEX(i + 1, j, k)] - p[INDEX(i - 1, j, k)]) / 3 / h;
				v[INDEX(i, j, k)] -= (p[INDEX(i, j + 1, k)] - p[INDEX(i, j - 1, k)]) / 3 / h;
				w[INDEX(i, j, k)] -= (p[INDEX(i, j, k + 1)] - p[INDEX(i, j, k - 1)]) / 3 / h;
			}
		}
	}
	set_boundaries(1, u); set_boundaries(2, v);
}


void Fire::vorticity_confinement(float dt)
{
	int i, j, k, ijk;
	float* curlx = u0, * curly = v0, * curlz = w0, * curl = T0;	
	float dt0 = dt * vc_eps;
	float x, y, z;


	for (k = 1; k < N; k++) {
		for (j = 1; j < N; j++) {
			for (i = 1; i < N; i++) {
				ijk = INDEX(i, j, k);
				
				x = curlx[ijk] = (w[INDEX(i, j + 1, k)] - w[INDEX(i, j - 1, k)]) * 0.5f -
					(v[INDEX(i, j, k + 1)] - v[INDEX(i, j, k - 1)]) * 0.5f;

				
				y = curly[ijk] = (u[INDEX(i, j, k + 1)] - u[INDEX(i, j, k - 1)]) * 0.5f -
					(w[INDEX(i + 1, j, k)] - w[INDEX(i - 1, j, k)]) * 0.5f;

				
				z = curlz[ijk] = (v[INDEX(i + 1, j, k)] - v[INDEX(i - 1, j, k)]) * 0.5f -
					(u[INDEX(i, j + 1, k)] - u[INDEX(i, j - 1, k)]) * 0.5f;

				
				curl[ijk] = sqrtf(x * x + y * y + z * z);
			}
		}
	}

	for (k = 1; k < N; k++) {
		for (j = 1; j < N; j++) {
			for (i = 1; i < N; i++) {
				ijk = INDEX(i, j, k);
				float Nx = (curl[INDEX(i + 1, j, k)] - curl[INDEX(i - 1, j, k)]) * 0.5f;
				float Ny = (curl[INDEX(i, j + 1, k)] - curl[INDEX(i, j - 1, k)]) * 0.5f;
				float Nz = (curl[INDEX(i, j, k + 1)] - curl[INDEX(i, j, k - 1)]) * 0.5f;
				float len1 = 1.0f / (sqrtf(Nx * Nx + Ny * Ny + Nz * Nz) + 0.0000001f);
				Nx *= len1;
				Ny *= len1;
				Nz *= len1;
				u[ijk] += (Ny * curlz[ijk] - Nz * curly[ijk]) * dt0;
				v[ijk] += (Nz * curlx[ijk] - Nx * curlz[ijk]) * dt0;
				w[ijk] += (Nx * curly[ijk] - Ny * curlx[ijk]) * dt0;
			}
		}
	}
}


void Fire::vel_step(float dt)
{
	add_source(su, u, dt);
	add_source(sv, v, dt);
	add_source(sw, w, dt);
	add_buoyancy(dt);
	vorticity_confinement(dt);


	SWAP(u0, u); SWAP(v0, v); SWAP(w0, w);
	diffuse(1, u0, u, viscosity, dt);
	diffuse(2, v0, v, viscosity, dt);
	diffuse(3, w0, w, viscosity, dt);
	set_equilibrium();


	SWAP(u0, u); SWAP(v0, v); SWAP(w0, w);
	advect(1, u0, u, u0, v0, w0, dt);
	advect(2, v0, v, u0, v0, w0, dt);
	advect(3, w0, w, u0, v0, w0, dt);
	set_equilibrium();
	//std::cout << "Vel Step Done" << std::endl;
}

void Fire::dens_step(float dt)
{
	add_source(sd, d, dt);
	SWAP(d0, d);
	diffuse(0, d0, d, diffusion, dt);

	SWAP(d0, d);
	advect(0, d0, d, u, v, w, dt);

}

void Fire::dens_temp_step(float dt)
{
	add_source(sd, d, dt);
	add_source(sT, T, dt);
	SWAP(d0, d);
	diffuse(0, d0, d, diffusion, dt);
	SWAP(d0, d);
	SWAP(T0, T);
	advect_cool(0, d0, d, T0, T, u, v, w, dt);
	//std::cout << "DensTemp Step Done" << std::endl;
}

void Fire::step(float dt)
{
	vel_step(dt);
	dens_temp_step(dt);
}


void Fire::clear_buffer(float* x)
{
	for (int i = 0; i < SIZE; i++) {
		x[i] = 0.0f;
	}
}

void Fire::clear_sources(void)
{
	for (int i = 0; i < SIZE; i++) {
		sd[i] = su[i] = sv[i] = 0.0f;
	}
}