// basic file operations
#include <thread>
#include <fstream>
#include <iostream>
#include <cmath>
using namespace std;
void foo(int t, double dx, double dy, double dt, double dz, ofstream& file, int
zmin, int zmax, double***& Ucurr, double***& Uprev1, double***& Uprev2,
double***& pos, double extra_heat, double Ax, double Ay, double Az, double
Ax_inner, double Ay_inner, double Az_inner, int Nx, int Ny, int Nz)
{
	//For the data storage I just continously overwrite the same 2 arrays, I just therefore repeat the same block of code twice, and save to each alternately (t%2==0).
	//This isnt best practice, if I was to do it again I'd probably combine Uprev1 and Uprev2 and assign to each alternately using a similar modulo thing within the block.
	if (t % 2 == 0)
	{
		for (int i = 1;i <= Nx - 1;i++)
		{
		for (int j = 1;j <= Ny - 1;j++)
		{
		for (int k = zmin/dz;k <zmax/dz;k++)
		{
			Ucurr[i][j][k] = 0;
			//this horrible block calculates the spatial double derivative of the temperature according to a formula derived from Taylor expansion, in the case of the inner region
			//it adds an amount of heat equivalent to the Joule heating effect
			//ax,ay,az includes a time step dt, so the time derivative is also numerically integrated forward Euler fashion 
			if (pos[i][j][k] == 1)
			{
			Ucurr[i][j][k] = Uprev2[i][j][k] + Ax * (Uprev2[i +
			1][j][k] + Uprev2[i - 1][j][k] - 2 * Uprev2[i][j][k]) + Ay
			* (Uprev2[i][j + 1][k] + Uprev2[i][j - 1][k] - 2 * Uprev2
			[i][j][k]) + Az * (Uprev2[i][j][k + 1] + Uprev2[i][j][k -
			1] - 2 * Uprev2[i][j][k]);
			}
			else if (pos[i][j][k] == 2)
			{
			Ucurr[i][j][k] = Uprev2[i][j][k] + Ax_inner *
			(Uprev2[i + 1][j][k] + Uprev2[i - 1][j][k] - 2 * Uprev2[i]
			[j][k]) + Ay_inner * (Uprev2[i][j + 1][k] + Uprev2[i][j -
			1][k] - 2 * Uprev2[i][j][k]) + Az_inner * (Uprev2[i][j][k +
			1] + Uprev2[i][j][k - 1] - 2 * Uprev2[i][j][k]) +
			extra_heat;
			}
			Uprev1[i][j][k] = Ucurr[i][j][k];
			//only saves every 10 seconds, plenty for the type of plotting being done
			if (t % 10000 == 0)
			{
			file << t * dt << "," << i * dx << "," << j * dy << ","
			<< k * dz << "," << Ucurr[i][j][k] << endl;
			}
		}
		}
		}
	}
	//same as previous block
	else
	{
		for (int i = 1;i <= Nx - 1;i++)
		{
		for (int j = 1;j <= Ny - 1;j++)
		{
		for (int k = zmin/dz;k <(zmax/dz);k++)
		{
			Ucurr[i][j][k] = 0;
			if (pos[i][j][k] == 1)
			{
			Ucurr[i][j][k] = Uprev1[i][j][k] + Ax * (Uprev1[i + 1]
			[j][k] + Uprev1[i - 1][j][k] - 2 * Uprev1[i][j][k]) + Ay *
			(Uprev1[i][j + 1][k] + Uprev1[i][j - 1][k] - 2 * Uprev1[i]
			[j][k]) + Az * (Uprev1[i][j][k + 1] + Uprev1[i][j][k - 1] -
			2 * Uprev1[i][j][k]);
			}
			else if (pos[i][j][k] == 2)
			{
			Ucurr[i][j][k] = Uprev1[i][j][k] + Ax_inner * (Uprev1[i
			+ 1][j][k] + Uprev1[i - 1][j][k] - 2 * Uprev1[i][j][k]) +
			Ay_inner * (Uprev1[i][j + 1][k] + Uprev1[i][j - 1][k] - 2 *
			Uprev1[i][j][k]) + Az_inner * (Uprev1[i][j][k + 1] +
			Uprev1[i][j][k - 1] - 2 * Uprev1[i][j][k]) + extra_heat;
			}
			Uprev2[i][j][k]=Ucurr[i][j][k];
			if (t % 10000==0)
			{
			file << t * dt << "," << i * dx << "," << j * dy << "," << k * dz << "," << Ucurr[i][j][k] << endl;
			}
		}
		}
		}
	}
}
int main()
{
int x_max = 250;
int y_max = 250;
int z_max = 200;
int t_max = 9000;
double dx = 5;
double dy = 5;
double dz = 2;
double dt = 0.05;
int Nx = int(x_max / dx) + 1;
int Ny = int(y_max / dy) + 1;
int Nz = int(z_max / dz) + 1;
int Nt = int(t_max / dt) + 1;
double*** Uprev1;
double*** Uprev2;
double*** Ucurr;
double*** pos;
Uprev1 = new double** [Nx];
Uprev2 = new double** [Nx];
Ucurr = new double** [Nx];
pos = new double** [Nx];
//Allocating the column space in heap dynamically
for (int i = 0; i < Nx; i++) {
Uprev1[i] = new double* [Ny];
Uprev2[i] = new double* [Ny];
Ucurr[i] = new double* [Ny];
pos[i] = new double* [Ny];
for (int j = 0; j < Ny; j++)
{
Uprev1[i][j] = new double[Nz];
Uprev2[i][j] = new double[Nz];
Ucurr[i][j] = new double[Nz];
pos[i][j] = new double[Nz];
for (int k = 0; k < Nx; k++)
{
Uprev1[i][j][k] = 0;
Uprev2[i][j][k] = 0;
Ucurr[i][j][k] = 0;
pos[i][j][k] = 0;
}
}
}
int inner_radius = 5;
int outer_radius = 120;
double alpha = 0.85;
double alpha_inner = 2.9;
double extra_heat = 8.27 * dt;
double Ax = (alpha * dt) / pow(dx, 2);
double Ay = (alpha * dt) / pow(dy, 2);
double Az = (alpha * dt) / pow(dz, 2);
double Ax_inner = (alpha_inner * dt) / pow(dx, 2);
double Ay_inner = (alpha_inner * dt) / pow(dy, 2);
double Az_inner = (alpha_inner * dt) / pow(dz, 2);
cout << "beepbopp";
//fix each point as heating element, sand, or free space
for (int i = 0;i < Nx;i++)
{
for (int j = 0;j <Ny;j++)
{
for (int k = 0;k <Nz ;k++)
{
if (pow((i * dx - x_max / 2), 2) + pow((j * dy - y_max / 2), 2)
<= pow(outer_radius, 2))
{
pos[i][j][k] = 1;
}
if ((pow((i * dx - x_max / 2), 2) + pow((j * dy - y_max / 2),
2) <= pow(inner_radius, 2)) && (k*dz >= z_max/2))
{
pos[i][j][k] = 2;
}
if (pow((i * dx - x_max / 2), 2) + pow((j * dy - y_max / 2), 2)
>= pow(outer_radius, 2))
{
pos[i][j][k] = 0;
}
}
}
}
//save to 4 files, as otherwise the files would get too large for convenience
//ended up reducing number of timesteps, so it wasnt that necessary in the end
ofstream file1;
file1.open("BucketData1.csv");
ofstream file2;
file2.open("BucketData2.csv");
ofstream file3;
file3.open("BucketData3.csv");
ofstream file4;
file4.open("BucketData4.csv");
file1 << "t,x,y,z,temperature" << endl;
file2 << "t,x,y,z,temperature" << endl;
file3 << "t,x,y,z,temperature" << endl;
file4 << "t,x,y,z,temperature" << endl;
//the method of numerical integration used never writes to the same point at the same time, so run the code on 4-threads to increase throughput
for (int t = 1; t < Nt; t++)
{
cout << "t = " << t * dt;
std::thread t1(foo, t, dx, dy, dt, dz, std::ref(file1), 1, 50, std::ref
(Ucurr), std::ref(Uprev1), std::ref(Uprev2), std::ref(pos),
extra_heat, Ax, Ay, Az, Ax_inner, Ay_inner, Az_inner, Nx, Ny, Nz);
std::thread t2(foo, t, dx, dy, dt, dz, std::ref(file2), 50, 100,
std::ref(Ucurr), std::ref(Uprev1), std::ref(Uprev2), std::ref(pos),
extra_heat, Ax, Ay, Az, Ax_inner, Ay_inner, Az_inner, Nx, Ny, Nz);
std::thread t3(foo, t, dx, dy, dt, dz, std::ref(file3), 100,150,
std::ref(Ucurr), std::ref(Uprev1), std::ref(Uprev2), std::ref(pos),
extra_heat, Ax, Ay, Az, Ax_inner, Ay_inner, Az_inner, Nx, Ny, Nz);
std::thread t4(foo, t, dx, dy, dt, dz, std::ref(file4), 150, 198,
std::ref(Ucurr), std::ref(Uprev1), std::ref(Uprev2), std::ref(pos),
extra_heat, Ax, Ay, Az, Ax_inner, Ay_inner, Az_inner, Nx, Ny, Nz);
t1.join();
t2.join();
t3.join();
t4.join();
}
file1.close();
file2.close();
file3.close();
file4.close();
return 0;
}
