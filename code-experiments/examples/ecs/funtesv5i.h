/* --------------------------------------------------------- */
/* --- File: funtesv5i.h  -------- Author: ACMO          --- */
/* --------------------------------------------------------- */
/*   
    ECS for non-linear function minimization. 
    Copyright 2020 ACMO
    e-mail: alexandre.cesar .AT. ufma.br
    SOURCE: 
        https://github.com/cavalcantigor/ecs
*/

/********* CONSTANTES E FUNCOES COMUNS ************/

static const double PI = 3.1415926535897932384626433832795;
static const double schw_ans = 418.9828872724336861210758797824382781982;

double DistAbso(double x1[], double x2[], int n)
{
	double dist, d;
	int i;

	dist = 0.0F;

	for (i = 0; i < n; i++)
	{
		dist += abs(x1[i] - x2[i]);
	}
	return (dist / n);
}

double DistEucl(double x1[], double x2[], int n)
{
	double dist = 0.0, d;
	int i;

	dist = 0;

	for (i = 0; i < n; i++)
	{
		d = x1[i] - x2[i];
		dist += d * d;
	}

	return (dist);
}

/**************** FUNCOES COM DUAS VARIAVEIS ***********/
/*      fEaso(x1,x2)=-cos(x1)�cos(x2)�exp(-((x1-pi)^2+(x2-pi)^2));
        -100<=x(i)<=100, i=1:2.
        f(x1,x2)=-1; (x1,x2)=(pi,pi).
        fea= -cos(x(1)) * cos(x(2)) * exp(-(    (x(1)-pi)^2 +     (x(2)-pi)^2));
*/
double fEaso(double x[], int n)
{
	register int i;
	double sum;

	sum = -cos(x[0]) * cos(x[1]) * exp(-(pow((x[0] - PI), 2) + pow((x[1] - PI), 2)));
	FuncaoTeste.numAval++;
	return sum;
}
/*      fGold(x1,x2)=[1+(x1+x2+1)^2�(19-14�x1+3�x1^2-14�x2+6�x1�x2+3�x2^2)]�[30+(2�x1-3�x2)^2�(18-32�x1+12�x1^2+48�x2-36�x1�x2+27�x2^2)];
        -2<=x(i)<=2, i=1:2.
        f(x[1],x[2])=3; (x[1],x[2])=(0,-1).
*/
double fGold(double x[], int n)
{
	register int i;
	double sum;
	sum = (1 + pow((x[0] + x[1] + 1), 2) * (19 - 14 * x[0] + 3 * pow(x[0], 2) - 14 * x[1] + 6 * x[0] * x[1] + 3 * pow(x[1], 2))) * (30 + pow((2 * x[0] - 3 * x[1]), 2) * (18 - 32 * x[0] + 12 * pow(x[0], 2) + 48 * x[1] - 36 * x[0] * x[1] + 27 * pow(x[1], 2)));
	FuncaoTeste.numAval++;
	return sum;
}

/*
   De Jong's Test function 5 (DUAS VARI�VEIS)
 * Xi Range: [-65.536, 65.536]
 */
double jongf5(double x[], int n)
{
	static int a[2][30] = {
		{-32, -16, 0, 16, 32, -32, -16, 0, 16, 32, -32, -16, 0, 16, 32, -32, -16, 0, 16, 32, -32, -16, 0, 16, 32, -32, -16, 0, 16, 32},
		{-32, -32, -32, -32, -32, -16, -16, -16, -16, -16, 0, 0, 0, 0, 0, 16, 16, 16, 16, 16, 32, 32, 32, 32, 32, 0, 0, 0, 0, 0}};
	int i, j;
	double r1, r2;

	for (j = 1, r2 = 0.002; j <= 25; j++)
	{
		for (i = 0, r1 = j; i < 2; i++)
			r1 += pow(x[i] - a[i][j], 6);
		r2 += 1 / r1;
	}
	FuncaoTeste.numAval++;
	return (1 / r2);
}

/******************** FUNCOES TEST DE JONG GENERALIZADAS ***************************/
/************************* Esfera - Funcao 1 de De Jong ***********************/
double sphere(double x[], int n)
{
	register int i;
	double sum;
	for (i = 0, sum = 100.; i < n; i++)
		sum += ((x[i] - 1.0) * (x[i] - 1.0));
	FuncaoTeste.numAval++;
	return sum;
}

double jongf1(double x[], int n)
{
	register int i;
	double sum;

	for (i = 0, sum = 0.0; i < n; i++)
		sum += (x[i] * x[i]);
	FuncaoTeste.numAval++;
	return (sum);
}

/*********************** Rosenbrock - Funcao 2 de De Jong *********************/
double rosenbrock(double x[], int n)
{
	int i;

	double d = 0;
	for (i = 0; i < n - 1; i++)
	{
		double p = x[i + 1] - x[i] * x[i];
		double q = x[i] - 1;
		d += 100 * p * p + q * q;
	}
	FuncaoTeste.numAval++;
	return d;
}

/*
 * De Jong's Test function 3
 * STEP  Xi Range: [-5.12, 5.12]
 */
double jongf3(double x[], int n)
{
	double result;
	double tmp;
	int i;

	for (i = 0, result = 0.0; i < n; i++)
	{
		tmp = floor(x[i]);
		if (tmp >= 0)
			result += tmp;
		else
			result -= tmp; // esta correto !!!
	}
	FuncaoTeste.numAval++;
	return (result);
}

/*
   De Jong's Test function 4
 * Xi Range: [-1.28, 1.28] * FALTA RUIDO GAUSSIANO
 */
double jongf4(double x[], int n)
{
	double result, ruido;
	int i;

	ruido = 0.0;
	for (i = 0, result = 0.0; i < n; i++)
		result += (i * pow(x[i], 4) + ruido);
	FuncaoTeste.numAval++;
	return (result);
}

/***************************
	HARTMAN (6,4)
****************************/

double hartman6(double *x, int n)
{
	static double a[6][4] = {{10.00, 0.05, 3.00, 17.00},
							 {3.00, 10.00, 3.50, 8.00},
							 {17.00, 17.00, 1.70, 0.05},
							 {3.50, 0.10, 10.00, 10.00},
							 {1.70, 8.00, 17.00, 0.10},
							 {8.00, 14.00, 8.00, 14.00}};

	static double p[6][4] = {{0.1312, 0.2329, 0.2348, 0.4047},
							 {0.1696, 0.4135, 0.1451, 0.8828},
							 {0.5569, 0.8307, 0.3522, 0.8732},
							 {0.0124, 0.3736, 0.2883, 0.5743},
							 {0.8283, 0.1004, 0.3047, 0.1091},
							 {0.5886, 0.9991, 0.6650, 0.0381}};

	static double c[] = {1.0, 1.2, 3.0, 3.2};
	double som1 = 0, som2 = 0;
	int i, j;

	for (i = 0; i < 4; i++)
	{
		som1 = 0;
		for (j = 0; j < n; j++)
		{
			som1 = som1 + a[j][i] * (x[j] - p[j][i]) * (x[j] - p[j][i]);
		}
		som2 = som2 + c[i] * exp(-som1);
	}

	FuncaoTeste.numAval++;
	return (-som2);
}

/********************************************************

                * 1st ICEO test functions

*********************************************************/
/*
 * Problem 3: Shekel's Foxholes (4 variaveis)
 */

double shekel(double *x, int n)
{
	int PARAM = 10; // 10, 7, 5
	static double a[10][4] = {{4., 4., 4., 4.},
							  {1., 1., 1., 1.},
							  {8., 8., 8., 8.},
							  {6., 6., 6., 6.},
							  {3., 7., 3., 7.},
							  {2., 9., 2., 9.},
							  {5., 5., 3., 3.},
							  {8., 1., 8., 1.},
							  {6., 2., 6., 2.},
							  {7., 3.6, 7., 3.6}};

	static double c[] = {0.1, 0.2, 0.2, 0.4, 0.4,
						 0.6, 0.3, 0.7, 0.5, 0.5};

	register int i, j;
	double sp, h, result = 0.0;

	for (i = 0; i < PARAM; i++)
	{ //
		for (j = 0, sp = 0.0; j < n; j++)
		{
			sp += (x[j] - a[i][j]) * (x[j] - a[i][j]);
		}
		result += 1.0 / (sp + c[i]);
	}
	FuncaoTeste.numAval++;
	return (-result);
}

/*
 * Problem 5: The Langerman's functions
 *
 * f11(x)=-sum(c(i)�(exp(-1/pi�sum((x-A(i))^2))�cos(pi�sum((x-A(i))^2)))), 
 * 	FUNCAO COM 1.5 E COM 0.1 EM C(1) DAO DIFERENTES MINIMOS
 */

double langerman0(double *x, int n) /* Langerman's function */
{
	static double a[30][10] = {
		{9.681, 0.667, 4.783, 9.095, 3.517, 9.325, 6.544, 0.211, 5.122, 2.020},
		{9.400, 2.041, 3.788, 7.931, 2.882, 2.672, 3.568, 1.284, 7.033, 7.374},
		{8.025, 9.152, 5.114, 7.621, 4.564, 4.711, 2.996, 6.126, 0.734, 4.982},
		{2.196, 0.415, 5.649, 6.979, 9.510, 9.166, 6.304, 6.054, 9.377, 1.426},
		{8.074, 8.777, 3.467, 1.863, 6.708, 6.349, 4.534, 0.276, 7.633, 1.567},
		{7.650, 5.658, 0.720, 2.764, 3.278, 5.283, 7.474, 6.274, 1.409, 8.208},
		{1.256, 3.605, 8.623, 6.905, 0.584, 8.133, 6.071, 6.888, 4.187, 5.448},
		{8.314, 2.261, 4.224, 1.781, 4.124, 0.932, 8.129, 8.658, 1.208, 5.762},
		{0.226, 8.858, 1.420, 0.945, 1.622, 4.698, 6.228, 9.096, 0.972, 7.637},
		{7.305, 2.228, 1.242, 5.928, 9.133, 1.826, 4.060, 5.204, 8.713, 8.247},
		{0.652, 7.027, 0.508, 4.876, 8.807, 4.632, 5.808, 6.937, 3.291, 7.016},
		{2.699, 3.516, 5.874, 4.119, 4.461, 7.496, 8.817, 0.690, 6.593, 9.789},
		{8.327, 3.897, 2.017, 9.570, 9.825, 1.150, 1.395, 3.885, 6.354, 0.109},
		{2.132, 7.006, 7.136, 2.641, 1.882, 5.943, 7.273, 7.691, 2.880, 0.564},
		{4.707, 5.579, 4.080, 0.581, 9.698, 8.542, 8.077, 8.515, 9.231, 4.670},
		{8.304, 7.559, 8.567, 0.322, 7.128, 8.392, 1.472, 8.524, 2.277, 7.826},
		{8.632, 4.409, 4.832, 5.768, 7.050, 6.715, 1.711, 4.323, 4.405, 4.591},
		{4.887, 9.112, 0.170, 8.967, 9.693, 9.867, 7.508, 7.770, 8.382, 6.740},
		{2.440, 6.686, 4.299, 1.007, 7.008, 1.427, 9.398, 8.480, 9.950, 1.675},
		{6.306, 8.583, 6.084, 1.138, 4.350, 3.134, 7.853, 6.061, 7.457, 2.258},
		{0.652, 2.343, 1.370, 0.821, 1.310, 1.063, 0.689, 8.819, 8.833, 9.070},
		{5.558, 1.272, 5.756, 9.857, 2.279, 2.764, 1.284, 1.677, 1.244, 1.234},
		{3.352, 7.549, 9.817, 9.437, 8.687, 4.167, 2.570, 6.540, 0.228, 0.027},
		{8.798, 0.880, 2.370, 0.168, 1.701, 3.680, 1.231, 2.390, 2.499, 0.064},
		{1.460, 8.057, 1.336, 7.217, 7.914, 3.615, 9.981, 9.198, 5.292, 1.224},
		{0.432, 8.645, 8.774, 0.249, 8.081, 7.461, 4.416, 0.652, 4.002, 4.644},
		{0.679, 2.800, 5.523, 3.049, 2.968, 7.225, 6.730, 4.199, 9.614, 9.229},
		{4.263, 1.074, 7.286, 5.599, 8.291, 5.200, 9.214, 8.272, 4.398, 4.506},
		{9.496, 4.830, 3.150, 8.270, 5.079, 1.231, 5.731, 9.494, 1.883, 9.732},
		{4.138, 2.562, 2.532, 9.661, 5.611, 5.500, 6.886, 2.341, 9.699, 6.500}};

	static double c[] = {
		0.806,
		0.517,
		0.1,
		0.908,
		0.965,
		0.669,
		0.524,
		0.902,
		0.531,
		0.876,
		0.462,
		0.491,
		0.463,
		0.714,
		0.352,
		0.869,
		0.813,
		0.811,
		0.828,
		0.964,
		0.789,
		0.360,
		0.369,
		0.992,
		0.332,
		0.817,
		0.632,
		0.883,
		0.608,
		0.326};

	register int i;
	double Sum, dist;
	Sum = 0.0F;
	for (i = 0; i < 5; i++)
	{
		dist = DistEucl(x, a[i], n);
		Sum -= c[i] * (exp(-dist / PI) * cos(PI * dist));
	}
	FuncaoTeste.numAval++;
	return (Sum);
}

double langerman1(double *x, int n) /* Langerman's function */
{
	static double a[30][10] = {
		{9.681, 0.667, 4.783, 9.095, 3.517, 9.325, 6.544, 0.211, 5.122, 2.020},
		{9.400, 2.041, 3.788, 7.931, 2.882, 2.672, 3.568, 1.284, 7.033, 7.374},
		{8.025, 9.152, 5.114, 7.621, 4.564, 4.711, 2.996, 6.126, 0.734, 4.982},
		{2.196, 0.415, 5.649, 6.979, 9.510, 9.166, 6.304, 6.054, 9.377, 1.426},
		{8.074, 8.777, 3.467, 1.863, 6.708, 6.349, 4.534, 0.276, 7.633, 1.567},
		{7.650, 5.658, 0.720, 2.764, 3.278, 5.283, 7.474, 6.274, 1.409, 8.208},
		{1.256, 3.605, 8.623, 6.905, 0.584, 8.133, 6.071, 6.888, 4.187, 5.448},
		{8.314, 2.261, 4.224, 1.781, 4.124, 0.932, 8.129, 8.658, 1.208, 5.762},
		{0.226, 8.858, 1.420, 0.945, 1.622, 4.698, 6.228, 9.096, 0.972, 7.637},
		{7.305, 2.228, 1.242, 5.928, 9.133, 1.826, 4.060, 5.204, 8.713, 8.247},
		{0.652, 7.027, 0.508, 4.876, 8.807, 4.632, 5.808, 6.937, 3.291, 7.016},
		{2.699, 3.516, 5.874, 4.119, 4.461, 7.496, 8.817, 0.690, 6.593, 9.789},
		{8.327, 3.897, 2.017, 9.570, 9.825, 1.150, 1.395, 3.885, 6.354, 0.109},
		{2.132, 7.006, 7.136, 2.641, 1.882, 5.943, 7.273, 7.691, 2.880, 0.564},
		{4.707, 5.579, 4.080, 0.581, 9.698, 8.542, 8.077, 8.515, 9.231, 4.670},
		{8.304, 7.559, 8.567, 0.322, 7.128, 8.392, 1.472, 8.524, 2.277, 7.826},
		{8.632, 4.409, 4.832, 5.768, 7.050, 6.715, 1.711, 4.323, 4.405, 4.591},
		{4.887, 9.112, 0.170, 8.967, 9.693, 9.867, 7.508, 7.770, 8.382, 6.740},
		{2.440, 6.686, 4.299, 1.007, 7.008, 1.427, 9.398, 8.480, 9.950, 1.675},
		{6.306, 8.583, 6.084, 1.138, 4.350, 3.134, 7.853, 6.061, 7.457, 2.258},
		{0.652, 2.343, 1.370, 0.821, 1.310, 1.063, 0.689, 8.819, 8.833, 9.070},
		{5.558, 1.272, 5.756, 9.857, 2.279, 2.764, 1.284, 1.677, 1.244, 1.234},
		{3.352, 7.549, 9.817, 9.437, 8.687, 4.167, 2.570, 6.540, 0.228, 0.027},
		{8.798, 0.880, 2.370, 0.168, 1.701, 3.680, 1.231, 2.390, 2.499, 0.064},
		{1.460, 8.057, 1.336, 7.217, 7.914, 3.615, 9.981, 9.198, 5.292, 1.224},
		{0.432, 8.645, 8.774, 0.249, 8.081, 7.461, 4.416, 0.652, 4.002, 4.644},
		{0.679, 2.800, 5.523, 3.049, 2.968, 7.225, 6.730, 4.199, 9.614, 9.229},
		{4.263, 1.074, 7.286, 5.599, 8.291, 5.200, 9.214, 8.272, 4.398, 4.506},
		{9.496, 4.830, 3.150, 8.270, 5.079, 1.231, 5.731, 9.494, 1.883, 9.732},
		{4.138, 2.562, 2.532, 9.661, 5.611, 5.500, 6.886, 2.341, 9.699, 6.500}};

	// OBS c[2] = 1.5 difere de c[2] = 1.4
	static double c[] = {
		0.806,
		0.517,
		1.5,
		0.908,
		0.965,
		0.669,
		0.524,
		0.902,
		0.531,
		0.876,
		0.462,
		0.491,
		0.463,
		0.714,
		0.352,
		0.869,
		0.813,
		0.811,
		0.828,
		0.964,
		0.789,
		0.360,
		0.369,
		0.992,
		0.332,
		0.817,
		0.632,
		0.883,
		0.608,
		0.326};

	register int i;
	double Sum, dist;
	Sum = 0;
	// n <= 10
	for (i = 0; i < 5; i++)
	{
		dist = DistEucl(x, a[i], n);
		Sum -= c[i] * (exp(-dist / PI) * cos(PI * dist));
	}
	FuncaoTeste.numAval++;
	return (Sum);
}

/*      ICEO Michalewitz  - segunda implementacao MESMA DO BILCHEV
        f12(x)=-sum(sin(x(i))�(sin(i�x(i)^2/pi))^(2�m)), i=1:n, m=10;
        0<=x(i)<=pi.
        global minimum
        f(x)=-4.687 (n=5); x(i)=???, i=1:n.
        f(x)=-9.66 (n=10); x(i)=???, i=1:n.
*/

double michalewicz(double x[], int n) /* Michalewitz */
{
	double u, m = 10.0F;
	int i;

	u = 0;

	for (i = 0; i < n; i++)
		u = u + sin(x[i]) * pow(sin((i + 1) * x[i] * x[i] / PI), 2.0 * m);

	FuncaoTeste.numAval++;
	return (-u);
}

/********************************************************

	* OUTRAS FUNCOES GENERALIZADAS FAMOSAS

*********************************************************/
/*************
        schwefel
**************/

double schwefel7(double x[], int n)
{
	int i;
	double d;
	for (i = 0, d = 0; i < n; i++)
		d -= x[i] * sin(sqrt(fabs(x[i])));

	FuncaoTeste.numAval++;
	return (d + schw_ans * (double)n);
}

/*************
        rastrigin
**************/

double rastrigin(double x[], int n)
{
	int i;
	double d;
	for (i = 0, d = 0; i < n; i++)
		d += x[i] * x[i] - 10 * cos(2 * PI * x[i]);
	FuncaoTeste.numAval++;
	return (10 * n + d);
}

/*************
        griewank
**************/

double griewank(double x[], int n)
{
	int i;
	double d, e;
	for (i = 0, d = e = 1; i < n; i++)
	{
		d += x[i] * x[i] / 4000;
		e *= cos(x[i] / sqrt(i + 1));
	}
	FuncaoTeste.numAval++;
	return (d - e);
}

// F4: Ridge [-64 64) : 0.0
double ridge(double *x, int n)
{
	int i, j;
	double d = 0;
	for (i = 0; i < n; i++)
	{
		double p = 0;
		for (j = 0; j <= i; j++)
			p += x[j];
		d += p * p;
	}
	FuncaoTeste.numAval++;
	return d;
}

/***********************

                FUNCOES GENERALIZADAS POUCO CONHECIDAS

*****************************/

/* 1.9 Ackley's Path function 10
Ackley's Path [Ack87] is a widely used multimodal test function.
function definition
f10(x)=-a�exp(-b�sqrt(1/n�sum(x(i)^2)))-exp(1/n�sum(cos(c�x(i))))+a+exp(1);
a=20; b=0.2; c=2�pi; i=1:n;
-32.768<=x(i)<=32.768.
global minimum
f(x)=0; x(i)=0, i=1:n.
*****************************/

double Ackley(double x[], int n)
{
	register int i;
	double sum;

	double a = 20.0F, b = 0.2F, c = 2.0F * PI;
	double sum1, sum2;

	for (i = 0, sum1 = 0, sum2 = 0; i < n; i++)
	{
		sum1 += pow(x[i], 2);
		sum2 += cos(c * x[i]);
	}
	sum = a + exp(1.0F) - a * exp(-b * sqrt(1.0F / n * sum1)) - exp(1.0F / n * sum2);
	//sum = (20 + exp(1.0F) + (-20)*exp((-0.2)*sqrt(sum1/n)) - exp(sum2/n));

	FuncaoTeste.numAval++;
	return sum;
}

double zakharov(double *x, int n)
{
	double som1 = 0.0F, som2 = 0.0F;
	int j;

	for (j = 0; j < n; j++)
	{
		som1 = som1 + (x[j] * x[j]);
		som2 = som2 + (0.5 * (j + 1) * x[j]);
	}

	FuncaoTeste.numAval++;

	return (som1 + som2 * som2 + pow(som2, 4));
}

/*
 * Unary test function
 * Problem 0: The Unary test function
 * Xi Range: {0,1}
 * Optimial Value: 0
 */
double unaria(double x[], int n)
{
	register int i;
	double Sum;
	Sum = n;
	for (i = 0; i < n; i++)
	{
		Sum -= x[i];
	}
	FuncaoTeste.numAval++;
	return (Sum);
}

/*
 * BUMP test function
 * Xi Range: ???
 * Optimial Value: ???
 */

/* RESTRICAO  */
void constraints_bump(double *x, double *con1, double *con2, int NVARS)
{
	int ii;

	*con1 = 1.0;
	*con2 = 0.0;
	for (ii = 0; ii < NVARS; ii++)
	{
		*con1 = *con1 * x[ii];
		*con2 = *con2 + x[ii];
	}
	return;
}

/* FUNCAO  */

double bump(double *x, int NVARS)
{
	int ii;
	double sumc4, prodc2, sumsq, obj;

	double resp, ress, pena = NVARS;

	sumc4 = 0.0;
	prodc2 = 1.0;
	sumsq = 0.0;
	for (ii = 0; ii < NVARS; ii++)
	{
		sumc4 = sumc4 + cos(x[ii]) * cos(x[ii]) * cos(x[ii]) * cos(x[ii]);
		prodc2 = prodc2 * cos(x[ii]) * cos(x[ii]);
		sumsq = sumsq + (double)(ii + 1) * x[ii] * x[ii];
	}
	sumsq = sqrt(sumsq);
	obj = (sumc4 - 2.0 * prodc2) / sumsq;
	if (obj < 0.0)
		obj = -obj;
	//ACMO penaliza se fora de restri��o
	constraints_bump(x, &resp, &ress, NVARS);
	obj = obj + pena * ((resp <= 0.75) + (ress >= 15 * NVARS / 2));
	FuncaoTeste.numAval++;
	return (obj);
}

/**********
1.3 Rotated hyper-ellipsoid function
An extension of the axis parallel hyper-ellipsoid is Schwefel's function1.2. With respect to the coordinate axes, this function produces rotated hyper-ellipsoids. It is continuos, convex and unimodal.

function definition
f1b(x)=sum(sum(x(j)^2), j=1:i), i=1:n;
-65.536<=x(i)<=65.536.

global minimum
f(x)=0; x(i)=0, i=1:n.
*************************/
double rotated(double x[], int n)
{
	double sum = 0;
	int i, j;
	for (i = 0; i < n; i++)
		for (j = 0; j < i; j++)
			sum = sum + pow(x[j], 2);

	FuncaoTeste.numAval++;
	return (sum);
}

/* ***************	ESTRUTURA DAS FUNCOES TESTE ***************** */

enum
{
	fea = 0,
	gol,
	dj5,
	she,
	lan0,
	mic,
	sph,
	dj1,
	ros,
	dj3,
	dj4,
	sch,
	ras,
	gri,
	ack,
	har,
	lan1,
	zak,
	rid,
	uno,
	rot,
	bmp
} enumFunc;

static double (*funccod[])(double[], int n) =
	{fEaso, fGold, jongf5,
	 shekel, langerman0, michalewicz,
	 sphere, jongf1, rosenbrock, jongf3, jongf4,
	 schwefel7, rastrigin, griewank, Ackley,
	 hartman6, langerman1, zakharov,
	 ridge, unaria, rotated, bump};

struct DadosFuncoes
{
	char nom[5];
	double inf;
	double sup;
	double opt[3];
} FuncoesTeste[] = {{"Fea", -100.000F, 100.000F, -1.0F},
					{"Gol", -2.0F, 2.0F, 3.0F},
					{"DjF5", -65.536F, 65.536F, 1.0F},
					{"She", -10.0F, 10.0F, -10.53641F, -10.40294F, -10.1532F}, // depende de PARAM
					{"Lan0", 0.0F, 10.0F, -0.96F},
					{"Mic", 0.0F, PI, -4.68765817908815, -9.66015167547059F},
					{"Sph", -5.12F, 5.12F, 100.0F},
					{"DjF1", -5.12F, 5.12F, 0.0F},
					{"Ros", -30.0F, 30.0F, 0.0F},
					{"DjF3", -5.12F, 5.12F, 0.0F},
					{"DjF4", -1.28F, 1.28F, 0.0F},
					{"Sch", -100.0F, 100.0F, 0.0F},
					{"Ras", -5.12, 5.12, 0.0F},
					{"Gri", -600.0, 600.0, 0.0F},
					{"Ack", -32.0F, 32.0F, 0.0F},
					{"Har", 0.0F, 1.0F, -3.32237F},
					{"Lan1", 0.0F, 10.0F, -1.40F},
					{"Zak", -5.0F, 10.0F, 0.0F},
					{"Rid", -64.0F, 64.0F, 0.0F},
					{"Uno", 0.0F, 1.0F, 0.0F},
					{"Rot", -65.536F, 65.536F, 0.0F},
					{"Bump", -10.0F, 10.0F, 0.0F}};
