#include "stdafx.h"
#include <iostream>
#include <inttypes.h>
#include <math.h>

constexpr uint8_t NETLIST_PORT_MAX = 4;
constexpr uint8_t NETLIST_PARAM_MAX = 5;
constexpr double MNA_JACOBIAN_STEP = 1e-6;

class utils_t {
public:
	template <typename T> static void copy_vect(T *, T*, uint16_t);
	template <typename T> static void copy_matrix(T**, T**, uint16_t, uint16_t);
	template <typename T> static void zero_vect(T *, uint16_t);
	template <typename T> static void zero_matrix(T**, uint16_t, uint16_t);
	template <typename T> static T mult_matrix_single_row_col(T**, T**, uint16_t, uint16_t, uint16_t);
	template <typename T> static void mult_matrix(T**, T**, T**, uint16_t, uint16_t, uint16_t);
	template <typename T> static void mult_mat_vect(T**, T*, T*, uint16_t, uint16_t);
};
class netlist_t {
public:
	uint8_t nodes;
	uint8_t components;
	enum port_order_t : uint8_t;
	enum data_order_t : uint8_t;
	enum component_order_t : uint8_t;
	struct data_t;
	struct port_t;
	data_t *value;
	port_t *port;
	component_order_t *component_list;
	netlist_t(uint8_t, uint8_t);
	~netlist_t(void);
};
class component_t {
public:
	static void insert(double **, double(**x)(netlist_t *, uint8_t, uint8_t, double *, double *), double(**y)(netlist_t *, uint8_t, uint8_t, double *, double *), uint8_t *, netlist_t *, uint8_t, uint8_t *, double *);
	static void insert_from_memory(double(**x)(netlist_t *, uint8_t, uint8_t, double *, double *), uint8_t);
	static void insert_copy_s(double(**x)(netlist_t *, uint8_t, uint8_t, double *, double *), uint8_t);
	static double copy_s(netlist_t *, uint8_t, uint8_t, double *, double*);
	static double const_val(netlist_t *, uint8_t, uint8_t, double *, double *);
	static const uint8_t v_mna_add = 1;
	static double v(netlist_t *, uint8_t, uint8_t, double *, double*);
	static const uint8_t resistor_mna_add = 0;
	static const uint8_t capacitor_mna_add = 0;
	static const uint8_t inductor_mna_add = 0;
	static const uint8_t diode_r_mna_add = 1;
	static double diode_vd(netlist_t *, uint8_t, uint8_t, double *, double*);
	static double diode_r_current(netlist_t *, uint8_t, uint8_t, double *, double *);

	//fet
	static double fet_current(netlist_t *, uint8_t, uint8_t, double *, double *);
};
class mna_t {
public:
	uint8_t nodes;
	uint8_t hidden;
	uint8_t size;
	uint8_t components;
	netlist_t *source;
	double **g;
	double **r;
	double **jn;
	double **jni;
	double **jni_old;
	double *z;
	double *s;
	double *s_old;
	double *s_diff;
	double(**y)(netlist_t *ref, uint8_t netlist_pos, uint8_t call_pos, double *input, double *memory);
	double(**x)(netlist_t *ref, uint8_t netlist_pos, uint8_t call_pos, double *input, double *memory);
	netlist_t *y_source;
	double *f_memory;
	uint8_t *f_netlist_pos;
	mna_t(netlist_t *);
	~mna_t(void);
	void generate(void);
	void update_z(void);
	double eval_row(uint8_t);
	double eval_row(uint8_t, double *);
	double* eval(double *);
	void iterate_once(void);
	void iterate_n(uint16_t);
	double** jacobian(void);
};
class math_t {
public:
	void sub(double *, double *, double *, uint8_t);
	void add(double *, double *, double *, uint8_t);
	double** invert(double **, double **, uint8_t);
private:
	int LUPDecompose(double **A, int N, double Tol, int *P);
	void LUPSolve(double **A, int *P, double *b, int N, double *x);
	void LUPInvert(double **A, int *P, int N, double **IA);
	double LUPDeterminant(double **A, int *P, int N);
}math;

template <typename T> void utils_t::copy_vect(T *from, T *to, uint16_t length) {
	while (length) {
		length--;
		to[length] = from[length];
	}
}
template <typename T> void utils_t::copy_matrix(T** from, T** to, uint16_t i, uint16_t j) {
	while (i) {
		i--;
		copy_vect(from[i], to[i], j);
	}
}
template <typename T> void utils_t::zero_vect(T *vect, uint16_t length) {
	while (length) {
		length--;
		vect[length] = 0;
	}
}
template <typename T> void utils_t::zero_matrix(T** M, uint16_t i, uint16_t j) {
	while (i) {
		i--;
		zero_vect(M[i], j);
	}
}
template <typename T> T utils_t::mult_matrix_single_row_col(T** M, T** N, uint16_t i, uint16_t j, uint16_t len) {
	T acc = 0;
	while (len) {
		len--;
		acc += M[i][len] * N[len][j];
	}
	return acc;
}
template <typename T> void utils_t::mult_matrix(T** M, T** N, T** OUT, uint16_t i1, uint16_t j2, uint16_t len) {
	uint16_t x;
	while (i1){
		i1--;
		x = j2;
		while (x){
			x--;
			OUT[i1][x] = mult_matrix_single_row_col(M, N, i1, x, len);
		}
	}
}
template <typename T> void utils_t::mult_mat_vect(T** M, T* N, T* OUT, uint16_t i_max, uint16_t j_max) {
	uint16_t j;
	while (i_max){
		i_max--;
		OUT[i_max] = (T)0;
		j = j_max;
		while (j){
			j--;
			OUT[i_max] += M[i_max][j] * N[j];
		}
	}
}



enum netlist_t::port_order_t : uint8_t {
	//2 port
	node_main = 0,
	node_dest,
	//opamp
	vp = 0,
	vn,
	vcc,
	vdd,
	vo,
	//jfet

	//bjt
	vc = 0,
	vb,
	ve,
	//mosfet_
	fet_n_drain = 0,
	fet_n_gate,
	fet_n_source
};
enum netlist_t::data_order_t : uint8_t{
	//voltage
	voltage = 0,
	//current
	current = 0,
	//resistor
	resistance = 0,
	//capacitor
	capacitance = 0,
	//inductor
	inductance = 0,
	//op amp
	a = 0,
	gain = 0,
	open_loop_gain = 0,
	//diode
	is = 0,
	vt,
	//jfet

	//bjt

	//mosfet
	fet_w = 0,
	fet_l,
	fet_kp,
	fet_vth,
	fet_lambda
};
enum netlist_t::component_order_t : uint8_t{
	voltage_src,
	current_src,
	resistor,
	inductor,
	capacitor,
	opamp_1,
	opamp_2,
	opamp_3,
	diode_r,
	diode_s,
	diode_z,
	jfet_n,
	jfet_p,
	bjt_npn,
	bjt_pnp,
	fet_n,
	fet_p
};

struct netlist_t::data_t {
	double item[NETLIST_PARAM_MAX];
};
struct netlist_t::port_t {
	uint8_t item[NETLIST_PORT_MAX];
};

netlist_t::netlist_t(uint8_t amnt_nodes, uint8_t amnt_comp) {
	nodes = amnt_nodes;
	components = amnt_comp;
	if (amnt_nodes) {
		value = new data_t[amnt_comp];
		port = new port_t[amnt_comp];
		component_list = new component_order_t[amnt_comp];
	}
}
netlist_t::~netlist_t() {
	uint8_t pos = nodes;
	if (nodes) {
		delete[] value;
		delete[] port;
		delete[] component_list;
	}
}



void component_t::insert(double **g, double(**x)(netlist_t *, uint8_t, uint8_t, double *, double *), double (**y)(netlist_t *, uint8_t, uint8_t, double *, double *), uint8_t *f_netlist_pos, netlist_t *netlist, uint8_t pos, uint8_t *hidden, double *f_memory) {
	switch (netlist->component_list[pos]){
	case(netlist_t::component_order_t::voltage_src):
		if (netlist->port[pos].item[netlist_t::port_order_t::node_main]) {
			g[netlist->port[pos].item[netlist_t::port_order_t::node_main] - 1][*hidden - 1] = 1;
			g[*hidden - 1][netlist->port[pos].item[netlist_t::port_order_t::node_main] - 1] = 1;
		}
		if (netlist->port[pos].item[netlist_t::port_order_t::node_dest]) {
			g[netlist->port[pos].item[netlist_t::port_order_t::node_dest] - 1][*hidden - 1] = -1;
			g[*hidden - 1][netlist->port[pos].item[netlist_t::port_order_t::node_dest] - 1] = -1;
		}
		if (netlist->port[pos].item[netlist_t::port_order_t::node_main] || netlist->port[pos].item[netlist_t::port_order_t::node_dest]) {
			x[*hidden - 1] = component_t::v;
			y[*hidden - 1] = component_t::copy_s;
			f_netlist_pos[*hidden - 1] = pos;
			(*hidden)--;
		}
		break;
	case(netlist_t::component_order_t::current_src):
		//remember to += memory (current) 
		break;
	case(netlist_t::component_order_t::resistor):
		if (netlist->port[pos].item[netlist_t::port_order_t::node_main]) {
			g[netlist->port[pos].item[netlist_t::port_order_t::node_main] - 1][netlist->port[pos].item[netlist_t::port_order_t::node_main] - 1] = 1 / netlist->value[pos].item[netlist_t::data_order_t::resistance];
			if (netlist->port[pos].item[netlist_t::port_order_t::node_dest]) {
				g[netlist->port[pos].item[netlist_t::port_order_t::node_main] - 1][netlist->port[pos].item[netlist_t::port_order_t::node_dest] - 1] = -1 / netlist->value[pos].item[netlist_t::data_order_t::resistance];
			}
		}
		if (netlist->port[pos].item[netlist_t::port_order_t::node_dest]) {
			g[netlist->port[pos].item[netlist_t::port_order_t::node_dest] - 1][netlist->port[pos].item[netlist_t::port_order_t::node_dest] - 1] = 1 / netlist->value[pos].item[netlist_t::data_order_t::resistance];
			if (netlist->port[pos].item[netlist_t::port_order_t::node_main]) {
				g[netlist->port[pos].item[netlist_t::port_order_t::node_dest] - 1][netlist->port[pos].item[netlist_t::port_order_t::node_main] - 1] = -1 / netlist->value[pos].item[netlist_t::data_order_t::resistance];
			}
		}
		break;
	case(netlist_t::component_order_t::capacitor):

		break;
	case(netlist_t::component_order_t::inductor):

		break;
	case(netlist_t::component_order_t::diode_r):
		if (netlist->port[pos].item[netlist_t::port_order_t::node_main]) {
			g[netlist->port[pos].item[netlist_t::port_order_t::node_main] - 1][*hidden - 1] = 1;
			g[*hidden - 1][netlist->port[pos].item[netlist_t::port_order_t::node_main] - 1] = 1;
		}
		if (netlist->port[pos].item[netlist_t::port_order_t::node_dest]) {
			g[netlist->port[pos].item[netlist_t::port_order_t::node_dest] - 1][*hidden - 1] = -1;
			g[*hidden - 1][netlist->port[pos].item[netlist_t::port_order_t::node_dest] - 1] = -1;
		}
		if (netlist->port[pos].item[netlist_t::port_order_t::node_main] || netlist->port[pos].item[netlist_t::port_order_t::node_dest]) {
			//x[*hidden - 1] = component_t::diode_vd;
			x[*hidden - 1] = component_t::copy_s;
			y[*hidden - 1] = component_t::diode_r_current;
			f_netlist_pos[*hidden - 1] = pos;
			(*hidden)--;
		}
		break;
	case (netlist_t::component_order_t::fet_n):
		if (netlist->port[pos].item[netlist_t::port_order_t::fet_n_drain]) {
			g[netlist->port[pos].item[netlist_t::port_order_t::fet_n_drain] - 1][*hidden - 1] = 1;
			g[*hidden - 1][netlist->port[pos].item[netlist_t::port_order_t::fet_n_drain] - 1] = 1;
		}
		if (netlist->port[pos].item[netlist_t::port_order_t::fet_n_source]) {
			g[netlist->port[pos].item[netlist_t::port_order_t::fet_n_source] - 1][*hidden - 1] = -1;
			g[*hidden - 1][netlist->port[pos].item[netlist_t::port_order_t::fet_n_source] - 1] = -1;
		}
		if (netlist->port[pos].item[netlist_t::port_order_t::fet_n_drain] || netlist->port[pos].item[netlist_t::port_order_t::fet_n_source]) {
			x[*hidden - 1] = component_t::copy_s;
			y[*hidden - 1] = component_t::fet_current;
			f_netlist_pos[*hidden - 1] = pos;
			(*hidden)--;
		}
		break;
	default:
		break;
	}
}
void component_t::insert_from_memory(double(**f)(netlist_t *, uint8_t, uint8_t, double *, double *), uint8_t pos) {
	f[pos] = component_t::const_val;
}
void component_t::insert_copy_s(double(**f)(netlist_t *, uint8_t,  uint8_t, double *, double *), uint8_t pos) {
	f[pos] = component_t::copy_s;
}
double component_t::copy_s(netlist_t *netlist, uint8_t netlist_pos, uint8_t call_pos, double *input, double *memory) {
	return input[call_pos];
}
double component_t::diode_vd(netlist_t *netlist, uint8_t netlist_pos, uint8_t call_pos, double *input, double *memory) {
	return input[netlist->port[netlist_pos].item[netlist_t::port_order_t::node_main] - 1] - input[netlist->port[netlist_pos].item[netlist_t::port_order_t::node_dest] - 1];
}
double component_t::diode_r_current(netlist_t *netlist, uint8_t netlist_pos, uint8_t call_pos, double *input, double *memory) {
	return (netlist->value[netlist_pos].item[netlist_t::data_order_t::is])*(exp((input[netlist->port[netlist_pos].item[netlist_t::port_order_t::node_main] - 1] - input[netlist->port[netlist_pos].item[netlist_t::port_order_t::node_dest] - 1]) / netlist->value[netlist_pos].item[netlist_t::data_order_t::vt]) - 1);
}
double component_t::v(netlist_t *netlist, uint8_t netlist_pos, uint8_t call_pos, double *input, double *memory) {
	return netlist->value[netlist_pos].item[netlist_t::data_order_t::voltage];
}
double component_t::const_val(netlist_t *netlist, uint8_t netlist_pos, uint8_t call_pos, double *input, double *memory) {
	return memory[call_pos];
}
//fet
double component_t::fet_current(netlist_t *netlist, uint8_t netlist_pos, uint8_t call_pos, double *variables, double *memory) {
	double vds = (netlist_t::port_order_t::fet_n_drain ? variables[netlist->port[netlist_pos].item[netlist_t::port_order_t::fet_n_drain] - 1] : 0) 
		- (netlist_t::port_order_t::fet_n_source ? variables[netlist->port[netlist_pos].item[netlist_t::port_order_t::fet_n_source] - 1] : 0);
	double vgs = (netlist_t::port_order_t::fet_n_gate ? variables[netlist->port[netlist_pos].item[netlist_t::port_order_t::fet_n_gate] - 1] : 0) 
		- (netlist_t::port_order_t::fet_n_source ? variables[netlist->port[netlist_pos].item[netlist_t::port_order_t::fet_n_source] - 1] : 0);
	double w = netlist->value[netlist_pos].item[netlist_t::data_order_t::fet_w];
	double l = netlist->value[netlist_pos].item[netlist_t::data_order_t::fet_l];
	double kp = netlist->value[netlist_pos].item[netlist_t::data_order_t::fet_kp];
	double vth = netlist->value[netlist_pos].item[netlist_t::data_order_t::fet_vth];
	double lambda = netlist->value[netlist_pos].item[netlist_t::data_order_t::fet_lambda];
	if (vth > vgs) {
		return 0;
	}
	else if (vds > (vgs - vth)) {
		return (w / l)*(kp / 2)*pow(vgs - vth, 2)*(1 + lambda * vds);
	}
	else {
		return (w / l)*(kp / 2)*(2*(vgs-vth)*vds - vds*vds)*(1 + lambda * vds);
	}
}



mna_t::mna_t(netlist_t *src) {
	uint8_t pos = src->components;
	source = src;
	nodes = src->nodes;
	size = src->nodes;
	components = src->components;
	if (size) {
		while (pos) {
			pos--;
			switch (src->component_list[pos]){
			case (netlist_t::component_order_t::voltage_src):
				size += component_t::v_mna_add;
				break;
			case (netlist_t::component_order_t::resistor):
				size += component_t::resistor_mna_add;
				break;
			case (netlist_t::component_order_t::capacitor):
				size += component_t::capacitor_mna_add;
				break;
			case (netlist_t::component_order_t::inductor):
				size += component_t::inductor_mna_add;
				break;
			case (netlist_t::component_order_t::diode_r):
				size += component_t::diode_r_mna_add;
				break;
			default:
				break;
			}
		}
		hidden = size;
		z = new double[size];
		s = new double[size];
		s_old = new double[size];
		s_diff = new double[size];
		typedef double(*comp_function)(netlist_t *, uint8_t, uint8_t, double *, double *);
		x = new comp_function[size];
		y = new comp_function[size];
		f_memory = new double[size];
		f_netlist_pos = new uint8_t[size];
		g = new double*[size];
		pos = size;
		while (pos) {
			pos--;
			g[pos] = new double[size];
		}
		r = new double*[size];
		pos = size;
		while (pos) {
			pos--;
			r[pos] = new double[size];
		}
		jn = new double*[size];
		pos = size;
		while (pos) {
			pos--;
			jn[pos] = new double[size];
		}
		jni = new double*[size];
		pos = size;
		while (pos) {
			pos--;
			jni[pos] = new double[size];
		}
		jni_old = new double*[size];
		pos = size;
		while (pos) {
			pos--;
			jni_old[pos] = new double[size];
		}
		utils_t::zero_vect(x, size);
		utils_t::zero_vect(y, size);
		utils_t::zero_vect(s, size);
		utils_t::zero_vect(f_memory, size);
		utils_t::zero_matrix(g, size, size);
		utils_t::zero_matrix(r, size, size);
		utils_t::zero_matrix(jni, size, size);
	}
}
mna_t::~mna_t() {
	uint8_t pos = size;
	if (size) {
		delete[] z;
		delete[] s;
		delete[] s_old;
		delete[] s_diff;
		delete[] x;
		delete[] y;
		delete[] f_memory;
		delete[] f_netlist_pos;
		while (pos) {
			pos--;
			delete[] g[pos];
		}
		delete[] g;
		pos = size;
		while (pos) {
			pos--;
			delete[] r[pos];
		}
		delete[] r;
		pos = size;
		while (pos) {
			pos--;
			delete[] jn[pos];
		}
		delete[] jn;
		pos = size;
		while (pos) {
			pos--;
			delete[] jni[pos];
		}
		delete[] jni;
		pos = size;
		while (pos) {
			pos--;
			delete[] jni_old[pos];
		}
		delete[] jni_old;
	}
}
void mna_t::generate() {
	uint8_t i = size;
	uint8_t pos = components;
	while (pos){
		pos--;
		component_t::insert(g, x, y, f_netlist_pos, source, pos, &hidden, f_memory);
	}
	pos = nodes;
	while (pos){
		pos--;
		f_netlist_pos[pos] = pos;
		component_t::insert_from_memory(x, pos);
		component_t::insert_copy_s(y,pos);
	}
	math.invert(g, r, size);
}
void mna_t::update_z() {
	uint8_t i = size;
	uint8_t j;
	while (i) {
		i--;
		z[i] = 0;	//limpando para usar +=
		j = size;
		while (j) {
			j--;
			z[j] += r[i][j] * x[j](source, f_netlist_pos[i], i, s, f_memory);
		}
	}
}
double mna_t::eval_row(uint8_t i) {
	double acc;
	uint8_t aux = size;
	acc = y[i](source, f_netlist_pos[i], i, s, f_memory);
	while (aux) {
		aux--;
		acc -= r[i][aux] * x[aux](source, f_netlist_pos[aux], aux, s, f_memory);
	}
	return acc;
}
double mna_t::eval_row(uint8_t i, double *input) {
	double acc;
	uint8_t aux = size;
	acc = y[i](source, f_netlist_pos[i], i, input, f_memory);
	while (aux) {
		aux--;
		acc -= r[i][aux] * x[aux](source, f_netlist_pos[aux], aux, input, f_memory);
	}
	return acc;
}
double* mna_t::eval(double* output) {
	uint8_t i_max = size;
	while (i_max){
		i_max--;
		output[i_max] = eval_row(i_max);
	}
	return output;
}
double** mna_t::jacobian() {
	uint8_t i = size;
	uint8_t j;
	double *s_temp;
	s_temp = new double[size];
	utils_t::copy_vect(s, s_temp, size);
	while (i) {
		i--;
		j = size;
		while (j) {
			j--;
			s_temp[j] += MNA_JACOBIAN_STEP;
			jn[i][j] = (eval_row(i, s_temp) - eval_row(i, s)) / MNA_JACOBIAN_STEP;
			s_temp[j] -= MNA_JACOBIAN_STEP; 
		}
	}
	math.invert(jn, jni, size);
	delete[] s_temp;
	return jn;
}
void mna_t::iterate_once() {
	double *s_new = new double[size];
	double *z_new = new double[size];
	eval(z);
	jacobian();
	utils_t::mult_mat_vect(jni, z, z_new, (uint16_t)size, (uint16_t)size);
	math.sub(s, z_new, s_new, size);
	utils_t::copy_vect(s_new, s, size);
	delete[] s_new;
	delete[] z_new;
}
void mna_t::iterate_n(uint16_t n) {
	while (n){
		n--;
		iterate_once();
	}
}



void math_t::sub(double *a, double *b, double *out, uint8_t length) {
	while (length) {
		length--;
		out[length] = a[length] - b[length];
	}
}
void math_t::add(double *a, double *b, double *out, uint8_t length) {
	while (length) {
		length--;
		out[length] = a[length] + b[length];
	}
}
int math_t::LUPDecompose(double **A, int N, double Tol, int *P) {

	int i, j, k, imax;
	double maxA, *ptr, absA;

	for (i = 0; i <= N; i++)
		P[i] = i; //Unit permutation matrix, P[N] initialized with N

	for (i = 0; i < N; i++) {
		maxA = 0.0;
		imax = i;

		for (k = i; k < N; k++)
			if ((absA = fabs(A[k][i])) > maxA) {
				maxA = absA;
				imax = k;
			}

		if (maxA < Tol) return 0; //failure, matrix is degenerate

		if (imax != i) {
			//pivoting P
			j = P[i];
			P[i] = P[imax];
			P[imax] = j;

			//pivoting rows of A
			ptr = A[i];
			A[i] = A[imax];
			A[imax] = ptr;

			//counting pivots starting from N (for determinant)
			P[N]++;
		}

		for (j = i + 1; j < N; j++) {
			A[j][i] /= A[i][i];

			for (k = i + 1; k < N; k++)
				A[j][k] -= A[j][i] * A[i][k];
		}
	}

	return 1;  //decomposition done 
}
void math_t::LUPSolve(double **A, int *P, double *b, int N, double *x) {

	for (int i = 0; i < N; i++) {
		x[i] = b[P[i]];

		for (int k = 0; k < i; k++)
			x[i] -= A[i][k] * x[k];
	}

	for (int i = N - 1; i >= 0; i--) {
		for (int k = i + 1; k < N; k++)
			x[i] -= A[i][k] * x[k];

		x[i] = x[i] / A[i][i];
	}
}
void math_t::LUPInvert(double **A, int *P, int N, double **IA) {

	for (int j = 0; j < N; j++) {
		for (int i = 0; i < N; i++) {
			if (P[i] == j)
				IA[i][j] = 1.0;
			else
				IA[i][j] = 0.0;

			for (int k = 0; k < i; k++)
				IA[i][j] -= A[i][k] * IA[k][j];
		}

		for (int i = N - 1; i >= 0; i--) {
			for (int k = i + 1; k < N; k++)
				IA[i][j] -= A[i][k] * IA[k][j];

			IA[i][j] = IA[i][j] / A[i][i];
		}
	}
}
double math_t::LUPDeterminant(double **A, int *P, int N) {

	double det = A[0][0];

	for (int i = 1; i < N; i++)
		det *= A[i][i];

	if ((P[N] - N) % 2 == 0)
		return det;
	else
		return -det;
}
double** math_t::invert(double **in, double** out, uint8_t size) {
	int *P = new int[size];
	double ** A;
	A = new double*[size];
	uint8_t pos = size;
	while (pos) {
		pos--;
		A[pos] = new double[size];
	}
	utils_t::copy_matrix(in, A, size, size);
	math.LUPDecompose(A, size, 0.001, P);
	math.LUPInvert(A, P, size, out);
	pos = size;
	while (pos) {
		pos--;
		delete[] A[pos];
	}
	delete[] A;
	return out;
}

using namespace std;
int main() {
	cout << "inicio" << endl;
	netlist_t net(3, 5);
	net.component_list[0] = netlist_t::component_order_t::resistor;
	net.component_list[1] = netlist_t::component_order_t::resistor;
	net.component_list[2] = netlist_t::component_order_t::voltage_src;
	net.component_list[3] = netlist_t::component_order_t::diode_r;
	net.component_list[4] = netlist_t::component_order_t::diode_r;

	cout << "componentes adicionados" << endl;

	net.port[0].item[netlist_t::port_order_t::node_main] = 1;
	net.port[0].item[netlist_t::port_order_t::node_dest] = 2;

	net.port[1].item[netlist_t::port_order_t::node_main] = 3;
	net.port[1].item[netlist_t::port_order_t::node_dest] = 0;

	net.port[2].item[netlist_t::port_order_t::node_main] = 1;
	net.port[2].item[netlist_t::port_order_t::node_dest] = 0;

	net.port[3].item[netlist_t::port_order_t::node_main] = 2;
	net.port[3].item[netlist_t::port_order_t::node_dest] = 3;

	net.port[4].item[netlist_t::port_order_t::node_main] = 2;
	net.port[4].item[netlist_t::port_order_t::node_dest] = 3;

	cout << "portas definidas" << endl;

	net.value[0].item[netlist_t::data_order_t::resistance] = 100;
	net.value[1].item[netlist_t::data_order_t::resistance] = 100;
	net.value[2].item[netlist_t::data_order_t::voltage] = 5;
	net.value[3].item[netlist_t::data_order_t::is] = 1e-12;
	net.value[3].item[netlist_t::data_order_t::vt] = 0.025;
	net.value[4].item[netlist_t::data_order_t::is] = 1e-12;
	net.value[4].item[netlist_t::data_order_t::vt] = 0.025;

	cout << "valores definidos" << endl;

	mna_t solver(&net);
	solver.generate();
	solver.iterate_n(2000);
	cout << "Va =     " << solver.s[0] << endl << "Vb =     " << solver.s[1] << endl << "Vc =     " << solver.s[2] << endl << "Is =     " << solver.s[3] << endl << "Vdiodo = " << solver.s[4] << endl;
	return 0;
}
