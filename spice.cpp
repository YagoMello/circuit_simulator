#include <stdio.h>
#include <inttypes.h>
#include <math.h>

constexpr uint8_t NETLIST_PORT_MAX = 5;
constexpr uint8_t NETLIST_PARAM_MAX = 5;
constexpr double MNA_JACOBIAN_STEP = 1e-6;

class utils_t {
public:
	template <typename T> void copy_vect(T *, T*, uint16_t);
	template <typename T> void copy_matrix(T**, T**, uint16_t , uint16_t);
	template <typename T> void zero_vect(T *, uint16_t);
	template <typename T> void zero_matrix(T**, uint16_t, uint16_t);
	template <typename T> T mult_matrix_single_row_col(T**, T**, uint16_t, uint16_t, uint16_t);
	template <typename T> void mult_matrix(T**, T**, T**, uint16_t, uint16_t, uint16_t);
}utils;
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
}
template <typename T> void utils_t::mult_matrix(T** M, T** N, T** OUT, uint16_t i1, uint16_t j2, uint16_t len) {
	uint16_t x;
	while (i){
		i1--;
		x = j2;
		while (x){
			x--;
			OUT[i1][x] = mult_matrix_single_row_col(M, N, i1, x, len);
		}
	}
}


class netlist_t{
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
enum netlist_t::port_order_t {
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
	ve
};
enum netlist_t::data_order_t {
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
	vs,
	vt,
	//jfet

	//bjt

	//mosfet
	m_w = 0,
	m_l,
	m_lambda,
	m_kp
};
enum netlist_t::component_order_t {
	voltage,
	current,
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
	mosfet_n,
	mosfet_p
};
struct netlist_t::data_t {
	double item[NETLIST_PARAM_MAX];
};
struct netlist_t::port_t {
	uint8_t item[NETLIST_PORT_MAX];
};

class component_t {
public:
	void insert(double **, double(**x)(netlist_t *, uint8_t, double *, double *), double(**y)(netlist_t *, uint8_t, double *, double *), uint8_t *, netlist_t *, uint8_t, uint8_t *);
	class generic_op_t {
		double copy_z(netlist_t *, uint8_t, double *);
		double const_val(netlist_t *, uint8_t, double *);
	}generic_op;
	class voltage_t {
		double voltage(netlist_t *, uint8_t , double *);
	}voltage;
	class current_t {

	}current;
	class resistor_t {
	public:
		const uint8_t mna_add = 0;
	}resistor;
	class capacitor_t {
	public:
		const uint8_t mna_add = 0;
	}capacitor;
	class inductor_t {
	public:
		const uint8_t mna_add = 0;
	}inductor;
	class diode_t {
	public:
		const uint8_t mna_add = 1;
		double current(netlist_t *, uint8_t , double *);
	}diode;
}component;
void component_t::insert(double **g, double(**x)(netlist_t *, uint8_t, double *, double *), double (**y)(netlist_t *, uint8_t, double *, double *), uint8_t *f_netlist_pos, netlist_t *source, uint8_t pos, uint8_t *hidden) {
	switch (source->component_list[pos]){
	case(netlist_t::component_order_t::voltage):
		if (netlist_t::port_order_t::node_main) {
			g[source->port[pos].item[netlist_t::port_order_t::node_main] - 1][*hidden] = 1;
			g[*hidden][source->port[pos].item[netlist_t::port_order_t::node_main] - 1] = 1;
		}
		if (netlist_t::port_order_t::node_dest) {
			g[source->port[pos].item[netlist_t::port_order_t::node_dest] - 1][*hidden] = -1;
			g[*hidden][source->port[pos].item[netlist_t::port_order_t::node_dest] - 1] = -1;
		}
		if (netlist_t::port_order_t::node_main | netlist_t::port_order_t::node_dest) {
			y[*hidden] = voltage.voltage;
			f_netlist_pos[*hidden] = pos;
			hidden++;
		}
		break;
	case(netlist_t::component_order_t::current):

		break;
	case(netlist_t::component_order_t::resistor):
		if (netlist_t::port_order_t::node_main) {
			g[source->port[pos].item[netlist_t::port_order_t::node_main] - 1][source->port[pos].item[netlist_t::port_order_t::node_main] - 1] = 1 / source->value[pos].item[netlist_t::data_order_t::resistance];
			if (netlist_t::port_order_t::node_dest) {
				g[source->port[pos].item[netlist_t::port_order_t::node_main] - 1][source->port[pos].item[netlist_t::port_order_t::node_dest] - 1] = -1 / source->value[pos].item[netlist_t::data_order_t::resistance];
			}
			y[source->port[pos].item[netlist_t::port_order_t::node_main] - 1] = generic_op.copy_z;
			f_netlist_pos[source->port[pos].item[netlist_t::port_order_t::node_main] - 1] = source->port[pos].item[netlist_t::port_order_t::node_main] - 1;
		}
		if (netlist_t::port_order_t::node_dest) {
			g[source->port[pos].item[netlist_t::port_order_t::node_dest] - 1][source->port[pos].item[netlist_t::port_order_t::node_dest] - 1] = 1 / source->value[pos].item[netlist_t::data_order_t::resistance];
			if (netlist_t::port_order_t::node_main) {
				g[source->port[pos].item[netlist_t::port_order_t::node_dest] - 1][source->port[pos].item[netlist_t::port_order_t::node_main] - 1] = -1 / source->value[pos].item[netlist_t::data_order_t::resistance];
			}
			y[source->port[pos].item[netlist_t::port_order_t::node_dest] - 1] = generic_op.copy_z;
			f_netlist_pos[source->port[pos].item[netlist_t::port_order_t::node_dest] - 1] = source->port[pos].item[netlist_t::port_order_t::node_dest] - 1;
		}
		break;
	case(netlist_t::component_order_t::capacitor):

		break;
	case(netlist_t::component_order_t::inductor):

		break;
	case(netlist_t::component_order_t::diode_r):
		if (netlist_t::port_order_t::node_main) {
			g[source->port[pos].item[netlist_t::port_order_t::node_main] - 1][*hidden] = 1;
			g[*hidden][source->port[pos].item[netlist_t::port_order_t::node_main] - 1] = 1;
		}
		if (netlist_t::port_order_t::node_dest) {
			g[source->port[pos].item[netlist_t::port_order_t::node_dest] - 1][*hidden] = -1;
			g[*hidden][source->port[pos].item[netlist_t::port_order_t::node_dest] - 1] = -1;
		}
		if (netlist_t::port_order_t::node_main || netlist_t::port_order_t::node_dest) {
			y[*hidden] = diode.current;
			f_netlist_pos[*hidden] = pos;
			hidden++;
		}
		break;
	default:
		break;
	}
}
double component_t::generic_op_t::copy_z(netlist_t *source, uint8_t pos, double *input) {
	return input[pos];
}
double component_t::diode_t::current(netlist_t *source, uint8_t pos, double *input) {
	return (source->value[pos].item[netlist_t::data_order_t::is])*(exp((input[source->port[pos].item[netlist_t::port_order_t::node_main] - 1] - input[source->port[pos].item[netlist_t::port_order_t::node_dest] - 1]) / source->value[pos].item[netlist_t::data_order_t::vt]) - 1);
}
double component_t::voltage_t::voltage(netlist_t *source, uint8_t pos, double *input) {
	return source->value[pos].item[netlist_t::data_order_t::voltage];
}
double component_t::generic_op_t::const_val(netlist_t *source, uint8_t pos, double *input) {
	return mna.f_memory;
}

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
	double (**y)(netlist_t *ref, uint8_t netlist_pos, double *input, double *memory);
	double (**x)(netlist_t *ref, uint8_t netlist_pos, double *input, double *memory);
	netlist_t *y_source;
	double *f_memory;
	uint8_t *f_netlist_pos;
	mna_t(netlist_t *, component_t *);
	~mna_t(void);
	void generate(void);//
	void update(void);
	void update_z(void);
	double eval_row(uint8_t);
	double eval_row(uint8_t, double *);
	void iterate(void);
	double** jacobian(void);
};
mna_t::mna_t(netlist_t *src, component_t *db) {
	uint8_t pos = src->components;
	source = src;
	hidden = src->nodes;
	nodes = src->nodes;
	size = src->nodes;
	components = src->components;
	if (size) {
		while (pos) {
			pos--;
			switch (src->component_list[pos]){
			case (netlist_t::component_order_t::resistor):
				size += db->resistor.mna_add;
				break;
			case (netlist_t::component_order_t::capacitor):
				size += db->capacitor.mna_add;
				break;
			case (netlist_t::component_order_t::inductor):
				size += db->inductor.mna_add;
				break;
			case (netlist_t::component_order_t::diode_r):
				size += db->diode.mna_add;
				break;
			default:
				break;
			}
		}
		z = new double[size];
		s = new double[size];
		s_old = new double[size];
		s_diff = new double[size];
		x = new (double(**)(netlist_t *, uint8_t, double *, double *))[size];
		y = new (double(**)(netlist_t *, uint8_t, double *, double *))[size];
		f_memory = new double[size];
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
		utils.zero_vect(x, size);
		utils.zero_vect(y, size);
		utils.zero_vect(s, size);
		utils.zero_vect(f_memory, size);
		utils.zero_matrix(g, size, size);
		utils.zero_matrix(r, size, size);
		utils.zero_matrix(jni, size, size);
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
	uint8_t j;
	double x_temp;
	double y_temp;
	while (pos){
		pos--;
		component.insert(g, x, y, f_netlist_pos, source, pos, &hidden);
	}
	math.invert(g, r, size);
	while (i) {
		i--;
		s[i] = 0;	//limpando para usar +=
		y_temp = y[i](source, f_netlist_pos[i], z, f_memory);
		j = size;
		while (j) {
			j--;
			s_old[j] += y_temp - r[i][j] * x[j](source, f_netlist_pos[i], z, f_memory);
		}
	}
	i = size;
	while (i) {
		i--;
		s[i] = 0;	//limpando para usar +=
		y_temp = y[i](source, f_netlist_pos[i], z, f_memory);
		j = size;
		while (j) {
			j--;
			s[j] += y_temp - r[i][j] * x[j](source, f_netlist_pos[i], z, f_memory);
		}
	}
}
void mna_t::update() {
	uint8_t i = size;
	uint8_t j;
	double x_temp;
	double y_temp;
	utils.copy_vect(s, s_old, size);
	while (i) {
		i--;
		s[i] = 0;	//limpando para usar +=
		y_temp = y[i](source, f_netlist_pos[i], s, f_memory);
		j = size;
		while (j) {
			j--;
			s[j] += y_temp - r[i][j] * x[j](source, f_netlist_pos[i], s, f_memory);
		}
	}
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
			z[j] += r[i][j] * x[j](source, f_netlist_pos[i], s, f_memory);
		}
	}
}
double mna_t::eval_row(uint8_t i) {
	double acc;
	uint8_t aux = size;
	acc = y[i](source, f_netlist_pos[i], s, f_memory);
	while (aux) {
		aux--;
		acc -= r[i][aux] * y[aux](source, f_netlist_pos[aux], s, f_memory);
	}
}
double mna_t::eval_row(uint8_t i, double *input) {
	double acc;
	uint8_t aux = size;
	acc = y[i](source, f_netlist_pos[i], input, f_memory);
	while (aux) {
		aux--;
		acc -= r[i][aux] * y[aux](source, f_netlist_pos[aux], input, f_memory);
	}
}
double** mna_t::jacobian() {
	uint8_t i = size;
	uint8_t j;
	utils.copy_vect(jni, jni_old, size);
	double *s_temp;
	s_temp = new double[size];
	while (i) {
		i--;
		j = size;
		while (j) {
			j--;////editar aqi
			jn[i][j] = (y[i](source, f_netlist_pos[i], s_temp, f_memory) - y[i](source, f_netlist_pos[i], s, f_memory)) / MNA_JACOBIAN_STEP;
		}
	}
	math.invert(jn, jni, size);
	delete[] s_temp;
	return jn;
}
void mna_t::iterate() {

	z
}

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
	utils.copy_matrix(in, A, size, size);
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
