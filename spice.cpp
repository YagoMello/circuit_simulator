#include <stdio.h>
#include <inttypes.h>
#include <math.h>

constexpr uint8_t NETLIST_PORT_MAX = 5;
constexpr uint8_t NETLIST_PARAM_MAX = 5;

class utils_t {
public:
	template <typename T> void copy_vect(T *, T*, uint16_t);
	template <typename T> void copy_matrix(T**, T**, uint16_t , uint16_t);
	template <typename T> void zero_vect(T *, uint16_t);
	template <typename T> void zero_matrix(T**, uint16_t, uint16_t);
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
	void insert(double **, double *, double(**y)(netlist_t *, uint8_t *, double *), netlist_t *, uint8_t *, netlist_t *, uint8_t , uint8_t *);
	class voltage_t {
		double voltage(double **arg) {
			
		}
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
		double current(netlist_t *, uint8_t *, double *);
	}diode;
}component;
void component_t::insert(double **g, double *x, double (**y)(netlist_t *, uint8_t *, double *), netlist_t *y_src, uint8_t *y_pos, netlist_t *source, uint8_t pos, uint8_t *hidden) {
	switch (source->component_list[pos]){
	case(netlist_t::component_order_t::voltage):
		if (netlist_t::port_order_t::node_main) {
			g[source->port[pos].item[netlist_t::port_order_t::node_main - 1]][*hidden] = 1;
			g[*hidden][source->port[pos].item[netlist_t::port_order_t::node_main - 1]] = 1;
		}
		if (netlist_t::port_order_t::node_dest) {
			g[source->port[pos].item[netlist_t::port_order_t::node_dest - 1]][*hidden] = -1;
			g[*hidden][source->port[pos].item[netlist_t::port_order_t::node_dest - 1]] = -1;
		}
		if (netlist_t::port_order_t::node_main | netlist_t::port_order_t::node_dest) {
			y[*hidden] = voltage.voltage;
		}
		hidden++;
		break;
	case(netlist_t::component_order_t::current):

		break;
	case(netlist_t::component_order_t::resistor):
		if (netlist_t::port_order_t::node_main) {
			g[source->port[pos].item[netlist_t::port_order_t::node_main - 1]][source->port[pos].item[netlist_t::port_order_t::node_main - 1]] = 1 / source->value[pos].item[netlist_t::data_order_t::resistance];
			if (netlist_t::port_order_t::node_dest) {
				g[source->port[pos].item[netlist_t::port_order_t::node_main - 1]][source->port[pos].item[netlist_t::port_order_t::node_dest - 1]] = -1 / source->value[pos].item[netlist_t::data_order_t::resistance];
			}
		}
		if (netlist_t::port_order_t::node_dest) {
			g[source->port[pos].item[netlist_t::port_order_t::node_dest - 1]][source->port[pos].item[netlist_t::port_order_t::node_dest - 1]] = 1 / source->value[pos].item[netlist_t::data_order_t::resistance];
			if (netlist_t::port_order_t::node_main) {
				g[source->port[pos].item[netlist_t::port_order_t::node_dest - 1]][source->port[pos].item[netlist_t::port_order_t::node_main - 1]] = -1 / source->value[pos].item[netlist_t::data_order_t::resistance];
			}
		}
		break;
	case(netlist_t::component_order_t::capacitor):

		break;
	case(netlist_t::component_order_t::inductor):

		break;
	case(netlist_t::component_order_t::diode_r):
		if (netlist_t::port_order_t::node_main) {
			g[source->port[pos].item[netlist_t::port_order_t::node_main - 1]][*hidden] = 1;
			g[*hidden][source->port[pos].item[netlist_t::port_order_t::node_main - 1]] = 1;
		}
		if (netlist_t::port_order_t::node_dest) {
			g[source->port[pos].item[netlist_t::port_order_t::node_dest - 1]][*hidden] = -1;
			g[*hidden][source->port[pos].item[netlist_t::port_order_t::node_dest - 1]] = -1;
		}
		if (netlist_t::port_order_t::node_main | netlist_t::port_order_t::node_dest) {
			y[*hidden] = diode.current;
			
		}
		hidden++;
		break;
	default:
		break;
	}
}
double component_t::diode_t::current(netlist_t *source, uint8_t *pos, double *state) {
	return (source->value[*pos].item[netlist_t::data_order_t::is])*(exp((state[source->port[*pos].item[netlist_t::port_order_t::node_main] - 1] - state[source->port[*pos].item[netlist_t::port_order_t::node_dest] - 1]) / source->value[*pos].item[netlist_t::data_order_t::vt]) - 1);
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
	double *x;
	double *s;
	double (**y)(netlist_t *, uint8_t *, double *);
	netlist_t *y_source;
	uint8_t *y_pos;
	mna_t(netlist_t *, component_t *);
	~mna_t(void);
	void generate(void);//
	void update(void);
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
		x = new double[size];
		s = new double[size];
		y = new (double(**)(netlist_t *, uint8_t *, double *))[size];
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
		utils.zero_vect(x, size);
		utils.zero_vect(y, size);
		utils.zero_vect(s, size);
		utils.zero_matrix(g, size, size);
		utils.zero_matrix(r, size, size);
	}
}
mna_t::~mna_t() {
	uint8_t pos = size;
	if (size) {
		delete[] x;
		delete[] s;
		delete[] y;
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
	}
}
void mna_t::generate() {
	uint8_t i = size;
	uint8_t pos = components;
	uint8_t j;
	double y_temp;

	while (pos){
		pos--;
		component.insert(g, x, y, y_source, y_pos,source, pos, &hidden);
	}

	math.invert(g, r, size);
	while (i) {
		i--;
		s[i] = 0;	//limpando para usar +=
		x[i] = 0; //cond inicial
		y_temp = y[i](source, &i, x);
		j = size;
		while (j) {
			j--;
			s[j] += x[j] - r[i][j] * y_temp;//tuderrado
		}
	}
}
void mna_t::update() {
	uint8_t i = size;
	uint8_t j;
	double y_temp;
	while (i) {
		i--;
		j = size;
		s[i] = 0;	//limpando para usar +=
		y_temp = y[i](source, &i, x);
		while (j) {
			j--;
			s[j] += x[j] - r[i][j] * y_temp;
		}
	}
}

class math_t {
public:
	double diff(double(*func)(netlist_t::data_t*));
	void invert(double **, double **, uint8_t);
private:
	int LUPDecompose(double **A, int N, double Tol, int *P);
	void LUPSolve(double **A, int *P, double *b, int N, double *x);
	void LUPInvert(double **A, int *P, int N, double **IA);
	double LUPDeterminant(double **A, int *P, int N);
}math;
double math_t::diff(double (*func)(netlist_t::data_t*)) {
	
	return 0.0;
}
/* INPUT: A - array of pointers to rows of a square matrix having dimension N
 *        Tol - small tolerance number to detect failure when the matrix is near degenerate
 * OUTPUT: Matrix A is changed, it contains both matrices L-E and U as A=(L-E)+U such that P*A=L*U.
 *        The permutation matrix is not stored as a matrix, but in an integer vector P of size N+1
 *        containing column indexes where the permutation matrix has "1". The last element P[N]=S+N,
 *        where S is the number of row exchanges needed for determinant computation, det(P)=(-1)^S
 */
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
/* INPUT: A,P filled in LUPDecompose; b - rhs vector; N - dimension
 * OUTPUT: x - solution vector of A*x=b
 */
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
/* INPUT: A,P filled in LUPDecompose; N - dimension
 * OUTPUT: IA is the inverse of the initial matrix
 */
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
/* INPUT: A,P filled in LUPDecompose; N - dimension.
 * OUTPUT: Function returns the determinant of the initial matrix
 */
double math_t::LUPDeterminant(double **A, int *P, int N) {

	double det = A[0][0];

	for (int i = 1; i < N; i++)
		det *= A[i][i];

	if ((P[N] - N) % 2 == 0)
		return det;
	else
		return -det;
}
void math_t::invert(double **in, double** out, uint8_t size) {
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
}
