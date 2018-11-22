
#include <nspireio/console.hpp>
#include <inttypes.h>
#include <math.h>
#include <cstdlib>
#include <string.h>


#ifndef TRUE
#define TRUE 1
#endif // TRUE

#ifndef FALSE
#define FALSE 0
#endif // FALSE

nio::console csl;

constexpr uint8_t NETLIST_CONNECTION_MAX = 4;
constexpr uint8_t NETLIST_PARAM_MAX = 5;
constexpr uint8_t NETLIST_ALIAS_SIZE = 16;
constexpr uint8_t COMPONENT_CODE_MAX_LENGTH = 16;
constexpr double NOD_DEF_ERROR_VECTOR_MAX_LENGTH = 1e-10;
constexpr double NOD_DEF_JACOBIAN_STEP = 1e-10;
constexpr double NOD_DEF_JNI_TOL = 1e-10;
constexpr double NOD_R_TOL = 1e-14;
double GMIN;
//NODE_REF
//2 port
constexpr uint8_t positive = 0;
constexpr uint8_t negative = 1;
//voltage_source trick
constexpr uint8_t hidden = 2;
//current_source
constexpr uint8_t from = 0;
constexpr uint8_t to = 1;
//mosfet
constexpr uint8_t fet_source = 0;
constexpr uint8_t fet_drain = 1;
constexpr uint8_t fet_gate = 2;

//VALUE_REF
//voltage_source
constexpr uint8_t voltage = 0;
//current_source
constexpr uint8_t current = 0;
//resistor
constexpr uint8_t resistance = 0;
//diode
constexpr uint8_t diode_is = 0;
constexpr uint8_t diode_vt = 1;
//mosfet
constexpr uint8_t fet_w = 0;
constexpr uint8_t fet_l = 1;
constexpr uint8_t fet_kp = 2;
constexpr uint8_t fet_vt = 3;
constexpr uint8_t fet_lambda = 4;

//TYPE_REF
constexpr uint16_t resistor = 0;
constexpr uint16_t voltage_source = 1;
constexpr uint16_t current_source = 2;
constexpr uint16_t capacitor = 3;
constexpr uint16_t inductor = 4;
constexpr uint16_t diode = 5;
constexpr uint16_t bjt_npn = 6;
constexpr uint16_t bjt_pnp = 7;
constexpr uint16_t fet_n = 8;
constexpr uint16_t fet_p = 9;

//NAME_REF
const char* resistor_name = "Resistor";
const char* voltage_source_name = "Voltage Source";
const char* current_source_name = "Current Source";
const char* capacitor_name = "Capacitor";
const char* inductor_name = "Inductor";
const char* diode_name = "Diode";
const char* bjt_npn_name = "BJT NPN";
const char* bjt_pnp_name = "BJT PNP";
const char* fet_n_name = "N channel MOSFET";
const char* fet_p_name = "P channel MOSFET";


//CODE_REF
const char* resistor_code = "r";
const char* voltage_source_code = "vs";
const char* current_source_code = "cs";
const char* capacitor_code = "cap";
const char* inductor_code = "ind";
const char* diode_code = "d";
const char* bjt_npn_code = "bjt npn";
const char* bjt_pnp_code = "bjt pnp";
const char* fet_n_code = "nmos";
const char* fet_p_code = "pmos";

class utils_nspire_t{
public:
    static char* double_to_ascii(double, char*);
    static char* safe_gets(char*);
    static void print_double(double);
    static void print_uint8(uint8_t);
};
char* utils_nspire_t::double_to_ascii(double number, char *str) {
    //copied from androider @ stackoverflow
    // handle special cases
    const double PRECISION = 0.00000000000001;
    const int MAX_NUMBER_STRING_SIZE = 32;
    double n = number;
    char *s = str;
    if (isnan(n)) {
        strcpy(s, "nan");
    } else if (isinf(n)) {
        strcpy(s, "inf");
    } else if (n == 0.0) {
        strcpy(s, "0");
    } else {
        int digit, m, m1;
        char *c = s;
        int neg = (n < 0);
        if (neg)
            n = -n;
        // calculate magnitude
        m = log10(n);
        int useExp = (m >= 14 || (neg && m >= 9) || m <= -9);
        if (neg)
            *(c++) = '-';
        // set up for scientific notation
        if (useExp) {
            if (m < 0)
               m -= 1.0;
            n = n / pow(10.0, m);
            m1 = m;
            m = 0;
        }
        if (m < 1.0) {
            m = 0;
        }
        // convert the number
        while (n > PRECISION || m >= 0) {
            double weight = pow(10.0, m);
            if (weight > 0 && !isinf(weight)) {
                digit = floor(n / weight);
                n -= (digit * weight);
                *(c++) = '0' + digit;
            }
            if (m == 0 && n > 0)
                *(c++) = '.';
            m--;
        }
        if (useExp) {
            // convert the exponent
            int i, j;
            *(c++) = 'e';
            if (m1 > 0) {
                *(c++) = '+';
            } else {
                *(c++) = '-';
                m1 = -m1;
            }
            m = 0;
            while (m1 > 0) {
                *(c++) = '0' + m1 % 10;
                m1 /= 10;
                m++;
            }
            c -= m;
            for (i = 0, j = m-1; i<j; i++, j--) {
                // swap without temporary
                c[i] ^= c[j];
                c[j] ^= c[i];
                c[i] ^= c[j];
            }
            c += m;
        }
        *(c) = '\0';
    }
    return s;
}
char* utils_nspire_t::safe_gets(char *out){
    do{
        csl >> out;
    }while(!out[0]);
    return out;
}
void utils_nspire_t::print_double(double num){
    char str[32];
    csl << double_to_ascii(num, str);
}
void utils_nspire_t::print_uint8(uint8_t value){
    csl << (int)value;
}

class utils_t {
public:
	template <typename T> static void copy_vect(T *, T*, uint16_t);
	template <typename T> static void copy_matrix(T**, T**, uint16_t, uint16_t);
	template <typename T> static void zero_vect(T *, uint16_t);
	template <typename T> static void zero_matrix(T**, uint16_t, uint16_t);
	template <typename T> static T mult_matrix_single_row_col(T**, T**, uint16_t, uint16_t, uint16_t);
	template <typename T> static void mult_mat_mat(T**, T**, T**, uint16_t, uint16_t, uint16_t);
	template <typename T> static void mult_mat_vect(T**, T*, T*, uint16_t, uint16_t);
	template <typename T> static void mult_vect_num(T*, T, T*, uint16_t);
	template <typename T> static double norm(T*, uint16_t);
};
class netlist_t {
public:
    struct row_t;
	uint8_t nodes;
	uint8_t components;
	row_t *row;
	netlist_t(uint8_t nod, uint8_t comp);
	~netlist_t(void);
};
class component_t {
public:
	static void insert_all(double **G, double *X, netlist_t *netlist, uint8_t netlist_pos, uint8_t *vs_pos);
	static void nop_f(uint8_t netlist_pos, netlist_t* netlist, double *data, double *target);
	static void voltage_source_f(uint8_t netlist_pos, netlist_t* netlist, double *data, double *target);
	static void current_source_f(uint8_t netlist_pos, netlist_t* netlist, double *data, double *target);
	static void diode_f(uint8_t netlist_pos, netlist_t* netlist, double *data, double *target);
    static void fet_n_f(uint8_t netlist_pos, netlist_t* netlist, double *data, double *target);
    static void fet_p_f(uint8_t netlist_pos, netlist_t* netlist, double *data, double *target);
    static double voltage_source_cf(uint8_t netlist_pos, netlist_t* netlist, double *data);
    static double current_source_cf(uint8_t netlist_pos, netlist_t* netlist, double *data);
    static double resistor_cf(uint8_t netlist_pos, netlist_t* netlist, double *data);
    static double diode_cf(uint8_t netlist_pos, netlist_t* netlist, double *data);
    static double fet_n_cf(uint8_t netlist_pos, netlist_t* netlist, double *data);
    static double fet_p_cf(uint8_t netlist_pos, netlist_t* netlist, double *data);
};
class nod_t{
public:
    uint8_t nodes;//how many nodes
	uint8_t dim;//nodal arrays length
	uint8_t components;//netlist size
	netlist_t *netlist;//source of a lot of data
	double LFS0;
	double *X;
	double *X_temp;
    double *S0;
    double *FS0;
    double *S0_temp;
    double *S1;
    double *FS1;
    double *S1_best;
    double *FS1_best;
    double *T1;
    double *NT1;
    double *U1;
    double *update_array_temp;
    double **G;//conductance
    double **R;//resistance
    double **JN;//Jacobian
    double **JNI;//inverse Jacobian
    double **updated_jacobian(double);//finds the numerical Jacobian
    int *inverse_aux_P;
    double **inverse_aux_A;
    void update_array(double *, double *);
    double eval_row(uint8_t, double *, double *);
    double* eval_f(double *, double *, double *);
	uint32_t iterate(uint32_t, double, double, double, double);
    nod_t(netlist_t *source);
    ~nod_t(void);
};
class math_t {
public:
	static void sub(double *, double *, double *, uint8_t);
	static void add(double *, double *, double *, uint8_t);
	static double** invert(double **, double **, uint8_t, double, double **, int *);
private:
	static int LUPDecompose(double **A, int N, double Tol, int *P);
	static void LUPSolve(double **A, int *P, double *b, int N, double *x);
	static void LUPInvert(double **A, int *P, int N, double **IA);
	static double LUPDeterminant(double **A, int *P, int N);
};

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
template <typename T> void utils_t::mult_mat_mat(T** M, T** N, T** OUT, uint16_t i1, uint16_t j2, uint16_t len) {
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
template <typename T> void utils_t::mult_vect_num(T* vect, T num, T* out, uint16_t len){
	while (len){
		len--;
		out[len] = vect[len]*num;
	}
}

template <typename T> double utils_t::norm(T* vect, uint16_t len){
    double acc = 0;
    while(len){
        len--;
        acc+=pow(vect[len],2);
    }
    return sqrt(acc);
}

struct netlist_t::row_t{
    uint8_t type;
    double value[NETLIST_PARAM_MAX];
    uint8_t node[NETLIST_CONNECTION_MAX];
    void (*update)(uint8_t netlist_pos, netlist_t* netlist, double *data, double *target);
    double (*current)(uint8_t netlist_pos, netlist_t* netlist, double *data);
    char alias[NETLIST_ALIAS_SIZE];
};

netlist_t::netlist_t(uint8_t amnt_nodes, uint8_t amnt_comp) {
	nodes = amnt_nodes - 1;
	components = amnt_comp;
	if (amnt_nodes) {
		row = new row_t[amnt_comp];
	}
}
netlist_t::~netlist_t() {
	if (nodes) {
		delete[] row;
	}
}


void component_t::insert_all(double **G, double *X, netlist_t *netlist, uint8_t netlist_pos, uint8_t *vs_pos) {
    while(netlist_pos){
        netlist_pos--;
        switch (netlist->row[netlist_pos].type){
        case(voltage_source):
            if(netlist->row[netlist_pos].node[positive]){
                G[netlist->row[netlist_pos].node[positive] - 1][*vs_pos - 1] = 1;
                G[*vs_pos - 1][netlist->row[netlist_pos].node[positive] - 1] = 1;
            }
            if(netlist->row[netlist_pos].node[negative]){
                G[netlist->row[netlist_pos].node[positive] - 1][*vs_pos - 1] = -1;
                G[*vs_pos - 1][netlist->row[netlist_pos].node[positive] - 1] = -1;
            }
            netlist->row[netlist_pos].node[hidden] = *vs_pos;
            netlist->row[netlist_pos].update = voltage_source_f;
            netlist->row[netlist_pos].current = voltage_source_cf;
            (*vs_pos)--;
            break;
        case(current_source):
            if (netlist->row[netlist_pos].node[positive]) {
                G[netlist->row[netlist_pos].node[positive] - 1][netlist->row[netlist_pos].node[positive] - 1] += GMIN;
                if (netlist->row[netlist_pos].node[negative]) {
                    G[netlist->row[netlist_pos].node[positive] - 1][netlist->row[netlist_pos].node[negative] - 1] += -GMIN;
                }
            }
            if (netlist->row[netlist_pos].node[negative]) {
                G[netlist->row[netlist_pos].node[negative] - 1][netlist->row[netlist_pos].node[negative] - 1] += GMIN;
                if (netlist->row[netlist_pos].node[positive]) {
                    G[netlist->row[netlist_pos].node[negative] - 1][netlist->row[netlist_pos].node[positive] - 1] += -GMIN;
                }
            }
            netlist->row[netlist_pos].update = current_source_f;
            netlist->row[netlist_pos].current = current_source_cf;
            break;
        case(resistor):
            if (netlist->row[netlist_pos].node[positive]) {
                G[netlist->row[netlist_pos].node[positive] - 1][netlist->row[netlist_pos].node[positive] - 1] += 1 / netlist->row[netlist_pos].value[resistance];
                if (netlist->row[netlist_pos].node[negative]) {
                    G[netlist->row[netlist_pos].node[positive] - 1][netlist->row[netlist_pos].node[negative] - 1] += -1 / netlist->row[netlist_pos].value[resistance];
                }
            }
            if (netlist->row[netlist_pos].node[negative]) {
                G[netlist->row[netlist_pos].node[negative] - 1][netlist->row[netlist_pos].node[negative] - 1] += 1 / netlist->row[netlist_pos].value[resistance];
                if (netlist->row[netlist_pos].node[positive]) {
                    G[netlist->row[netlist_pos].node[negative] - 1][netlist->row[netlist_pos].node[positive] - 1] += -1 / netlist->row[netlist_pos].value[resistance];
                }
            }
            netlist->row[netlist_pos].update = nop_f;
            netlist->row[netlist_pos].current = resistor_cf;
            break;
        case(capacitor):

            break;
        case(inductor):

            break;
        case(diode):
            if (netlist->row[netlist_pos].node[positive]) {
                G[netlist->row[netlist_pos].node[positive] - 1][netlist->row[netlist_pos].node[positive] - 1] += GMIN;
                if (netlist->row[netlist_pos].node[negative]) {
                    G[netlist->row[netlist_pos].node[positive] - 1][netlist->row[netlist_pos].node[negative] - 1] += -GMIN;
                }
            }
            if (netlist->row[netlist_pos].node[negative]) {
                G[netlist->row[netlist_pos].node[negative] - 1][netlist->row[netlist_pos].node[negative] - 1] += GMIN;
                if (netlist->row[netlist_pos].node[positive]) {
                    G[netlist->row[netlist_pos].node[negative] - 1][netlist->row[netlist_pos].node[positive] - 1] += -GMIN;
                }
            }
            netlist->row[netlist_pos].update = diode_f;
            netlist->row[netlist_pos].current = diode_cf;
            break;
        case (fet_n):
            if (netlist->row[netlist_pos].node[fet_source]) {
                G[netlist->row[netlist_pos].node[fet_source] - 1][netlist->row[netlist_pos].node[fet_source] - 1] += GMIN;
                if (netlist->row[netlist_pos].node[fet_drain]) {
                    G[netlist->row[netlist_pos].node[fet_source] - 1][netlist->row[netlist_pos].node[fet_drain] - 1] += -GMIN;
                }
                if (netlist->row[netlist_pos].node[fet_gate]) {
                    G[netlist->row[netlist_pos].node[fet_source] - 1][netlist->row[netlist_pos].node[fet_gate] - 1] += -GMIN;
                }
            }
            if (netlist->row[netlist_pos].node[fet_drain]) {
                G[netlist->row[netlist_pos].node[fet_drain] - 1][netlist->row[netlist_pos].node[fet_drain] - 1] += GMIN;
                if (netlist->row[netlist_pos].node[fet_source]) {
                    G[netlist->row[netlist_pos].node[fet_drain] - 1][netlist->row[netlist_pos].node[fet_source] - 1] += -GMIN;
                }
                if (netlist->row[netlist_pos].node[fet_gate]) {
                    G[netlist->row[netlist_pos].node[fet_drain] - 1][netlist->row[netlist_pos].node[fet_gate] - 1] += -GMIN;
                }
            }
            if (netlist->row[netlist_pos].node[fet_gate]) {
                G[netlist->row[netlist_pos].node[fet_gate] - 1][netlist->row[netlist_pos].node[fet_gate] - 1] += GMIN;
                if (netlist->row[netlist_pos].node[fet_source]) {
                    G[netlist->row[netlist_pos].node[fet_gate] - 1][netlist->row[netlist_pos].node[fet_source] - 1] += -GMIN;
                }
                if (netlist->row[netlist_pos].node[fet_drain]) {
                    G[netlist->row[netlist_pos].node[fet_gate] - 1][netlist->row[netlist_pos].node[fet_drain] - 1] += -GMIN;
                }
            }
            netlist->row[netlist_pos].update = fet_n_f;
            netlist->row[netlist_pos].current = fet_n_cf;
            break;
        case(fet_p):
            if (netlist->row[netlist_pos].node[fet_source]) {
                G[netlist->row[netlist_pos].node[fet_source] - 1][netlist->row[netlist_pos].node[fet_source] - 1] += GMIN;
                if (netlist->row[netlist_pos].node[fet_drain]) {
                    G[netlist->row[netlist_pos].node[fet_source] - 1][netlist->row[netlist_pos].node[fet_drain] - 1] += -GMIN;
                }
            }
            if (netlist->row[netlist_pos].node[fet_drain]) {
                G[netlist->row[netlist_pos].node[fet_drain] - 1][netlist->row[netlist_pos].node[fet_drain] - 1] += GMIN;
                if (netlist->row[netlist_pos].node[fet_source]) {
                    G[netlist->row[netlist_pos].node[fet_drain] - 1][netlist->row[netlist_pos].node[fet_source] - 1] += -GMIN;
                }
            }
            netlist->row[netlist_pos].update = fet_p_f;
            netlist->row[netlist_pos].current = fet_p_cf;
            break;
        default:
            break;
        }
    }
}

void component_t::nop_f(uint8_t netlist_pos, netlist_t* netlist, double *data, double *target){
}
void component_t::voltage_source_f(uint8_t netlist_pos, netlist_t* netlist, double *data, double *target){
    target[netlist->row[netlist_pos].node[hidden] - 1] = netlist->row[netlist_pos].value[voltage];
}
void component_t::current_source_f(uint8_t netlist_pos, netlist_t* netlist, double *data, double *target){
    if(netlist->row[netlist_pos].node[from]){
        target[netlist->row[netlist_pos].node[from] - 1] -= netlist->row[netlist_pos].value[current];
    }
    if(netlist->row[netlist_pos].node[to]){
        target[netlist->row[netlist_pos].node[to] - 1] += netlist->row[netlist_pos].value[current];
    }
}
void component_t::diode_f(uint8_t netlist_pos, netlist_t* netlist, double *data, double *target){
    double d_current;
    double d_is = netlist->row[netlist_pos].value[diode_is];
    double d_vt = netlist->row[netlist_pos].value[diode_vt];
    double d_vd = 0;
    if(netlist->row[netlist_pos].node[positive]){
        d_vd += data[netlist->row[netlist_pos].node[positive] - 1];
    }
    if(netlist->row[netlist_pos].node[negative]){
        d_vd -= data[netlist->row[netlist_pos].node[negative] - 1];
    }
    d_current = d_is * (exp(d_vd / d_vt) - 1);
    if(netlist->row[netlist_pos].node[positive]){
        target[netlist->row[netlist_pos].node[positive] - 1] -= d_current;
    }
    if(netlist->row[netlist_pos].node[negative]){
        target[netlist->row[netlist_pos].node[negative] - 1] += d_current;
    }
}
void component_t::fet_n_f(uint8_t netlist_pos, netlist_t* netlist, double *data, double *target){
    double vs, vd, vg, vds, vgs, id, w, l, kp, vt, lambda;
    w = netlist->row[netlist_pos].value[fet_w];
    l = netlist->row[netlist_pos].value[fet_l];
    kp = netlist->row[netlist_pos].value[fet_kp];
    vt = netlist->row[netlist_pos].value[fet_vt];
    lambda = netlist->row[netlist_pos].value[fet_lambda];
    if(netlist->row[netlist_pos].node[fet_source]){
        vs = data[netlist->row[netlist_pos].node[fet_source] - 1];
    }
    else{
        vs = 0;
    }
    if(netlist->row[netlist_pos].node[fet_drain]){
        vd = data[netlist->row[netlist_pos].node[fet_drain] - 1];
    }
    else{
        vd = 0;
    }
    if(netlist->row[netlist_pos].node[fet_gate]){
        vg = data[netlist->row[netlist_pos].node[fet_gate] - 1];
    }
    else{
        vg = 0;
    }
    vgs = vg-vs;
    vds = vd-vs;
    if(vgs <= vt){
        id = 0;
    }
    else if(vd<(vg-vt)){
        id = (kp/2)*(w/l)*(2*(vgs-vt)*(vds)-pow(vds, 2));
    }
    else{
        id = (kp/2)*(w/l)*pow(vgs-vt, 2)*(1+lambda*vds);
    }
    if(netlist->row[netlist_pos].node[fet_source]){
        target[netlist->row[netlist_pos].node[fet_source] - 1] += id;
    }
    if(netlist->row[netlist_pos].node[fet_drain]){
        target[netlist->row[netlist_pos].node[fet_drain] - 1] -= id;
    }
}
void component_t::fet_p_f(uint8_t netlist_pos, netlist_t* netlist, double *data, double *target){
    double vs, vd, vg, vds, vgs, id, w, l, kp, vt, lambda;
    w = netlist->row[netlist_pos].value[fet_w];
    l = netlist->row[netlist_pos].value[fet_l];
    kp = netlist->row[netlist_pos].value[fet_kp];
    vt = netlist->row[netlist_pos].value[fet_vt];
    lambda = netlist->row[netlist_pos].value[fet_lambda];
    if(netlist->row[netlist_pos].node[fet_source]){
        vs = data[netlist->row[netlist_pos].node[fet_source] - 1];
    }
    else{
        vs = 0;
    }
    if(netlist->row[netlist_pos].node[fet_drain]){
        vd = data[netlist->row[netlist_pos].node[fet_drain] - 1];
    }
    else{
        vd = 0;
    }
    if(netlist->row[netlist_pos].node[fet_gate]){
        vg = data[netlist->row[netlist_pos].node[fet_gate] - 1];
    }
    else{
        vg = 0;
    }
    vgs = vg-vs;
    vds = vd-vs;
    if(vgs >= vt){
        id = 0;
    }
    else if(vd>(vg-vt)){
        id = (kp/2)*(w/l)*(2*(vgs-vt)*(vds)-pow(vds, 2));
    }
    else{
        id = (kp/2)*(w/l)*pow(vgs-vt, 2)*(1+lambda*vds);
    }
    if(netlist->row[netlist_pos].node[fet_source]){
        target[netlist->row[netlist_pos].node[fet_source] - 1] -= id;
    }
    if(netlist->row[netlist_pos].node[fet_drain]){
        target[netlist->row[netlist_pos].node[fet_drain] - 1] += id;
    }
}

double component_t::voltage_source_cf(uint8_t netlist_pos, netlist_t* netlist, double *data){
    return data[netlist->row[netlist_pos].node[hidden] - 1];
}
double component_t::current_source_cf(uint8_t netlist_pos, netlist_t* netlist, double *data){
    return netlist->row[netlist_pos].value[current];
}
double component_t::resistor_cf(uint8_t netlist_pos, netlist_t* netlist, double *data){
    double r_vp = 0;
    double r_vn = 0;
    double r_r = netlist->row[netlist_pos].value[resistance];
    if(netlist->row[netlist_pos].node[positive]){
        r_vp = data[netlist->row[netlist_pos].node[positive] - 1];
    }
    if(netlist->row[netlist_pos].node[negative]){
        r_vn = data[netlist->row[netlist_pos].node[negative] - 1];
    }
    return (r_vp - r_vn)/r_r;
}
double component_t::diode_cf(uint8_t netlist_pos, netlist_t* netlist, double *data){
    double d_is = netlist->row[netlist_pos].value[diode_is];
    double d_vt = netlist->row[netlist_pos].value[diode_vt];
    double d_vd = 0;
    if(netlist->row[netlist_pos].node[positive]){
        d_vd += data[netlist->row[netlist_pos].node[positive] - 1];
    }
    if(netlist->row[netlist_pos].node[negative]){
        d_vd -= data[netlist->row[netlist_pos].node[negative] - 1];
    }
    return d_is * (exp(d_vd / d_vt) - 1);
}
double component_t::fet_n_cf(uint8_t netlist_pos, netlist_t* netlist, double *data){
    double vs, vd, vg, vds, vgs, w, l, kp, vt, lambda;
    w = netlist->row[netlist_pos].value[fet_w];
    l = netlist->row[netlist_pos].value[fet_l];
    kp = netlist->row[netlist_pos].value[fet_kp];
    vt = netlist->row[netlist_pos].value[fet_vt];
    lambda = netlist->row[netlist_pos].value[fet_lambda];
    if(netlist->row[netlist_pos].node[fet_source]){
        vs = data[netlist->row[netlist_pos].node[fet_source] - 1];
    }
    else{
        vs = 0;
    }
    if(netlist->row[netlist_pos].node[fet_drain]){
        vd = data[netlist->row[netlist_pos].node[fet_drain] - 1];
    }
    else{
        vd = 0;
    }
    if(netlist->row[netlist_pos].node[fet_gate]){
        vg = data[netlist->row[netlist_pos].node[fet_gate] - 1];
    }
    else{
        vg = 0;
    }
    vgs = vg-vs;
    vds = vd-vs;
    if(vgs <= vt){
        return 0;
    }
    else if(vds<(vgs-vt)){
        return (kp/2)*(w/l)*(2*(vgs-vt)*(vds)-pow(vds, 2));
    }
    else{
        return (kp/2)*(w/l)*pow(vgs-vt, 2)*(1+lambda*vds);
    }
}
double component_t::fet_p_cf(uint8_t netlist_pos, netlist_t* netlist, double *data){
    double vs, vd, vg, vds, vgs, w, l, kp, vt, lambda;
    w = netlist->row[netlist_pos].value[fet_w];
    l = netlist->row[netlist_pos].value[fet_l];
    kp = netlist->row[netlist_pos].value[fet_kp];
    vt = netlist->row[netlist_pos].value[fet_vt];
    lambda = netlist->row[netlist_pos].value[fet_lambda];
    if(netlist->row[netlist_pos].node[fet_source]){
        vs = data[netlist->row[netlist_pos].node[fet_source] - 1];
    }
    else{
        vs = 0;
    }
    if(netlist->row[netlist_pos].node[fet_drain]){
        vd = data[netlist->row[netlist_pos].node[fet_drain] - 1];
    }
    else{
        vd = 0;
    }
    if(netlist->row[netlist_pos].node[fet_gate]){
        vg = data[netlist->row[netlist_pos].node[fet_gate] - 1];
    }
    else{
        vg = 0;
    }
    vgs = vg-vs;
    vds = vd-vs;
    if(vgs >= vt){
        return 0;
    }
    else if(vds>(vgs-vt)){
        return -(kp/2)*(w/l)*(2*(vgs-vt)*(vds)-pow(vds, 2));
    }
    else{
        return -(kp/2)*(w/l)*pow(vgs-vt, 2)*(1+lambda*vds);
    }
}



nod_t::nod_t(netlist_t *src) {
	uint8_t pos = src->components;
	uint8_t voltages;
	netlist = src;
	nodes = src->nodes;
	dim = src->nodes;
	components = src->components;
	if (dim) {
		while (pos) {
			pos--;
            switch (src->row[pos].type){
			case (voltage_source):
				dim++;
				break;
			default:
				break;
			}
		}
		X = new double[dim];
		X_temp = new double[dim];
        S0 = new double[dim];
        FS0 = new double[dim];
        S0_temp = new double[dim];
        S1 = new double[dim];
        FS1 = new double[dim];
        S1_best = new double[dim];
        FS1_best = new double[dim];
        T1 = new double[dim];
        NT1 = new double[dim];
        U1 = new double[dim];
        update_array_temp = new double[dim];
        G = new double*[dim];
        R = new double*[dim];
        JN = new double*[dim];
        JNI = new double*[dim];
        inverse_aux_P = new int[dim+1];
        inverse_aux_A = new double*[dim];
		pos = dim;
		while (pos) {
			pos--;
			G[pos] = new double[dim];
		}
		pos = dim;
		while (pos) {
			pos--;
			R[pos] = new double[dim];
		}
		pos = dim;
		while (pos) {
			pos--;
			JN[pos] = new double[dim];
		}
		pos = dim;
		while (pos) {
			pos--;
			JNI[pos] = new double[dim];
		}
		pos = dim;
		while (pos) {
			pos--;
			inverse_aux_A[pos] = new double[dim];
		}
		utils_t::zero_vect(X, dim);
		utils_t::zero_vect(S0, dim);
		utils_t::zero_matrix(G, dim, dim);
		utils_t::zero_matrix(R, dim, dim);
		//utils_t::zero_matrix(JNI, dim, dim);
		voltages = dim;
		component_t::insert_all(G, X, netlist, components, &voltages);
		math_t::invert(G, R, dim, NOD_R_TOL, inverse_aux_A, inverse_aux_P);
	}
}
nod_t::~nod_t() {
    delete[] X;
	delete[] X_temp;
    delete[] S0;
    delete[] FS0;
    delete[] S0_temp;
    delete[] S1;
    delete[] FS1;
    delete[] S1_best;
    delete[] FS1_best;
    delete[] T1;
    delete[] NT1;
    delete[] U1;
    delete[] update_array_temp;
    delete[] inverse_aux_P;
	uint8_t pos = dim;
	if (dim) {
		while (pos) {
			pos--;
			delete[] G[pos];
		}
		delete[] G;
		pos = dim;
		while (pos) {
			pos--;
			delete[] R[pos];
		}
		delete[] R;
		pos = dim;
		while (pos) {
			pos--;
			delete[] JN[pos];
		}
		delete[] JN;
		pos = dim;
		while (pos) {
			pos--;
			delete[] JNI[pos];
		}
		delete[] JNI;
		pos = dim;
		while (pos) {
			pos--;
			delete[] inverse_aux_A[pos];
		}
		delete[] inverse_aux_A;
	}
}

double** nod_t::updated_jacobian(double jacobian_step) {
	uint8_t i = dim;
	uint8_t j;
	update_array(S0, X);
	utils_t::copy_vect(S0, S0_temp, dim);
	while (i) {
		i--;
		j = dim;
		while (j) {
			j--;
			S0_temp[j] += jacobian_step;
			update_array(S0_temp, X_temp);
			JN[i][j] = (eval_row(i, S0_temp, X_temp) - eval_row(i, S0, X)) / jacobian_step;
			S0_temp[j] -= jacobian_step;
		}
	}
	return JN;
}
void nod_t::update_array(double *data, double *target){
    uint8_t pos = components;
    utils_t::zero_vect(target, dim);
    while(pos){
        pos--;
        netlist->row[pos].update(pos, netlist, data, target);
    }
}

double nod_t::eval_row(uint8_t i, double *yy, double *xx) {
	double acc;
	uint8_t aux = dim;
	acc = yy[i];
	while (aux) {
		aux--;
		acc -= R[i][aux] * xx[aux];
	}
	return acc;
}
double* nod_t::eval_f(double *yy, double *xx, double *output) {
	uint8_t i_max = dim;
	update_array(yy, xx);
	while (i_max){
		i_max--;
		output[i_max] = eval_row(i_max, yy, xx);
	}
	return output;
}
uint32_t nod_t::iterate(uint32_t k, double sub_step, double error_max, double jacobian_step, double jacobian_inverse_tol) {
    uint8_t unlimited = !k;
    uint8_t not_done;
    uint8_t first_sub_step;
    double LFS1_best;
    double LT1;
    double LFS1;
    double lambda;
    eval_f(S0, X, FS0);
    LFS0 = utils_t::norm(FS0, dim);
    while((LFS0 > error_max) && (k || unlimited)){
        if(k){
            k--;
        }
        lambda = 1;
        not_done = TRUE;
        //while (not_done){
            /* T1 = Jn^-1*F(S0);
             * S1 = S0 - Lam * T1;
             * if(len(F(S1)) > len(F(S0))){
             *     if(Lam*len(T1) > STEP){
             *         Lam /= 10;
             *     }
             *     else{
             *         S0 -= STEP*T1/len(T1);
             *         not_done = FALSE;
             *     }
             * }
             * else{
             *     S0 = S1;
             *     not_done = FALSE;
             * }
             */

            /* "static double FS0"
             * "static double LFS0"
             * JNI = Jn^-1;
             * T1 = JNI*FS0;
             * LT1 = len(T1);
             * NT1 = T1/LT1;
             * U1 = Lam*T1;
             * S1 = S0 - U1;
             * FS1 = F(S1);
             * LFS1 = len(FS1);
             * if(LFS1 > LFS0){
             *     if(U1 > STEP){
             *         Lam /= 10;
             *     }
             *     else{
             *         S0 -= STEP*NT1;
             *         not_done = FALSE;
             *     }
             * }
             * else{
             *     S0 = S1;
             *     FS0  = FS1;
             *     LFS0 = LFS1;
             *     not_done = FALSE;
             * }
             */

            /* "static double FS0"
             * "static double LFS0"
             * "FIRST = TRUE"
             * "BEST_S1";
             * "BEST_LFS1";
             * JNI = Jn^-1;  //<<<BEFORE WHILE
             * T1 = JNI*FS0; //<<<BEFORE WHILE
             * LT1 = len(T1);
             * NT1 = T1/LT1;
             * U1 = Lam*T1;
             * S1 = S0 - U1;
             * FS1 = F(S1);
             * LFS1 = len(FS1);
             * if(LFS1 >= LFS0){
             *     if(U1 >= STEP){
             *         Lam /= 10;
             *         if((BEST_LFS1 > LFS1) || FIRST){
             *             FIRST = FALSE;
             *             BEST_S1 = S1;
             *             BEST_LFS1 = LFS1;
             *         }
             *     }
             *     else{
             *         S0 = BEST_S1;
             *         not_done = FALSE;
             *     }
             * }
             * else{
             *     S0 = S1;
             *     FS0  = FS1;
             *     LFS0 = LFS1;
             *     not_done = FALSE;
             * }
             */
        //}
        first_sub_step = TRUE;
        math_t::invert(updated_jacobian(jacobian_step), JNI, dim, jacobian_inverse_tol, inverse_aux_A, inverse_aux_P);
        utils_t::mult_mat_vect(JNI, FS0, T1, dim, dim);
        LT1 = utils_t::norm(T1, dim);
        utils_t::mult_vect_num(T1, 1/LT1, NT1, dim);
        while(not_done){
            utils_t::mult_vect_num(T1, lambda, U1, dim);
            math_t::sub(S0, U1, S1, dim);
            eval_f(S1, X, FS1);
            LFS1 = utils_t::norm(FS1, dim);
            if(LFS1 > LFS0){
                if(lambda*LT1 >= sub_step){
                    lambda /= 10;
                    if((LFS1_best > LFS1) || first_sub_step){
                        first_sub_step = FALSE;
                        utils_t::copy_vect(S1, S1_best, dim);
                        utils_t::copy_vect(FS1, FS1_best, dim);
                        LFS1_best = LFS1;
                    }
                }
                else{
                    not_done = FALSE;
                    if(first_sub_step){
                        utils_t::copy_vect(S1, S0, dim);
                        utils_t::copy_vect(FS1, FS0, dim);
                        LFS0 = LFS1;
                    }
                    else{
                        utils_t::copy_vect(S1_best, S0, dim);
                        utils_t::copy_vect(FS1_best, FS0, dim);
                        LFS0 = LFS1_best;
                    }
                }
            }
            else{
                not_done = FALSE;
                utils_t::copy_vect(S1, S0, dim);
                utils_t::copy_vect(FS1, FS0, dim);
                LFS0 = LFS1;
            }
        }
    }
    return k;
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
			(P[N])++;
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
double** math_t::invert(double **in, double** out, uint8_t dim, double tol, double **A, int *P) {
	utils_t::copy_matrix(in, A, dim, dim);
	math_t::LUPDecompose(A, dim, tol, P);
	math_t::LUPInvert(A, P, dim, out);
	return out;
}

class spice_t{
public:
    uint32_t iterations;
    spice_t();
    ~spice_t();
    netlist_t *netlist;
    nod_t *solver;
    void simulate(void);
    void request_components(void);
    void small_signal_analysis(void);
    void show_result(void);
    void show_voltages(void);
    void show_currents(void);
    void show_statistics(void);
    static void show_init_screen(void);
};
spice_t::spice_t(){
    uint16_t nodes;
    uint16_t components;
    char str[32];
    csl << "Nodes: ";
    nodes = atoi(utils_nspire_t::safe_gets(str));
    csl << "Components: ";
    components = atoi(utils_nspire_t::safe_gets(str));
    csl << "GMIN (0 = default): ";
    gmin = atoi(utils_nspire_t::safe_gets(str));;
    gmin = gmin ? gmin : 10e-12;
    netlist = new netlist_t(nodes, components);
    request_components();
    solver = new nod_t(netlist);
}
spice_t::~spice_t(){
    delete netlist;
    delete solver;
}
void spice_t::show_init_screen(void){
    csl << "BSpice Alpha 0.2" << nio::endl;
    csl << "GND = Node 0 = 0" << nio::endl;
    csl << "Component code:" << nio::endl;
    csl << resistor_name << ": " << resistor_code << nio::endl;
    csl << voltage_source_name << ": " << voltage_source_code << nio::endl;
    csl << current_source_name << ": " << current_source_code << nio::endl;
    csl << capacitor_name << ": " << capacitor_code << nio::endl;
    csl << inductor_name << ": " << inductor_code << nio::endl;
    csl << diode_name << ": " << diode_code << nio::endl;
    csl << bjt_npn_name << ": " << bjt_npn_code << nio::endl;
    csl << bjt_pnp_name << ": " << bjt_pnp_code << nio::endl;
    csl << fet_n_name << ": " << fet_n_code << nio::endl;
    csl << fet_p_name << ": " << fet_p_code << nio::endl;
}
void spice_t::simulate(){
    char initial_guess;
    uint16_t pos = netlist->nodes;
    double step;
    double close_enough;
    double jacobian_step;
    double jacobian_tolerance;
    csl << "Iterations (0 = unlimited): ";
    iterations = atoi(utils_nspire_t::safe_gets(str));
    csl << "Bad step max length [V] (0 = default): ";
    step = atof(utils_nspire_t::safe_gets(str));
    step = step ? step : 1e-3;
    csl << "Error vector max length [V] (0 = default): ";
    close_enough = atof(utils_nspire_t::safe_gets(str));
    close_enough = close_enough ? close_enough : NOD_DEF_ERROR_VECTOR_MAX_LENGTH;
    csl << "Numeric derivative step [V] (0 = default): ";
    jacobian_step = atof(utils_nspire_t::safe_gets(str));
    jacobian_step = jacobian_step ? jacobian_step : NOD_DEF_JACOBIAN_STEP;
    csl << "Jacobian inverse tolerance (0 = default): ";
    jacobian_tolerance = atof(utils_nspire_t::safe_gets(str));
    jacobian_tolerance = jacobian_tolerance ? jacobian_tolerance : NOD_DEF_JNI_TOL;
    csl << "Initial guess? (y/n): ";
    csl >> initial_guess;
    while(pos && (initial_guess != 'n')){
        csl << "Node " << pos << ": ";
        pos--;
        solver->S0[pos] = atof(utils_nspire_t::safe_gets(str));
    }
    iterations -= solver->iterate(iterations, step, close_enough, jacobian_step, jacobian_tolerance);
}
void spice_t::request_components(){
    char code[COMPONENT_CODE_MAX_LENGTH];
    uint8_t confirmation;
    uint8_t pos = netlist->components;
    char str[32];
    int gambiarra;
    while(pos){
        csl << "----- Component " << (unsigned int)pos << " -----" << nio::endl;
        pos--;
        csl << "Component: ";
        utils_nspire_t::safe_gets(code); =
        csl << "Alias: ";
        utils_nspire_t::safe_gets(netlist->row[pos].alias);
        if(!strcmp(code, resistor_code)){
            netlist->row[pos].type = resistor;
            csl << "Node +: ";
            netlist->row[pos].node[positive] = atoi(utils_nspire_t::safe_gets(str));
            csl << "Node -: ";
            netlist->row[pos].node[negative] = atoi(utils_nspire_t::safe_gets(str));
            csl << "Resistance: ";
            netlist->row[pos].value[resistance] = atof(utils_nspire_t::safe_gets(str));
        }
        else if(!strcmp(code, voltage_source_code)){
            netlist->row[pos].type = voltage_source;
            csl << "Node +: ";
            netlist->row[pos].node[positive] = atoi(utils_nspire_t::safe_gets(str));
            csl << "Node -: ";
            netlist->row[pos].node[negative] = atoi(utils_nspire_t::safe_gets(str));
            csl << "Voltage: ";
            netlist->row[pos].value[voltage] = atof(utils_nspire_t::safe_gets(str));
        }
        else if(!strcmp(code, current_source_code)){
            netlist->row[pos].type = current_source;
            csl << "Current from: ";
            netlist->row[pos].node[from] = atoi(utils_nspire_t::safe_gets(str));
            csl << "to: ";
            netlist->row[pos].node[to] = atoi(utils_nspire_t::safe_gets(str));
            csl << "Current: ";
            netlist->row[pos].value[current] = atof(utils_nspire_t::safe_gets(str));
        }
        else if(!strcmp(code, diode_code)){
            netlist->row[pos].type = diode;
            csl << "Node +: ";
            netlist->row[pos].node[from] = atoi(utils_nspire_t::safe_gets(str));
            csl << "Node -: ";
            netlist->row[pos].node[to] = atoi(utils_nspire_t::safe_gets(str));
            csl << "Is: ";
            netlist->row[pos].value[diode_is] = atof(utils_nspire_t::safe_gets(str));
            csl << "Vt: ";
            netlist->row[pos].value[diode_vt] = atof(utils_nspire_t::safe_gets(str));
        }
        else if(!strcmp(code, fet_n_code)){
            netlist->row[pos].type = fet_n;
            csl << "Drain node: ";
            netlist->row[pos].node[fet_drain] = atoi(utils_nspire_t::safe_gets(str));
            csl << "Source node: ";
            netlist->row[pos].node[fet_source] = atoi(utils_nspire_t::safe_gets(str));
            csl << "Gate node: ";
            netlist->row[pos].node[fet_gate] = atoi(utils_nspire_t::safe_gets(str));
            csl << "W: ";
            netlist->row[pos].value[fet_w] = atof(utils_nspire_t::safe_gets(str));
            csl << "L: ";
            netlist->row[pos].value[fet_l] = atof(utils_nspire_t::safe_gets(str));
            csl << "Kp: ";
            netlist->row[pos].value[fet_kp] = atof(utils_nspire_t::safe_gets(str));
            csl << "Lambda: ";
            netlist->row[pos].value[fet_lambda] = atof(utils_nspire_t::safe_gets(str));
            csl << "Vth: ";
            netlist->row[pos].value[fet_vt] = atof(utils_nspire_t::safe_gets(str));
        }
        else if(!strcmp(code, fet_p_code)){
            netlist->row[pos].type = fet_p;
            csl << "Drain node: ";
            netlist->row[pos].node[fet_drain] = atoi(utils_nspire_t::safe_gets(str));
            csl << "Source node: ";
            netlist->row[pos].node[fet_source] = atoi(utils_nspire_t::safe_gets(str));
            csl << "Gate node: ";
            netlist->row[pos].node[fet_gate] = atoi(utils_nspire_t::safe_gets(str));
            csl << "W: ";
            netlist->row[pos].value[fet_w] = atof(utils_nspire_t::safe_gets(str));
            csl << "L: ";
            netlist->row[pos].value[fet_l] = atof(utils_nspire_t::safe_gets(str));
            csl << "Kp: ";
            netlist->row[pos].value[fet_kp] = atof(utils_nspire_t::safe_gets(str));
            csl << "Lambda: ";
            netlist->row[pos].value[fet_lambda] = atof(utils_nspire_t::safe_gets(str));
            csl << "Vth: ";
            netlist->row[pos].value[fet_vt] = atof(utils_nspire_t::safe_gets(str));
        }
        csl << "confirm? (y/n): ";
        csl >> confirmation;
        if(confirmation == 'n'){
            pos++;
        }
    }
    csl << "Press any key to continue... ";
    wait_key_pressed();
}
void spice_t::small_signal_analysis(void){
    uint8_t pos = netlist->components;
    char str[32];
    uint8_t line = 0;
    while(pos){
        pos--;
        if(line >= 28){
            csl << "Press any key to continue... ";
            wait_key_pressed();
            line = 0;
        }
        switch(netlist->row[pos].type){
        case(fet_p):
            //don't insert break here
        case(fet_n):
            double vs, vd, vg, vds, vgs, w, l, kp, vt, lambda, id;
            w = netlist->row[pos].value[fet_w];
            l = netlist->row[pos].value[fet_l];
            kp = netlist->row[pos].value[fet_kp];
            vt = netlist->row[pos].value[fet_vt];
            lambda = netlist->row[pos].value[fet_lambda];
            id = netlist->row[pos].current(pos, netlist, solver->S0);
            if(netlist->row[pos].node[fet_source]){
                vs = solver->S0[netlist->row[pos].node[fet_source] - 1];
            }
            else{
                vs = 0;
            }
            if(netlist->row[pos].node[fet_drain]){
                vd = solver->S0[netlist->row[pos].node[fet_drain] - 1];
            }
            else{
                vd = 0;
            }
            if(netlist->row[pos].node[fet_gate]){
                vg = solver->S0[netlist->row[pos].node[fet_gate] - 1];
            }
            else{
                vg = 0;
            }
            vgs = vg-vs;
            vds = vd-vs;
            if(vds >= vgs-vt){
                csl << netlist->row[pos].alias << " - pentodo - Rds: " << utils_nspire_t::double_to_ascii((1+lambda*vds)/(lambda*id), str) << nio::endl;
                csl << netlist->row[pos].alias << " - pentodo - Gm: " << utils_nspire_t::double_to_ascii(2*id/(vgs-vt), str) << nio::endl;
            }
            else{
                csl << netlist->row[pos].alias << " - triodo - Rds: " << utils_nspire_t::double_to_ascii(1/((w/l)*kp*((1+2*lambda*vds)*(vgs-vt) - vds*(1+1.5*lambda*vds))), str) << nio::endl;
                csl << netlist->row[pos].alias << " - triodo - Gm: " << utils_nspire_t::double_to_ascii((w/l)*kp*vds*(1+lambda*vds), str) << nio::endl;
            }
            line += 2;
        default:
            break;
        }
    }
    csl << "Press any key to continue... ";
    wait_key_pressed();
}
void spice_t::show_result(){
    uint8_t pos = netlist->nodes;
    uint8_t line = 0;
    while (pos){
        pos--;
        if(line == 29){
            csl << "Press any key to continue... ";
            wait_key_pressed();
            line = 0;
        }
        csl << "Node " << (int)(pos+1) << ": " << utils_nspire_t::double_to_ascii(solver->S0[pos], str) << "V" << nio::endl;
        line++;
    }
    csl << "Press any key to continue... ";
    wait_key_pressed();
}
void spice_t::show_voltages(){
    uint16_t pos = netlist->components;
    uint8_t line = 0;
    while(pos){
        pos--;
        if(line >= 28){
            csl << "Press any key to continue... ";
            wait_key_pressed();
            line = 0;
        }
        switch(netlist->row[pos].type){
        case(current_source):
            csl << netlist->row[pos].alias << ": " << utils_nspire_t::double_to_ascii(solver->S0[netlist->row[pos].node[to] - 1] - solver->S0[netlist->row[pos].node[from] - 1], str) << "V" <<nio::endl;
            line++;
            break;
        case(fet_n):
            csl << netlist->row[pos].alias << " - Vds: " << utils_nspire_t::double_to_ascii(solver->S0[netlist->row[pos].node[fet_drain] - 1] - solver->S0[netlist->row[pos].node[fet_source] - 1], str) << "V" <<nio::endl;
            csl << netlist->row[pos].alias << " - Vgs: " << utils_nspire_t::double_to_ascii(solver->S0[netlist->row[pos].node[fet_gate] - 1] - solver->S0[netlist->row[pos].node[fet_source] - 1], str) << "V" <<nio::endl;
            line += 2;
            break;
        case(fet_p):
            csl << netlist->row[pos].alias << " - Vds: " << utils_nspire_t::double_to_ascii(solver->S0[netlist->row[pos].node[fet_drain] - 1] - solver->S0[netlist->row[pos].node[fet_source] - 1], str) << "V" <<nio::endl;
            csl << netlist->row[pos].alias << " - Vgs: " << utils_nspire_t::double_to_ascii(solver->S0[netlist->row[pos].node[fet_gate] - 1] - solver->S0[netlist->row[pos].node[fet_source] - 1], str) << "V" <<nio::endl;
            line += 2;
            break;
        default:
            csl << netlist->row[pos].alias << ": " << utils_nspire_t::double_to_ascii(solver->S0[netlist->row[pos].node[positive] - 1] - solver->S0[netlist->row[pos].node[negative] - 1], str) << "V" <<nio::endl;
            line++;
            break;
        }
    }
    csl << "Press any key to continue... ";
    wait_key_pressed();
}
void spice_t::show_currents(){
    uint16_t pos = netlist->components;
    uint8_t line = 0;
    while(pos){
        pos--;
        if(line >= 29){
            csl << "Press any key to continue... ";
            wait_key_pressed();
            line = 0;
        }
        switch(netlist->row[pos].type){
        case(current_source):
            csl << netlist->row[pos].alias << ": " << utils_nspire_t::double_to_ascii(netlist->row[pos].current(pos, netlist, solver->S0), str) << "A" <<nio::endl;
            break;
        case(fet_n):
            csl << netlist->row[pos].alias << " - Id: " << utils_nspire_t::double_to_ascii(netlist->row[pos].current(pos, netlist, solver->S0), str) << "A" <<nio::endl;
            break;
        case(fet_p):
            csl << netlist->row[pos].alias << " - Id: " << utils_nspire_t::double_to_ascii(netlist->row[pos].current(pos, netlist, solver->S0), str) << "A" <<nio::endl;
            break;
        default:
            csl << netlist->row[pos].alias << ": " << utils_nspire_t::double_to_ascii(netlist->row[pos].current(pos, netlist, solver->S0), str) << "A" <<nio::endl;
            break;
        }
        line++;
    }
    csl << "Press any key to continue... ";
    wait_key_pressed();
}
void spice_t::show_statistics(){
    csl << "LFS0: " << utils_nspire_t::double_to_ascii(solver->LFS0, str) << nio::endl;
    csl << "Iterations: " << (int)iterations << nio::endl;
}

int main() {
    spice_t::show_init_screen();
    csl << "-----" << nio::endl;
	spice_t spice;
	spice.simulate();
	csl << "-----" << nio::endl;
	spice.show_statistics();
	csl << "-----" << nio::endl;
    spice.show_result();
    csl << "-----" << nio::endl;
    spice.show_voltages();
    csl << "-----" << nio::endl;
    spice.show_currents();
    csl << "-----" << nio::endl;
    spice.small_signal_analysis();
    csl << "Press any key to return to OS... ";
    wait_key_pressed();
	return 0;
}

