/*
William Zhang

Last Updated: March 2, 2021

Implementation of the Hydra protocol described in "Hydra: Succinct Fully
Pipelineable Interactive Arguments of Knowledge" by Zhang and Xia.

Polynomial and arithmetic logic based on Thaler's implementation described in
"Practical Verified Computation with Streaming Interactive Proofs" by Cormode,
Mitzenmacher, and Thaler.

This version is specific to testing with the SHA-256 boolean circuit.
*/

#include <bits/stdc++.h>
#include <omp.h>

#define MASK 4294967295 //2^32-1
#define PRIME 2305843009213693951 //2^61-1

//fast log2 (g++)
#define FAST_LOG2(x) (sizeof(unsigned long)*8 - 1 - __builtin_clzl((unsigned long)(x)))
#define FAST_LOG2_UP(x) ({ unsigned long log = FAST_LOG2(x); ((x) - (1 << log)) ? log + 1 : log; })

using namespace std;

typedef struct gate{
  int type; //-1 passthru, 0 and, 1 xor, 2, inv
  uint64_t val; //value of the gate
  int in1; //left parent
  int in2; //right parent
  vector<uint64_t> wire; //vector wire contribution
  vector<uint64_t> beta; //vector beta contribution
} gate;

//circuit represented as 2d vector of gates
typedef vector<vector<gate>> circ;
circ c;

//fast modular arithmetic for %PRIME
uint64_t mod(uint64_t x){
  uint64_t ret = (x >> 61) + (x & PRIME);
  return ret >= PRIME ? ret - PRIME : ret;
}

//fast modular multiplication for %PRIME
inline uint64_t modMult(uint64_t x, uint64_t y){
  uint64_t hi_x = x >> 32;
  uint64_t hi_y = y >> 32;
  uint64_t low_x = x & MASK;
  uint64_t low_y = y & MASK;
  uint64_t piece1 = mod((hi_x * hi_y) << 3);
  uint64_t z = (hi_x * low_y + hi_y * low_x);
  uint64_t hi_z = z >> 32;
  uint64_t low_z = z & MASK;
  uint64_t piece2 = mod((hi_z << 3) + mod((low_z << 32)));
  uint64_t piece3 = mod(low_x * low_y);
  uint64_t ret = mod(piece1 + piece2 + piece3);
  return ret;
}

//fast modular exponentiation for %PRIME
uint64_t modPow(uint64_t b, uint64_t e){
  uint64_t ret;
  if(e == 1){
    return b;
  }
  if(e == 0){
    return 1;
  }
  if((e & 1) == 0){
    ret = modPow(b, (e >> 1));
    return modMult(ret, ret);
  }
  else{
    return modMult(modPow(b, e-1), b);
  }
}

//extended euclidian algorithm
void eucl(uint64_t u, uint64_t &u1, uint64_t &u2, uint64_t &u3){
  u1 = 1;
  u2 = 0;
  u3 = u;
  uint64_t v1 = 0;
  uint64_t v2 = 1;
  uint64_t v3 = PRIME;
  uint64_t q;
  uint64_t t1;
  uint64_t t2;
  uint64_t t3;
  do{
    q = u3 / v3;
    t1 = mod(u1+PRIME - modMult(q, v1));
    t2 = mod(u2+PRIME - modMult(q, v2));
    t3 = mod(u3+PRIME - modMult(q, v3));
    u1 = v1;
    u2 = v2;
    u3 = v3;
    v1 = t1;
    v2 = t2;
    v3 = t3;
  } while(v3 != 0 && v3 != PRIME);
}

//modular multiplicative inverse
uint64_t inv(uint64_t a)
{
  uint64_t u1;
  uint64_t u2;
  uint64_t u3;
  eucl(a, u1, u2, u3);
  if(u3 == 1){
    return mod(u1);
  }
  else{
    return 0;
  }
}

//fast modular division for %PRIME
uint64_t modDiv(uint64_t a, uint64_t b){
  uint64_t i = inv(b);
  return modMult(i, a);
}

//lagrange interpolation of v
uint64_t interp(uint64_t v, vector<uint64_t> &r, uint64_t n, int r_i){
  uint64_t x = v;
  uint64_t c = 1;
  for(uint64_t i = 0; i < n; ++i){
    if(x & 1){
      c = modMult(c, r[r_i + i]);
    }
    else{
      c = modMult(c, 1+PRIME - r[r_i + i]);
    }
    x >>= 1;
  }
  return c;
}

//beta polynomial at k
uint64_t beta(int mi, vector<uint64_t> &z, uint64_t k){
  uint64_t ret = 1;
  uint64_t x = k;
  for(int i = 0; i < mi; ++i){
    if(x & 1){
      ret = modMult(ret, z[i]);
    }
    else{
      ret = modMult(ret, 1+PRIME - z[i]);
    }
    x >>= 1;
  }
  return ret;
}

//evaluate w_i polynomial
uint64_t w_i(int mi, int ni, int i, vector<uint64_t> &r, int r_i){
  uint64_t ret = 0;
  for(uint64_t k = 0; k < ni; ++k){
    ret = mod(ret + modMult(c[i][k].val, interp(k, r, mi, r_i)));
  }
  return ret;
}

//extrapolate polynomial v to r
uint64_t extrap(vector<uint64_t> &v, uint64_t n, uint64_t r){
  uint64_t ret = 0;
  uint64_t mult = 1;
  for(uint64_t i = 0; i < n; ++i){
    mult = 1;
    for(uint64_t j = 0; j < n; ++j){
      if(i > j){
        mult = (modMult(modMult(mult, mod(r - j+PRIME)), inv(i - j)));
      }
      if(i < j){
        mult = modMult(modMult(mult, mod(r - j+PRIME)), inv(mod(i+PRIME - j)));
      }
    }
    ret = mod(ret + modMult(mult, v[i]));
  }
  return ret;
}

uint64_t log2_max, gates_max;

bool HYDRA;

//main sumcheck protocol
pair<uint64_t, pair<uint64_t, uint64_t>> sumcheck(int lvl, vector<uint64_t> &z, uint64_t ri, pair<vector<uint64_t>, vector<uint64_t>> &challenges, int tmp_index){
  
  cout << "start" << endl;

  int mi = FAST_LOG2_UP(c[lvl].size());
  int ni = c[lvl].size();
  int mip1 = FAST_LOG2_UP(c[lvl-1].size());
  int nip1 = c[lvl-1].size();

  int num_var = mi + 2 * mip1;
  int num_effective = modPow(2, mip1);
  uint64_t alpha_p, beta_p;

  vector<uint64_t> r;

  if(!HYDRA){
    for(int i = 0; i < num_var; ++i){
      r.push_back(rand());
    }
  }
  else{
    for(int i = 0; i < mi; ++i){
      r.push_back(rand());
    }
    r.insert(r.end(), challenges.first.begin(), challenges.first.end());
    r.insert(r.end(), challenges.second.begin(), challenges.second.end());
  }

  cout << "here1" << endl;

  for(uint64_t k = 0; k < ni; ++k){
    c[lvl][k].beta[tmp_index] = beta(mi, z, k);
    c[lvl][k].wire[tmp_index] = 1;
  }

  vector<vector<uint64_t>> poly;
  for(int k = 0; k < num_var; ++k){
    poly.push_back({0, 0, 0});
  }

  vector<uint64_t> val;
  for(int k = 0; k < num_effective; ++k){
    val.push_back(k < nip1 ? c[lvl-1][k].val : 0);
  }

  cout << "here" << endl;

  uint64_t beta = 0;
  uint64_t k_j;
  uint64_t wire = 0;
  uint64_t w1 = 0, w2 = 0;
  uint64_t index;
  
  for(int j = 0; j < num_var; ++j){
    if(j > mi && (j != mi + mip1)){
      for(int k = 0; k < num_effective; ++k){
        index = k >> 1;
        if(k & 1){
          val[index] = mod(val[index] + modMult(val[k], r[j-1]));
        }
        else{
          val[index] = modMult(val[k], 1+PRIME - r[j-1]);
        }
      }
      num_effective >>= 1;
    }
    if(j == mi + mip1){
      num_effective = modPow(2, mip1);
      for(int k = 0; k < num_effective; ++k){
        val[k] = k < nip1 ? c[lvl-1][k].val : 0;
      }
      w1 = w_i(mip1, nip1, lvl-1, r, mi);
    }
    if(j == mi){
      beta = c[lvl][0].beta[tmp_index];
    }
    for(int k = 0; k < ni; ++k){
      if(j < mi){
        w1 = c[lvl-1][c[lvl][k].in1].val;
        w2 = c[lvl-1][c[lvl][k].in2].val;
      }
      else if(j < mi + mip1){
        w2 = c[lvl-1][c[lvl][k].in2].val;
      }
      if(j < mi){
        k_j = k >> j;
        if(k_j & 1){
          c[lvl][k].beta[tmp_index] = modMult(c[lvl][k].beta[tmp_index], inv(z[j]));
        }
        else{
          c[lvl][k].beta[tmp_index] = modMult(c[lvl][k].beta[tmp_index], inv(1+PRIME - z[j]));
        }
      }

      for(int m = 0; m <= 2; ++m){
        if(j < mi){
          if(m == 0){
            beta = modMult(c[lvl][k].beta[tmp_index], modMult(1+PRIME - z[j], 1+PRIME - m));
          }
          else if(m == 1){
            beta = modMult(c[lvl][k].beta[tmp_index], modMult(z[j], m));
          }
          else{
            beta = modMult(c[lvl][k].beta[tmp_index], mod(modMult(1+PRIME - z[j], 1+PRIME - m) + modMult(z[j], m)));
          }
          if(m == 2){
            c[lvl][k].beta[tmp_index] = modMult(c[lvl][k].beta[tmp_index], mod(modMult(1+PRIME - z[j], 1+PRIME - r[j]) + modMult(z[j], r[j])));
          }
          if(k_j & 1){
            if(m == 0){
              continue;
            }
            else{
              wire = modMult(c[lvl][k].wire[tmp_index], m);
            }
            if(m == 2){
              c[lvl][k].wire[tmp_index] = modMult(c[lvl][k].wire[tmp_index], r[j]);
            }
          }
          else{
            if(m==1){
              continue;
            }
            else{
              wire = modMult(c[lvl][k].wire[tmp_index], 1+PRIME - m);
            }
            if(m == 2){
              c[lvl][k].wire[tmp_index] = modMult(c[lvl][k].wire[tmp_index], 1+PRIME - r[j]);
            }
          }
        }
        if(j >= mi && j < mi+mip1){
          index = c[lvl][k].in1 >> (j - mi);
          if(index & 1){
            if(m == 0){
              continue;
            }
            else if(m == 1){
              w1 = val[index];
            }
            else{
              w1 = mod(modMult(val[index], m) + modMult(val[index-1], 1+PRIME - m));
            }
            wire = modMult(c[lvl][k].wire[tmp_index], m);
            if(m == 2){
              c[lvl][k].wire[tmp_index] = modMult(c[lvl][k].wire[tmp_index], r[j]);
            }
          }
          else{
            if(m == 1){
              continue;
            }
            else if(m == 0){
              w1 = val[index];
            }
            else{
              w1 = mod(modMult(val[index+1], m) + modMult(val[index], 1+PRIME - m));
            }
            wire = modMult(c[lvl][k].wire[tmp_index], 1+PRIME - m);
            if(m == 2){
              c[lvl][k].wire[tmp_index] = modMult(c[lvl][k].wire[tmp_index], 1+PRIME - r[j]);
            }
          }
        }
        else if(j >= mi+mip1){
          index = c[lvl][k].in2 >> (j - mi - mip1);
          if(index & 1){
            if(m == 0){
              continue;
            }
            else if(m == 1){
              w2 = val[index];
            }
            else{
              w2 = mod(modMult(val[index], m) + modMult(val[index-1], 1+PRIME - m));
            }
            wire = modMult(c[lvl][k].wire[tmp_index], m);
            if(m == 2){
              c[lvl][k].wire[tmp_index] = modMult(c[lvl][k].wire[tmp_index], r[j]);
            }
          }
          else{
            if(m == 1){
              continue;
            }
            else if(m == 0){
              w2 = val[index];
            }
            else{
              w2 = mod(modMult(val[index+1], m) + modMult(val[index], 1+PRIME - m));
            }
            wire = modMult(c[lvl][k].wire[tmp_index], 1+PRIME - m);
            if(m == 2){
              c[lvl][k].wire[tmp_index] = modMult(c[lvl][k].wire[tmp_index], 1+PRIME - r[j]);
            }
          }
        }

        if(c[lvl][k].type == -1){
          poly[j][m] = mod(poly[j][m] + mod(modMult(beta, modMult(wire, w1))));
        }
        else if(c[lvl][k].type == 0){
          poly[j][m] = mod(poly[j][m] + mod(modMult(beta, modMult(wire, modMult(w1, w2)))));
        }
        else if(c[lvl][k].type == 1){
          poly[j][m] = mod(poly[j][m] + mod(modMult(beta, modMult(wire, mod(w1 + w2)))));
        }
        else if(c[lvl][k].type == 2){
          poly[j][m] = mod(poly[j][m] + mod(modMult(beta, modMult(wire, mod(1+PRIME - w1)))));
        }
      }
    }
  }

  uint64_t check = 0;
  for(int j = 1; j < num_var; j++) {
    check = extrap(poly[j-1], 3, r[j-1]);
    if((check != mod(poly[j][0]+poly[j][1])) && (check != (mod(poly[j][0]+poly[j][1])+PRIME)) && (check != mod(poly[j][0]+poly[j][1])-PRIME)){
      cout << "PROTOCOL FAILED: intermediate sumcheck message is incorrect" << endl;
    }
  }
  w2 = w_i(mip1, nip1, lvl-1, r, mi+mip1);

  if(!HYDRA){
    if((lvl != c.size()-1) && (mod(poly[0][0]+poly[0][1]) != ri) && (mod(poly[0][0]+poly[0][1]) != ri-PRIME) && (mod(poly[0][0]+poly[0][1]) != ri+PRIME)){
      cout << "PROTOCOL FAILED: claimed value is incorrect" << endl;
    }

    uint64_t test_wire_passthru = w1;
    uint64_t test_wire_and = modMult(w1, w2);
    uint64_t test_wire_xor = mod(w1 + w2);
    uint64_t test_wire_inv = mod(1+PRIME - w1);

    uint64_t fz = 0;
    for(int k = 0; k < ni; ++k){
      if(c[lvl][k].type == -1){
        fz = mod(fz + modMult(beta, modMult(c[lvl][k].wire[tmp_index], test_wire_passthru)));
      }
      if(c[lvl][k].type == 0){
        fz = mod(fz + modMult(beta, modMult(c[lvl][k].wire[tmp_index], test_wire_and)));
      }
      else if(c[lvl][k].type == 1){
        fz = mod(fz + modMult(beta, modMult(c[lvl][k].wire[tmp_index], test_wire_xor)));
      }
      if(c[lvl][k].type == 2){
        fz = mod(fz + modMult(beta, modMult(c[lvl][k].wire[tmp_index], test_wire_inv)));
      }
    }
    check = extrap(poly[num_var-1], 3, r[num_var-1]);
    if((check != fz) && (check != fz + PRIME) && (check != fz-PRIME)){
      cout << "PROTOCOL FAILED: final check is incorrect" << endl;
    }

    vector<uint64_t> lpoly;
    vector<uint64_t> point(mip1);
    vector<uint64_t> vec(2);

    for(int k = 0; k < mip1+1; k++){
      for(int j = 0; j < mip1; j++){
        vec[0]=r[mi+j];
        vec[1]=r[mi+mip1+j];
        point[j] = extrap(vec, 2, k);
      }
      lpoly.push_back(w_i(mip1, nip1, lvl-1, point, 0)); 
    }

    if((w1 != lpoly[0]) && (w1 != lpoly[0] + PRIME) && (w1 != lpoly[0]-PRIME)){
      cout << "PROTOCOL FAILED: new claim is incorrect" << endl;
    }
    if((w2 != lpoly[1]) && (w2 != lpoly[1] + PRIME) && (w2 != lpoly[1]-PRIME)){
      cout << "PROTOCOL FAILED: new claim is incorrect" << endl;
    }

    uint64_t t = mod(0);
    for(int j = 0; j < mip1; j++)
    {
      vec[0]=r[mi+j];
      vec[1]=r[mi+mip1+j];
      z[j] = extrap(vec, 2, t);
    }

    uint64_t ret = extrap(lpoly, mip1+1, t);
    return {ret, {0, 0}};
  }
  return {mod(poly[0][0]+poly[0][1]), {w1, w2}};
}

typedef struct info{
  int root;
  vector<int> lvls;
} info;

const int NUM_GATES = 116246;
const int NUM_WIRES = 116758;
int WIDTH = 512;
const int DEPTH = 5383;
info memo[NUM_WIRES]; //key: id; value: lvl, id

int get_lvl(info &i, int lvl){
  return i.lvls[lvl - i.root];
}

void op(int a_id, int b_id, int x_id, int op){
  auto a = memo[a_id], b = memo[b_id];
  if(c.size() == max(a.root, b.root)+1){
    c.push_back(vector<gate>());
  }
  if(max(a.root, b.root) > min(a.root + a.lvls.size()-1, b.root + b.lvls.size()-1)){
    if(a.root > b.root){
      for(int i = b.root + b.lvls.size(); i <= a.root; ++i){
        c[i].push_back({-1, false, b.lvls.back()});
        b.lvls.push_back(c[i].size()-1);
      }
    }
    else{
      for(int i = a.root + a.lvls.size(); i <= b.root; ++i){
        c[i].push_back({-1, false, a.lvls.back()});
        a.lvls.push_back(c[i].size()-1);
      }
    }
  }
  c[max(a.root, b.root)+1].push_back({op, false, get_lvl(a, max(a.root, b.root)), get_lvl(b, max(a.root, b.root))});
  memo[x_id] = {max(a.root, b.root)+1, vector<int>(1, c[max(a.root, b.root)+1].size()-1)};
}

void create_circuit(){
  ifstream fin("raw.circuit");
  int tmp;
  fin >> tmp >> tmp >> tmp >> tmp >> tmp;

  vector<gate> input_level;
  for(int i = 0; i < WIDTH; ++i){
    gate g;
    g.val = rand()%2;
    input_level.push_back(g);
    memo[i] = {0, vector<int>(1, i)};
  }
  c.push_back(input_level);

  int num_in, num_out, a, b, x;
  string type;
  for(int i = 0; i < NUM_WIRES; ++i){
    fin >> num_in >> num_out;
    if(num_in == 1){
      fin >> a;
    }
    else{
      fin >> a >> b;
    }
    fin >> x >> type;

    if(type == "AND"){
      op(a, b, x, 0);
    }
    else if(type == "XOR"){
      op(a, b, x, 0);
    }
    else if(type == "INV"){
      op(a, a, x, 2);
    }
  }
}

void add_zeroed(){
  int n = log2_max+1;
  cout << n << endl;
  vector<uint64_t> zeroed(n, 0);
  for(int i = 0; i < c.size(); ++i){
    for(int j = 0; j < c[i].size(); ++j){
      c[i][j].wire = zeroed;
      c[i][j].beta = zeroed;
    }
  }
}

//evaluates the circuit
void eval_circuit(){
  for(int i = 1; i < c.size(); ++i){
    log2_max = max(log2_max, FAST_LOG2_UP(c[i].size()));
    for(int j = 0; j < c[i].size(); ++j){
      if(c[i][j].type == -1){
        c[i][j].val = c[i-1][c[i][j].in1].val;
      }
      else if(c[i][j].type == 0){
        c[i][j].val = c[i-1][c[i][j].in1].val & c[i-1][c[i][j].in2].val;
      }
      else if(c[i][j].type == 1){
        c[i][j].val = c[i-1][c[i][j].in1].val ^ c[i-1][c[i][j].in2].val;
      }
      else if(c[i][j].type == 2){
        c[i][j].val = 1 - c[i-1][c[i][j].in1].val;
      }
    }
  }
}

//the f polynomial, {a, b} => f_i(x) = ax+b
vector<pair<uint64_t, uint64_t>> f;
//interpolation points for queries
vector<vector<pair<uint64_t, vector<uint64_t>>>> query_interp_points;
//interpolation points for challenges
vector<vector<pair<pair<uint64_t, uint64_t>, pair<vector<uint64_t>, vector<uint64_t>>>>> challenge_interp_points;

//interpolate the query points
vector<uint64_t> univar_interp_query(vector<pair<uint64_t, vector<uint64_t>>> &points){
  vector<uint64_t> coeffs(points.size());
  vector<uint64_t> basis(points.size());
  basis[0] = 1;
  vector<uint64_t> ddif;
  for(auto &point : points){
    ddif.push_back(point.second.back());
  }

  #pragma omp parallel for
  for(int i = 0; i < points.size(); ++i){
    for(int j = 0; j < points.size(); ++j){
      coeffs[j] = mod(coeffs[j] + modMult(ddif[0], basis[j]));
    }
    if(i < points.size()-1){
      uint64_t neg_x = PRIME - points[i].first;
      for(int j = i+1; j > 0; --j){
        basis[j] = mod(basis[j-1] + modMult(neg_x, basis[j]));
      }
      basis[0] = modMult(basis[0], neg_x);
      for(int j = 0; j+i+1 < points.size()+1; ++j){
        uint64_t num = mod(ddif[j+1]+PRIME - ddif[j]);
        uint64_t den = mod(points[j+i+1].first+PRIME - points[j].first);
        ddif[j] = modMult(num, inv(den));
      }
    }
  }
  return coeffs;
}

//interpolate the challenge points
vector<uint64_t> univar_interp_challenge(vector<pair<pair<uint64_t, uint64_t>, pair<vector<uint64_t>, vector<uint64_t>>>> &points){
  vector<uint64_t> coeffs(points.size());
  vector<uint64_t> basis(points.size());
  basis[0] = 1;
  vector<uint64_t> ddif;
  for(auto &point : points){
    ddif.push_back(point.second.first.back());
  }

  #pragma omp parallel for
  for(int i = 0; i < points.size(); ++i){
    for(int j = 0; j < points.size(); ++j){
      coeffs[j] = mod(coeffs[j] + modMult(ddif[0], basis[j]));
    }
    if(i < points.size()-1){
      uint64_t neg_x = PRIME - points[i].first.first;
      for(int j = i+1; j > 0; --j){
        basis[j] = mod(basis[j-1] + modMult(neg_x, basis[j]));
      }
      basis[0] = modMult(basis[0], neg_x);
      for(int j = 0; j+i+1 < points.size()+1; ++j){
        uint64_t num = mod(ddif[j+1]+PRIME - ddif[j]);
        uint64_t den = mod(points[j+i+1].first.first+PRIME - points[j].first.first);
        ddif[j] = modMult(num, inv(den));
      }
    }
  }
  return coeffs;
} 

double LATENCY;

int main(){
  cout << "INPUT: LATENCY(MS) HYDRA" << endl;
  cin >> LATENCY >> HYDRA;

  if(HYDRA){
    omp_set_max_active_levels(2);
    omp_set_num_threads(omp_get_max_threads());
  }

  cout << "creating circuit..." << endl;
  auto start = chrono::steady_clock::now();
  create_circuit();
  auto end = chrono::steady_clock::now();
  cout << "creation time: " << chrono::duration <double, milli> (end-start).count()/1000 << " s\n" << endl;

  cout << "evaluating circuit..." << endl;
  start = chrono::steady_clock::now();
  eval_circuit();
  add_zeroed();
  end = chrono::steady_clock::now();
  cout << "done, output: " << c.back()[0].val << endl;
  cout << "evaluation time: " << chrono::duration <double, milli> (end-start).count()/1000 << " s\n" << endl;

  if(!HYDRA) cout << "proving/verifying..." << endl;
  else cout << "generating interpolation points (offline)..." << endl;
  start = chrono::steady_clock::now();

  if(!HYDRA){
    uint64_t ri;
    vector<uint64_t> z; 
    for(int i = 0; i < log2_max; ++i){
      z.push_back(rand());
    }
    vector<uint64_t> empty_vec;
    pair<vector<uint64_t>, vector<uint64_t>> empty_pair = {empty_vec, empty_vec};
    for(int i = c.size()-1; i > 0; --i){
      ri = sumcheck(i, z, ri, empty_pair, 0).first;
    }
  }
  else {
    //offline interpolation point generation
    for(int i = 0; i < log2_max; ++i){
      f.push_back({rand(), rand()});
    }
    for(int i = 0; i < c.size(); ++i){
      vector<pair<uint64_t, vector<uint64_t>>> query_layer_points;
      for(int j = 0; j < log2_max+1; ++j){
        vector<uint64_t> point;
        uint64_t x = rand();
        for(int k = 0; k < log2_max; ++k){
          point.push_back(mod(modMult(f[k].first, x) + f[k].second));
        }
        query_layer_points.push_back({x, point});
      }
      query_interp_points.push_back(query_layer_points);

      vector<pair<pair<uint64_t, uint64_t>, pair<vector<uint64_t>, vector<uint64_t>>>> challenge_layer_points;
      for(int j = 0; j < log2_max+1; ++j){
        vector<uint64_t> point1, point2;
        uint64_t x1 = rand(), x2 = rand();
        for(int k = 0; k < log2_max; ++k){
          point1.push_back(mod(modMult(f[k].first, x1) + f[k].second));
          point2.push_back(mod(modMult(f[k].first, x2) + f[k].second));
        }
        challenge_layer_points.push_back({{x1, x2}, {point1, point2}});
      }
      challenge_interp_points.push_back(challenge_layer_points);
    }
    end = chrono::steady_clock::now();
    cout << "generation time: " << chrono::duration <double, milli> (end-start).count()/1000 << " s\n" << endl;

    cout << "proving/verifying..." << endl;
    start = chrono::steady_clock::now();
    
    //first commit all w_i polynomials (proofs are opened and verified here,
    //should be after the protocol but in terms of timing it does not make a
    //difference)
    #pragma omp parallel for
    for (int i = c.size()-1; i > 0; --i){
      string cmd = "cd ../../deps/poly-commit/target/release && ./hydra-cmd ";
      cmd += to_string(c[i].size()) + " ";
      for(int j = 0; j < c[i].size(); ++j){
        cmd += to_string(c[i][j].val) + " ";
      }
      cmd += to_string(log2_max+1) + " ";
      cmd += to_string(log2_max) + " ";
      for(int j = 0; j < log2_max+1; ++j) {
        for(int k = 0; k < log2_max; ++k){
          cmd += to_string(query_interp_points[i][j].second[k]) + " ";
        }
      }
      int err = system(cmd.c_str());
      if(err){
        cout << "PROTOCOL FAILED: polynomial commitment error" << endl;
      }
    }

    uint64_t ri;
    vector<uint64_t> z; 
  
    //engage in parallel interactive sumcheck with GKR
    //#pragma omp parallel for
    for(int i = c.size()-1; i > 0; --i){
      //#pragma omp parallel for
      for(int j = 0; j < log2_max+1; ++j){
        cout << i << " " << j << endl;
        pair<vector<uint64_t>, vector<uint64_t>> challenge_pair = {challenge_interp_points[i][j].second.first, challenge_interp_points[i][j].second.second};
        auto ret = sumcheck(i, query_interp_points[i][j].second, 0, challenge_pair, j);
        query_interp_points[i][j].second.push_back(ret.first);
        challenge_interp_points[i][j].second.first.push_back(ret.second.first);
        challenge_interp_points[i][j].second.second.push_back(ret.second.second);
      }
    }
    
    //consistency check
    #pragma omp parallel for
    for(int i = c.size()-1; i > 1; --i){
      auto poly_coeffs_query = univar_interp_query(query_interp_points[i-1]);
      auto poly_coeffs_challenge = univar_interp_challenge(challenge_interp_points[i]);
      if(!(poly_coeffs_query == poly_coeffs_challenge)){
        cout << "PROTOCOL FAILED: query/challenge polynomial interpolation error" << endl;
      }
    }
  }
  
  end = chrono::steady_clock::now();
  cout << "done." << endl;

  auto time = chrono::duration <double, milli> (end-start).count()/1000;

  //we add the artificial latency here
  if(!HYDRA){
    //note that the traditional GKR latency does in fact depend on the depth
    //since it is not parallelized, instead it is a sequential process
    time += double(DEPTH) * 2.0 * double(log2_max) * LATENCY/1000;
  }
  else{
    //note that hydra latency does not depend on the depth since it is
    //parallelized (batched messages all at once), so the latency is only a
    //factor of a single sumcheck protocol instance
    time += double(log2_max) * 2.0 * LATENCY/1000;
  }
  cout << fixed;
  cout << "protocol time (with latency): " << time << " s\n" << endl;
}
