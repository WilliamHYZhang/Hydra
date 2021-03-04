/*
William Zhang

Last Updated: March 2, 2021

Implementation of the subcircuit protocol described in "Hydra: Succinct Fully
Pipelineable Interactive Arguments of Knowledge" by Zhang and Xia.

Polynomial and arithmetic logic based on Thaler's implementation described in
"Practical Verified Computation with Streaming Interactive Proofs" by Cormode,
Mitzenmacher, and Thaler.

This version is specific to testing with the SHA-256 boolean circuit.
*/

#include <bits/stdc++.h>

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
  uint64_t wire; //scalar wire contribution
  uint64_t beta; //scalar beta contribution
} gate;

//circuit represented as 2d vector of gates
typedef vector<vector<gate>> circ;
circ c;

//subcircuits represented as a vector of circuits
vector<circ> subcirc;

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
uint64_t w_i(int mi, int ni, int i, vector<uint64_t> &r, int r_i, int b){
  uint64_t ret = 0;
  for(uint64_t k = 0; k < ni; ++k){
    ret = mod(ret + modMult(subcirc[b][i][k].val, interp(k, r, mi, r_i)));
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

//main sumcheck protocol
uint64_t sumcheck(int lvl, vector<uint64_t> &z, uint64_t ri, int b){
  int mi = FAST_LOG2_UP(subcirc[b][lvl].size());
  int ni = subcirc[b][lvl].size();
  int mip1 = FAST_LOG2_UP(subcirc[b][lvl-1].size());
  int nip1 = subcirc[b][lvl-1].size();

  int num_var = mi + 2 * mip1;
  int num_effective = modPow(2, mip1);
  uint64_t alpha_p, beta_p;

  vector<uint64_t> r;
  for(int i = 0; i < num_var; ++i){
    r.push_back(rand());
  }

  for(uint64_t k = 0; k < ni; ++k){
    subcirc[b][lvl][k].beta = beta(mi, z, k);
    subcirc[b][lvl][k].wire = 1;
  }

  vector<vector<uint64_t>> poly;
  for(int k = 0; k < num_var; ++k){
    poly.push_back({0, 0, 0});
  }

  vector<uint64_t> val;
  for(int k = 0; k < num_effective; ++k){
    val.push_back(k < nip1 ? subcirc[b][lvl-1][k].val : 0);
  }

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
        val[k] = k < nip1 ? subcirc[b][lvl-1][k].val : 0;
      }
      w1 = w_i(mip1, nip1, lvl-1, r, mi, b);
    }
    if(j == mi){
      beta = subcirc[b][lvl][0].beta;
    }
    for(int k = 0; k < ni; ++k){
      if(j < mi){
        w1 = subcirc[b][lvl-1][subcirc[b][lvl][k].in1].val;
        w2 = subcirc[b][lvl-1][subcirc[b][lvl][k].in2].val;
      }
      else if(j < mi + mip1){
        w2 = subcirc[b][lvl-1][subcirc[b][lvl][k].in2].val;
      }
      if(j < mi){
        k_j = k >> j;
        if(k_j & 1){
          subcirc[b][lvl][k].beta = modMult(subcirc[b][lvl][k].beta, inv(z[j]));
        }
        else{
          subcirc[b][lvl][k].beta = modMult(subcirc[b][lvl][k].beta, inv(1+PRIME - z[j]));
        }
      }

      for(int m = 0; m <= 2; ++m){
        if(j < mi){
          if(m == 0){
            beta = modMult(subcirc[b][lvl][k].beta, modMult(1+PRIME - z[j], 1+PRIME - m));
          }
          else if(m == 1){
            beta = modMult(subcirc[b][lvl][k].beta, modMult(z[j], m));
          }
          else{
            beta = modMult(subcirc[b][lvl][k].beta, mod(modMult(1+PRIME - z[j], 1+PRIME - m) + modMult(z[j], m)));
          }
          if(m == 2){
            subcirc[b][lvl][k].beta = modMult(subcirc[b][lvl][k].beta, mod(modMult(1+PRIME - z[j], 1+PRIME - r[j]) + modMult(z[j], r[j])));
          }
          if(k_j & 1){
            if(m == 0){
              continue;
            }
            else{
              wire = modMult(subcirc[b][lvl][k].wire, m);
            }
            if(m == 2){
              subcirc[b][lvl][k].wire = modMult(subcirc[b][lvl][k].wire, r[j]);
            }
          }
          else{
            if(m==1){
              continue;
            }
            else{
              wire = modMult(subcirc[b][lvl][k].wire, 1+PRIME - m);
            }
            if(m == 2){
              subcirc[b][lvl][k].wire = modMult(subcirc[b][lvl][k].wire, 1+PRIME - r[j]);
            }
          }
        }
        if(j >= mi && j < mi+mip1){
          index = subcirc[b][lvl][k].in1 >> (j - mi);
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
            wire = modMult(subcirc[b][lvl][k].wire, m);
            if(m == 2){
              subcirc[b][lvl][k].wire = modMult(subcirc[b][lvl][k].wire, r[j]);
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
            wire = modMult(subcirc[b][lvl][k].wire, 1+PRIME - m);
            if(m == 2){
              subcirc[b][lvl][k].wire = modMult(subcirc[b][lvl][k].wire, 1+PRIME - r[j]);
            }
          }
        }
        else if(j >= mi+mip1){
          index = subcirc[b][lvl][k].in2 >> (j - mi - mip1);
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
            wire = modMult(subcirc[b][lvl][k].wire, m);
            if(m == 2){
              subcirc[b][lvl][k].wire = modMult(subcirc[b][lvl][k].wire, r[j]);
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
            wire = modMult(subcirc[b][lvl][k].wire, 1+PRIME - m);
            if(m == 2){
              subcirc[b][lvl][k].wire = modMult(subcirc[b][lvl][k].wire, 1+PRIME - r[j]);
            }
          }
        }

        if(subcirc[b][lvl][k].type == -1){
          poly[j][m] = mod(poly[j][m] + mod(modMult(beta, modMult(wire, w1))));
        }
        else if(subcirc[b][lvl][k].type == 0){
          poly[j][m] = mod(poly[j][m] + mod(modMult(beta, modMult(wire, modMult(w1, w2)))));
        }
        else if(subcirc[b][lvl][k].type == 1){
          poly[j][m] = mod(poly[j][m] + mod(modMult(beta, modMult(wire, mod(w1 + w2)))));
        }
        else if(subcirc[b][lvl][k].type == 2){
          poly[j][m] = mod(poly[j][m] + mod(modMult(beta, modMult(wire, mod(1+PRIME - w1)))));
        }
       
      }
    }
  }
  if((lvl != subcirc[b].size()-1) && (mod(poly[0][0]+poly[0][1]) != ri) && (mod(poly[0][0]+poly[0][1]) != ri-PRIME) && (mod(poly[0][0]+poly[0][1]) != ri+PRIME)){
    cout << "PROTOCOL FAILED: claimed value is incorrect" << endl;
  }

  uint64_t check = 0;
  for(int j = 1; j < num_var; j++) {
    check = extrap(poly[j-1], 3, r[j-1]);
    if((check != mod(poly[j][0]+poly[j][1])) && (check != (mod(poly[j][0]+poly[j][1])+PRIME)) && (check != mod(poly[j][0]+poly[j][1])-PRIME)){
      cout << "PROTOCOL FAILED: intermediate sumcheck message is incorrect" << " b " << b << endl;
    }
  }
  w2 = w_i(mip1, nip1, lvl-1, r, mi+mip1, b);

    uint64_t test_wire_passthru = w1;
    uint64_t test_wire_and = modMult(w1, w2);
    uint64_t test_wire_xor = mod(w1 + w2);
    uint64_t test_wire_inv = mod(1+PRIME - w1);

    uint64_t fz = 0;
    for(int k = 0; k < ni; ++k){
      if(subcirc[b][lvl][k].type == -1){
        fz = mod(fz + modMult(beta, modMult(subcirc[b][lvl][k].wire, test_wire_passthru)));
      }
      if(subcirc[b][lvl][k].type == 0){
        fz = mod(fz + modMult(beta, modMult(subcirc[b][lvl][k].wire, test_wire_and)));
      }
      else if(subcirc[b][lvl][k].type == 1){
        fz = mod(fz + modMult(beta, modMult(subcirc[b][lvl][k].wire, test_wire_xor)));
      }
      if(subcirc[b][lvl][k].type == 2){
        fz = mod(fz + modMult(beta, modMult(subcirc[b][lvl][k].wire, test_wire_inv)));
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
    lpoly.push_back(w_i(mip1, nip1, lvl-1, point, 0, b)); 
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
  return ret;
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

//creates the circuit
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

int SUB_DEPTH;

//pipeline the circuit by splitting it up into subcircuits
void pipeline(){
  for(int i = 0; i < c.size(); i += SUB_DEPTH){
    vector<vector<gate>> sub(c.begin()+(i==0 ? 0 : i-1), i+SUB_DEPTH >= c.size() ? c.end() : c.begin()+i+SUB_DEPTH);
    subcirc.push_back(sub);
  }
}

bool SUBCIRCUIT;

int main(){
  cout << "INPUT: SUB_DEPTH SUBCIRCUIT" << endl;
  cin >> SUB_DEPTH >> SUBCIRCUIT;

  cout << "creating circuit..." << endl;
  auto start = chrono::steady_clock::now();
  create_circuit();
  auto end = chrono::steady_clock::now();
  cout << "creation time: " << chrono::duration <double, milli> (end-start).count()/1000 << " s\n" << endl;

  cout << "evaluating circuit..." << endl;
  start = chrono::steady_clock::now();
  eval_circuit();
  end = chrono::steady_clock::now();
  cout << "done, output: " << c.back()[0].val << endl;
  cout << "evaluation time: " << chrono::duration <double, milli> (end-start).count()/1000 << " s\n" << endl;

  cout << "pipelining circuit..." << endl;
  start = chrono::steady_clock::now();
  pipeline();
  end = chrono::steady_clock::now();
  cout << "pipeline time: " << chrono::duration <double, milli> (end-start).count()/1000 << " s\n" << endl;


  cout << "proving/verifying..." << endl;
  start = chrono::steady_clock::now();

  if(!SUBCIRCUIT){
    uint64_t ri;
    vector<uint64_t> z; 
    for(int i = 0; i < log2_max; ++i){
      z.push_back(rand());
    }
    for(int b = 0; b < subcirc.size(); ++b){
      for(int i = subcirc[b].size()-1; i > 0; --i){
        ri = sumcheck(i, z, ri, b);
      }
    }
  }
  else {
    //parallelized circuit proving
    #pragma omp parallel for
    for(int b = 0; b < subcirc.size(); ++b){
      uint64_t ri = 0;
      vector<uint64_t> z;
      for(int i = 0; i < log2_max; ++i){
        z.push_back(rand());
      }
      for(int i = subcirc[b].size()-1; i > 0; --i){
        ri = sumcheck(i, z, ri, b);
      }
      
      //commit w_in/w_out
      string cmd = "cd ../../deps/poly-commit/target/release && ./hydra-cmd ";
      cmd += to_string(subcirc[b][0].size()) + " ";
      for(int i = 0; i < subcirc[b][0].size(); ++i){
        cmd += to_string(subcirc[b][0][i].val) + " ";
      }
      cmd += "1 ";
      cmd += to_string(log2_max) + " ";
      for(int i = 0; i < log2_max; ++i){
        cmd += to_string(z[i]) + " ";
      }
      int err = system(cmd.c_str());
      if(err){
        cout << "PROTOCOL FAILED: polynomial commitment error" << endl;
      }
    }
  }
  
  end = chrono::steady_clock::now();
  cout << "done." << endl;
  cout << "protocol time: " << chrono::duration <double, milli> (end-start).count()/1000 << " s\n" << endl;
}
