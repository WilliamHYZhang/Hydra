#include <bits/stdc++.h>
using namespace std;

typedef struct info{
  int root;
  vector<int> lvls;
} info;

typedef struct gate{
  int type; // 0 NAND, -1 PASSTHRU, 1 CONSTANT
  bool val;
  int in1;
  int in2;
} gate;

vector<vector<gate>> c;

unordered_map<int, info> memo; //key: id; value: lvl, id

int get_lvl(info &i, int lvl){
  return i.lvls[lvl - i.root];
}

info nand(info &a, info &b){
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
  c[max(a.root, b.root)+1].push_back({0, false, get_lvl(a, max(a.root, b.root)), get_lvl(b, max(a.root, b.root))});
  info ret = {max(a.root, b.root)+1, vector<int>(1, c[max(a.root, b.root)+1].size()-1)};
  return ret;
}

info g1, g2, g3, g4; //intermediary gates

void add_gate(int op, int id, int id1 = 0, int id2 = 0){
  if(op != 0 && op != 15){
    memo[id1] = memo[id1], memo[id2] = memo[id2];
  }
  if(op == 0){ //F
    if(!c.size()){
      c.push_back(vector<gate>());
    }
    gate c_g = {1, false};
    c[0].push_back(c_g);
    memo[id] = {0, vector<int>(1, c[0].size()-1)};
  }
  else if(op == 2){ //(NOT B) AND A, ((B NAND B) NAND A) NAND ((B NAND B) NAND A)
    g1 = nand(memo[id2], memo[id2]); //B NAND B
    g2 = nand(g1, memo[id1]); //(B NAND B) NAND A
    memo[id] = nand(g2, g2); //((B NAND B) NAND A) NAND ((B NAND B) NAND A)
  }
  else if(op == 4){ //(NOT A) AND B, ((A NAND A) NAND B) NAND ((A NAND A) NAND B)
    g1 = nand(memo[id1], memo[id1]); //B NAND B
    g2 = nand(g1, memo[id2]); //(B NAND B) NAND A
    memo[id] = nand(g2, g2); //((B NAND B) NAND A) NAND ((B NAND B) NAND A)    
  }
  else if(op == 6){ //XOR, (A NAND (A NAND B)) NAND (B NAND (A NAND B))
    g1 = nand(memo[id1], memo[id2]); //A NAND B
    g2 = nand(memo[id1], g1); //A NAND (A NAND B)
    g3 = nand(memo[id2], g1); //B NAND (A NAND B)
    memo[id] = nand(g2, g3); //(A NAND (A NAND B)) NAND (B NAND (A NAND B))
  }
  else if(op == 7){ //NAND
    memo[id] = nand(memo[id1], memo[id2]);
  }
  else if(op == 8){ //AND, (A NAND B) NAND (A NAND B)
    g1 = nand(memo[id1], memo[id2]); //A NAND B
    memo[id] = nand(g1, g1); //(A NAND B) NAND (A NAND B)
  }
  else if(op == 9){ //XNOR ((A NAND A) NAND (B NAND B)) NAND (A NAND B)
    g1 = nand(memo[id1], memo[id1]); //A NAND A
    g2 = nand(memo[id2], memo[id2]); //B NAND B
    g3 = nand(g1, g2); //((A NAND A) NAND (B NAND B))
    g4 = nand(memo[id1], memo[id2]); //A NAND B
    memo[id] = nand(g3, g4); //((A NAND A) NAND (B NAND B)) NAND (A NAND B)
  }
  else if(op == 15){ //T
    if(!c.size()){
      c.push_back(vector<gate>());
    }
    gate c_g = {1, true};
    c[0].push_back(c_g);
    memo[id] = {0, vector<int>(1, c[0].size()-1)};    
  }
}


ifstream circuit_fin("../master/src/out"); //circuit output file from frigate compiler
//ifstream circuit_fin("sample"); //small test circuit file
string str, str_op;
int op, g, in1, in2;

int cnt = 0, cur = 1000000;

int main(){
  auto start = chrono::steady_clock::now();
  while(getline(circuit_fin, str)){
    istringstream ss(str);
    ss >> str_op >> g;
    if(str_op == "IN"){
      op = rand() % 2 * 15; //random T/F for input now
    }
    else if(str_op == "copy(6)"){
      op = 6;
    }
    else if(str_op == "OUT"){
      continue; //last layer of circuit already known to be output
    }
    else{
      op = stoi(str_op);
    }
    if(op == 0 || op == 15){
      add_gate(op, g);
    }
    else{
      ss >> in1 >> in2;
      add_gate(op, g, in1, in2);
    }

    if(cnt % cur == 0){
      cout << "circuit size: " << c.size() << " lines: " << cnt << endl;
    }
    ++cnt;
  }
  auto end = chrono::steady_clock::now();
  cout << "file io and creation time: " << chrono::duration <double, milli> (end-start).count()/1000 << " s" << endl;
}
