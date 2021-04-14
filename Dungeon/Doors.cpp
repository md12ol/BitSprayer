//
// Created by micha on 2021-04-09.
//

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <queue>
#include <climits>
#include <algorithm>
#include <random>

using namespace std;

#include "../BitSprayer/stat.h"
#include "../BitSprayer/bitspray.h"
#include "setu.h"
#define RNS 91207819

#define runs 30
#define mevs 2000
#define RI 20
#define popsize 100
#define tsize 7
#define MNM 3
#define verbose 1

#define states 12
//#define alp 2
//#define alt 2
#define Qz 2000
int Q[Qz];

#define BB 10 //Grid size

#define Rz 256
#define reqDiam 20
#define diamTests 7

void initalg();
void initpop();
vector<pair<int, int>> BFS(graph &G, int source, int *maxs);
int fitness(bitspray &A);
int diameter(graph &G);
void report(ostream &aus);
void matingevent();
bool getbit(int &val, int bits, int &psn);
void randomOrder(int size);
void fillEdges();
int approxDiameter(graph &G);
bool getbit(int &val, int &psn);
int validation(bitspray &A);
double reportbest(ostream &aus, int run);

bitspray pop[popsize];  //Population of bitsprayers
int fit[popsize];  //Fitness values
int dx[popsize];  //Sorting index

graph *initG;
graph *curG;
graph *D;
int doors[1000][2];
int posDoors;
vector<int> order;
bool over = false;

int main() {
  fstream stat;   //statistics reporting stream
  fstream best;   //best reporting stream
  fstream initGraph;
  fstream outGraph;
  fstream doorGraph;
  int numDoors;
  char fn[100];
  int g = 0;
  sprintf(fn, "Input/Graphs/graph%02d.dat", g);
  initGraph.open(fn, ios::in);
  initG = new graph(Rz);
  initG->read(initGraph);
  initGraph.close();

  sprintf(fn, "Input/Doors/doors%02d.dat", g);
  doorGraph.open(fn, ios::in);
  D = new graph(Rz);
  D->read(doorGraph);
  doorGraph.close();
  posDoors = D->edges();
  fillEdges();

  curG = new graph(Rz);

  order.reserve(posDoors);
  randomOrder(posDoors);
  best.open("Output/best.sda", ios::out);

//  sprintf(fn, "Output/newg%02d.dat", g);
//  outGraph.open(fn, ios::out);
  initalg();
  for (int run = 0; run < runs; run++) {
    cout << "Run=" << run << endl;
    sprintf(fn, "Output/run%02d.dat", run);
    stat.open(fn, ios::out);
    initpop();
    report(stat);
    for (int mev = 0; mev < mevs; mev++) {
      matingevent();
      if ((mev + 1) % RI == 0) {
        if (verbose == 1) {
          cout << run << " " << (mev + 1) / RI << " ";
        }
        report(stat);
      }
    }
    reportbest(best, run);
    stat.close();
  }
//  best.close();
  outGraph.close();
  best.close();

  delete initG;
  delete D;
  delete curG;
  return 0;
}

void randomOrder(int size) {
  for (int i = 0; i < size; i++) {
    order.push_back(i);
  }
  shuffle(order.begin(), order.end(), default_random_engine(RNS));
//  for (int i = 0; i < size; i++) {
//    cout << order[i] << " ";
//  }
//  cout << endl;
}

void initalg() {
  srand48(RNS);
  for (int i = 0; i < popsize; i++) {
    pop[i].create(states);
  }
}

void initpop() {
  cout << "Initial population fitness values" << std::endl;
  for (int i = 0; i < popsize; i++) {
    pop[i].randomize();
    fit[i] = fitness(pop[i]);
    cout << fit[i] << " ";
    dx[i] = i;
  }
  cout << endl;
}

void developQ(bitspray &A) {//unpack the queue
  int h, t;  //head and tail of queue

  over = false;
  for (t = 0; t < Qz; t++)Q[t] = 0;        //clear the queue
  A.reset(Q, h, t);                 //reset the self driving automata
  while (t < Qz - 2)A.next(Q, h, t, Qz);  //run the automata
}

bool necroticFilter() {
  int count = 0;
  for (int i = 0; i < Qz; i++) {
    if (Q[i] == 1) {
      count++;
    }
  }
  double ratio = (double) count / (double) Qz;
  if (ratio < 0.5 || ratio > 0.9) {
    return true;
  } else {
    return false;
  }
}

int fitness(bitspray &A) {
  int addedDoors = 0;
  int acpt = 0;
  int psn = 0;
  int idx = 0;
  curG->empty(Rz);
  curG->copy(*initG);
  developQ(A);
  if (necroticFilter()) {
    return posDoors;
  }

  while (approxDiameter(*curG) > reqDiam && getbit(acpt, psn)) {
    if (acpt == 1) {
      if (curG->add(doors[order[idx]][0], doors[order[idx]][1])) {
        addedDoors++;
      }
    }
    idx = (idx + 1) % posDoors;
  }

  if (over) {
    return posDoors;
  } else {
    return addedDoors;
  }
}

int validation(bitspray &A) {
  int addedDoors = 0;
  int acpt = 0;
  int psn = 0;
  int idx = 0;
  curG->empty(Rz);
  curG->copy(*initG);
  developQ(A);

  while (diameter(*curG) > reqDiam && getbit(acpt, psn)) {
    if (acpt == 1) {
      if (curG->add(doors[order[idx]][0], doors[order[idx]][1])) {
        addedDoors++;
      }
    }
    idx = (idx + 1) % posDoors;
  }

  if (over) {
    return -1;
  } else {
    return addedDoors;
  }
}

bool getbit(int &val, int &psn) {//get a number from the queue
  if (psn + 1 >= Qz) {
    over = true;
    return false; //safety first
  }
  val = Q[psn++];  //zero the value
//  for (int i = 0; i < bits; i++)
//    val = 2 * val + Q[psn++];  //assemble the integer
  return true;
}

vector<pair<int, int>> BFS(graph &G, int source) {
  int dist[G.size()];
  int path[G.size()];
  queue<int> qu;
  int node, numnbr;

  for (int i = 0; i < G.size(); i++) {
    dist[i] = INT_MAX;
    path[i] = -1;
  }
  dist[source] = 0;

  qu.push(source);
  int nbrs[G.size()];
  while (!qu.empty()) {
    node = qu.front();
    qu.pop();
    numnbr = G.Nbrs(node, nbrs);
    if (numnbr > 0) {
      for (int i = 0; i < numnbr; i++) {
        if (dist[nbrs[i]] > dist[node] + 1) {
          dist[nbrs[i]] = dist[node] + 1;
          path[nbrs[i]] = node;
          qu.push(nbrs[i]);
        }
      }
    }
  }

  vector<pair<int, int>> so;
  for (int i = 0; i < G.size(); i++) {
    so.push_back(make_pair(dist[i], i));
  }
  sort(so.begin(), so.end());
  return so;
}

int diameter(graph &G) {
  int maxs[G.size()];
  for (int i = 0; i < G.size(); i++) {
    maxs[i] = BFS(G, i).back().first;
  }
  int max = 0;
  for (int i = 0; i < G.size(); i++) {
    if (maxs[i] > max) {
      max = maxs[i];
    }
  }
  return max;
}

int approxDiameter(graph &G) {
  int maxs[diamTests];
  int from = 0;
  int nextF;
  vector<int> tested;
  vector<pair<int, int>> so;
  for (int i = 0; i < diamTests; i++) {
    maxs[i] = 0;
    tested.push_back(from);
    so = BFS(G, from);
    maxs[i] = so.back().first;
    nextF = so.size() - 1;
    do {
      from = so.at(nextF).second;
      nextF--;
    } while (*find(tested.begin(), tested.end(), from) == from);
  }
  int max = 0;
  for (int i = 0; i < diamTests; i++) {
    if (maxs[i] > max) {
      max = maxs[i];
    }
  }
  return max;
}

void report(ostream &aus) {
  dset D;
  D.add(fit, popsize);
  aus << D.Rmu() << " " << D.RCI95() << " " << D.Rsg() << " " << D.Rmin();
  aus << endl;
  if (verbose == 1) {
    cout << D.Rmu() << " " << D.RCI95() << " " << D.Rsg() << " " << D.Rmin();
    cout << endl;
  }
}

double reportbest(ostream &aus, int run) {
  int b = 0;
  for (int i = 1; i < popsize; i++) {
    if (fit[i] < fit[b]) b = i;
  }
  cout << "Best Approx Fitness: " << fit[b] << endl;
  fit[b] = validation(pop[b]);
  aus << "Run" << run << " fitness: " << fit[b] << endl;
  cout << "Best Fitness: " << fit[b] << endl;
//  render(run);
//  printGraph(run);
//  printBoard(run);
//  printDoors(run);
//  pop[b]->print(aus);
  return fit[b];
}

void matingevent() {
  int nm;

  Tselect(fit, dx, tsize, popsize);
  pop[dx[0]].copy(pop[dx[tsize - 2]]);
  pop[dx[1]].copy(pop[dx[tsize - 1]]);
  pop[dx[0]].tpc(pop[dx[1]]);
  nm = (int) lrand48() % MNM + 1;
  pop[dx[0]].mutate(nm);
  nm = (int) lrand48() % MNM + 1;
  pop[dx[1]].mutate(nm);
  fit[dx[0]] = fitness(pop[dx[0]]);
  fit[dx[1]] = fitness(pop[dx[1]]);
}

void fillEdges() {
  int psn = 0;
  int numNbrs;
  int nbrs[D->size()];

  for (int i = 0; i < D->size(); i++) {
    numNbrs = D->Nbrs(i, nbrs);
    for (int j = 0; j < numNbrs; j++) {
      if (i < nbrs[j]) {
        doors[psn][0] = i;
        doors[psn][1] = nbrs[j];
        psn++;
      }
    }
  }
}