#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <queue>
#include <climits>
#include <algorithm>
#include <random>
#include <filesystem>

#define RNS 91207819
#define runs 30
#define mevs 500
#define RI (long)(mevs/100)
#define popsize 64
#define tsize 7
#define MNM 2
#define verbose true

#define states 16
#define Qz 2000
int Q[Qz];

#define BB 10 //Grid size
#define Rz 512
#define reqDiam 35
#define diamTests 3

void initalg();
void initpop();
vector <pair<int, int>> BFS(graph &G, int source);
int fitness(bitspray &A);
int diameter(graph &G);
void report(ostream &aus);
void matingevent();
bool getnum(int &val, int bits, int &psn);
void randomOrder(int size);
void fillEdges();
int approxDiameter(graph &G);
bool getbit(int &val, int &psn);
int validation(bitspray &A);
double reportbest(ostream &aus, int run);
int bestDiam();
void cmdLineRun(int run, ostream &aus);
void cmdLineIntro(ostream &aus, int g);

bitspray pop[popsize];  //Population of bitsprayers
int fit[popsize];  //Fitness values
int dx[popsize];  //Sorting index

graph *initG;
graph *G;
graph *D;
int doors[1000][2];
int posDoors;
vector<int> order;
bool over = false;
int startD = -1;
int bestD = -1;