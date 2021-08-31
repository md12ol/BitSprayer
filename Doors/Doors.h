#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <queue>
#include <climits>
#include <algorithm>
#include <random>
#include <filesystem>
#include <bits/stdc++.h>

#define verbose false
#define RNS 91204658

// GA Parameters
#define runs 30
#define mevs 50000
#define RI (long)(mevs/100)
#define VI (long)(mevs/10)

// Parameters for tuning
int popsize;
#define tsize 7
int MNM;
int states;

// Board parameters
#define BB 10   //Grid size
#define Rz 512  //Number of rooms
#define Z 80    // Board size
int Rx[Rz], Ry[Rz], Dx[Rz], Dy[Rz];     // Room data

// Bitsprayer parameters
int Qz;
int *Q;

#define diamTests 15

// Evolution methods
void initAlg();
void initPop(int run);
void matingevent();
void tournament();

// Bitsprayer methods
bool getnum(int &val, int bits, int &psn);
void randomOrder(int size);
void fillEdges();
bool getbit(int &val, int &psn);
void getRooms(fstream &rmData);
void developQ(bitspray &a);

// Output methods
void cmdLineRun(int run, ostream &aus);
void cmdLineIntro(ostream &aus, int g);
void report(ostream &aus);
string reportBest(ostream &aus, int run, char *outLoc);
void printGraph(int run, char *outLoc);
pair<pair<int, int>, pair<int, int>> getOverlap(int r1, int r2);
void render(int run, char *outLoc);

// Fitness methods
pair<int, int> fitness(bitspray &a, bool actDiam);
bool necroticFilter();
pair<int, int> dsOfGFrA(bitspray &a);
void actualFitnessUpdate(ostream &aus);
vector <pair<int, int>> getBests();

// Diameter methods
int bestDiam();
int diameter(graph &graph);
int approxDiameter(graph &graph);
vector <pair<int, int>> bfs(graph &graph, int source);

// Debugging methods
vector <pair<int, int>> getEdges(int numEdges);
void diameterTest();
void fitness_check();

// Extra methods
string printPair(int idx);
void approxTest();
int old_fitness(bitspray &A);
int validation(bitspray &A);
int singleFitness(pair<int, int> vals);

bitspray *pop;  //Population of bitsprayers
int *dfit;
int *efit;
//pair<int, int> fit[popsize];  //Fitness values
//int sfit[popsize];
int *dx;  //Sorting index

graph *initG;
graph *G;
graph *D;
int doors[1000][2];
int posDoors;
vector<int> order;
bool over = false;
int startD = -1;
int bestD = -1;
int actBest;

#define diamMult 10
#define doorMult 1