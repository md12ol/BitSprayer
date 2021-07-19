#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <queue>
#include <climits>
#include <algorithm>
#include <random>
#include <filesystem>

#define verbose true
#define RNS 91207819

// GA Parameters
#define runs 30
#define mevs 50000
#define RI (long)(mevs/100)

// Parameters for tuning
#define popsize 50
#define tsize 7
#define MNM 3
#define states 12

// Board parameters
#define BB 10   //Grid size
#define Rz 512  //Number of rooms
#define Z 80    // Board size
int Rx[Rz], Ry[Rz], Dx[Rz], Dy[Rz];     // Room data

// Bitsprayer parameters
int Qz;
int *Q;

#define diamTests 10

void initalg();
void initpop(int run);
vector <pair<int, int>> BFS(graph &graph, int source);

pair<int, int> fitness(bitspray &A, bool actDiam);
int diameter(graph &graph);
void report(ostream &aus);
void matingevent();
bool getnum(int &val, int bits, int &psn);
void randomOrder(int size);
void fillEdges();
int approxDiameter(graph &graph);
bool getbit(int &val, int &psn);
double reportbest(ostream &aus, int run, char *outLoc);
int bestDiam();
void fitness_check();
void cmdLineRun(int run, ostream &aus);
void cmdLineIntro(ostream &aus, int g);
void printGraph(int run, char *outLoc);
void approxTest();
void getRooms(fstream &rmData);
void developQ(bitspray &A);
bool getbit(int &val, int &psn);
void diameterTest();
int old_fitness(bitspray &A);
int validation(bitspray &A);
bool necroticFilter();
int singleFitness(pair<int, int> vals);
string printPair(pair<int, int> vals);

bitspray pop[popsize];  //Population of bitsprayers
pair<int, int> fit[popsize];  //Fitness values
int sfit[popsize];
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

#define diamMult 10
#define doorMult 1