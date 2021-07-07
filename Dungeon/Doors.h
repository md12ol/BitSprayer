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
#define mevs 20
#define RI (long)(mevs/5)

// Parameters for tuning
#define popsize 12
#define tsize 7
#define MNM 2
#define states 16

// Board parameters
#define BB 10   //Grid size
#define Rz 512  //Number of rooms
#define Z 80    // Board size
int Rx[Rz], Ry[Rz], Dx[Rz], Dy[Rz];     // Room data

// Bitsprayer parameters
#define Qz 2000
int Q[Qz];

#define diamTests 3

void initalg();
void initpop(int run);
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
double reportbest(ostream &aus, int run, char *outLoc);
int bestDiam();
void cmdLineRun(int run, ostream &aus);
void cmdLineIntro(ostream &aus, int g);
void printGraph(int run, char *outLoc);
void approxTest();
void getRooms(fstream rmData);

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