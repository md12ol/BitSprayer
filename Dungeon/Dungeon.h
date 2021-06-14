#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <list>
#include <cmath>
#include <filesystem>
#include <iomanip>

#define RNS 91207819
#define runs 30
#define mevs 10000
#define RI (long)(mevs/100)
#define popsize 24
#define tsize 7
#define MNM 2
#define verbose 1

#define states 12
#define Qz 100000
int Q[Qz];

#define Z 80 //Board size
int Bd[Z][Z];

#define BB 10 //Grid size

#define Rz 512
int cnr;
int Rx[Rz], Ry[Rz], Dx[Rz], Dy[Rz];

#define Rbs 8
#define RecRs pow(2,Rbs)
vector<int> Rr;

void initalg();                //initialize the algorithm
void developQ(bitspray &A);         //unpack the queue
double fitness(bitspray &A);   //compute the fitness of an alternator
void initpop();                //initialize a population
void matingevent();            //run a mating event
void report(ostream &aus);     //report current summary statistics
void render(int run, char *outLoc);          //render a picture
double
reportbest(ostream &aus, int run, char *outLoc); //report current best creature
void printGraph(int run, char *outLoc);
void printBoard(int run, char *outLoc);
int getRoomNm(int x, int y);
void printDoors(int run, char *outLoc);
void cmdLineRun(int run, ostream &aus);

bitspray *pop[popsize];  //Population of bitsprayers
double fit[popsize];  //Fitness values
int dx[popsize];  //Sorting index

graph *G;
