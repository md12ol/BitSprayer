#include "../BitSprayer/stat.h"
#include "../BitSprayer/bitspray.h"
#include "setu.h"
#include "Doors.h"
using namespace std;

int main() {
    fstream stat;   //statistics reporting stream
    fstream best;   //best reporting stream
    fstream initGraph;
    fstream outGraph;
    fstream doorGraph;
    char fn[100];
    char *outLoc = new char[45];
    char *outRoot = new char[20];
    sprintf(outRoot, "./DoorOut/");

    for (int gNum = 1; gNum <= 3; gNum++) {
        sprintf(outLoc, "%sG%d Output - %02dS, %02dP, %dM/",
                outRoot, gNum, states, popsize, MNM);
        std::filesystem::create_directory(outLoc);

        sprintf(fn, "./Input/graph%d.dat", gNum);
        initGraph.open(fn, ios::in);
        initG = new graph(Rz);
        initG->read(initGraph);
        initGraph.close();

        sprintf(fn, "./Input/doors%d.dat", gNum);
        doorGraph.open(fn, ios::in);
        D = new graph(Rz);
        D->read(doorGraph);
        doorGraph.close();
        posDoors = D->edges();
        fillEdges();

        G = new graph(Rz);
        order.reserve(posDoors);
        randomOrder(posDoors);
        sprintf(fn, "%sbest.sda", outLoc);
        best.open(fn, ios::out);

        startD = diameter(*initG);
        bestD = bestDiam();
        cmdLineIntro(cout, gNum);

        //  sprintf(fn, "Output/newg%02d.dat", g);
        //  outGraph.open(fn, ios::out);
        initalg();

        for (int run = 0; run < runs; run++) {
            sprintf(fn, "%srun%02d.dat", outLoc, run);
            stat.open(fn, ios::out);
            if (verbose) cmdLineRun(run, cout);
            initpop();
            report(stat);
            for (int mev = 0; mev < mevs; mev++) {
                matingevent();
                if ((mev + 1) % RI == 0) {
                    if (verbose) {
                        cout << left << setw(5) << run;
                        cout << left << setw(4) << (mev + 1) / RI;
                    }
                    report(stat);
                }
            }
            reportbest(best, run);
            stat.close();
        }
        best.close();
        outGraph.close();
    }
    delete initG;
    delete D;
    delete G;
    return 0;
}

void cmdLineIntro(ostream &aus, int g) {
    aus << "Dungeon with Doors Generator." << endl;
    aus << "Starting Dungeon: " << g << endl;
    aus << "Starting Dungeon Diameter: " << startD << endl;
    aus << "Best Dungeon Diameter Possible: " << bestD << endl;
//  aus << "Check readme.dat for more information about parameters/output.";
//  aus << endl;
}

void cmdLineRun(int run, ostream &aus) {
    aus << endl << "Beginning Run " << run << " of " << runs - 1 << endl;
    aus << left << setw(5) << "Run";
    aus << left << setw(4) << "RI";
    aus << left << setw(10) << "Mean";
    aus << left << setw(12) << "95% CI";
    aus << left << setw(10) << "SD";
    aus << left << setw(8) << "Best";
    aus << endl;
    aus << left << setw(5) << run;
    aus << left << setw(4) << "0";
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
    cout << "Best D is " << bestD << std::endl;
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
    if (ratio < 0.1 || ratio > 0.9) {
        return true;
    } else {
        return false;
    }
}

int bestDiam() {
    G->empty(Rz);
    G->copy(*initG);
    for (int i = 0; i < posDoors; i++) {
        G->add(doors[i][0], doors[i][1]);
    }
    return diameter(*G);
}

int fitness(bitspray &A) {
    int addedDoors = 0;
    int acpt = 0;
    int psn = 0;
    int idx = 0;
    G->empty(Rz);
    G->copy(*initG);
    developQ(A);
    if (necroticFilter()) {
        return posDoors;
    }

    while (approxDiameter(*G) > bestD && getbit(acpt, psn)) {
        if (acpt == 1) {
            if (G->add(doors[order[idx]][0], doors[order[idx]][1])) {
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
    G->empty(Rz);
    G->copy(*initG);
    developQ(A);

    while (diameter(*G) > bestD && getbit(acpt, psn)) {
        if (acpt == 1) {
            if (G->add(doors[order[idx]][0], doors[order[idx]][1])) {
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
    int size = G.size();
    int dist[size];
    int path[size];
    queue<int> qu;
    int node, numnbr;

    for (int i = 0; i < size; i++) {
        dist[i] = INT_MAX;
        path[i] = -1;
    }
    dist[source] = 0;

    qu.push(source);
    int nbrs[size];
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
    for (int i = 0; i < size; i++) {
        so.push_back(make_pair(dist[i], i));
    }
    sort(so.begin(), so.end());
    return so;
}

int diameter(graph &G) {
    int size = G.size();
    int maxs[size];
    for (int i = 0; i < size; i++) {
        maxs[i] = BFS(G, i).back().first;
    }
    int max = 0;
    for (int i = 0; i < size; i++) {
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
    aus << left << setw(10) << D.Rmu();
    aus << left << setw(12) << D.RCI95();
    aus << left << setw(10) << D.Rsg();
    aus << left << setw(8) << D.Rmin();
    aus << endl;
    if (verbose) {
        cout << left << setw(10) << D.Rmu();
        cout << left << setw(12) << D.RCI95();
        cout << left << setw(10) << D.Rsg();
        cout << left << setw(8) << D.Rmin();
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