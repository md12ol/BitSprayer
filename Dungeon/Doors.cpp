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
    fstream roomData;
    char fn[100];
    char *outLoc = new char[45];
    char *outRoot = new char[20];
    sprintf(outRoot, "./DoorOut/");

    for (int gNum = 3; gNum <= 3; gNum++) {
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
        Qz = posDoors;
        Q = new int[Qz];

        sprintf(fn, "./Input/rooms%d.dat", gNum);
        roomData.open(fn, ios::in);
        getRooms(roomData);
        roomData.close();

        G = new graph(Rz);
        order.reserve(posDoors);
        randomOrder(posDoors);

        sprintf(fn, "%sbest.sda", outLoc);
        best.open(fn, ios::out);

        startD = diameter(*initG);
        bestD = bestDiam();
        cmdLineIntro(cout, gNum);

        //  outGraph.open(fn, ios::out);
        initalg();

        for (int run = 0; run < runs; run++) {
            sprintf(fn, "%srun%02d.dat", outLoc, run);
            stat.open(fn, ios::out);
            initpop(run);
            if (verbose) cmdLineRun(run, cout);
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
            reportbest(best, run, outLoc);

            stat.close();
        }
        best.close();
        outGraph.close();
//        cout << "Graph " << gNum << " with " << posDoors << " potential edges."
//             << endl;
//        cout << "Start: " << startD << endl;
//        cout << "Best: " << bestD << endl;
//        diameterTest();
        cout << endl;
    }
    delete initG;
    delete D;
    delete G;
    delete[] Q;
    return 0;
}

vector<pair<int, int>> getDoors() {
    int acpt = 0;
    int psn = 0;
    int idx = 0;
    int r1, r2;
    vector<pair<int, int>> el;
    vector<int> ds;

    while (getbit(acpt, psn)) {
        if (acpt == 1) {
            ds.push_back(idx);
        }
        idx++;
    }
    for (int d:ds) {
        r1 = doors[order[d]][0];
        r2 = doors[order[d]][1];
        el.emplace_back(make_pair(r1, r2));
    }
    return el;
}

vector<pair<int, int>> getEdges(int numEdges) {
    vector<pair<int, int>> edges;
    vector<int> ds;
    edges.reserve(numEdges);
    ds.reserve(numEdges);
    int rnd, r1, r2;
    for (int i = 0; i < numEdges; i++) {
        do {
            rnd = (int) lrand48() % posDoors;
        } while (*find(ds.begin(), ds.end(), rnd) == rnd);
        ds.push_back(rnd);
    }
    for (int d:ds) {
        r1 = doors[order[d]][0];
        r2 = doors[order[d]][1];
        edges.emplace_back(r1, r2);
    }
    return edges;
}

void diameterTest() {
    int apx, act;
    int numTests = 1000;
    int numEdges = (int) (0.80 * posDoors);
    int errors = 0;
    int bests = 0;
    int bestWEr = 0;
    srand48(RNS);
    G->empty(Rz);
    G->copy(*initG);
    vector<pair<int, int>> edges;
    cout << "Adding " << numEdges << " doors." << endl;
//    cout << "Apx" << '\t' << "Act" << endl;
    for (int idx = 0; idx < numTests; idx++) {
        edges = getEdges(numEdges);
        for (pair<int, int> e:edges) {
            G->add(e.first, e.second);
        }
        apx = approxDiameter(*G);
        act = diameter(*G);
        if (act == bestD) {
//            cout << "B\t";
            bests++;
        }
        if (apx != act || apx == bestD) {
            cout << "Test " << idx << '\t';
            cout << apx << '\t' << act;
            if (apx == act) {
                cout << "\tB";
            } else {
                errors++;
                cout << '\t' << apx - act;
            }
//            cout << '\t' << "Edges: ";
//            for (pair<int, int> e:edges) {
//                cout << "[" << e.first << ", " << e.second << "], ";
//            }
            cout << endl;
        }

        for (pair<int, int> e:edges) {
            G->del(e.first, e.second);
        }
    }
    cout << "Total errors: " << errors << '\t' << (float) 100 * errors /
        numTests <<
         '%' << endl;
    cout << "Total bests: " << bests << '\t' << (float) 100 * bests / numTests
         <<
         '%'
         << endl;
}

void getRooms(fstream &rmData) {
    char buf[1000];
    int xin, yin, dxin, dyin, k, rms;

    rmData.getline(buf, 999);
    rms = atoi(buf);

    for (int rm = 0; rm < rms; rm++) {
        rmData.getline(buf, 999);
        xin = atoi(buf);
        k = 0;
        while (buf[k] != '\t')k++;
        while (buf[k] == '\t')k++;
        yin = atoi(buf + k);
        while (buf[k] != '\t')k++;
        while (buf[k] == '\t')k++;
        dxin = atoi(buf + k);
        while (buf[k] != '\t')k++;
        while (buf[k] == '\t')k++;
        dyin = atoi(buf + k);

        Rx[rm] = xin;
        Ry[rm] = yin;
        Dx[rm] = dxin;
        Dy[rm] = dyin;
    }
}

void cmdLineIntro(ostream &aus, int g) {
    aus << "Dungeon with Doors Generator." << endl;
    aus << "Starting Dungeon: " << g << endl;
    aus << "Starting Dungeon Diameter: " << startD << endl;
    aus << "Best Dungeon Diameter Possible: " << bestD << endl;
    aus << "Maximum Doors to Add (old_fitness): " << posDoors << endl;
//  aus << "Check readme.dat for more information about parameters/output.";
//  aus << endl;
}

void cmdLineRun(int run, ostream &aus) {
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
}

void initalg() {
    srand48(RNS);
    for (int i = 0; i < popsize; i++) {
        pop[i].create(states);
    }
}

void initpop(int run) {
    cout << endl << "Beginning Run " << run << " of " << runs - 1 << endl;
    cout << "Initial population old_fitness values" << std::endl;
    for (int i = 0; i < popsize; i++) {
        do {
            pop[i].randomize();
            developQ(pop[i]);
        } while (necroticFilter());
        fit[i] = fitness(pop[i], false);
        sfit[i] = singleFitness(fit[i]);
        cout << "[" << fit[i].first << ", " << fit[i].second << "] ";
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
    if (ratio < 0.4 || ratio > 0.9) {
        return true;
    } else {
        return false;
    }
}

void fitness_check() {
    for (int i = 0; i < popsize; ++i) {
        fit[i] = fitness(pop[i], true);
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

pair<int, int> fitness(bitspray &A, bool actDiam) {
    int addedDoors = 0;
    int curDiam = startD;
    vector<pair<int, int>> el;
    G->empty(Rz);
    G->copy(*initG);
    developQ(A);

    if (necroticFilter()) {
        return make_pair(startD, posDoors);
    }

    el = getDoors();
    //  Add all the doors
    for (pair<int, int> p : el) {
        G->add(p.first, p.second);
    }
    addedDoors = (int) el.size();
    if (actDiam) {
        curDiam = diameter(*G);
    } else {
        curDiam = approxDiameter(*G);
    }

    return make_pair(curDiam, addedDoors);

//    if (curDiam == bestD) {
//        if (curDiam == diameter(*G)){
//            return el.size();
//        } else {
//            return el.size() + curDiam;
//        }
//        G->empty(Rz);
//        G->copy(*initG);
//        curDiam = startD;
//        auto ds = el.begin();
//        while (curDiam > bestD && idx < el.size()) {
//            G->add(ds[idx].first, ds[idx].second);
//            curDiam = approxDiameter(*G);
//            addedDoors++;
//            idx = idx + 1;
//        }
//        return addedDoors;
//    } else {
//        return posDoors;
//    }
}

int singleFitness(pair<int, int> vals) {
    int distFromBest = vals.first - bestD;
    int penalty = 0;
    if (distFromBest == 0) {
        penalty = 0;
    } else {
        penalty = (int) (10 * pow(2, distFromBest)) + 50;
    }
    if (penalty > posDoors) {
        penalty = posDoors;
    }
    return vals.second + penalty;
}

int old_fitness(bitspray &A) {
    int addedDoors = 0;
    int acpt = 0;
    int psn = 0;
    int idx = 0;
    G->empty(Rz);
    G->copy(*initG);
    developQ(A);
    if (necroticFilter()) {
        return posDoors + startD;
    }

    while (approxDiameter(*G) > bestD && getbit(acpt, psn)) {
        if (acpt == 1) {
            if (G->add(doors[order[idx]][0], doors[order[idx]][1])) {
                addedDoors++;
            }
        }
        idx = (idx + 1) % posDoors;
    }

    if (approxDiameter(*G) > bestD) {
        return posDoors + diameter(*G);
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
        return posDoors + diameter(*G);
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

vector<pair<int, int>> BFS(graph &graph, int source) {
    int size = graph.size();
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
        numnbr = graph.Nbrs(node, nbrs);
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
        so.emplace_back(make_pair(dist[i], i));
    }
    sort(so.begin(), so.end());
    return so;
}

int diameter(graph &graph) {
    int size = graph.size();
    int maxs[size];
    for (int i = 0; i < size; i++) {
        maxs[i] = BFS(graph, i).back().first;
    }
    int max = 0;
    for (int i = 0; i < size; i++) {
        if (maxs[i] > max) {
            max = maxs[i];
        }
    }
    return max;
}

int approxDiameter(graph &graph) {
    int maxs[diamTests];
    int from = 0;
    int nextF;
    vector<int> tested;
    vector<pair<int, int>> so;
    for (int i = 0; i < diamTests; i++) {
        maxs[i] = 0;
        tested.push_back(from);
        so = BFS(graph, from);
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
    dset dataset;
    int singIdx = 0;
    int sBest = posDoors + startD; // TODO: Fix
    int diamIdx = 0;
    int dBest = startD;
    int edgeIdx = 0;
    int eBest = posDoors;
    for (int i = 0; i < popsize; i++) {
        if (sfit[i] < sfit[singIdx]) {
            singIdx = i;
//            sBest = sfit[singIdx];
        }
    }
    for (int i = 0; i < popsize; i++) {
        if (fit[i].first < fit[diamIdx].first) {
            diamIdx = i;
//            dBest = fit[dBest].first;
        }
    }
    for (int i = 0; i < popsize; i++) {
        if (fit[i].second < fit[edgeIdx].second) {
            edgeIdx = i;
//            eBest = fit[edgeIdx].second;
        }
    }
    dataset.add(sfit, popsize);
    aus << left << setw(10) << dataset.Rmu();
    aus << left << setw(12) << dataset.RCI95();
    aus << left << setw(10) << dataset.Rsg();
    aus << left << setw(8) << dataset.Rmin();
//    aus << left << printPair(fit[singIdx]) << '\t';
//    aus << left
    aus << endl;
    if (verbose) {
        cout << left << setw(10) << dataset.Rmu();
        cout << left << setw(12) << dataset.RCI95();
        cout << left << setw(10) << dataset.Rsg();
        cout << left << setw(8) << dataset.Rmin();
        cout << left << printPair(fit[singIdx]) << '\t';
        cout << left << printPair(fit[diamIdx]) << '\t';
        cout << left << printPair(fit[edgeIdx]);
        cout << endl;
    }

}

string printPair(pair<int, int> vals) {
    string rtn = "[";
    rtn += to_string(vals.first);
    rtn += ", ";
    rtn += to_string(vals.second);
    rtn += "]";
    return rtn;
}

pair<int, int> dsOfGFrA(bitspray &A) {
    vector<pair<int, int>> el;
    G->empty(Rz);
    G->copy(*initG);
    developQ(A);

    el = getDoors();
    for (pair<int, int> p : el) {
        G->add(p.first, p.second);
    }

    int apx, act;
    apx = approxDiameter(*G);
    act = diameter(*G);
    return make_pair(apx, act);
}

void approxTest() {
    int bestApx = 0;
    int bestAct = 0;
    int diffs = 0;
    cout << "Beginning Approx Test" << endl;
    cout << "B?\t" << "Ds\t" << "App\t" << "Act" << endl;
    for (int i = 0; i < popsize; i++) {
        if (fit[i].first == bestD) {
            cout << "B\t";
            bestApx++;
        } else {
            cout << '\t';
        }
        cout << fit[i].second << '\t';
        pair<int, int> diams = dsOfGFrA(pop[i]);
        if (diams.first != diams.second) {
            diffs++;
        }
        if (fit[i].second == bestD && diams.first == diams.second) {
            bestAct++;
        }
        cout << diams.first << '\t' << diams.second << endl;
    }
    cout << "Errors: " << diffs << " " << (float) diffs * 100 / popsize << "%"
         << endl;
    cout << "Approx. Best: " << bestApx << " " << (float) bestApx * 100 /
        popsize <<
         "%" << endl;
    cout << "Actual Best: " << bestAct << " " << (float) bestAct * 100 /
        (float) bestApx << "%" << endl;
    cout << "End of Approx Test" << endl;
}

double reportbest(ostream &aus, int run, char *outLoc) {
    int b = 0;
    for (int i = 1; i < popsize; i++) {
        if (sfit[i] < sfit[b]) b = i;
    }
    cout << "Best Approx Fitness: " << sfit[b] << endl;
//    fit[b] = validation(pop[b]);
    approxTest();
    aus << "Run" << run << " fitness: " << sfit[b] << endl;
//    cout << "Best Fitness: " << fit[b] << endl;
//  render(run);
    printGraph(run, outLoc);
    pop[b].print(aus);
    return sfit[b];
}

void printGraph(int run, char *outLoc) {
    fstream graph;
    char fn[100];

    sprintf(fn, "%sgraph%02d.dat", outLoc, run);
    graph.open(fn, ios::out);
    G->resetV();
    G->write(graph);
    cout << "Graph printed." << endl;
    graph.close();
}

void matingevent() {
    int nm;

    Tselect(sfit, dx, tsize, popsize);
    pop[dx[0]].copy(pop[dx[tsize - 2]]);
    pop[dx[1]].copy(pop[dx[tsize - 1]]);
    pop[dx[0]].tpc(pop[dx[1]]);
    nm = (int) lrand48() % MNM + 1;
    pop[dx[0]].mutate(nm);
    nm = (int) lrand48() % MNM + 1;
    pop[dx[1]].mutate(nm);
    fit[dx[0]] = fitness(pop[dx[0]], false);
    sfit[dx[0]] = singleFitness(fit[dx[0]]);
    fit[dx[1]] = fitness(pop[dx[1]], false);
    sfit[dx[1]] = singleFitness(fit[dx[1]]);
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