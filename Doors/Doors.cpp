#include "BitSprayer/stat.h"
#include "BitSprayer/bitspray.h"
#include "Dungeon/setu.h"
#include "Dungeon/ps.h"
#include "Doors.h"
using namespace std;

int main(int argc, char *argv[]) {
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

    char *end;
    states = (int) strtol(argv[1], &end, 10);
    popsize = (int) strtol(argv[2], &end, 10);
    MNM = (int) strtol(argv[3], &end, 10);
    int gNum = (int) strtol(argv[4], &end, 10);

    pop = new bitspray[popsize];
    dfit = new int[popsize];
    efit = new int[popsize];
    dx = new int[popsize];

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
    G->copy(*initG);
    render(31, outLoc);

    sprintf(fn, "%sbest.sda", outLoc);
    best.open(fn, ios::out);

    startD = diameter(*initG);
    bestD = bestDiam();
    cmdLineIntro(cout, gNum);
    initAlg();

    for (int run = 0; run < runs; run++) {
        sprintf(fn, "%srun%02d.dat", outLoc, run);
        stat.open(fn, ios::out);
        initPop(run);
        if (verbose) cmdLineRun(run, cout);
        report(stat);
        for (int mev = 0; mev < mevs; mev++) {
            matingevent();
            if ((mev + 1) % RI == 0) {
                if ((mev + 1) % VI == 0) {
                    actualFitnessUpdate(stat);
                }
                stat << endl;
                if (verbose) {
                    cout << left << setw(5) << run;
                    cout << left << setw(4) << (mev + 1) / RI;
                }
                report(stat);
            }
        }
        reportBest(best, run, outLoc);
        stat.close();
    }
    best.close();
    outGraph.close();
    cout << endl;

    delete initG;
    delete D;
    delete G;
    delete[] Q;
    delete[] pop;
    delete[] dfit;
    delete[] efit;
    delete[] dx;
    delete[] outRoot;
    delete[] outLoc;
    for (int i = 0; i < popsize; i++) {
        pop[i].destroy();
    }
    return 0;
}

void getRooms(fstream &rmData) {
    char buf[1000];
    char *end;
    int xin, yin, dxin, dyin, k, rms;

    rmData.getline(buf, 999);
    rms = (int) strtol(buf, &end, 10);

    for (int rm = 0; rm < rms; rm++) {
        rmData.getline(buf, 999);
        xin = (int) strtol(buf, &end, 10);
        k = 0;
        while (buf[k] != '\t')k++;
        while (buf[k] == '\t')k++;
        yin = (int) strtol(buf + k, &end, 10);
        while (buf[k] != '\t')k++;
        while (buf[k] == '\t')k++;
        dxin = (int) strtol(buf + k, &end, 10);
        while (buf[k] != '\t')k++;
        while (buf[k] == '\t')k++;
        dyin = (int) strtol(buf + k, &end, 10);

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
    aus << "Best Dungeon Diameter Possible (Diam Fitness): " << bestD << endl;
    aus << "Maximum Doors to Add (Door Fitness): " << posDoors << endl;
//  aus << "Check readme.dat for more information about parameters/output.";
//  aus << endl;
}

void cmdLineRun(int run, ostream &aus) {
    aus << left << setw(5) << "Run";
    aus << left << setw(4) << "dRI";
    aus << left << setw(10) << "dMean";
    aus << left << setw(12) << "d95% CI";
    aus << left << setw(10) << "dSD";
    aus << left << setw(8) << "dBest";
    aus << left << setw(10) << "eMean";
    aus << left << setw(12) << "e95% CI";
    aus << left << setw(10) << "eSD";
    aus << left << setw(8) << "eBest";
    aus << left << "Best Overall";
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

void initAlg() {
    srand48(RNS);
    for (int i = 0; i < popsize; i++) {
        pop[i].create(states);
    }
}

void initPop(int run) {
    cout << endl << "Beginning Run " << run << " of " << runs - 1 << endl;
    cout << "Initial population fitness values" << std::endl;
    pair<int, int> fits;
    for (int i = 0; i < popsize; i++) {
        do {
            pop[i].randomize();
            developQ(pop[i]);
        } while (necroticFilter());
        fits = fitness(pop[i], true);
        dfit[i] = fits.first;
        efit[i] = fits.second;
        cout << "[" << dfit[i] << ", " << efit[i] << "] ";
        dx[i] = i;
    }
    vector<pair<int, int>> bests;
    bests.reserve(popsize);
    bests = getBests();
    actBest = bests.front().first;
    cout << endl;
}

void developQ(bitspray &a) {//unpack the queue
    int h, t;  //head and tail of queue

    over = false;
    for (t = 0; t < Qz; t++)Q[t] = 0;        //clear the queue
    a.reset(Q, h, t);                 //reset the self driving automata
    while (t < Qz - 2)a.next(Q, h, t, Qz);  //run the automata
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

int bestDiam() {
    G->empty(Rz);
    G->copy(*initG);
    for (int i = 0; i < posDoors; i++) {
        G->add(doors[i][0], doors[i][1]);
    }
    return diameter(*G);
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

void actualFitnessUpdate(ostream &aus) {
    pair<int, int> fits;
    dset ddat, edat;
    ddat.add(dfit, popsize);
    edat.add(efit, popsize);
    cout << "Fitness Update\t" << ddat.Rmu() << ", " << edat.Rmu();
    aus << "\tFitness Update\t" << ddat.Rmu() << ", " << edat.Rmu();
    for (int i = 0; i < popsize; i++) {
        fits = fitness(pop[i], true);
        dfit[i] = fits.first;
        efit[i] = fits.second;
    }
    ddat.clear();
    edat.clear();
    ddat.add(dfit, popsize);
    edat.add(efit, popsize);
    cout << " --> " << ddat.Rmu() << ", " << edat.Rmu() << endl;
    aus << "-->" << ddat.Rmu() << ", " << edat.Rmu();
    vector<pair<int, int>> bests = getBests();
    bests.reserve(popsize);
    if (actBest != bests.front().first) {
        actBest = bests.front().first;
    }
}

pair<int, int> fitness(bitspray &a, bool actDiam) {
    int addedDoors = 0;
    int curDiam = startD;
    vector<pair<int, int>> el;
    G->empty(Rz);
    G->copy(*initG);
    developQ(a);

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

vector<pair<int, int>> bfs(graph &graph, int source) {
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
    so.reserve(size);
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
        maxs[i] = bfs(graph, i).back().first;
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
    tested.reserve(diamTests);
    vector<pair<int, int>> so;
    so.reserve(diamTests);
    for (int &max : maxs) {
        max = 0;
        tested.push_back(from);
        so = bfs(graph, from);
        max = so.back().first;
        nextF = (int) so.size() - 1;
        do {
            from = so.at(nextF).second;
            nextF--;
        } while (*find(tested.begin(), tested.end(), from) == from);
    }
    int max = 0;
    for (int i : maxs) {
        if (i > max) {
            max = i;
        }
    }
    return max;
}

void tournament() {
    int sw, rp;
    for (int i = 0; i < tsize; i++) {
        rp = (int) (lrand48() % popsize); //select random member
        sw = dx[i];
        dx[i] = dx[rp];
        dx[rp] = sw;
    }

    do {
        rp = 0;
        for (int i = 0; i < tsize - 1; i++) {
            if (dfit[dx[i]] < dfit[dx[i + 1]]) { // diam dominates
                sw = dx[i];
                dx[i] = dx[i + 1];
                dx[i + 1] = sw; //swap
                rp = 1;  //set flag
            } else if (dfit[dx[i]] == dfit[dx[i + 1]]) { // equal diam
                if (efit[dx[i]] < efit[dx[i + 1]]) {
                    sw = dx[i];
                    dx[i] = dx[i + 1];
                    dx[i + 1] = sw; //swap
                    rp = 1;  //set flag
                }
            }
        }
    } while (rp == 1);
}

void matingevent() {
    int nm;
    pair<int, int> fits;

    do {
        tournament();
    } while (actBest == dx[0] || actBest == dx[1]);
    pop[dx[0]].copy(pop[dx[tsize - 2]]);
    pop[dx[1]].copy(pop[dx[tsize - 1]]);
    pop[dx[0]].tpc(pop[dx[1]]);
    nm = (int) lrand48() % MNM + 1;
    pop[dx[0]].mutate(nm);
    nm = (int) lrand48() % MNM + 1;
    pop[dx[1]].mutate(nm);
    fits = fitness(pop[dx[0]], false);
    dfit[dx[0]] = fits.first;
    efit[dx[0]] = fits.second;
    fits = fitness(pop[dx[1]], false);
    dfit[dx[1]] = fits.first;
    efit[dx[1]] = fits.second;
}

bool sortbysec(const pair<int, int> &a, const pair<int, int> &b) {
    return a.second < b.second;
}

vector<pair<int, int>> getBests() {
    int dBest = startD;
    vector<pair<int, int>> bests; // idx, edges
    bests.reserve(popsize);
    for (int i = 0; i < popsize; i++) {
        if (dfit[i] < dBest) {
            dBest = dfit[i];
        }
    }
    for (int i = 0; i < popsize; i++) {
        if (dfit[i] == dBest) {
            bests.emplace_back(make_pair(i, efit[i]));
        }
    }
    sort(bests.begin(), bests.end(), sortbysec);
    return bests;
}

void report(ostream &aus) {
    dset diamdata, edgedata;
    int bestIdx;
    vector<pair<int, int>> bests; // idx, edges
    bests.reserve(popsize);
    bests = getBests();
    bestIdx = bests.front().first;
    diamdata.add(dfit, popsize);
    edgedata.add(efit, popsize);
    aus << left << setw(10) << diamdata.Rmu();
    aus << left << setw(12) << diamdata.RCI95();
    aus << left << setw(10) << diamdata.Rsg();
    aus << left << setw(8) << diamdata.Rmin();
    aus << left << setw(10) << edgedata.Rmu();
    aus << left << setw(12) << edgedata.RCI95();
    aus << left << setw(10) << edgedata.Rsg();
    aus << left << setw(8) << edgedata.Rmin();
    aus << left << printPair(bestIdx);
    if (verbose) {
        cout << left << setw(10) << diamdata.Rmu();
        cout << left << setw(12) << diamdata.RCI95();
        cout << left << setw(10) << diamdata.Rsg();
        cout << left << setw(8) << diamdata.Rmin();
        cout << left << setw(10) << edgedata.Rmu();
        cout << left << setw(12) << edgedata.RCI95();
        cout << left << setw(10) << edgedata.Rsg();
        cout << left << setw(8) << edgedata.Rmin();
        cout << left << printPair(bestIdx);
        cout << endl;
    }
    pair<int, int> tfits = fitness(pop[bestIdx], true);
    if (bestIdx != actBest) {
        pair<int, int> bfits = fitness(pop[actBest], true);
        if (tfits.first <= bfits.first && tfits.second <= bfits.second) {
            actBest = bestIdx;
        }
    } else {
        dfit[bestIdx] = tfits.first;
        efit[bestIdx] = tfits.second;
    }
}

string reportBest(ostream &aus, int run, char *outLoc) {
    int b = 0;
    vector<pair<int, int>> bests;
    bests.reserve(popsize);
    bests = getBests();
    int bestIdx = bests.front().first;
    cout << "Best Approx Fitness: " << printPair(bestIdx) << endl;
    aus << "Run" << run << " Approx Fitness: " << printPair(bestIdx) << endl;
    pair<int, int> fits = fitness(pop[bestIdx], true);
    dfit[bestIdx] = fits.first;
    efit[bestIdx] = fits.second;
    cout << "Best Actual Fitness: " << printPair(bestIdx) << endl;
    aus << "Run" << run << " Actual Fitness: " << printPair(bestIdx) << endl;
    render(run, outLoc);
    printGraph(run, outLoc);
    pop[b].print(aus);
    return printPair(bestIdx);
}

void render(int run, char *outLoc) {
    char fn[60];
    sprintf(fn, "%sDungeon%02d.eps", outLoc, run);
    psDoc pic(fn, 0, 0, BB * Z + 2, BB * Z + 2);
    pic.setcol(225, 225, 225);
    for (int i = 0; i < Z + 1; i++) {
        pic.startpath();
        pic.moveto(0, i * BB);
        pic.lineto(Z * BB, i * BB);
        pic.stroke();
        pic.startpath();
        pic.moveto(i * BB, 0);
        pic.lineto(i * BB, Z * BB);
        pic.stroke();
    }

    for (int i = 0; i < Rz; i++) {//loop over rooms
        if ((Dx[i] == 1) || (Dy[i] == 1))pic.setcol(0, 0, 180); //blue corridors
        else pic.setcol(180, 180, 180); //light gray rooms
        if (i == 0)pic.setcol(180, 0, 0); //first room
        pic.startpath();
        pic.moveto(Rx[i] * BB, Ry[i] * BB);
        pic.lineto((Rx[i] + Dx[i]) * BB, Ry[i] * BB);
        pic.lineto((Rx[i] + Dx[i]) * BB, (Ry[i] + Dy[i]) * BB);
        pic.lineto(Rx[i] * BB, (Ry[i] + Dy[i]) * BB);
        pic.closepath();
        pic.fill();
    }

    pic.setcol(0, 0, 0);
    for (int i = 0; i < Rz; i++) { //bound the rooms
        if (Dx[i] > 1 && Dy[i] > 1) {
            pic.startpath();
            pic.moveto(Rx[i] * BB, Ry[i] * BB);
            pic.lineto((Rx[i] + Dx[i]) * BB, Ry[i] * BB);
            pic.lineto((Rx[i] + Dx[i]) * BB, (Ry[i] + Dy[i]) * BB);
            pic.lineto(Rx[i] * BB, (Ry[i] + Dy[i]) * BB);
            pic.closepath();
            pic.stroke();
        }
    }

    pic.setcol(0, 255, 0);
    int nnbrs = 0;
    int nbrs[Rz];
    pair<pair<int, int>, pair<int, int>> points;
    for (int from = 0; from < Rz; from++) {
        nnbrs = G->Nbrs(from, nbrs);
        for (int to = 0; to < nnbrs; to++) {
            points = getOverlap(from, nbrs[to]);
            pic.startpath();
            if (points.first.first == points.second.first) {
                switch (abs(points.first.second - points.second.second)) {
                    case 1:
                        pic.moveto(points.first.first * BB,
                                   points.first.second * BB
                                       + 2);
                        pic.lineto(points.second.first * BB,
                                   points.second.second *
                                       BB - 2);
                        break;
                    case 2:
                        pic.moveto(points.first.first * BB,
                                   points.first.second * BB
                                       + 5);
                        pic.lineto(points.second.first * BB,
                                   points.second.second *
                                       BB - 5);
                        break;
                    case 3:
                        pic.moveto(points.first.first * BB,
                                   points.first.second * BB
                                       + 10);
                        pic.lineto(points.second.first * BB,
                                   points.second.second *
                                       BB - 10);
                        break;
                    case 4:
                        pic.moveto(points.first.first * BB,
                                   points.first.second * BB
                                       + 15);
                        pic.lineto(points.second.first * BB,
                                   points.second.second *
                                       BB - 15);
                        break;
                    default:
                        break;
                }
            } else {
                switch (abs(points.first.first - points.second.first)) {
                    case 1:
                        pic.moveto(points.first.first * BB + 2,
                                   points.first.second *
                                       BB);
                        pic.lineto(points.second.first * BB - 2,
                                   points.second.second *
                                       BB);
                        break;
                    case 2:
                        pic.moveto(points.first.first * BB + 5,
                                   points.first.second *
                                       BB);
                        pic.lineto(points.second.first * BB - 5,
                                   points.second.second *
                                       BB);
                        break;
                    case 3:
                        pic.moveto(points.first.first * BB + 10,
                                   points.first.second *
                                       BB);
                        pic.lineto(points.second.first * BB - 10,
                                   points.second.second *
                                       BB);
                        break;
                    case 4:
                        pic.moveto(points.first.first * BB + 15,
                                   points.first.second *
                                       BB);
                        pic.lineto(points.second.first * BB - 15,
                                   points.second.second *
                                       BB);
                        break;
                }
            }
            pic.stroke();
        }
    }
    cout << "Render Completed." << endl;
}

pair<pair<int, int>, pair<int, int>> getOverlap(int r1, int r2) {
    int r1x1, r1x2, r1y1, r1y2, r2x1, r2x2, r2y1, r2y2;
    r1x1 = Rx[r1];
    r1x2 = Rx[r1] + Dx[r1];
    r1y1 = Ry[r1];
    r1y2 = Ry[r1] + Dy[r1];
    r2x1 = Rx[r2];
    r2x2 = Rx[r2] + Dx[r2];
    r2y1 = Ry[r2];
    r2y2 = Ry[r2] + Dy[r2];
    int sx, sy, ex, ey;
    pair<int, int> start;
    pair<int, int> end;
    if (r1x1 == r2x2) { // R2 to the left of R1
        sx = r1x1;
        ex = r1x1;
        sy = max(r1y1, r2y1);
        ey = min(r1y2, r2y2);
    } else if (r1x2 == r2x1) { // R2 to the right of R1
        sx = r1x2;
        ex = r1x2;
        sy = max(r1y1, r2y1);
        ey = min(r1y2, r2y2);
    } else if (r1y1 == r2y2) { // R2 above R1
        sx = max(r1x1, r2x1);
        ex = min(r1x2, r2x2);
        sy = r1y1;
        ey = r1y1;
    } else if (r1y2 == r2y1) { // R2 below R1
        sx = max(r1x1, r2x1);
        ex = min(r1x2, r2x2);
        sy = r1y2;
        ey = r1y2;
    }
    start = make_pair(sx, sy);
    end = make_pair(ex, ey);
    pair<pair<int, int>, pair<int, int>> rtn = make_pair(start, end);
    return rtn;
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

string printPair(int idx) {
    string rtn = "[";
    rtn += to_string(dfit[idx]);
    rtn += ", ";
    rtn += to_string(efit[idx]);
    rtn += "]";
    return rtn;
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

pair<int, int> dsOfGFrA(bitspray &a) {
    vector<pair<int, int>> el;
    G->empty(Rz);
    G->copy(*initG);
    developQ(a);

    el = getDoors();
    for (pair<int, int> p : el) {
        G->add(p.first, p.second);
    }

    int apx, act;
    apx = approxDiameter(*G);
    act = diameter(*G);
    return make_pair(apx, act);
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

void approxTest() {
    int bestApx = 0;
    int bestAct = 0;
    int diffs = 0;
    cout << "Beginning Approx Test" << endl;
    cout << "B?\t" << "Ds\t" << "App\t" << "Act" << endl;
//    for (int i = 0; i < popsize; i++) {
//        if (fit[i].first == bestD) {
//            cout << "B\t";
//            bestApx++;
//        } else {
//            cout << '\t';
//        }
//        cout << fit[i].second << '\t';
//        pair<int, int> diams = dsOfGFrA(pop[i]);
//        if (diams.first != diams.second) {
//            diffs++;
//        }
//        if (fit[i].second == bestD && diams.first == diams.second) {
//            bestAct++;
//        }
//        cout << diams.first << '\t' << diams.second << endl;
//    }
    cout << "Errors: " << diffs << " " << (float) diffs * 100 / popsize << "%"
         << endl;
    cout << "Approx. Best: " << bestApx << " " << (float) bestApx * 100 /
        popsize <<
         "%" << endl;
    cout << "Actual Best: " << bestAct << " " << (float) bestAct * 100 /
        (float) bestApx << "%" << endl;
    cout << "End of Approx Test" << endl;
}