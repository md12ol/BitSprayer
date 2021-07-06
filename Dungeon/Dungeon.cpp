#include "../BitSprayer/stat.h"
#include "../BitSprayer/bitspray.h"
#include "ps.h"
#include "setu.h"
#include "Dungeon.h"

using namespace std;

int main() {
    fstream stat;   //statistics reporting stream
    fstream best;   //best reporting stream
    char fn[100];
    char *outLoc = new char[45];
    char *outRoot = new char[20];
    sprintf(outRoot, "./Output/");
    double sum = 0.0;   //sum of best fitness from each run
    G = new graph(Rz);
    sprintf(outLoc, "%sOutput - %02dS, %03dP, %dM/",
            outRoot, states, popsize, MNM);
    std::filesystem::create_directory(outLoc);
    sprintf(fn, "%sGraphs/", outLoc);
    std::filesystem::create_directory(fn);
    sprintf(fn, "%sBoards/", outLoc);
    std::filesystem::create_directory(fn);
    sprintf(fn, "%sDoors/", outLoc);
    std::filesystem::create_directory(fn);
    sprintf(fn, "%sDungeons/", outLoc);
    std::filesystem::create_directory(fn);


    initalg();
    sprintf(fn, "%sbest.sda", outLoc);
    best.open(fn, ios::out);
    for (int run = 0; run < runs; run++) {
        sprintf(fn, "%srun%02d.dat", outLoc, run);
        stat.open(fn, ios::out);
        if (verbose == 1) cmdLineRun(run, cout);
        initpop();
        report(stat);
        for (int mev = 0; mev < mevs; mev++) {
            matingevent();
            if ((mev + 1) % RI == 0) {
                if (verbose == 1) {
                    cout << left << setw(5) << run;
                    cout << left << setw(4) << (mev + 1) / RI;
                }
                report(stat);
            }
        }
        stat.close();
        sum += reportbest(best, run, outLoc);
    }
    cout << "Mean fitness: " << sum / runs << endl;
    best.close();
    delete G;
    return (0);
}

//This routine fills the queue
void developQ(bitspray &A) {//unpack the queue
    int h, t;  //head and tail of queue

    for (t = 0; t < Qz; t++)Q[t] = 0;        //clear the queue
    A.reset(Q, h, t);                 //reset the self driving automata
    while (t < Qz - 2)A.next(Q, h, t, Qz);  //run the automata
}

void initalg() {
    srand48(RNS);
    for (int i = 0; i < popsize; i++) {
        pop[i] = new bitspray(states);
    }
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

void initpop() {
    for (int i = 0; i < popsize; i++) {
        pop[i]->randomize();
        fit[i] = fitness(*pop[i]);
        dx[i] = i;
    }
}

void clearboard() {
    for (int i = 0; i < Z; i++) {
        for (int j = 0; j < Z; j++) {
            Bd[i][j] = 0;
        }
    }
    G->empty(Rz);
}

void fill(int a, int b, int da, int db) {
    for (int i = a; i < a + da; i++) {
        for (int j = b; j < b + db; j++) {
            Bd[i][j] = 1;
        }
    }
}

bool available(int a, int b, int da, int db) {
    if (a < 0 || b < 0 || a + da > Z || b + db > Z) return false; //off board
    for (int i = a; i < a + da; i++) {
        for (int j = b; j < b + db; j++) {
            if (Bd[i][j] == 1) return false; //room full
        }
    }
    return true;
}

bool getnum(int &val, int bits, int &psn) {//get a number from the queue
    if (psn + bits >= Qz) return false; //safety first
    val = 0;  //zero the value
    for (int i = 0; i < bits; i++)
        val = 2 * val + Q[psn++];  //assemble the integer
    return true;
}

double fitness(bitspray &A) {
    int psn;    //position within queue bit buffer
    int rmn;    //current active room number
    int side;   //which side? 4
    int typ;    //room type 0=corridor, 1-7=room
    int slide;  //How far down the wall, 16
    int ttlz;   //total size
    int a, b, da, db;  //New room descriptor - ULC and side lengths
    int ttl;      //filled space of board
    double fitv;  //fitness

    clearboard();
    developQ(A);
    Rx[0] = Z / 2 - 2;  //x-value of initial room
    Ry[0] = Z / 2 - 2;  //y-value
    Dx[0] = 4;          //width
    Dy[0] = 4;          //height
    fill(Rx[0], Ry[0], Dx[0], Dy[0]);
    Rr.push_back(0); //RECENT
    cnr = 1;
    psn = 0;
    while (cnr < Rz && getnum(rmn, Rbs, psn)) {
//    rmn = rmn % cnr;  //get an already created room
        //RECENT
        if (cnr > RecRs) {
            int extra = cnr - (int) RecRs;
            rmn = extra + (rmn % (int) RecRs);
        } else {
            rmn = rmn % cnr;
        }
        if (!getnum(side, 2, psn)) break;
        if (!getnum(typ, 3, psn)) break;
        if (!getnum(slide, 4, psn)) break;
        if (!getnum(ttlz, 4, psn)) break;

        if (typ == 0) {  //hallway
            if (side < 2) { //horizontal hallway
                da = ttlz + 1;
                db = 1;
            } else {  //vertical hallway
                db = ttlz + 1;
                da = 1;
            }
        } else {  //room
            da = ttlz % 4 + 1;
            db = (ttlz / 4) % 4 + 1;
        }

        switch (side) {
            case 0: //right side
                a = Rx[rmn] + Dx[rmn];
                b = Ry[rmn] + slide % Dy[rmn];
                break;
            case 1: //left side
                a = Rx[rmn] - da;
                b = Ry[rmn] + slide % Dy[rmn];
                break;
            case 2: //bottom side
                b = Ry[rmn] + Dy[rmn];
                a = Rx[rmn] + slide % Dx[rmn];
                break;
            case 3: //top side
                b = Ry[rmn] - db;
                a = Rx[rmn] + slide % Dx[rmn];
                break;
            default:
                break;
        } //end switch00

        if (available(a, b, da, db)) {  //room available
            fill(a, b, da, db);
            Rx[cnr] = a;
            Ry[cnr] = b;
            Dx[cnr] = da;
            Dy[cnr] = db;
            G->add(rmn, cnr);
            cnr++;
        }

    }
    //calculate the fitness
    a = Z;
    b = Z;
    da = 0;
    db = 0;
    ttl = 0;
    for (int i = 0; i < Z; i++) {
        for (int j = 0; j < Z; j++) {
            if (Bd[i][j] == 1) { //full square
                if (i < a)a = i;  //left side
                if (j < b)b = j;  //top side
                if (i > da)da = i;  //right side
                if (j > db) db = j; //bottom side
                ttl++;
            }
        }
    }
    fitv = ((double) ttl);
    fitv *= ttl;
    fitv /= ((double) (da - a));
    fitv /= ((double) (db - b));

    return (fitv);
}

void matingevent() {
    int nm;

    tselect(fit, dx, tsize, popsize);
    pop[dx[0]]->copy(*pop[dx[tsize - 2]]);
    pop[dx[1]]->copy(*pop[dx[tsize - 1]]);
    pop[dx[0]]->tpc(*pop[dx[1]]);
    nm = (int) lrand48() % MNM + 1;
    pop[dx[0]]->mutate(nm);
    nm = (int) lrand48() % MNM + 1;
    pop[dx[1]]->mutate(nm);
    fit[dx[0]] = fitness(*pop[dx[0]]);
    fit[dx[1]] = fitness(*pop[dx[1]]);
}

void report(ostream &aus) {
    dset D;
    D.add(fit, popsize);
    aus << left << setw(10) << D.Rmu();
    aus << left << setw(12) << D.RCI95();
    aus << left << setw(10) << D.Rsg();
    aus << left << setw(8) << D.Rmax();
    aus << endl;
    if (verbose == 1) {
        cout << left << setw(10) << D.Rmu();
        cout << left << setw(12) << D.RCI95();
        cout << left << setw(10) << D.Rsg();
        cout << left << setw(8) << D.Rmax();
        cout << endl;
    }
}

void render(int run, char *outLoc) {
    char fn[60];
    sprintf(fn, "%sDungeons/Dungeon%02d.eps", outLoc, run);
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

    for (int i = 0; i < cnr; i++) {//loop over rooms
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
    for (int i = 0; i < cnr; i++) { //bound the rooms
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
}

double reportbest(ostream &aus, int run, char *outLoc) {
    int b = 0;
    for (int i = 1; i < popsize; i++) {
        if (fit[i] > fit[b]) b = i;
    }
    fit[b] = fitness(*pop[b]);
    aus << "Run" << run << " fitness: " << fit[b] << endl;
    if (verbose == 1) cout << "Best Fitness: " << fit[b] << endl;
    render(run, outLoc);
    printGraph(run, outLoc);
    printBoard(run, outLoc);
    printDoors(run, outLoc);
    pop[b]->print(aus);
    return fit[b];
}

void printGraph(int run, char *outLoc) {
    fstream graph;
    char fn[100];

    sprintf(fn, "%sGraphs/graph%02d.dat", outLoc, run);
    graph.open(fn, ios::out);
    G->resetV();
    G->write(graph);
    cout << "Graph printed." << endl;
    graph.close();
}

void printBoard(int run, char *outLoc) {
    fstream board;
    char fn[100];

    sprintf(fn, "%sBoards/board%02d.dat", outLoc, run);
    board.open(fn, ios::out);
    for (int i = 0; i < Z; i++) {
        for (int j = 0; j < Z; j++) {
            board << Bd[i][j] << " ";
        }
        board << endl;
    }
    cout << "Board printed." << endl;
    board.close();
}

void printDoors(int run, char *outLoc) {
    graph *D = new graph(Rz);
    D->empty(cnr);
    int Bx1, Bx2, By1, By2, Orm; //boarder xs and ys

    for (int rm = 0; rm < cnr; rm++) { //for each room
        Bx1 = Rx[rm];
        Bx2 = Rx[rm] + Dx[rm];
        By1 = Ry[rm];
        By2 = Ry[rm] + Dy[rm];
        //for each cell in the room
        for (int x = Bx1; x < Bx2; x++) {
            for (int y = By1; y < By2; y++) {
                if (x == Bx1 && x > 0 && Bd[x - 1][y] == 1) { //to the left
                    Orm = getRoomNm(x - 1, y);
                    if (Orm != -1 && !G->edgeP(rm, Orm)) {
                        D->add(rm, Orm);
                    }
                }
                if (x == Bx2 && x < Z - 1 && Bd[x + 1][y] == 1) { //to the right
                    Orm = getRoomNm(x + 1, y);
                    if (Orm != -1 && !G->edgeP(rm, Orm)) {
                        D->add(rm, Orm);
                    }
                }
                if (y == By1 && y > 0 && Bd[x][y - 1] == 1) { //above
                    Orm = getRoomNm(x, y - 1);
                    if (Orm != -1 && !G->edgeP(rm, Orm)) {
                        D->add(rm, Orm);
                    }
                }
                if (y == By2 && y < Z - 1 && Bd[x][y + 1] == 1) { //below
                    Orm = getRoomNm(x, y + 1);
                    if (Orm != -1 && !G->edgeP(rm, Orm)) {
                        D->add(rm, Orm);
                    }
                }
            }
        }
    }
    fstream doors;
    char fn[100];
    sprintf(fn, "%sDoors/doors%02d.dat", outLoc, run);
    doors.open(fn, ios::out);
    D->write(doors);
    cout << "Doors printed." << endl;
    doors.close();
}

int getRoomNm(int x, int y) {
    int Bx1, Bx2, By1, By2; //boarder xs and ys

    for (int rm = 0; rm < cnr; rm++) { //for each room
        Bx1 = Rx[rm];
        Bx2 = Rx[rm] + Dx[rm];
        By1 = Ry[rm];
        By2 = Ry[rm] + Dy[rm];
        if (Bx1 <= x && x <= Bx2 && By1 <= y
            && y <= By2) { //if x, y in room cells
            return rm;
        }
    }
    return -1;
}