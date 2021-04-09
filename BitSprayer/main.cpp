//
// Created by micha on 2021-03-27.
//

#include <iostream>
#include "bitspray.h"

using namespace std;

int main() {

  auto *b = new bitspray(5, 5);
  b->randomize();
  b->reset();
//  b->print(cout);
  for (int i = 0; i < 20; ++i) {
//    b->getBuf(cout);
    cout << b->bit() << endl;
  }

}
