#include <cstdio>
#include <vector>

#include "pairing_heap.h"

using cluster_approx::PairingHeap;
using std::vector;

int main() {
  vector<PairingHeap<double, int>::ItemHandle> buffer;
  PairingHeap<double, int> heap(&buffer);
  PairingHeap<double, int>::ItemHandle handle1 = heap.insert(10.0, 1);
  PairingHeap<double, int>::ItemHandle handle2 = heap.insert(20.0, 2);
  PairingHeap<double, int>::ItemHandle handle0 = heap.insert(0.0, 0);
  PairingHeap<double, int>::ItemHandle handle3 = heap.insert(30.0, 3);

  double val;
  int payload;
  heap.delete_min(&val, &payload);
  printf("%lf %d\n", val, payload);

  heap.delete_min(&val, &payload);
  printf("%lf %d\n", val, payload);

  heap.add_to_heap(5.0);

  heap.decrease_key(handle3, 35.0, 1.0);

  PairingHeap<double, int> heap2(&buffer);
  PairingHeap<double, int>::ItemHandle handle4 = heap.insert(4.0, 4);
  PairingHeap<double, int>::ItemHandle handle5 = heap.insert(50.0, 5);

  PairingHeap<double, int> heap3 = PairingHeap<double, int>::meld(&heap,
                                                                  &heap2);

  heap3.delete_min(&val, &payload);
  printf("%lf %d\n", val, payload);

  heap3.delete_min(&val, &payload);
  printf("%lf %d\n", val, payload);

  heap3.delete_min(&val, &payload);
  printf("%lf %d\n", val, payload);

  heap3.delete_min(&val, &payload);
  printf("%lf %d\n", val, payload);
  return 0;
};
