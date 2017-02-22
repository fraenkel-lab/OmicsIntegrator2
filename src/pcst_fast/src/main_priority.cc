#include <iostream>

#include "priority_queue.h"

using cluster_approx::PriorityQueue;

using namespace std;

int main() {
  PriorityQueue<double, int> queue;
  queue.insert(10.0, 1);
  queue.insert(20.0, 2);
  queue.insert(0.0, 0);
  queue.insert(30.0, 3);

  double val;
  int index;

  queue.delete_min(&val, &index);
  cout << val << " " << index << endl;
  queue.delete_min(&val, &index);
  cout << val << " " << index << endl;

  queue.decrease_key(5.0, 3);

  queue.delete_min(&val, &index);
  cout << val << " " << index << endl;

  queue.insert(40.0, 4);
  queue.delete_element(2);

  queue.delete_min(&val, &index);
  cout << val << " " << index << endl;

  return 0;
}
