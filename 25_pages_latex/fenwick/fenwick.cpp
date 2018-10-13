#include <vector>
using namespace std;

struct Fenwick {
  vector<long long> tree;

  //Constructs the tree initialized to all 0
  Fenwick(int n) 
  { tree.resize(n, 0); }

  //Constructs the tree initialized to the given values
  Fenwick(const vector<long long>& val)
  {
    tree.resize(val.size(), 0);
    for(int i=0; i<val.size(); i++)
      increase(i, val[i]);
  }

  //Increase item i by delta d
  void increase(int i, long long d) 
  { for(; i < tree.size(); i |= i+1) tree[i] += d; }

  //get sum of the range [0, i)
  //sum(0) = 0
  long long sum(int i)
  {
    long long s = 0;
    for(; i > 0; i &= i-1) s += tree[i-1];
    return s;
  }

  //get sum of the range [left, right)
  long long sum(int left, int right)
  { return sum(right) - sum(left); }
};
