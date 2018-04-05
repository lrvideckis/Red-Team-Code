// #include <bits/stdc++.h>
#include <vector>
#include <iostream>
#include <complex>
#include <algorithm>
using namespace std;
typedef long long ll;

int maxDist = 0, maxLength = 0;

void mult(vector<int> &a, vector<int> &b, vector<int> &c) {
    c.clear();
    for(int degree = 0; degree <= maxLength; ++degree) {
        int sum = 0;
        for(int i = 0; i <= degree; ++i) {
            sum += a[degree-i] * b[i];
        }
        c.push_back(sum);
    }
}

typedef complex<double> cd;
typedef vector<cd> vcd;

vcd fft(const vcd &as) {
    int n = as.size();
    int k = 0;
    while ((1 << k) < n) k++;
    vector<int> rev(n);
    rev[0] = 0;
    int high1 = -1;
    for (int i = 1; i < n; i++) {
        if ((i & (i - 1)) == 0)
            high1++;
        rev[i] = rev[i ^ (1 << high1)];
        rev[i] |= (1 << (k - high1 - 1));
    }
    vcd roots(n);
    for (int i = 0; i < n; i++) {
        double alpha = 2 * M_PI * i / n;
        roots[i] = cd(cos(alpha), sin(alpha));
    }

    vcd cur(n);
    for (int i = 0; i < n; i++)
        cur[i] = as[rev[i]];

    for (int len = 1; len < n; len <<= 1) {
        vcd ncur(n);
        int rstep = roots.size() / (len * 2);
        for (int pdest = 0; pdest < n;) {
            int p1 = pdest;
            for (int i = 0; i < len; i++) {
                cd val = roots[i * rstep] * cur[p1 + len];
                ncur[pdest] = cur[p1] + val;
                ncur[pdest + len] = cur[p1] - val;
                pdest++, p1++;
            }
            pdest += len;
        }
        cur.swap(ncur);
    }
    return cur;
}

void reverseFFT(int n, vector<complex<double> > &coef, vector<complex<double> > &roots) {
    roots.clear();
    roots = fft(coef);
    reverse(roots.begin()+1,roots.end());
    for(complex<double> &r : roots) r = r.real() / n;
}

void multFFT(vector<int> &a, vector<int> &b, vector<int> &c) {
    int n = a.size() + b.size() - 1;
    n = (ll)(ceil(log2(n)));
    n = (ll)pow(2, n);
    vector<complex<double> > coefA, coefB, rootsA, rootsB;
    for(int i = 0; i < n; ++i) {
        if(i < a.size()) {
            coefA.push_back(complex<double>(a[i],0));
        } else {
            coefA.push_back(complex<double>(0,0));
        }
        if(i < b.size()) {
            coefB.push_back(complex<double>(b[i],0));
        } else {
            coefB.push_back(complex<double>(0,0));
        }
    }
    rootsA = fft(coefA);
    rootsB = fft(coefB);

    vector<complex<double> > prod(n);
    for(int i = 0; i < n; ++i) {
        prod.at(i) = rootsA[i] * rootsB[i];
    }
    vector<complex<double> > answer;
       
    reverseFFT(n,prod,answer);
       
    //c.clear();
    for(auto &x : answer) {
        c.push_back(round(x.real()));
    }
    c.resize(a.size() + b.size() - 1);
}
