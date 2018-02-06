#include <bits/stdc++.h>
using namespace std;
typedef long long ll;

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
       
    c.clear();
    for(auto &x : answer) {
        c.push_back(round(x.real()));
    }
    c.resize(a.size() + b.size() - 1);
}

int32_t main() {ios::sync_with_stdio(false);cin.tie(0);cout.tie(0);
    int n;
    cin >> n;
    vector<int> pos(50001,0), neg(50001,0);
    int zeroCount = 0, posCount = 0, negCount = 0;
    vector<int> arr;
    while(n--) {
        int temp;
        cin >> temp;
        arr.push_back(temp);
        if(temp == 0) {
            zeroCount++;
        } else if(temp > 0) {
            posCount++;
            pos[temp]++;
        } else {
            negCount++;
            neg[-temp]++;
        }
    }
    int counter = zeroCount*(zeroCount-1)*(zeroCount-2);//0 + 0 = 0
    int i = 0;
    for(int x : pos) {
        counter += 2*x*(x-1)*zeroCount;//pos + 0 = pos and 0 + pos = pos
        if(i+i < pos.size()) {
            counter -= pos[i+i] * x * x;
            counter += pos[i+i] * x * (x-1);
        }
        i++;
    }
    i=0;
    for(int x : neg) {
        counter += 2*x*(x-1)*zeroCount;//neg + 0 = neg and 0 + neg = neg
        if(i+i < neg.size()) {
            counter -= neg[i+i] * x * x;
            counter += neg[i+i] * x * (x-1);
        }
        i++;
    }
    vector<int> product;
    multFFT(pos, pos, product);
    i = 0;
    for(int i = 0; i < pos.size(); ++i) {
        counter += pos[i] * product[i];
    }
    multFFT(neg, neg, product);
    for(int i = 0; i < neg.size(); ++i) {
        counter += product[i] * neg[i];
    }
    if(zeroCount > 0) {
        for(int i = 0; i < pos.size(); ++i) {
            counter += zeroCount * pos[i]*neg[i]*2;//pos + neg = 0 && neg + pos = 0 
        }
    }
    //pos + neg = pos && pos + neg = neg:
    int posNegCount = 0;
    multFFT(pos,neg,product);
    for(int i = 0; i < neg.size(); ++i) {
        posNegCount += product[i] * neg[i];
        posNegCount += product[i] * pos[i];
    }
    multFFT(neg,pos,product);
    for(int i = 0; i < neg.size(); ++i) {
        posNegCount += product[i] * neg[i];
        posNegCount += product[i] * pos[i];
    }
    cout << counter+posNegCount << '\n';
    return 0;
}










