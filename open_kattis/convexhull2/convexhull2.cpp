#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <map>
#include <set>
#include <limits>
#include <stack>
#include <random>

using namespace std;

constexpr int Colinear = -1, NoIntersect = 0, Intersect = 1;
constexpr int CW = 2, CCW = 3;
constexpr int Inside = 4, Outside = 5, OnEdge = 6;

const double EPSILON = 0.0001;

struct point{
    long double x, y;
    point(long double x_=0, long double y_=0) : x(x_), y(y_){}

    // Only < operator is unusual behavior
    // the rest can probably all be memorized
    bool operator <(const point& other) const
    {return (x < other.x ? true : (x == other.x && y < other.y));}

        bool operator == (const point& other) const
    {
        return abs(other.x - x) < EPSILON && abs(other.y - y) < EPSILON;
    }

    point operator *(const long double& d) const
    {
        return point(x*d, y*d);
    }

    point operator /(const long double& d) const
    {
        return point(x/d, y/d);
    }

    //Add other operators as needed
    point operator-(const point& other) const {
    	return point(x - other.x, y - other.y);
    }
    point operator+(const point& other) const {
    	return point(x + other.x, y + other.y);
    }
};

struct segment {point p1, p2;};

//Dot product ab.ac
double dot(const point& a, const point& b, const point& c){
    point AB = b - a;
    point AC = c - a;
    return AB.x*AC.x + AB.y*AC.y;
}

//Cross product AB X AC
double cross(const point& A, const point& B, const point& C){
    point AB = B - A, AC = C - A;
    return(AB.x * AC.y - AB.y * AC.x);
}

int orientation(const point& p, const point& q, const point& r){
    double val = cross(p, q, r);
    if(abs(val) < EPSILON) return Colinear;
    return (val > 0) ? CCW : CW;
}

double sqmag(const point& p1) { return p1.x*p1.x + p1.y*p1.y; }
double mag(const point& p1) { return sqrt(sqmag(p1)); }

vector<point> sortPoints(vector<point> points, bool keepColinear = false) {
  if(points.size() < 2)
    return points;

  if(points.size() == 2) {
    if(!keepColinear && points[0] == points[1]) { return {points[0]}; }
    else { return points; }
  }

  point lowestPoint = points[0];

  for(int i=0; i<points.size(); i++)
    if(points[i].y < lowestPoint.y ||
      (abs(points[i].y - lowestPoint.y) < EPSILON && points[i].x < lowestPoint.x))
        lowestPoint = points[i];

  point horiz = lowestPoint + point(1, 0);
  sort(points.begin(), points.end(),
  [&](const point& l, const point& r) {
    if(r == lowestPoint || r == l) return false;
    if(l == lowestPoint) return true;

    long double scoreL = dot(lowestPoint, horiz, l) / mag(l-lowestPoint);
    long double scoreR = dot(lowestPoint, horiz, r) / mag(r-lowestPoint);

    return scoreL > scoreR;
  });

  size_t firstDiff = 1;
  while(firstDiff < points.size() && points[0] == points[firstDiff]) firstDiff++;

  if(firstDiff == points.size()) {
      return keepColinear ? points : vector<point>{points[0]};
  }

  size_t firstColinearEnd = firstDiff+1;
  while(firstColinearEnd < points.size() && orientation(points[0], points[firstDiff], points[firstColinearEnd]) == Colinear)
    firstColinearEnd++;

  sort(points.begin() + 1, points.begin() + firstColinearEnd, 
       [&](const point& l, const point& r) {return mag(l - points[0]) < mag(r - points[0]);});

  if(firstColinearEnd == points.size()) {
      return keepColinear ? points : vector<point>{points[0], points[points.size()-1]};
  }

  firstDiff = points.size() - 1;
  while(firstDiff >= firstColinearEnd && points[0] == points[firstDiff]) firstDiff--;

  size_t lastColinearStart = firstDiff - 1;
  while(lastColinearStart >= firstColinearEnd && orientation(points[0], points[firstDiff], points[lastColinearStart]) == Colinear)
    lastColinearStart--;

  lastColinearStart++;
  sort(points.begin() + lastColinearStart, points.end(), 
       [&](const point& l, const point& r) {return mag(l - points[0]) > mag(r - points[0]);}); 

  return points;
}

int main() {
	int n;
	cin >> n;
	string ignore;
	vector<point> pts;
	for(int i = 0; i < n; i++) {
		point pt;
		cin >> pt.x >> pt.y >> ignore;
		if(ignore == "Y") {
			pts.push_back(pt);
		}
	}

    vector<point> hull = sortPoints(pts);

	int mpi = 0;
	for(int i = 0; i < hull.size(); i++) {
    if(hull[i].x < hull[mpi].x ||
      (abs(hull[i].x - hull[mpi].x) < EPSILON && hull[i].y < hull[mpi].y))
      mpi = i;
	}
	cout << fixed << setprecision(0) << hull.size() << endl;
	for(int i = mpi; i < hull.size(); i++) {
		cout << hull[i].x << " " << hull[i].y << endl;
	}
	for(int i = 0; i < mpi; i++) {
		cout << hull[i].x << " " << hull[i].y << endl;
	}
}