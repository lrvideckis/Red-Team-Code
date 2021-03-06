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

bool onSegmentPossible(const point& p, const segment& s){
    bool x = (abs(s.p1.x - s.p2.x) < EPSILON && abs(p.x - s.p2.x) < EPSILON) || (p.x <= max(s.p1.x, s.p2.x) && p.x >= min(s.p1.x, s.p2.x));
    bool y = (abs(s.p1.y - s.p2.y) < EPSILON && abs(p.y - s.p2.y) < EPSILON) || (p.y <= max(s.p1.y, s.p2.y) && p.y >= min(s.p1.y, s.p2.y));
    return x && y;
}

//If they do not intersect, the result is an empty vector
//If they intersect at exactly 1 point, the result contains that point
//If they overlap for non-0 distance, the left and right points of that intersection are returned
vector<point> intersect(const segment& s1, const segment& s2){
    point a = s1.p1, b = s1.p2, c = s2.p1, d = s2.p2;

    if(orientation(a, b, c) == Colinear && orientation(a, b, d) == Colinear &&
        orientation(c, d, a) == Colinear && orientation(c, d, b) == Colinear){
        point min_s1 = min(a, b), max_s1 = max(a, b);
        point min_s2 = min(c, d), max_s2 = max(c, d);
        if(max_s1 < min_s2 || max_s2 < min_s1) return {};

        point start = max(min_s1, min_s2), end = min(max_s1, max_s2);
        if(start == end) return {start};
        else return {min(start, end), max(start, end)};
    }

    long double a1 = b.y - a.y, a2 = d.y - c.y;
    long double b1 = a.x - b.x, b2 = c.x - d.x;
    long double c1 = a1*a.x + b1*a.y, c2 = a2*c.x + b2*c.y;
    long double det = a1*b2 - a2*b1;
    if(abs(det) > EPSILON){
        point inter((b2*c1 - b1*c2)/det, (a1*c2 - a2*c1)/det);
        if(onSegmentPossible(inter, s1) && onSegmentPossible(inter, s2)) return {inter};
    }
    return {};
}

double sqmag(const point& p1) { return p1.x*p1.x + p1.y*p1.y; }
double mag(const point& p1) { return sqrt(sqmag(p1)); }
double sproject(const point& a, const point& b) { return dot(point(0, 0), a, b)/mag(b); }
point vproject(const point& a, const point& b) { return b * sproject(a, b) / mag(b); }

bool straddle(const segment& s1, const segment& s2) {
    long double cross1 = cross(s1.p1, s1.p2, s2.p1);
    long double cross2 = cross(s1.p1, s1.p2, s2.p2);

    if((cross1 > 0 && cross2 > 0) || (cross1 < 0 && cross2 < 0)) return false;
    if(abs(cross1) < EPSILON && abs(cross2) < EPSILON && orientation(s1.p2, s2.p1, s2.p2) != Colinear) return false;
    return true;
}

long double linePointDist(const segment& s, const point& p, bool isSegment=false) {
    double segLen = sqmag(s.p1 - s.p2);
    if(segLen < EPSILON) return mag(p - s.p1);
    double proj = dot(s.p1, s.p2, p) / segLen;
    if(isSegment) { proj = min(1.0, max(0.0, proj)); }

    point projected = s.p1 + (s.p2 - s.p1) * proj;
    return mag(p - projected);
}

long double polyArea(const vector<point>& points){
    long double result = 0;
    for(int i=0, j=1; i<points.size(); i++, j=(j+1)%points.size()) {
        result += points[i].x * points[j].y;
        result -= points[i].y * points[j].x;
    }
    return result/2;
}

//Returns Inside, Outside, or OnEdge
int pointInPoly(const vector<point>& poly, const point& p){
    bool inside = false;

    long double maxX = numeric_limits<long double>::lowest();
    for(const point& p : poly) maxX = max(maxX, p.x);

    point outside(maxX+1, p.y);
    vector<point> intersection;

    for(int i=0, j = poly.size()-1; i < poly.size(); i++, j=i-1) { 
        if(p == poly[i] || p == poly[j]) return OnEdge;
        if(orientation(p, poly[i], poly[j]) == Colinear &&
           onSegmentPossible(p, segment{poly[i], poly[j]})) return OnEdge;

        intersection = intersect(segment{p, outside}, segment{poly[i], poly[j]}); 
        if(intersection.size() == 1) {
            if(poly[i] == intersection[0] && poly[j].y <= p.y) continue;
            if(poly[j] == intersection[0] && poly[i].y <= p.y) continue;
            inside = !inside;
        }
    }   
    return (inside ? Inside : Outside);    
}

vector<point> convexHull(vector<point> points, bool keepColinear = false) {
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

  size_t i = firstColinearEnd-1;
  size_t m = 0;
  if(!keepColinear) { points[++m] = points[i++]; }
  else { m = i++; }

  points[++m] = points[i++];

  for(; i < lastColinearStart+1; ++i) {
    auto orient = orientation(points[m-1], points[m], points[i]);

    while(orient == CW || (!keepColinear && orient == Colinear)) {
        m--;
        orient = orientation(points[m-1], points[m], points[i]);
    }

    points[++m] = points[i];
  }

  //Optionally add all of colinear chunk at end
  if(keepColinear) {
      while(i < points.size())
        points[++m] = points[i++];
  }

  points.resize(m+1);
  return points;
}