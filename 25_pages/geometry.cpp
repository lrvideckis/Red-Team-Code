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

//Constant values to be returned
constexpr int Colinear = -1, NoIntersect = 0, Intersect = 1;
constexpr int CW = 2, CCW = 3;
constexpr int Inside = 4, Outside = 5, OnEdge = 6;

//Epsilon for all double comparisons
const double EPSILON = 0.0001;

struct point
{
    long double x, y;

    point(long double x_=0, long double y_=0) : x(x_), y(y_){}

    // Only < operator is unusual behavior
    // the rest can probably all be memorized
    bool operator <(const point& other) const
    {
        return (x < other.x ? true : (x == other.x && y < other.y));
    }

    bool operator == (const point& other) const
    {
        return abs(other.x - x) < EPSILON && abs(other.y - y) < EPSILON;
    }

    bool operator != (const point& other) const
    {
        return abs(other.x - x) > EPSILON || abs(other.y - y) > EPSILON;
    }

    point operator *(const long double& d) const
    {
        return point(x*d, y*d);
    }

    point operator /(const long double& d) const
    {
        return point(x/d, y/d);
    }

    point operator +(const point& r) const
    {
        return point(x+r.x, y+r.y);
    }

    point operator -(const point& r) const
    {
        return point(x-r.x, y-r.y);
    }

    point operator -() const
    {
        return point(-x, -y);
    }
};

ostream& operator << (ostream& out, const point& p)
{
    return out << "(" << p.x << ", " << p.y << ")";
}

//Container for line segment
struct segment
{
    point p1, p2;
};

//Dot product ab.ac
double dot(const point& a, const point& b, const point& c)
{
    point AB = b - a;
    point AC = c - a;
    return AB.x*AC.x + AB.y*AC.y;
}

//Cross product
//AB X AC
double cross(const point& A, const point& B, const point& C)
{
    point AB = B - A, AC = C - A;
    return(AB.x * AC.y - AB.y * AC.x);
}

//Finds orientation of triplet of points p, q, r
//Returns Colinear, CW, or CCW
int orientation(const point& p, const point& q, const point& r)
{
    double val = cross(p, q, r);
    if(abs(val) < EPSILON) return Colinear;
    return (val > 0) ? CCW : CW;
}

//Checks if point p is possibly on the segment s
//but doesn't guarantee it is
bool onSegment(const point& p, const segment& s)
{
    bool x = (abs(s.p1.x - s.p2.x) < EPSILON && abs(p.x - s.p2.x) < EPSILON) || (p.x <= max(s.p1.x, s.p2.x) && p.x >= min(s.p1.x, s.p2.x));
    bool y = (abs(s.p1.y - s.p2.y) < EPSILON && abs(p.y - s.p2.y) < EPSILON) || (p.y <= max(s.p1.y, s.p2.y) && p.y >= min(s.p1.y, s.p2.y));
    return x && y;
}

//Returns of list of intersection points between segments s1, and s2
//If they do not intersect, the result is an empty vector
//If they intersect at exactly 1 point, the result contains that point
//If they overlap for non-0 distance, the left and right points of that intersection
//  are returned
vector<point> intersect(const segment& s1, const segment& s2)
{
/*
    cout << "Intersect:" << endl;
    cout << s1.p1 << " -> " << s1.p2 << endl;
    cout << s2.p1 << " -> " << s2.p2 << endl;
*/
    point a = s1.p1, b = s1.p2, c = s2.p1, d = s2.p2;

    if(orientation(a, b, c) == Colinear && orientation(a, b, d) == Colinear &&
        orientation(c, d, a) == Colinear && orientation(c, d, b) == Colinear)
    {
        point min_s1 = min(a, b), max_s1 = max(a, b);
        point min_s2 = min(c, d), max_s2 = max(c, d);
/*
        cout << "Colinear" << endl;
        cout << min_s1 << " -> " << max_s1 << endl;
        cout << min_s2 << " -> " << max_s2 << endl;
*/
        if(max_s1 < min_s2 || max_s2 < min_s1) return {};

        point start = max(min_s1, min_s2), end = min(max_s1, max_s2);
        if(start == end)
            return {start};
        else 
            return {min(start, end), max(start, end)};
    }

    long double a1 = b.y - a.y, a2 = d.y - c.y;
    long double b1 = a.x - b.x, b2 = c.x - d.x;
    long double c1 = a1*a.x + b1*a.y, c2 = a2*c.x + b2*c.y;
    long double det = a1*b2 - a2*b1;
    if(abs(det) > EPSILON)
    {
        point inter((b2*c1 - b1*c2)/det, (a1*c2 - a2*c1)/det);
        //cout << "Checking " << inter << " vs segments" << endl;
        //cout << onSegment(inter, s1) << " " << onSegment(inter, s2) << endl;
        if(onSegment(inter, s1) && onSegment(inter, s2))
            return {inter};
    }
    return {};
}

//Squared magnitude of point vector
double sqmag(const point& p1)
{
    return p1.x*p1.x + p1.y*p1.y;
}

//Magnitude of point vector
double mag(const point& p1)
{
    return sqrt(sqmag(p1));
}

//Scalar projection of vector a onto vector b
double sproject(const point& a, const point& b)
{
    return dot(point(0, 0), a, b)/mag(b);
}

//Vector projection of vector a onto vector b
point vproject(const point& a, const point& b)
{
    return b * sproject(a, b) / mag(b);
}

//Checks if two segments straddle each other
bool straddle(const segment& s1, const segment& s2)
{
    long double cross1 = cross(s1.p1, s1.p2, s2.p1);
    long double cross2 = cross(s1.p1, s1.p2, s2.p2);

    if((cross1 > 0 && cross2 > 0) || 
       (cross1 < 0 && cross2 < 0)) return false;

    if(abs(cross1) < EPSILON && abs(cross2) < EPSILON &&
       orientation(s1.p2, s2.p1, s2.p2) != Colinear)
       return false;
    
    return true;
}

//Returns distance from line (or segment) to point
long double linePointDist(const segment& s, const point& p, bool isSegment=false)
{
    double segLen = sqmag(s.p1 - s.p2);
    if(segLen < EPSILON)
        return mag(p - s.p1);

    //Project point onto line
    double proj = dot(s.p1, s.p2, p) / segLen;

    //If segment, bound to range of line 0, 1
    if(isSegment)
    {
        proj = min(1.0, max(0.0, proj));
    }

    point projected = s.p1 + (s.p2 - s.p1) * proj;
    return mag(p - projected);
}

//Returns positive area if points are counterclockwise,
//negative area if clockwise
long double polyArea(const vector<point>& points)
{
    long double result = 0;
    for(int i=0, j=1; i<points.size(); i++, j=(j+1)%points.size())
    {
        result += points[i].x * points[j].y;
        result -= points[i].y * points[j].x;
    }
    return result/2;
}

//Checks if point p is inside the polygon
//Returns Inside, Outside, or OnEdge
int pointInPoly(const vector<point>& poly, const point& p)
{
    //cout << "Point in poly " << p.x << " " << p.y << endl;
    bool inside = false;

    long double maxX = numeric_limits<long double>::lowest();
    for(const point& p : poly)
        maxX = max(maxX, p.x);

    //Create point definitely outside polygon
    point outside(maxX+1, p.y);

    vector<point> intersection;

    for(int i=0, j = poly.size()-1; i < poly.size(); i++, j=i-1)
    { 
        if(p == poly[i] || p == poly[j]) return OnEdge;
        if(orientation(p, poly[i], poly[j]) == Colinear &&
           onSegment(p, segment{poly[i], poly[j]})) return OnEdge;

        intersection = intersect(segment{p, outside}, segment{poly[i], poly[j]});
        //cout << intersection.size() << " intersections with " << poly[i].x << ", " << poly[i].y << " -> " << poly[j].x << ", " << poly[j].y << endl;
        if(intersection.size() == 1)
        {
            if(poly[i] == intersection[0] && poly[j].y <= p.y) continue;
            if(poly[j] == intersection[0] && poly[i].y <= p.y) continue;
            //cout << intersection[0].x << " " << intersection[0].y << endl;
            inside = !inside;
        }
    }   

    //cout << "Is inside? " << inside << endl;
    return (inside ? Inside : Outside);    
}

//Computes the convex hull of a set of points
//Using the graham scan algorithm
vector<point> convexHull(vector<point> points, bool keepColinear = false) {
  //Trivial cases; 0, 1, or 2 points
  if(points.size() < 2)
    return points;

  if(points.size() == 2)
  {
    if(!keepColinear && points[0] == points[1])
    {
        return {points[0]};
    }
    else
    {
        return points;
    }
  }

  point lowestPoint = points[0];

  //Don't just use point < operator because that checks
  //x then y; we need y then x
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

  //Sort the colinear chunks of points at the start and end of the list

  //Find first point going from start+1 that is different from first point
  size_t firstDiff = 1;
  while(firstDiff < points.size() && points[0] == points[firstDiff]) firstDiff++;

  //Check for degenerate case: all points the same
  if(firstDiff == points.size())
  {
      return keepColinear ? points : vector<point>{points[0]};
  }

  //Find first point that is not colinear with point 0 and point firstDiff
  size_t firstColinearEnd = firstDiff+1;
  while(firstColinearEnd < points.size() && orientation(points[0], points[firstDiff], points[firstColinearEnd]) == Colinear)
    firstColinearEnd++;

  //Sort points between 0 and firstColinearEnd by distance from point 0
  sort(points.begin() + 1, points.begin() + firstColinearEnd, 
       [&](const point& l, const point& r) {return mag(l - points[0]) < mag(r - points[0]);});

  //Check for degenerate case: all points colinear
  if(firstColinearEnd == points.size())
  {
      return keepColinear ? points : vector<point>{points[0], points[points.size()-1]};
  }

  //Find first point going backward that is not the same as point 0
  firstDiff = points.size() - 1;
  while(firstDiff >= firstColinearEnd && points[0] == points[firstDiff]) firstDiff--;

  //Find first point going backward that is not colinear with point 0 and point from above
  size_t lastColinearStart = firstDiff - 1;
  while(lastColinearStart >= firstColinearEnd && orientation(points[0], points[firstDiff], points[lastColinearStart]) == Colinear)
    lastColinearStart--;

  //Sort points between there and end by distance from point 0
  lastColinearStart++;
  sort(points.begin() + lastColinearStart, points.end(), 
       [&](const point& l, const point& r) {return mag(l - points[0]) > mag(r - points[0]);}); 

  /*
  cerr << "Hull sort" << endl;
  for(point& p : points)
      cerr << p << " ";
  cerr << endl;
  */

  //Either add all of colinear chunk at front, or just the last point in it
  size_t i = firstColinearEnd-1;
  size_t m = 0;
  if(!keepColinear)
  {
      points[++m] = points[i++];
  }
  else
  {
      m = i++;
  }

  points[++m] = points[i++];

  //Process all points up to the chunk of colinear points at the end using Graham Scan
  for(; i < lastColinearStart+1; ++i) {
    auto orient = orientation(points[m-1], points[m], points[i]);

    while(orient == CW || (!keepColinear && orient == Colinear)) {
        m--;
        orient = orientation(points[m-1], points[m], points[i]);
    }

    points[++m] = points[i];
  }

  //Optionally add all of colinear chunk at end
  if(keepColinear)
  {
      while(i < points.size())
        points[++m] = points[i++];
  }

  points.resize(m+1);
  return points;
}


#include <random>
#include <array>
int main()
{
    mt19937_64 reng;
    uniform_int_distribution<> pt_dist(-10, 10);
    uniform_int_distribution<> size_dist(1, 10);

    segment s{{-7, 3}, {-2.9, -1}};
    vector<tuple<point, bool, double>> tests {
        //Line segment points
        make_tuple(point(-7, 3), false, 0),
        make_tuple(point(-2.9, -1), false, 0),
        make_tuple(point(-7, 3), true, 0),
        make_tuple(point(-2.9, -1), true, 0),

        //Points projecting between segment
        make_tuple(point(-8.1, -1), false, 3.63128),
        make_tuple(point(-2.4, 3.7), false, 3.71334),
        make_tuple(point(-4.69455, 0.75078), false, 0),
        make_tuple(point(-8.1, -1), true, 3.63128),
        make_tuple(point(-2.4, 3.7), true, 3.71334),
        make_tuple(point(-4.69455, 0.75078), true, 0),

        //Points projecting outside of segment
        make_tuple(point(-10.8, 3.5), true, 3.83275),
        make_tuple(point(-7.5, 6.5), true, 3.53553),
        make_tuple(point(-3, -4.7), true, 3.70135),
        make_tuple(point(1.2, -0.8), true, 4.10488),
        make_tuple(point(-10.8, 3.5), false, 2.29574),
        make_tuple(point(-7.5, 6.5), false, 2.15607),
        make_tuple(point(-3, -4.7), false, 2.71823),
        make_tuple(point(1.2, -0.8), false, 3.00628),

        //Points on line outside of segment
        make_tuple(point(0.82454, -4.63369), false, 0),
        make_tuple(point(-10.39683, 6.31398), false, 0),
        make_tuple(point(0.82454, -4.63369), true, 5.20345),
        make_tuple(point(-10.39683, 6.31398), true, 4.74561)
    };

    for(const auto& t : tests)
    {
        double dist = linePointDist(s, get<0>(t), get<1>(t));
        if(abs(dist-get<2>(t)) > EPSILON)
        {
            cerr << "Incorrect dist: " << get<0>(t) << " " << get<1>(t) << " Expected " << get<2>(t) << " got " << " -> " << dist << endl;
        }
    }

    return 0;

    uint64_t i = 0;
    while(true)
    {
        cout << i++ << endl;
        vector<point> pts(size_dist(reng));
        //cerr << endl << "Points: " << endl;
        for(point& p : pts)
        {  
            p.x = pt_dist(reng);
            p.y = pt_dist(reng);

            //cerr << p << " ";
        }
        //cerr << endl;

        vector<point> hull = convexHull(pts);

/*
        cerr << "Hull: " << endl;
        for(point& p : hull)
            cerr << p << " ";
        cerr << endl;
*/
        set<point> hullMap(hull.begin(), hull.end());
        for(point& p : pts)
        {
            int in = pointInPoly(hull, p);

            //Look for colinear points in polygon
            bool hasColinear = false;
            auto back = hull.end()-2, mid = hull.end()-1, front = hull.begin();
            if(hull.size() > 2)
            {
                for(; front != hull.end(); back = mid, mid = front++)
                {
                    if(orientation(*back, *mid, *front) == Colinear)
                    {
                        hasColinear = true;
                        break;
                    }
                }
            }

            //Look for point inside hull that is part of the hull or point outside hull
            if(hasColinear || in == Outside || (in == Inside && hullMap.count(p)))
            {
                cerr << "--------------------------------------------" << endl;
                cerr << "Points:" << endl;
                cerr << "\t" << pts.size() << endl;
                for(const point& p : pts)
                    cerr << "\t" << p << endl;
                cerr << "\nHull:" << endl;
                cerr << "\t" << hull.size() << endl;
                for(const point& p : hull)
                    cerr << "\t" << p << endl;
                cerr << endl;

                if(hasColinear)
                    cerr << "Colinear: " << *back << " | " << *mid << " | " << *front << endl;
                else
                    cerr << "Wrong: " << p << " : " << in << endl;
            }
        }
    }
}