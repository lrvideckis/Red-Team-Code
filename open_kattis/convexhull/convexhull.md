# Convex Hull
https://open.kattis.com/problems/convexhull

Tags: Geometry, Convex Hull

## Problem Summary
You are given a set of points, output the convex hull in counter-clockwise order, starting with any point

## Solution
This problem requires any convex hull algorithm faster than O(n^2). The Graham Scan algorithm is O(n) after the points are sorted, bringing the complexity to O(nlgn) in this point. Points are sorted with a special sort function that arranges them in order going counter clockwise from the lower left most point, and then the points at the beginning and end of the list are sorted by distance away from that origin point. Special care must be taken to handle edge cases correctly:

    * 1 Point
    * All points the same
    * 2 different points
    * Duplicates of any point
    * Duplicates of specifically the lower-left point
    * Multiple points colinear on the left or right after sorting counter clockwise