# Convex Hull
https://open.kattis.com/problems/convexhull2

Tags: Geometry, Convex Hull

## Problem Summary
You are given a set of points, each one is indicated to be either on the convex hull or not. Output the number of points on the hull and the points themselves in counterclockwise order, starting with the furthest left, then furthest down.

## Solution
This problem uses the same sort algorithm as the Graham Scan convex hull. First, the lower-left point is found (minimum y, then minimum x) and points are sorted according to their angle relative to that point. Then the group of colinear points at the beginning and end of the sorted list must be sorted again, but by their distance from the lower left point. After points are sorted, the start point for the output must be found, and then points are outputted in order following that point.