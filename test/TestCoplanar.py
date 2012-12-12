#!/usr/bin/env python

"""
TODO: Finish.
"""

def grid(side, n):
        """
        Returns a grid with specified side length and number of points,
        centred at the origin and with z=0.
        """
        if n < 2:
                raise TypeError('n must be >= 2')
        coords = []
        for yi in range(0,n):
                y = float(yi) / (n-1) * side - (side/2.0)
                for xi in range(0,n):
                        x = float(xi) / (n-1) * side - (side/2.0)
                        coords.append([x, y, 0.0])
        return coords

