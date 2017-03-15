"""
Generate random data for testing.
"""

import numpy as np 
from numpy.random import random
from numpy.linalg import norm

class Walker(object):

    def __init__(self, mins, maxes, delta=0.1, attractor_weight=()):
        self.delta = delta
        self.mins = np.array(mins, np.float)
        self.maxes = np.array(maxes, np.float)
        self.offset = self.maxes - self.mins
        self.side = norm(self.offset)
        self.attract = [(np.array(a, np.float), wt) for (a, wt) in attractor_weight]
        r = random((2,))
        start = self.mins + r * (self.maxes - self.mins)
        self.points = [start]

    def next(self, epsilon=1e-3):
        last = self.points[-1]
        last_offset = last - self.mins
        #print "last_offset", last_offset
        shift = random((2,)) - 0.5
        shift = random() * shift / (norm(shift) + epsilon)
        #print "initial shift", shift
        for (attract, wt) in self.attract:
            diff = attract - last
            #print attract, "diff", attract
            ndiff = norm(diff)
            if abs(ndiff) > epsilon:
                diff = diff / ndiff
                attraction = self.side / (self.side + ndiff)
                #print "attraction", attraction
                rattract = random() * wt * attraction * diff
                shift += rattract
                #print "rattract", rattract, shift
        #print "final shift", shift
        initial_offset = next_offset = last_offset + self.delta * shift
        offset = self.offset
        #next_offset = next_offset % self.offset
        too_big = (next_offset - offset)
        #print "too big before", too_big
        too_big[np.where(too_big < 0)] = 0
        #print "too_big after", too_big
        next_offset -= 2 * too_big
        #print "too_big", too_big, initial_offset, "fixed", next_offset
        too_small = next_offset.copy()
        too_small[np.where(too_small > 0)] = 0
        #print "too_small adjusted", too_small
        next_offset -= 2 * too_small
        #print "too small", too_small, initial_offset, "fixed", next_offset
        assert np.all(next_offset >= 0), repr((initial_offset, next_offset))
        assert np.all(next_offset <= offset), repr((initial_offset, next_offset))
        next =  self.mins + next_offset
        self.points.append(next)
        return next

def visual_test():
    from jp_svg_canvas import canvas, cartesian_svg
    import time
    mins = (minx, miny) = (-2, -5)
    maxes = (maxx, maxy) = (8, 5)
    C = cartesian_svg.doodle(minx, miny, maxx, maxy)
    C.show()
    attractor_weight = [((0,3), 0.2), ((5,0), 0.1)]
    delta = 0.5
    W = Walker(mins, maxes, delta, attractor_weight)
    for ((x,y), r) in attractor_weight:
        C.circle(None, x, y, r, "pink")
    C.flush()
    C.axes()
    last = None
    for i in range(1000):
        next = (x, y) = W.next()
        if last is not None:
            C.line(None, x, y, last[0], last[1], "magenta")
        C.circle(None, x, y, 0.05, "red")
        last = next
        C.flush()

def random_point(mins, maxes):
    r = random((2,))
    return mins + r * (maxes - mins)

def test_tsv(filename="../examples/test_data.tsv"):
    outfile = open(filename, "w")
    mins_slide = np.array((0,0))
    maxes_slide = np.array((40, 40))
    mins_tsne = np.array((-10,-10))
    maxes_tsne = np.array((10,10))
    delta_slide = 2
    delta_tsne = 0.5
    a_s = 0.4
    a_t = 0.4
    header = "SLIDE\tSLIDE_x\tSLIDE_y\tTSNE_x\tTSNE_y\tINDICATOR\n"
    outfile.write(header)
    count = 0
    for indicator in range(15):
        indicator_s = "ind" + repr(indicator)
        nslides = 1 + int(2 * random())
        for slide in range(nslides):
            slide_s = "slide" + repr(indicator) + "_" + repr(slide)
            s_attractors = [(random_point(mins_slide, maxes_slide), a_s) for i in range(2)]
            t_attractors = [(random_point(mins_tsne, maxes_tsne), a_t) for i in range(2)]
            #print "s_attractors", s_attractors
            #print "t_attractors", t_attractors
            Ws = Walker(mins_slide, maxes_slide, delta_slide, s_attractors)
            Wt = Walker(mins_tsne, maxes_tsne, delta_tsne, t_attractors)
            npoints = 10 + int(random() * 500)
            for i in range(npoints):
                (sx, sy) = Ws.next()
                (tx, ty) = Wt.next()
                data = map(str, [slide_s, sx, sy, tx, ty, indicator_s])
                line = ("\t".join(data)) + "\n"
                outfile.write(line)
                count += 1
    outfile.close()
    print "wrote", count, "to", filename

def test0():
    mins = [-1, -10]
    maxes = [9, 0]
    delta = 1
    attractor_weight = [((5,3), 0.5)]
    W = Walker(mins, maxes, delta, attractor_weight)
    print W.points, "offset", W.offset, W.side
    for i in range(10):
        print W.next()
    print W.points

if __name__ == "__main__":
    #test0()
    test_tsv()
