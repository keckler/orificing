function [center, nw, ne, w, e, sw, se] = convertToReducedMatrix(map, center, nw, ne, w, e, sw, se)

center = map(center);
nw = map(nw);
ne = map(ne);
w = map(w);
e = map(e);
sw = map(sw);
se = map(se);