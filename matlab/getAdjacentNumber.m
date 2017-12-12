function [center, nw, ne, w, e, sw, se] = getAdjacentNumber(r, c, ncols, nwr, ner, wr, er, swr, ser, nwc, nec, wc, ec, swc, sec)

center = (r-1)*ncols+c;
nw = (nwr-1)*ncols+nwc;
ne = (ner-1)*ncols+nec;
w = (wr-1)*ncols+wc;
e = (er-1)*ncols+ec;
sw = (swr-1)*ncols+swc;
se = (ser-1)*ncols+sec;