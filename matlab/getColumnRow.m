function [nwr, nwc, ner, nec, wr, wc, er, ec, swr, swc, ser, sec] = getColumnRow(r, c)
        
nwr = r-1; nwc = c; %northwest row, column
ner = r-1; nec = c+1; %northeast row, column
wr = r; wc = c-1; %west row, column
er = r; ec = c+1; %east row, column
swr = r+1; swc = c-1; %southwest row, column
ser = r+1; sec = c; %southeast row, column