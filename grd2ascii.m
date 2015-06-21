ccc
filename='VestaShapeDLR64_pixel.grd';
ri=ncread(filename,'z')';
outfilename='VestaShapeDLR64_pixel.txt';
WriteASCIIGrid(outfilename,ri)