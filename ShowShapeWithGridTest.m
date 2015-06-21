ccc

load VestaHASTALAVESTAshape_sh720.mat % loads lmcosi_shape

Resolution=1; 
GridResolution=10;
MaxLatitude=81;
MaxDegreeTopo=100;

ShowShapeWithGrid(lmcosi_shape,MaxDegreeTopo,Resolution,GridResolution,MaxLatitude)