ccc

A=[1 2 0; 1 2 1; 1 2 2; 0 2 1; 0 2 2];
 
B=[1 2 0; 1 2 1; 0 2 1; 1 2 2; 0 2 2];
 
[c,ia,ib]=intersect(A,B,'rows');

A(ia,:)-B

A(ib,:)-B

B(ia,:)-A

B(ib,:)-A

[junk, P] = ismember(A,B,'rows');