clear;

% 4 states
syms x1 x2 dx1 dx2
x = [x1;x2;dx1;dx2];

% 1 known input and 1 unknown input
syms u;
u = u;
w = [];

% 3 unknown parameters 
syms p1 p2 p3 p4
p =[p1; p2; p3; p4];

% 2 outputs
h = [x3];

% initial conditions
ics  = []; 
known_ics = [0,0,0];

% dynamic equations
f = [- p1*x1 + p2 + u;
- p3* + p4;
- (p1+p3)*x3 + (p4*x1+p2*x2)];

save('LLW1987','x','p','u','w','h','f','ics','known_ics');
