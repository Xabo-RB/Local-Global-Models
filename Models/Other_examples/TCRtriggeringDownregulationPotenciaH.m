syms S T A
x = [S; T; A];

% 1 output
h = S/(lambda+1) + ((T+A)*lambda/(lambda+1));

% one input
syms u1;
u = [];

% 4 unknown parameters 
syms lambda phi s ki h
p =[lambda; phi; s; ki; h];

% 1 output
h = S/(lambda+1) + ((T+A)*lambda/(lambda+1));

% dynamic equations
f = [-lambda * phi * (S-T) + s*(1-S);
    phi * (S - T) + s*(1-T) - k * (T^h)*(L^h);
    k *(T^h)*(L^h) - ki*A];

% initial conditions
ics  = []; 
known_ics = [0,0];

save('C2M','x','p','h','f','u','ics','known_ics');







