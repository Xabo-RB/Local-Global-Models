%% MODELO D (h = unknown)

syms T A
x = [T; A];

% one input
syms u1;
u = [];

% 4 unknown parameters 
syms phi s ki hh k L
p =[s; ki; hh; k; L];

% 1 output
h = T + A;

% dynamic equations
f = [s*(1-T) - k * (T^hh)*(L^hh);
    k *(T^hh)*(L^hh) - ki*A];


% initial conditions
ics  = []; 
known_ics = [0,0];

save('TCRtrig','x','p','h','f','u','ics','known_ics');



%% MODELO D (h = cte)

syms T A
x = [T; A];

% one input
syms u1;
u = [];

% 4 unknown parameters 
syms phi s ki hh k L
p =[s; ki; hh; k; L];

% 1 output
h = T + A;

% dynamic equations
f = [s*(1-T) - k * (T^hh)*(L^hh);
    k *(T^hh)*(L^hh) - ki*A];


% initial conditions
ics  = []; 
known_ics = [0,0];

save('TCRtrig','x','p','h','f','u','ics','known_ics');





