%% MODELO E (h = cte)

syms S T A
x = [S; T; A];

% one input
syms u1;
u = [];

% 4 unknown parameters 
syms lambda phi s ki hh k L
p =[lambda; phi; s; ki; hh; k; L];

% 1 output
h = S/(lambda+1) + ((T+A)*lambda/(lambda+1));

% dynamic equations
f = [-lambda * phi * (S-T) + s*(1-S);
    phi * (S - T) + s*(1-T) - k * (T^hh)*(L^hh);
    k *(T^hh)*(L^hh) - ki*A];

% 1 outpu
h = S/(lambda+1) + ((T+A)*lambda/(lambda+1));

% initial conditions
ics  = []; 
known_ics = [0,0,0];

save('TCRtrig','x','p','h','f','u','ics','known_ics');



%% MODELO E (h unknown)

syms S T A
x = [S; T; A];

% one input
syms u1;
u = [];

% 4 unknown parameters 
syms lambda phi s ki k L
p =[lambda; phi; s; ki; k; L];


% dynamic equations
f = [-lambda * phi * (S-T) + s*(1-S);
    phi * (S - T) + s*(1-T) - k * (T^5)*(L^5);
    k *(T^5)*(L^5) - ki*A];

% 1 output
h = S/(lambda+1) + ((T+A)*lambda/(lambda+1));

% initial conditions
ics  = []; 
known_ics = [0,0,0];

save('TCRtrig','x','p','h','f','u','ics','known_ics');





