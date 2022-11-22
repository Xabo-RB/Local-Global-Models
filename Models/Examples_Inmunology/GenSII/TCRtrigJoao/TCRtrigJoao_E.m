function model = TCRtrigJoao_E()
% HIV Model with Constant and Time Varying Parameters
% Miao H, Xia X, Perelson AS, Wu H. 
% "On identifiability of nonlinear ODE models and applications in viral dynamics." 
% SIAM review 53.1 (2011): 3-39.

%COMO SI EL PARÁMETRO DESCONOCIDO FUESE UNA ENTRADA CONOCIDA U

    % Symbolic variables
	syms lambda phi s ki hh k L
    syms S T A
	syms S0 T0 A0

    % Parameters
	model.sym.p = [lambda; phi; s; ki; hh; k; L];
    %model.sym.p = [phi; s; ki; hh; k; L];

    % State variables
	model.sym.x = [S T A];

    % Control vectors (g)
	model.sym.g = [0
                   0
                   0];

    % Autonomous dynamics (f)
	model.sym.xdot = [-lambda * phi * (S-T) + s*(1-S)
					phi * (S - T) + s*(1-T) - k * (T^hh)*(L^hh)
					k *(T^hh)*(L^hh) - ki*A];

    % Initial conditions
	model.sym.x0 = [S0;T0;A0];

    % Observables    
	model.sym.y = [S/(lambda+1) + ((T+A)*lambda/(lambda+1))];
end