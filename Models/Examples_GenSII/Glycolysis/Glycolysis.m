function model = Glycolysis()
    % Glycolysis provides the GenSSI implementation of the model
    % of the glycolysis metabolic pathway as introduced by
    % 
    %    Bartl et al. (2010). Just-in-time activation of a glycolysis
    %              inspired metabolic network - solution with a dynamic 
    %              optimization approach. Proc. 55nd International 
    %              Scientific Colloquium. Ilmenau, Germany.

    % Symbolic variables
    syms k1 k2 k3 k4 kM x10 x20 x30 x40 x50
    syms x1 x2 x3 x4 x5

    % Parameters
    model.sym.p = [k1;k2;k3;k4;kM;x10;x20;x30;x40;x50];

    % State variables
    model.sym.x = [x1;x2;x3;x4;x5];
    
    % Control vectors (g)
    model.sym.g = [-(k1*x1)/(kM+x1),               0,               0,               0
                    (k1*x1)/(kM+x1),-(k2*x2)/(kM+x2),               0,               0
                                  0, (k2*x2)/(kM+x2),-(k3*x3)/(kM+x3),               0
                                  0, (k2*x2)/(kM+x2), (k3*x3)/(kM+x3),-(k4*x4)/(kM+x4)
                                  0,               0,               0, (k4*x4)/(kM+x4)];

    % Autonomous dynamics (f)
    model.sym.xdot = [0
                      0
                      0
                      0
                      0];

    % Initial conditions
    model.sym.x0 = [x10;x20;x30;x40;x50];
    
    % Observables
    model.sym.y = [x1;x2;x3;x4;x5];
end



