function model = HighDimNonLin()
    % HighDimNonLin provides the GenSSI implementation of the model by
    % 
    %    Saccomani et al. (2005). Examples of testing global 
    %    identifiability of biological and biomedical models with the
    %    DAISY software, Computers in Biology and Medicine, 40, 402-407.
    %
    % which is denoted as "high-dimensional nonlinear model".

    % Symbolic variables
    syms x04 x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 x15 x16 x17 x18 x19 x20
    syms x01 x02 x03 x04 x05 x06 x07 x08 x09 x010 x011 x012 x013 x014 x015 x016 x017 x018 x019 x020
    syms vm km p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13 p14 p15 p16 p17 p18 p19 p20

    % Parameters
    % (Note: As the initial conditions are not part of the parameter
    % vector, they are assumed to be known.)
    model.sym.p = [vm;km;p1;p2;p3;p4;p5;p6;p7;p8;p9;p10;p11;p12;p13;p14;p15;p16;p17;p18;p19;p20];

    % State variables
    model.sym.x = [x1;x2;x3;x4;x5;x6;x7;x8;x9;x10;x11;x12;x13;x14;x15;x16;x17;x18;x19;x20];

    % Control vectors (g)
    model.sym.g = [1
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0];

    % Autonomous dynamics (f)
    model.sym.xdot=[-vm*x1/(km+x1)-p1*x1
                    p1*x1-p2*x2
                    p2*x2-p3*x3
                    p3*x3-p4*x4
                    p4*x4-p5*x5
                    p5*x5-p6*x6
                    p6*x6-p7*x7
                    p7*x7-p8*x8
                    p8*x8-p9*x9
                    p9*x9-p10*x10
                    p10*x10-p11*x11
                    p11*x11-p12*x12
                    p12*x12-p13*x13
                    p13*x13-p14*x14
                    p14*x14-p15*x15
                    p15*x15-p16*x16
                    p16*x16-p17*x17
                    p17*x17-p18*x18
                    p18*x18-p19*x19
                    p19*x19-p20*x20];

    % Initial conditions
    model.sym.x0=[x01;x02;x03;x04;x05;x06;x07;x08;x09;x010;x011;x012;x013;x014;x015;x016;x017;x018;x019;x020];

    % Observables
    model.sym.y = model.sym.x;
end

