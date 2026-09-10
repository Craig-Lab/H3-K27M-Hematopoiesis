%--------------------------------------------------------------------------
% File that creates the Transiting HSCs with MEP Bypass model. 
% It stores it in a mat-file in the 'model' file.
%--------------------------------------------------------------------------
clear all;

syms k1 k2 p1 p2 p3 p4 p5 p6 p7 p8 p9 d3 d4 d5 d6 d7 d8 d9 r12 r21 ...
     HSC1 HSC2 MPP CMP GMP MEP EPCD71p EPCD71n MLP ...
    
% states:
x = [HSC1; HSC2; MPP; CMP; GMP; MEP; EPCD71p; EPCD71n; MLP]; 

% outputs:
h = x;  

% no input:
u = [];

% parameters:
p = [k1; k2; p1; p2; p3; p4; p5; p6; p7; p8; p9; ...
        d3; d4; d5; d6; d7; d8; d9; r12; r21];  
    
% dynamic equations:    
f = [(k1-p1-r12)*HSC1 + r21*HSC2;
     (k2-p2-r21)*HSC2 + r12*HSC1;
     p1*HSC1 + p2*HSC2 - (d3+p3+p6+p9)*MPP;
     p3*MPP - (d4+p4+p5)*CMP;
     p4*CMP - d5*GMP;
     p5*CMP - (d6+p7)*MEP;
     p8*EPCD71p - d9*EPCD71n;
     p9*MPP + p7*MEP - (d8+p8)*EPCD71p;
     p6*MPP - d7*MLP;]; 
    
% initial conditions:    
ics = [16578.9; 4950.22; 7959.44; 34212.1; 68763.4; 8402.0; 76082.2; 57766.0; 11298.4]; % H3-K27M initial conditions

% which initial conditions are known:
known_ics = [1,1,1,1,1,1,1,1,1];

save('../models/TransitingHSC_MEPBypass','x','h','u','p','f','ics','known_ics');
