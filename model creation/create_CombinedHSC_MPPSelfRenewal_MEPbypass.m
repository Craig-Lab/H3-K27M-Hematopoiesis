%--------------------------------------------------------------------------
% File that creates the Combined HSCs with MPP Self-Renewal & MEP Bypass model. 
% It stores it in a mat-file in the 'model' file.
%--------------------------------------------------------------------------
clear all;

syms k1 p1 p3 p4 p5 p6 p7 p8 p9 k3 d4 d5 d6 d7 d8 d9 ...
     HSC MPP CMP GMP MEP EPCD71p EPCD71n MLP ...
    
% states:
x = [HSC; MPP; CMP; GMP; MEP; EPCD71p; EPCD71n; MLP]; 

% outputs:
h = x;  

% no input:
u = [];

% parameters:
p = [k1; p1; p3; p4; p5; p6; p7; p8; p9; ...
        k3; d4; d5; d6; d7; d8; d9];  
    
% dynamic equations:    
f = [(k1-p1)*HSC;
     p1*HSC - (p3+p6-k3+p9)*MPP;
     p3*MPP - (d4+p4+p5)*CMP;
     p4*CMP - d5*GMP;
     p5*CMP - (d6+p7)*MEP;
     p8*EPCD71p - d9*EPCD71n;
     p9*MPP + p7*MEP - (d8+p8)*EPCD71p;
     p6*MPP - d7*MLP;]; 
    
% initial conditions:    
ics = [7536.11; 5395.67; 20271.3; 115966.0; 5997.0; 66362.4; 57536.0; 32338.0]; % H3-WT initial conditions
% which initial conditions are known:
known_ics = [1,1,1,1,1,1,1,1];

save('../models/CombinedHSC_MPPSelfRenewal_MEPbypass','x','h','u','p','f','ics','known_ics');
