function dxdt = ml_rhs(t,x,CM,gCa,gK,gL,VCa,VK,VL,V1,V2,V3,V4,phi,I,gamma_)
VA = x(1);
NA = x(2);
VB = x(3);
NB = x(4);
MinfA = (1+tanh((VA-V1)/V2))/2;
NinfA = (1+tanh((VA-V3)/V4))/2;
tauNA = 1/phi/cosh((VA-V3)/V4/2);
MinfB = (1+tanh((VB-V1)/V2))/2;
NinfB = (1+tanh((VB-V3)/V4))/2;
tauNB = 1/phi/cosh((VB-V3)/V4/2);
dVAdt = (-gL*(VA-VL) - gCa*MinfA*(VA-VCa) - gK*NA*(VA-VK) + I)/CM ...
    + gamma_*(VB-VA);
dNAdt = (NinfA-NA)/tauNA;
dVBdt = (-gL*(VB-VL) - gCa*MinfB*(VB-VCa) - gK*NB*(VB-VK) + I)/CM ...
    + gamma_*(VA-VB);
dNBdt = (NinfB-NB)/tauNB;
dxdt = [dVAdt; dNAdt; dVBdt; dNBdt];