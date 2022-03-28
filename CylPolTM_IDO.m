function[alphaTM]=CylPolTM_IDO(params, aCyl,k0,epsilon_r1,epsilon_r2);
%
% A funcrion that computes the electric polarizability of thin dielectic
% cylinder in TM illomination (E field along the cylinder axis)
omega_mu0=params.k0*params.c*params.mu0; %k0*3e8*1.256637061435917e-06;

n1=params.n_out;
n2=params.n_in;
k0n1a=k0*n1*aCyl;
k0n2a=k0*n2*aCyl;
j0k0n2a=besselj(0,k0n2a);
j1k0n2a=besselj(1,k0n2a);
Dnum=-besselj(1,k0n1a)*j0k0n2a*n1+j1k0n2a*besselj(0,k0n1a)*n2;
Den=-besselh(1,1,k0n1a)*j0k0n2a*n1+j1k0n2a*besselh(0,1,k0n1a)*n2;
minus_b0TM=Dnum/Den;
alphaTM=4*minus_b0TM/omega_mu0;  
end