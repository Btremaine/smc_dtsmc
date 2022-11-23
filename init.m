% init.m
% dtsmc tutorial for 2nd order system

clear, close all

% MKS units
kv= 0.05;
ks= 0.02;
M= 5;
Ts = 5E-3;
d=1.0;   % multiplier on plant mass

% c.t. plant, states {x, x_dot}%
%   M*dx_dot/dt =  U - ks*x - kv*x_dot
Ac = [  0      1;                  
     -ks/M -kv/M;];
Bc = [ 0;
      1/M ];
Cc = [1 0];
Dc= 0;
% check plat step response
sysc= ss(Ac,Bc,Cc,Dc);
step(sysc)

% discretize
sysd = c2d(sysc,Ts);

% matricies for sliding mode controller\
A= sysd.A;
B= sysd.B;
C= sysd.C;
D= sysd.D;

%% Sliding mode parameters ==============================================
dmin = -0.05;
dmax = +0.2;

% initialize defaults
d1 = 0;
del_d1 = 0;
d2 = 0;
del_d2 = 1;
alpha = 0;
beta= 0;
c1= [1; 1];
c2= [1; 0];

% design parameters s0>0 epsi>0, :: h(s) = 1         for |s|>=s0
%                                          |s|/s0    for |s|<s0
epsi = 0.1;          % RD1 0.0002;
s0   = 1.0;         % 1.0;

method = 'gao';   % <==== set method desired, not case sensitive

switch upper(method)
    case 'RD1'
      c1 = [0.25; 
            1.0];
      d1= 0.5*c1'*B*(dmax + dmin);
      del_d1 = 0.5*abs(c1'*B)*(dmax - dmin);
      % check constraints
      assert(abs(c1'*B)> eps,'c1^T*B is zero !')
    case 'RD2'  % RD2 not debugged 11/22/22
      c2 = 10*[ 0.25; 
            -0.25*B(1)/B(2)];
      d2= 0.5*c2'*A*B*(dmax + dmin);
      del_d2 = 0.5*abs(c2'*A*B)*(dmax - dmin);
      % check constraints
      assert(abs(c2'*B)<eps,' c2^T*B not zerou !')
      assert(abs(c2'*A*B)>eps,' c2^T*A*B is zero !')
      % check eigenvalues
      Ac2 = A - B*inv(c2'*A*B)*c2'*A^2;
      disp('eig of Ac2')
      disp(eig(Ac))
    case 'GAO'
      alpha = 0.002;
      beta= 0.001;
      c1 = [0.25; 
            1.0];
      assert(abs(c1'*B)> eps,'c1^T*B is zero !')
    otherwise
        disp('Unknown method')
end

%%
% LQR design follows:
% augment plant with new state, integral of error {x, x_dot, xa_dot}
Aa = [  0      1       0;                  
      -ks/M -kv/M     0;
       -1     0       0];
Ba = [ 0;
      1/M; 
      0];
Ca = [1 0 0];
Da= 0;

% Use lqr to calculate gain for all 3 states
sys= ss(Aa,Ba,Ca,Da);
Q= 0.95*ones(3);
R= 0.05;
K = lqr(sys, Q, R);

% 2nd order system
Q2= 0.95*ones(2);
R2= 0.05;
K2 = lqr( ss(A,B,C,D), Q2, R2);
kf1= K2(1);
kf2= K2(2);

