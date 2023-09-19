%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ----------------------------------------------------------------------- %
%                  HIGH-ORDER IMMERSED INTERFACE METHOD                   %
%            H. EScamilla Puc,  R. Itza Balam & M. Uh Zapata              %
%                                Sept 2023                                %
% ----------------------------------------------------------------------- %
%  It solves the one-dimensional Poisson equation:                        %
%                           u_xx  = f                                     %
%  knowing jump conditions: [u], [u_x], [f], [fx], and [fxx] at the       %
%  interface and Dirichlet boundary conditions.                           %
% ----------------------------------------------------------------------- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROBLEM PARAMETERS & FUNCTIONS

%--------------------------------
% Discretization
xI   = 0.0;               % Initial of the domain: (xI,xF)
xF   = 1.0;               % Final   of the domain: (xI,xF)
alf  = 0.4;               % Location of the interface
Mvec = [10,20,40,80,160]; % Number of sub-divitions (vector) 
%--------------------------------
% Method
b = 1/12;                 % b=1/12 (4th-order), b=0 (2nd-order)
%--------------------------------
% Functions
fun_uL    = @(x)         sin(pi*x);
fun_uxL   = @(x)      pi*cos(pi*x);
fun_fL    = @(x) -(pi^2)*sin(pi*x);
fun_fxL   = @(x) -(pi^3)*cos(pi*x);
fun_fxxL  = @(x)  (pi^4)*sin(pi*x);

fun_uR    = @(x)         cos(pi*x);
fun_uxR   = @(x)     -pi*sin(pi*x);
fun_fR    = @(x) -(pi^2)*cos(pi*x);
fun_fxR   = @(x)  (pi^3)*sin(pi*x);
fun_fxxR  = @(x)  (pi^4)*cos(pi*x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WORKSPACE

M = length(Mvec);
hvec     = zeros(M,1);
NormE1   = zeros(M,1);
NormE2   = zeros(M,1);
OrderOr1 = zeros(M,1);
OrderOr2 = zeros(M,1);

for s=1:M
    %----------------------------------------------------------------------
    % DISCRETIZATION
    n = Mvec(s);
    %---------------------
    % Points and step size
    x  = linspace(xI,xF,n+1)';
    h  = x(2)-x(1); 
    h2 = h*h;
    %---------------------
    % Find I: the interval where alpha is
    for i=1:n
        if  (x(i)<=alf) && (alf<x(i+1))
            I = i;
            break;
        end
    end
    %----------------------------------------------------------------------
    % MATRIX & RIGHT-HAND SIDE
    A1 = zeros(n,1);
    A2 = zeros(n+1,1);
    A3 = zeros(n,1);
    rhs= zeros(n+1,1);
    %---------------------
    f(1:I)     = fun_fL(x(1:I));
    f(I+1:n+1) = fun_fR(x(I+1:n+1));    
    %---------------------
    % Boundary
    A2(1)    = 1/h2;
    rhs(1)   = fun_uL(x(1))/h2; 
    A2(n+1)  = 1/h2;
    rhs(n+1) = fun_uR(x(n+1))/h2;
    %---------------------
    % Regular points: Equation (22) in Escamilla et al. 2023
    A1(1:n-1) =  1/h2;
    A3(2:n)   =  1/h2;
    A2(2:n)   = -2/h2;
    rhs(2:n)  =  b*f(3:n+1)+(1-2*b)*f(2:n)+b*f(1:n-1);
    %---------------------
    % Irregular points: Equation (26) in Escamilla et al. 2023
    hL   = x(I)   - alf;
    hR   = x(I+1) - alf;
    %-----  
    uJ   = fun_uR(alf)   - fun_uL(alf);    % [u]
    uxJ  = fun_uxR(alf)  - fun_uxL(alf);   % [ux]
    fJ   = fun_fR(alf)   - fun_fL(alf);    % [f]
    fxJ  = fun_fxR(alf)  - fun_fxL(alf);   % [fx]
    fxxJ = fun_fxxR(alf) - fun_fxxL(alf);  % [fxx]
    %-----
    CI   =  (1/h2)*(uJ + hR*uxJ + 0.5*(hR^2)*fJ) ...
            - b*(fJ + hR*(1-2*hR^2/h2)*fxJ + 0.5*hR^2*(1-hR^2/h2)*fxxJ);
    CIp1 = -(1/h2)*(uJ + hL*uxJ + 0.5*(hL^2)*fJ) ...
            + b*(fJ + hL*(1-2*hL^2/h2)*fxJ + 0.5*hL^2*(1-hL^2/h2)*fxxJ);
    %-----    
    rhs(I)   = rhs(I)   + CI;
    rhs(I+1) = rhs(I+1) + CIp1;
    %----------------------------------------------------------------------
    % LINEAR SYSTEM SOLUTION BY THOMAS
    U = zeros(n+1,1);
    for i=1:n
        A2(i+1)=A2(i+1)-A3(i)*A1(i)/A2(i);
        rhs(i+1)=rhs(i+1)-rhs(i)*A1(i)/A2(i);
    end
    U(n+1)=rhs(n+1)/A2(n+1);
    for i=n:-1:1
        U(i)=(rhs(i)-U(i+1)*A3(i))/A2(i);
    end       
    %----------------------------------------------------------------------
    % ERRORS & ORDER
    %---------------------
    % Exact solution
    Uexact(1:I)     = fun_uL(x(1:I));
    Uexact(I+1:n+1) = fun_uR(x(I+1:n+1));
    %---------------------
    % Norm errors
    Err       = abs(Uexact'-U);
    NormE1(s) = norm(Err,inf);
    NormE2(s) = sqrt(h*sum(Err.^2));
    %---------------------
    % Estimated order
    hvec(s) = h;
    if s~=1
        div = log((hvec(s-1))/(hvec(s)));
        OrderOr1(s) = log(NormE1(s-1)/NormE1(s))/div;
        OrderOr2(s) = log(NormE2(s-1)/NormE2(s))/div;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DISPLAY & PLOT

disp('   N   | Max-norm Order | L2-norm  Order')
disp('----------------------------------------')
fprintf(' %5i | %.2e  ---- | %.2e  ---- \n',...
    Mvec(1),NormE1(1),NormE2(1))
for i=2:M
    fprintf(' %5i | %.2e  %.2f | %.2e  %.2f \n', ...
    Mvec(i),NormE1(i),OrderOr1(i),NormE2(i),OrderOr2(i))
end

figure
hold on
plot(x,U,'o')
plot([x(1:I);alf],[Uexact(1:I)';fun_uL(alf)],'-k','LineWidth',2)
plot([alf;x(I+1:end)],[fun_uR(alf);Uexact(I+1:end)'],'-k','LineWidth',2)
plot([alf,alf],[min(U),max(U)],'--r')
legend('Numerical','Analytical')
xlabel('x')
ylabel('u')
title(sprintf('Solution (N=%d)',n))

figure
hold on
plot(x,Err,'-k','LineWidth',2)
plot([alf,alf],[min(Err),max(Err)],'--r')
xlabel('x')
ylabel('|U-u|')
title(sprintf('Absolute error (N=%d)',n))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
