%% 4 step Adams-Bashforth scheme to solve ODEs of control systems
function [xsolu,usolu,unsolu]=odeab4(RHS,k,t,x0,ut0,par)
% calculation of the fixed time step
if length(t)>1
    dt=t(2)-t(1);
end

% memory allocation for solution
xsolu=nan(length(x0),length(t));
usolu=nan(size(ut0,1),length(t));
unsolu=nan(size(ut0,1),length(t));

% solution of the differential equation
for kt=1:length(t)
    % subsequent time step
    tnew=t(kt);
    
    % evaluation of initial conditions
    if kt==1
        xnew=x0;
        utnew=ut0;
    else
        % evaluation of the right-hand side of the ODE
        rhs1=RHS(told,xold,udelay,par);
        
        % calculation of the solution
        if kt==2
            xnew = xold+dt*rhs1;
        elseif kt==3
            xnew = xold+dt*(3*rhs1-rhs2)/2;
        elseif kt==4
            xnew = xold+dt*(23*rhs1-16*rhs2+5*rhs3)/12;
        elseif kt>4
            xnew = xold+dt*(55*rhs1-59*rhs2+37*rhs3-9*rhs4)/24;
        end
        
        % calculation of the input history
        utnew=[utold(:,2:end),uold];
    end
    
    % given a controller get the input, otherwise rely on input history only
    if ~isempty(k)
        [unew,unnew] = k(tnew,xnew,utnew,par);
    else
        unew = nan(size(ut0,1),1);
        unnew = nan(size(unew));
    end
    
    % store right-hand side expressions for the multi step method
    if kt>3
        rhs4=rhs3;
    end
    if kt>2
        rhs3=rhs2;
    end
    if kt>1
        rhs2=rhs1;
    end
    
    % store solution
    xsolu(:,kt)=xnew;
    usolu(:,kt)=unew;
    unsolu(:,kt)=unnew;

    % update
    told=tnew;
    xold=xnew;
    uold=unew;
    utold=utnew;

    % get input
    if isempty(utold)   % delay-free case
        udelay=uold;
    else                % delayed input at the start of history
        udelay=utold(:,1);
    end
end