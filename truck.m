%% TISSf-CBFs for connected automated trucks with input delays and disturbance
clear; clc;

%% Problem setup
% Delay
par.tau=0.5;        % delay [s]
% par.tau=0;        % no delay

% Type of system, predictor and safety filter
% simple or high-fidelity dynamics (without or with first-order lag)
par.model='fidel';	% 'simple', 'fidel'
% no predictor, ideal predictor with known future dynamics or approximate predictor
par.pred='approx';	% 'none', 'ideal', 'approx'
% controller based on CBF, ISSf-CBF or TISSf-CBF
par.filter='TISSf';	% 'CBF', 'ISSf', 'TISSf'

% Parameters
par.A=0.4;          % control gain [1/s]
par.B=0.5;          % control gain [1/s]
par.Dst=5;          % range policy parameter [m]
par.kappa=0.5;      % range policy parameter [1/s]
par.vmax=20;        % range policy parameter [m/s]
par.Dsf=3;          % CBF coefficient [m]
par.T=1/par.kappa;  % CBF coefficient [s]
par.sigma0=1;       % TISSf parameter [s^3/m]
par.lambda=0.3;     % TISSf parameter [1/m]
par.xi=0.25;        % first-order lag [s]

% Simulation settings
t0=0;               % start time [s]
tend=20;            % end time [s]
dt=0.01;            % fixed time step [s]
t=t0:dt:tend;       % time [s]

% Initial conditions
v0=15;                      % truck velocity
vL0=v0;                     % leader's velocity
D0=v0/par.kappa+par.Dst;    % distance
switch par.filter	% (T)ISSf increases the equilibrium distance, start closer to equilibrium
    case 'ISSf'
        D0=D0+2.5;
    case 'TISSf'
        D0=D0+2.5;
end
switch par.model
    case 'simple'
        x0=[D0;v0;vL0];     % state
    case 'fidel'
        a0=0;               % acceleration
        x0=[D0;v0;vL0;a0];	% state
end

% Initial input history excluding headpoint (empty for delay-free system)
uhist=@(t)0*t;              % input before controller is turned on
hist=t0-par.tau:dt:t0-dt;   % history
ut0=uhist(hist);            % initial input history

% Plot settings
tmin=0; tmax=10;
Dmin=0; Dmax=40;
vmin=0; vmax=15;
amin=-10; amax=2;
hmin=-5; hmax=10;
umin=-10; umax=2;
dmin=-6; dmax=6;
purple=[170,0,170]/256;
orange=[255,170,0]/256;
black=[0,0,0];
gray=[0.5,0.5,0.5];
blue=[0,0,1];
darkgreen=[0,170,0]/256;
darkred=[230,0,0]/256;

%% Simulation
% Simulate system
switch par.model
    case 'simple'
        [x,u,un]=odeab4(@rhs,@k,t,x0,ut0,par);
    case 'fidel'
        [x,u,un]=odeab4(@RHS,@k,t,x0,ut0,par);
end

% Evaluate solution
h=CBF(x,par);
[D,v,vL]=states(x);
a=gradient(v)./gradient(t);
aL=accL(t);

% Calculate disturbance
[d,dhat]=disturbance(t,x,u,ut0,par);

%% Plot results
figure; clf;

% Distance
subplot(2,3,1); hold on; box on;
plot(t,D,'Color',purple,'LineWidth',2,'DisplayName','$D$');
PlotFinalize({'time, $t$ (s)','distance, $D$ (m)'},[tmin,tmax,Dmin,Dmax]);

% Velocity
subplot(2,3,2); hold on; box on;
plot(t,vL,'Color',darkred,'LineWidth',2,'DisplayName','$v_{\rm L}$');
plot(t,v,'Color',purple,'LineWidth',2,'DisplayName','$v$');
PlotFinalize({'time, $t$ (s)','velocity, $v$ (m/s)'},[tmin,tmax,vmin,vmax]);

% Acceleration
subplot(2,3,3); hold on; box on;
plot(t,aL,'Color',darkred,'LineWidth',2,'DisplayName','$a_{\rm L}$');
plot(t,a,'Color',purple,'LineWidth',2,'DisplayName','$a$');
PlotFinalize({'time, $t$ (s)','acceleration, $a$ (m/s$^2$)'},[tmin,tmax,amin,amax]);

% Distance-velocity plot
subplot(2,3,4); hold on; box on;
Dplot=0:0.1:50;
Vdes=min(par.kappa*(Dplot-par.Dst),par.vmax);   % range policy
Vsafe=(Dplot-par.Dsf)/par.T;                    % safe set boundary
plot(Dplot,Vdes,'Color',blue,'LineWidth',2,'DisplayName','$V$')
plot(Dplot,Vsafe,'Color',black,'LineWidth',2,'DisplayName','$\partial S$')
plot(D,v,'Color',purple,'LineWidth',2,'DisplayName','$v$')
PlotFinalize({'distance, $D$ (m)','velocity, $v$ (m/s)'},[0,50,0,par.vmax]);

% CBF
% subplot(2,3,4); hold on; box on;
% plot([tmin,tmax],[0,0],'k','LineWidth',1,'HandleVisibility','off');
% plot(t,h,'Color',darkgreen,'LineWidth',2,'DisplayName','$h$');
% PlotFinalize({'time, $t$ (s)','CBF, $h$ (m)'},[tmin,tmax,hmin,hmax]);

% Disturbance
subplot(2,3,5); hold on; box on;
plot(t,dhat,'Color',orange,'LineWidth',2,'DisplayName','$\hat{d}$');
plot(t,d,'Color',darkgreen,'LineWidth',2,'DisplayName','$d$');
PlotFinalize({'time, $t$ (s)','disturbance, $d$ (m/s$^2$)'},[tmin,tmax,dmin,dmax]);

% Control input
subplot(2,3,6); hold on; box on;
plot(t,un,'Color',orange,'LineWidth',2,'DisplayName','$u_{\rm n}$');
plot(t,u,'Color',darkgreen,'LineWidth',2,'DisplayName','$u$');
PlotFinalize({'time, $t$ (s)','input, $u$ (m/s$^2$)'},[tmin,tmax,umin,umax]);

%% Save results
filename=[par.model,'_',par.pred,'_',par.filter,'_tau',num2str(par.tau,'%3.2f')];
filename=strrep(filename,'.','');
% saveas(gcf,[filename,'.fig']);
% saveas(gcf,[filename,'.svg']);

%% Functions for dynamics
% Right hand side
function dxdt = rhs(t,x,u,par)
    dxdt = f(t,x,par)+mvp(g(t,x,par),u);
end
function dxdt = RHS(t,x,u,par)
    dxdt = F(t,x,par)+mvp(G(t,x,par),u);
end

% Lead vehicle's acceleration
function aL = accL(t)
    aL = -10*(t-3).*(3<=t & t<=4)...
         -10*(4<t & t<=4.5)...
         +(10*(t-4.5)-10).*(4.5<t & t<=5.5);
end

% States
function [D,v,vL,a,xc] = states(x)
    D = x(1,:);
    v = x(2,:);
    vL = x(3,:);
    if size(x,1)==4
        a = x(4,:);
    else
        a = nan(size(D));
    end
    xc = x(1:3,:);  % states for control design
end

% Model for control design - double integrator car-following
function f = f(t,x,~)
    [~,v,vL] = states(x);
    aL = accL(t);
    f = [vL-v;
         zeros(size(v));
         aL];
end
function g = g(~,x,~)
    g = repmat(...
          [0;
           1;
           0],...
          1,1,size(x,2));
end

% Model for high-fidelity simulation - added first-order lag
function F = F(t,x,par)
    [~,v,vL,a] = states(x);
    aL = accL(t);
    F = [vL-v;
         a;
         aL;
         -a/par.xi];
end
function G = G(~,x,par)
    G = repmat(...
          [0;
           0;
           0;
           1/par.xi],...
          1,1,size(x,2));
end

%% Functions for control
% Nominal controller
function un = kd(~,x,~,par)
    [D,v,vL] = states(x);
    un = par.A*(min(par.kappa*(D-par.Dst),par.vmax)-v) + ...
         par.B*(min(vL,par.vmax)-v);
end

% CBF
function [h,gradh] = CBF(x,par)
    [D,v,vL] = states(x);
    h = D-par.Dsf-par.T*v;
    gradh = [ones(size(D));
             -par.T*ones(size(v));
             zeros(size(vL))];
end

% Controller
function [u,un,h] = k(t,x,ut,par)
    % predictor for delay compensation
    [xp,tp] = pred(t,x,ut,par);
   
    % nominal controller
%     un = kd(t,x,ut,par);  % commented lines allow for safety filtering
                            % instead of using a nominal controller
    un = kd(tp,xp,[],par);
    
    % safety filter
    [h,gradh] = CBF(xp,par);
%     Lfh = dot(gradh,f(tp,xp),1);
    Lgh = vmp(gradh,g(tp,xp));
%     Phi = (Lfh+dot(Lgh,un,1)+par.alpha*h)./dot(Lgh,Lgh,1);
    switch par.filter
        case 'CBF'      % CBF filter
            sigma = 0;
        case 'ISSf'     % ISSf filter
            sigma = par.sigma0;
        case 'TISSf'	% TISSf filter
            sigma = par.sigma0*exp(-par.lambda*h);
    end
%     u = un+max(0,-Phi+sigma).*Lgh;
    u = un+sigma.*Lgh;
end

% Predictor
function [xp,tp] = pred(t,x,ut,par)
    % use only states available for control design
    [~,~,~,~,xc] = states(x);
    % no prediction for zero delay
    if par.tau==0
        xp=xc;
        tp=t;
    % prediction when there is delay
    else
        dt = par.tau/size(ut,2);    % time step
        hor = t+(0:dt:par.tau);     % prediction horizon
        switch par.pred
            case 'none'     % no prediction
                xp = xc;
                tp = t;
            case 'ideal'    % ideal prediction with known future dynamics
                xpred = odeab4(@rhs,[],hor,xc,ut,par);
                xp = xpred(:,end);
                tp = t+par.tau;
            case 'approx'   % approximate prediction with dynamics at t
                xpred = odeab4(@(~,x,u,par)rhs(t,x,u,par),[],hor,xc,ut,par);
                xp = xpred(:,end);
                tp = t;
            case 'fidel'	% ground truth prediction with high-fidelity model for disturbance analysis
                xpred = odeab4(@RHS,[],hor,x,ut,par);
                [~,~,~,~,xp] = states(xpred(:,end));
                tp = t+par.tau;
        end
    end
end

% Disturbance
function [d,dhat] = disturbance(t,x,u,ut0,par)
    [~,v,~] = states(x);
    a = gradient(v)./gradient(t);
    if par.tau==0
        % disturbance in model
        d = a-u;
        % disturbance including prediction errors (same for zero delay)
        dhat = a-u;
    else
        % disturbance in model
        utau = interp1(t,u,t-par.tau,'linear','extrap');
        utau(t-par.tau<0) = ut0;
        d = a-utau;
        % ideal control input with ground truth info (no prediction error)
        uid=nan(size(u));
        ut=[ut0,u];
        pred0=par.pred;
        switch par.model	% change predictor to one with ground truth
            case 'simple'
                par.pred='ideal';
            case 'fidel'
                par.pred='fidel';
        end
        for kt=1:size(x,2)  % calculate ideal input at each state
            uid(:,kt)=k(t(kt),x(:,kt),ut(:,kt:kt+size(ut0,2)),par);
        end
        par.pred = pred0;	% change predictor back
        % disturbance including prediction errors
        uidtau = interp1(t,uid,t-par.tau,'linear','extrap');
        uidtau(t-par.tau<0) = ut0;
        dhat = a-uidtau;
    end
end

%% Functions for vectorized computation
% Matrix-vector product (for multiple matrices and column vectors)
function p = mvp(a,b)
    p = nan(size(a,1),size(b,2));
    for kt=1:size(b,2)
        p(:,kt) = a(:,:,kt)*b(:,kt);
    end
end

% Vector-matrix product (for multiple row vectors and matrices)
function p = vmp(a,b)
    p = nan(size(b,2),size(a,2));
    for kt=1:size(a,2)
        p(:,kt) = (a(:,kt).'*b(:,:,kt)).';
    end
end

%% Functions for plotting
% Finalize plot with axis limits, labels, legends, title
function PlotFinalize(axislabels,axislimits)
    axis(axislimits);
    pbaspect([1,1,1]);
    xlabel(axislabels{1},'Interpreter','latex');
    ylabel(axislabels{2},'Interpreter','latex');
    if length(axislabels)>2
        zlabel(axislabels{3},'Interpreter','latex');
    end
    set(gca,'TickLabelInterpreter','latex','FontSize',12);
    legend('Location','best','Interpreter','latex','FontSize',14);
    if isempty(get(get(gca,'Legend'),'String'))
        legend off;
    end
end