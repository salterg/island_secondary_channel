clear all
%%
runid= 'mississippi_scenario.mat'
load('initial_condition_island_archive.mat')
hotstart=1;
%numerical
N=20000; %number of cells
deta_thresh=1*10^-2;
watlev_accuracy=10^-9;
r=1.000; %stretching factor for grid
timesteps=1000000;
saveinterval=500;
dt_initial=10000;
dt=dt_initial;

%parameters
L_star=10;
offset_star=.3; %only for cold start simulation
S_factor=2;
SLR_factor=1;
ds_factor=1;
b_star=.05; %reference height for channel rouse profile
% switches
momentum_fix=1;%1=include head loss correction at upstream boundary when iterating for WL (default), 0= no head loss correction
variable_chezy=0; %0=const friction coefficient (default), 1= calculated friction coefficient
dtver=1;%1=dynamic time step based on amount of sedimentation,0=fixed time step

%field parameters based on chadwick et al. compilation
S_island=S_factor* 4.3000e-05;
H_M=21; %main channel depth
cf_channel=.005;
Q_channel=29000;
qs_channel=0.007362733;
B_channel=650;
u_starM=sqrt(cf_channel)*(Q_channel/(H_M*B_channel)); %main channel shear velocity m/s
SLR=SLR_factor*0.0170; %m/yr, use 5/1000 if not specified in chadwick et al.
ds=ds_factor*.3/1000; %grain size (m)

%other parameters
L=H_M*L_star;
lambda=.4;%porosity
R=1*1.65; %specific density
rho=1000;
g=9.81;
timeyear=31557600;

exner_prefactor=1; %used in flow expansion calculation
threshold_factor=1; %multiplies threshold of motion
r_multiplier=7; %ratio of near bed and average concentraiton
cf_island=.01; %cf used in island when assuming constant chezy coeff

%make the grid (r>1 gives tigher grid spacing near upstream island edge)
stretch=ones(1,N);
for i=2:N
    stretch(i)=r*stretch(i-1);
end
dx=stretch*L/sum(stretch);
dx_edge=.5*dx(1:end-1)+.5*dx(2:end);
x=zeros(1,N);
x(1)=dx(1)/2;
for i=2:N;
    x(i)=x(i-1)+dx_edge(i-1);
end

%dimensionalize
H_diff=S_island*L; %difference in water level of the two boundary conditions (m);
b=b_star*H_M;%ref height (m)
initial_offset=offset_star*H_M;

%% initial conditions
if hotstart==0
H_sea=0;
H_main=H_sea+H_diff;
initial_slope=0;
if initial_slope==1
    eta=zeros(1,N);
    eta(1)=H_main-initial_offset;
    for i=2:N;
        eta(i)=eta(i-1)-dx_edge(i-1)*S_island;
    end
else
eta=zeros(1,N)+H_main-initial_offset;%initial condition for eta
end
h=zeros(1,N);
S=zeros(1,N);
c=zeros(1,N);%sed concentration c(x) for each grain size class
E=zeros(1,N);%entrainment
cf=zeros(1,N);%friction coefficient cf(x)
end


C1=18;
C2=1;
viscosity=10^-6;
w_s=R*g*ds.^2./(C1*(1*viscosity)+(.75*C2*g*R.*ds.^3).^.5);%ds-specific settling velocity

cmean_method=2; %2 (default) corrects for log law of wall velocity profile when computing reference concenration c_b, 1 uses rouse directly
if cmean_method==2
    ksM=11*H_M/exp(.41./sqrt(cf_channel));
    q_loglaw=(u_starM./.41)*(log(30*H_M/ksM)*H_M-H_M+ksM/30);
    z_up=linspace(b,H_M,10000);
    qsz=log(30*z_up/ksM).*((H_M./z_up-1)./(H_M/b-1)).^(w_s/(.41*u_starM));
    c_b=1/((u_starM/(.41*qs_channel))*(log(30*b/ksM)*b-b+ksM/30+sum((.5*qsz(2:end)+.5*qsz(1:end-1)).*(z_up(2:end)-z_up(1:end-1)))));
elseif cmean_method==1
    cmean=B_channel*qs_channel/Q_channel;
    z=linspace(b,H_M,10000);%calculating profile just from H_M-h(1) to H_M
    cvert=((H_M./z-1)./(H_M/b-1)).^(w_s/(.41*u_starM));%rouse profile
    c_b=cmean/mean(cvert);
end

    

r_ratio=r_multiplier.*ones(1,length(x));

%% initial calculations

Re_p=sqrt(R*g.*ds).*ds/(viscosity);
tausc=.22*(sqrt(R*g.*ds).*ds./(1*viscosity)).^(-.6)+.06.*10.^(-7.7.*(sqrt(R*g.*ds).*ds/(1*viscosity)).^(-.6));
tausc=tausc*threshold_factor;

q_guess0=1;%used in iteration for q
q_guess=1.1;
q=0;

if variable_chezy==1
    ks=3*ds*ones(1,N);
     Cz=(1/.4)*log(11.*h./ks);
     cf=Cz.^(-2);
elseif variable_chezy==0
    cf=cf_island*ones(1,N);
end
cf_edge=[cf(1) .5*cf(1:N-1)+.5*cf(2:N) cf(N)];

t=zeros(1,timesteps);

%% saves
eta_save=zeros(timesteps/saveinterval,N);
h_save=zeros(timesteps/saveinterval,N);
q_save=zeros(timesteps/saveinterval,1);
c_save=zeros(timesteps/saveinterval,N);
E_save=zeros(timesteps/saveinterval,N);
sea_save=zeros(1,timesteps/saveinterval);
cf_save=zeros(timesteps/saveinterval,N);
t_save=[saveinterval:dt*saveinterval:timesteps*dt];
dt_save=zeros(timesteps/saveinterval,1);
%% big loop for updating model

for j=1:timesteps
    if j>1
    t(j)=t(j-1)+dt;
    end
    %% backwater equation with secant method for discharge

    S(1:N-1)=(eta(1:N-1)-eta(2:N))./dx_edge;
    S(N)=S(N-1);
    h(N)=H_sea-eta(N);
    if q==0|imag(q)~=0
        q_guess0=1*10^-6;
        q_guess=1.1*10^-6;
    end
    q=q_guess0;
    for i=1:N-1 %predictor corrector method to get h(x) from downstream to upstream
        if variable_chezy==1
        dh_pred=(S(N-i)-(((1/.4)*log(11*h(N+1-i)./ks(N+1-i))).^-2)*(q^2)/(g*h(N+1-i)^3))/(1-(q^2)/(g*h(N+1-i)^3));
        h(N-i)=h(N+1-i)-dx_edge(N-i)*dh_pred;
        dh_corr=(S(N-i)-(((1/.4)*log(11*h(N-i)./ks(N-i))).^-2)*(q^2)/(g*h(N-i)^3))/(1-(q^2)/(g*h(N-i)^3));
        h(N-i)=h(N+1-i)-.5*dx_edge(N-i)*(dh_pred+dh_corr);
        else
        dh_pred=(S(N-i)-cf_edge(N-i)*(q^2)/(g*h(N+1-i)^3))/(1-(q^2)/(g*h(N+1-i)^3));
        h(N-i)=h(N+1-i)-dx_edge(N-i)*dh_pred;
        dh_corr=(S(N-i)-cf_edge(N-i)*(q^2)/(g*h(N-i)^3))/(1-(q^2)/(g*h(N-i)^3));
        h(N-i)=h(N+1-i)-.5*dx_edge(N-i)*(dh_pred+dh_corr);
        end
    end
    if momentum_fix==1
        H_guess0=(eta(1)+h(1))+(q^2)/(2*g*h(1)^2);
    else
        H_guess0=(eta(1)+h(1));
    end
    
    q=q_guess;
    for i=1:N-1
        if variable_chezy==1
        dh_pred=(S(N-i)-(((1/.4)*log(11*h(N+1-i)./ks(N+1-i))).^-2)*(q^2)/(g*h(N+1-i)^3))/(1-(q^2)/(g*h(N+1-i)^3));
        h(N-i)=h(N+1-i)-dx_edge(N-i)*dh_pred;
        dh_corr=(S(N-i)-(((1/.4)*log(11*h(N-i)./ks(N-i))).^-2)*(q^2)/(g*h(N-i)^3))/(1-(q^2)/(g*h(N-i)^3));
        h(N-i)=h(N+1-i)-.5*dx_edge(N-i)*(dh_pred+dh_corr);
        else
        dh_pred=(S(N-i)-cf_edge(N-i)*(q^2)/(g*h(N+1-i)^3))/(1-(q^2)/(g*h(N+1-i)^3));
        h(N-i)=h(N+1-i)-dx_edge(N-i)*dh_pred;
        dh_corr=(S(N-i)-cf_edge(N-i)*(q^2)/(g*h(N-i)^3))/(1-(q^2)/(g*h(N-i)^3));
        h(N-i)=h(N+1-i)-.5*dx_edge(N-i)*(dh_pred+dh_corr);
        end
    end
    if momentum_fix==1
        H_guess=(eta(1)+h(1))+(q^2)/(2*g*h(1)^2);
    else
        H_guess=(eta(1)+h(1));
    end

    numit=0;
    while abs(H_guess-H_main)>watlev_accuracy %secant method iteration to find q such that WL at upstream boundary matches prescribed within specified tolerance
        %%
        numit=numit+1;
        H_guess1=H_guess;
        q_guess1=q_guess;
        %q_guess=min(sqrt(g*(H_main-eta(1))^3),max(0,q_guess-(H_guess-H_main)*(q_guess-q_guess0)/((H_guess-H_main)-(H_guess0-H_main))));
        q_guess=min(1.1*q_guess,max(q_guess/2,q_guess-(H_guess-H_main)*(q_guess-q_guess0)/((H_guess-H_main)-(H_guess0-H_main))));

        q=q_guess;
        for i=1:N-1
        if variable_chezy==1
        dh_pred=(S(N-i)-(((1/.4)*log(11*h(N+1-i)./ks(N+1-i))).^-2)*(q^2)/(g*h(N+1-i)^3))/(1-(q^2)/(g*h(N+1-i)^3));
        h(N-i)=h(N+1-i)-dx_edge(N-i)*dh_pred;
        dh_corr=(S(N-i)-(((1/.4)*log(11*h(N-i)./ks(N-i))).^-2)*(q^2)/(g*h(N-i)^3))/(1-(q^2)/(g*h(N-i)^3));
        h(N-i)=h(N+1-i)-.5*dx_edge(N-i)*(dh_pred+dh_corr);
        else
        dh_pred=(S(N-i)-cf_edge(N-i)*(q^2)/(g*h(N+1-i)^3))/(1-(q^2)/(g*h(N+1-i)^3));
        h(N-i)=h(N+1-i)-dx_edge(N-i)*dh_pred;
        dh_corr=(S(N-i)-cf_edge(N-i)*(q^2)/(g*h(N-i)^3))/(1-(q^2)/(g*h(N-i)^3));
        h(N-i)=h(N+1-i)-.5*dx_edge(N-i)*(dh_pred+dh_corr);
        end
        end
        if momentum_fix==1
            H_guess=(eta(1)+h(1))+(q^2)/(2*g*h(1)^2);
        else
            H_guess=(eta(1)+h(1));
        end
        q_guess0=q_guess1;
        H_guess0=H_guess1; 
        %%
        if numit>1000
            q=0;
            break
        end
    end
    if min(h)<=0|q<0|imag(q)~=0|max(abs(imag(q)))~=0|sum(isnan(h))>0|max(abs(imag(h)))>0
        q=0;
        h=zeros(1,N);
    end

    if variable_chezy==1
        ks=3*ds*ones(1,N);
         Cz=(1/.4)*log(11*h./ks);
         cf=Cz.^(-2);
         cf_edge=[cf(1) .5*cf(1:N-1)+.5*cf(2:N) cf(N)];
    end
    %% calculate boundary condition concentration
    if q==0
        c=zeros(1,N);
    else
        if momentum_fix==0
            z=linspace(max(b,H_M-h(1)),H_M,10000);%calculating profile just from H_M-h(1) to H_M
        else
            z=linspace(max(b,H_M-(H_main-eta(1))),H_M,10000);%calculating profile just from H_M-h(1) to H_M
        end
        cvert=((H_M./z-1)./(H_M/b-1)).^(w_s/(.41*u_starM));%rouse profile
        c(1)=c_b*mean(cvert); 
    end
    
    %% calculate entrainment
    tau_star=(cf.*(q./h).^2)./(R*g*ds);
    ks=11*h./exp(.41./sqrt(cf));
    E(tau_star>tausc)=(((.015*ds./max(ks(tau_star>tausc),.01*(h(tau_star>tausc)))).*((tau_star(tau_star>tausc)./tausc-1).^1.5)).*(sqrt(R*g*ds)*ds/(1*viscosity)).^(-.2));
    E(tau_star<=tausc)=0;
    if q==0
        E(:,:)=0;
    end
    
    %% calculate sediment concentration
    if q==0
        c(:,:)=0;
    else
    for i=2:N
    c(i)=c(i-1)-dx(i-1)*w_s*(r_ratio(i-1)*c(i-1)-E(i-1))/q;
    end
    end

    %% update bed and sea
    if q>0
    if dtver==1
    deta=exner_prefactor*dt*((r_ratio.*c-E).*w_s'/(1-lambda));%sedimentation/erosion
    dt=min(dt_initial,dt*deta_thresh./(max(abs(deta)./h(1))));
    deta=exner_prefactor*dt*((r_ratio.*c-E).*w_s'/(1-lambda));% sedimentation/erosion
    
    else
    deta=exner_prefactor*dt*((r_ratio.*c-E).*w_s'/(1-lambda));% sedimentation/erosion    

    end

    else
        dt=dt_initial;
     deta=zeros(size(c));
    end
    eta=eta+deta; %update bed elevation
    H_sea=H_sea+dt*SLR/timeyear;%update sea level
    H_main=H_sea+H_diff;%update channel water level

    
    %% save data
    if rem(j,saveinterval)==0
    eta_save(j/saveinterval,:)=eta;
    h_save(j/saveinterval,:)=h;
    q_save(j/saveinterval)=q;
    c_save(j/saveinterval,:,:)=c';
    E_save(j/saveinterval,:,:)=E';
    sea_save(j/saveinterval)=H_sea;
    cf_save(j/saveinterval,:)=cf;
    t_save(j/saveinterval)=t(j);
    dt_save(j/saveinterval)=dt;
    end
    if rem(j,100)==0
disp([num2str(j) ' ' sprintf('%.15g',q)]) 
    end
    if rem(j,200000)==0&&saveon==1
        save(runid)
    end
end

    h_normal=(cf.*(q.^2)./(g.*S)).^(1/3);
    h_critical=((q.^2)/g).^(1/3);
   q_star=q_save/((g^.5)*(H_M^1.5));
    t_star=t_save*w_s*c_b/(H_M*(1-lambda));
    eta_star=eta/H_M;
    x_star=x/L;


