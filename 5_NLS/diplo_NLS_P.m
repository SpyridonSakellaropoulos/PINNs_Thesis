%% Προσέγγιση NLS με RK4 για σύγκιση με Pytorch

%έχω θέσει ψ=u+iv

clear;
clc;
clf;
close all;
format short;

%% Σταθερές NLS εξίσωσης
p=0.5;                        %συντελεστής διασποράς
q=1.0;                        %συντελεστής μη γραμμικότητας

psi_zero = 1.0;               %παράμετρος σχήματος / εύρους του σολιτονίου
A = sqrt(2*p/q)*psi_zero;     %μέγιστο πλάτος σολιτονίου
c = 0.5;                      %σχετίζεται με την φάση/ταχύτητα του σολιτονίου

%% Ορισμός χρόνου και χώρου του προβλήματος
t_min = 0.0;           %αρχικός χρόνος
t_max = 5.0;           %τελικός χρόνος
dt = 0.001;           %βήμα χρόνου
t = t_min:dt:t_max;    %διαμερισμός χρόνου
nt = length(t);        %πλήθος σημείων χρόνου

L = 40.0;              %συνολικό πλάτος χώρου
L_min = -L/2;          %πλάτος χώρου ελάχιστο
L_max = L/2;           %πλάτος χώρου μέγιστο
dx = sqrt(2*p*dt)*1.5; %βήμα χώρου σεβούμενο το CFL dt<dx^2/(2p)
x = L_min:dx:L_max;    %διαμερισμός χώρου
nx = length(x);        %πλήθος σημείων χώρου


%% Αρχική συνθήκη σε πραγματικό και φανταστικό μέρος
psi_re0 = (A * sech(psi_zero*x) .* cos(c*x)).'; %πραγματικό μέρος Α.Σ.
psi_im0 = (A * sech(psi_zero*x) .* sin(c*x)).'; %φανταστικό μέρος Α.Σ.

max(psi_re0)
min(psi_re0)
max(psi_im0)
min(psi_im0)


%% Αρχικοποίηση πίνακα λύσης
psi_re = zeros(nx,nt); %αρχικοποίηση πραγματικού μέρους λύσης
psi_im = zeros(nx,nt); %αρχικοποίηση φανταστικού μέρους λύσης

psi_re(:,1) = psi_re0; %πραγματικό μέρος Α.Σ.
psi_im(:,1) = psi_im0; %φανταστικό μέρος Α.Σ.


%% Μέθοδος Runge Kutta 4
for i=1:nt-1
    %psi_re(:,i+1)=RK4re(psi_re(:,i),t,dx,dt,p,q,psi_im(:,i));
    %psi_im(:,i+1)=RK4im(psi_im(:,i),t,dx,dt,p,q,psi_re(:,i));

    [psi_re(:,i+1), psi_im(:,i+1)] = RK4_all(psi_re(:,i), psi_im(:,i), dx, dt, p, q);
end


%% Υπολογισμός απόλυτης τιμής λύσης
psi_abs = sqrt(psi_re.^2 + psi_im.^2);

max(psi_abs(:))
min(psi_abs(:))

writematrix(psi_abs, 'psi_NLS_p.csv');


%% Γραφική παράσταση
figure;
dntdraw=500;
for i=1:dntdraw:nt
    plot(x,psi_abs(:,i), 'c', 'DisplayName', 'Eksisosi NLS', 'LineWidth', 1);
    xlabel('x');
    ylabel('|ψ|');
    grid on;
    %legend('Box', 'on', 'Location','north');
    title(['Εξίσωση NLS t = ',num2str(t(i))]);
    axis([L_min,L_max,0, max(psi_abs(:))*1.1]);
    drawnow;
end

%% 3D γραφική παράσταση |ψ(x,t)| 
figure;
step_t = 500; % πάρε κάθε 500 χρονικά βήματα
psi_plot = psi_abs(:,1:step_t:end); % μικρότερη μήτρα
t_plot = t(1:step_t:end);

[X,T] = meshgrid(x, t_plot);
surf(X, T, psi_plot');  % transpose όπως πριν
shading interp;
colormap jet;
colorbar;
xlabel('x');
ylabel('t');
zlabel('|ψ(x,t)|');
title('3D απεικόνιση λύσης NLS');
view(45,30);
axis tight;
zlim([0, max(psi_plot(:))*1.1]);

% figure;
% plot(x, psi_re0, 'b', 'LineWidth', 1.5); hold on;
% plot(x, psi_im0, 'r', 'LineWidth', 1.5);
% xlabel('x');
% ylabel('\psi');
% legend('Re(\psi)', 'Im(\psi)');
% title('Αρχική συνθήκη NLS: ψ(x,0)');
% grid on;


%% Υποπρογράμματα Μεθόδων


function [u_ip1, v_ip1]=RK4_all(u, v, dx, dt, p, q)
    [uk1, vk1] = f(u             , v             , dx, p,q);
    [uk2, vk2] = f(u + 0.5*dt*uk1, v + 0.5*dt*vk1, dx, p,q);
    [uk3, vk3] = f(u + 0.5*dt*uk2, v + 0.5*dt*vk2, dx, p,q);
    [uk4, vk4] = f(u + dt*uk3    , v + dt*vk3    , dx, p,q);
    
    %τελικό βήμα RK4
    u_ip1 = u + dt/6 * (uk1 + 2*uk2 + 2*uk3 + uk4);
    v_ip1 = v + dt/6 * (vk1 + 2*vk2 + 2*vk3 + vk4);
end

function [du, dv] = f(u, v, dx, p, q)
    N=length(u); %μήκος των u και v
    
    %Πεπερασμένες διαφορές για u_xx
    up=[u(2:N);u(2)];
    %up=[u(2:N);u(1)];
    uu=[u(1:N)];
    um=[u(N-1);u(1:N-1)];
    %um=[u(N);u(1:N-1)];

    u_xx = (up-2*uu+um) / dx^2;

    %Πεπερασμένες διαφορές για v_xx
    vp=[v(2:N);v(2)];
    %vp=[v(2:N);v(1)];
    vv=[v(1:N)];
    vm=[v(N-1);v(1:N-1)];
    %vm=[v(N);v(1:N-1)];

    v_xx = (vp-2*vv+vm) / dx^2;

    
    %Υπολογισμών των v_t και u_t
    dv =  p * u_xx + q * (u.^2 + v.^2) .* u;
    du = -p * v_xx - q * (u.^2 + v.^2) .* v;
end

%% Τα βασανισμένα μου λάθη

% function z=RK4re(psi_re,t,dx,dt,p,q,psi_im)
%     k1 = fre(psi_re            , t         , dx, p,q,psi_im);
%     k2 = fre(psi_re + 0.5*dt*k1, t + 0.5*dt, dx, p,q,psi_im);
%     k3 = fre(psi_re + 0.5*dt*k2, t + 0.5*dt, dx, p,q,psi_im);
%     k4 = fre(psi_re + dt*k3    , t + dt    , dx, p,q,psi_im);
% 
%     %τελικό βήμα RK4
%     z = psi_re + dt/6 * (k1 + 2*k2 + 2*k3 + k4);
% end
% 
% function z=RK4im(psi_im,t,dx,dt,p,q,psi_re)
%     k1 = fim(psi_im            , t         , dx, p,q,psi_re);
%     k2 = fim(psi_im + 0.5*dt*k1, t + 0.5*dt, dx, p,q,psi_re);
%     k3 = fim(psi_im + 0.5*dt*k2, t + 0.5*dt, dx, p,q,psi_re);
%     k4 = fim(psi_im + dt*k3    , t + dt    , dx, p,q,psi_re);
% 
%     %τελικό βήμα RK4
%     z = psi_im + dt/6 * (k1 + 2*k2 + 2*k3 + k4);
% end

% function dv=fre(u, ~, dx, p,q,v)
%     N=length(u);
% 
%     %up=[u(2:N);u(2)];
%     up=[u(2:N);u(1)];
%     uu=[u(1:N)];
%     %um=[u(N-1);u(1:N-1)];
%     um=[u(N);u(1:N-1)];
% 
%     oros_1 = p .* (up-2*uu+um) / dx^2;
%     oros_2 = q .* (u.^2 + v.^2).*u;
% 
%     dv = oros_1 + oros_2;
% end
% 
% function du=fim(v, ~, dx, p,q,u)
%     N=length(v);
% 
%     %vp=[v(2:N);v(2)];
%     vp=[v(2:N);v(1)];
%     vv=[v(1:N)];
%     %vm=[v(N-1);v(1:N-1)];
%     vm=[v(N);v(1:N-1)];
% 
%     oros_1 = p .* (vp-2*vv+vm) / dx^2;
%     oros_2 = q .* (u.^2 + v.^2).*v;
% 
%     du = -oros_1 - oros_2;
% end
