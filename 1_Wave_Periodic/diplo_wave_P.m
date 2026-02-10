%Πεπερασμένες Διαφορές με RK4 για την ΜΔΕ Κύματος με Περιοδικές

clear;
clc;
clf;

%Διάφορα Global που θα χρειαστώ
%global c

%Διάφορες σταθερές που θα χρειαστώ
c=1; %σταθερά κύματος
A=1; %σταθερά πλάτους κύματος

%Ορισμός των παραμέτρων του προβλήματος
L = 2;
x0 = -L; % αρχική θέση
xend = L; % τελική θέση
t0 = 0; %αρχικός χρόνος
tend = 5; %τελικός χρόνος

%Απόσταση βήματος χώρου/χρόνου
dx = 0.01;
dt = 0.001;

%Διαμερισμός χώρου/χρόνου και πλήθος σημείων
x = x0:dx:xend;
nx = length(x);

t = t0:dt:tend;
nt = length(t);

%Κριτήριο CFL
CFL = c * dt / dx;

% Έλεγχος του κριτηρίου CFL
if CFL > 1
    error('Το κριτήριο CFL δεν ικανοποιείται: CFL = %f', CFL);
else 
    fprintf('Το κριτήριο CFL ικανοποιείται: CFL = %f\n', CFL);
end

fprintf('Αριθμός σημείων στο διάστημα x: %d\n', nx);
fprintf('Αριθμός σημείων στο διάστημα t: %d\n', nt);

%Αρχική Συνθήκη και Αρχική Ταχύτητα
u0=A*sin(pi/L*x).'; %διάνυσμα στήλης
v0=-A*pi/L*c*cos(pi/L*x).'; %διάνυσμα στήλης

%Αρχικοποίηση πίνακα λύσης
w=zeros(2*nx,nt);
w(:,1)=[u0;v0];


%-----------------------Μέθοδος Runge Kutta 4------------------------------

%Η λύση αποθηκεύεται ανά στήλη
for i=1:nt-1
    w(:,i+1)=RK4(w(:,i),t,dx,dt,c);
end

%Εξαγωγή λύσεων θέσης/ταχύτητας
u=w(1:nx,:);
v=w(nx+1:2*nx,:);

%save('wave_D_rk4.mat','x','t','u');
writematrix(u, 'u_wave_p.csv');
%writematrix(v, 'v_wave_p.csv');
%writematrix(x', 'x_wave_p.csv');
%writematrix(t', 't_wave_p.csv');


%-------------------Grafikes Parastaseis-----------------------------------

dntdraw=500;
for i=1:dntdraw:nt
    plot(x,u(:,i), 'c', 'DisplayName', 'Eksisosi Kimatos', 'LineWidth', 1);
    xlabel('x');
    ylabel('u');
    grid on;
    %legend('Box', 'on', 'Location','north');
    title(['Εξίσωση Κύματος t = ',num2str(t(i))]);
    axis([x0,xend,-1.1,1.1]);
    drawnow;
end


%------------------Ipoprogrammata Methodon---------------------------------

function z=RK4(w,t,dx,dt,c)
    k1 = f(w, t, dx, c);
    k2 = f(w + 0.5*dt*k1, t + 0.5*dt, dx, c);
    k3 = f(w + 0.5*dt*k2, t + 0.5*dt, dx, c);
    k4 = f(w + dt*k3, t + dt, dx, c);

    %τελικό βήμα RK4
    z = w + dt/6 * (k1 + 2*k2 + 2*k3 + k4);
end

function dw=f(w,t,dx,c)

    N=length(w)/2; %χωρισμός στήλης στην μέση
    u=w(1:N); %πρώτο μισό στήλης η θέση
    v=w(N+1:2*N); %δεύτερο μισό στήλης η ταχύτητα
    
    %du=[v(1:N)]; %το ίδιο με το από κάτω αφού είναι ήδη υπολογισμένο
    du=[v(1:N-1);v(1)]; %αποτέλεσμα υπολογισμού πρώτου μισού στήλης

    up=[u(2:N);u(2)];
    uu=[u(1:N)];
    um=[u(N-1);u(1:N-1)];

    dv=c^2.*(up-2*uu+um)/dx^2; %αποτελέσματα υπολογισμού δεύτερου μισού στήλης

    dw=[du;dv]; %αποτέλεσμα

end
