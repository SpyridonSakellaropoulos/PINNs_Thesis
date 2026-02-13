%Ψευδο-φασματική μέθοδος με RK4 για την KdV Κύματος με Περιοδικές

clear;
clc;
clf;

%Διάφορα Global που θα χρειαστώ
%global 

%Διάφορες σταθερές που θα χρειαστώ
c=1; %ταχύτητα σολιτονίου
p=6; %συντελεστής μη γραμμικότητας
q=1; %συντελεστής διασποράς

%Ορισμός των παραμέτρων του προβλήματος
L = 24;
x0 = -L/2; % αρχική θέση
xend = L/2; % τελική θέση
t0 = 0; %αρχικός χρόνος
tend = 5; %τελικός χρόνος
N = 124; %αριθμός διακριτοποίησης χώρου 64

%Απόσταση βήματος χώρου/χρόνου
dt = 0.0001;
dx = L/N;

nt = round(tend/dt); %αριθμός σημείων χρόνου

%Διαμερισμός χώρου και πλήθος σημείων
x = (-L/2:dx:L/2-dx)'; 
disp(length(x));

%φασματικά (Fourier) κύματα
k = [0:N/2-1 -N/2:-1]'*2*pi/L; 
k1 = 1i*k;     %τελεστής πρώτης παραγώγου
k2 = -k.^2;    %τελεστής δεύτερης παραγώγου
k3 = -1i*k.^3; %τελεστής τρίτης παραγώγου

% Αρχικοποίηση της αρχικής συνθήκης u0
u = c/2 * sech(sqrt(c)/2 * x).^2; % παράδειγμα αρχικής συνθήκης

%Αρχικοποίηση πίνακα λύσης
udata = u; 
tdata = 0;

%RK4
for nn = 1:nt                    
    du1 = real(-ifft(k3.*fft(u))) - real(3*ifft(k1.*fft(u.^2)));  
    v = u + 0.5*du1*dt;
    du2 = real(-ifft(k3.*fft(v))) - real(3*ifft(k1.*fft(v.^2)));  
    v = u + 0.5*du2*dt;
    du3 = real(-ifft(k3.*fft(v))) - real(3*ifft(k1.*fft(v.^2)));  
    v = u + du3*dt;
    du4 = real(-ifft(k3.*fft(v))) - real(3*ifft(k1.*fft(v.^2)));           
%     du1 = -ifft(k3.*fft(u)) + real(6*ifft(k1.*fft(u)).*u);  
%     v = u + 0.5*du1*dt;
%     du2 = -ifft(k3.*fft(v)) + real(6*ifft(k1.*fft(v)).*v);  
%     v = u + 0.5*du2*dt;
%     du3 = -ifft(k3.*fft(v)) + real(6*ifft(k1.*fft(v)).*v);  
%     v = u + du3*dt;
%     du4 = -ifft(k3.*fft(v)) + real(6*ifft(k1.*fft(v)).*v);
    u = u + (du1 + 2*du2 + 2*du3 + du4)*dt/6;
    if mod(nn, 20) == 0
       udata = [udata u]; 
       tdata = [tdata nn*dt];
    end
 end
 
u_abs = abs(udata');
disp(size(u_abs', 1)); %αριθμός γραμμών
disp(size(u_abs', 2)); %αριθμός στηλών
writematrix(u_abs', 'u_kdv_p.csv');

%γραφική παράσταση λύσης
waterfall(x, tdata, u_abs);           % solution plotting
colormap(jet(128)); view(10, 60)
xlabel('x')
ylabel('t')
zlabel('|u|')
axis([-L/2 L/2 0 max(tdata) min(u_abs(:)) max(u_abs(:))]); 
grid on;
ax = gca;
ax.XTick = (-L/2):2:(L/2);
ax.YTick = (0):1:(tend);
ax.ZTick = min(u_abs(:)):0.05:max(u_abs(:));
