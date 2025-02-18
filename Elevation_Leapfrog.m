close all
clear all

pathToHere         = mfilename('fullpath');
I = regexp(pathToHere,'\');
pathToHere = pathToHere(1:I(end));
pathToDataFileStorage = [pathToHere 'dataFiles/'];

storageStride  = 500;
dt             = 0.0001;
A  = dir([pathToDataFileStorage '*.txt']);
nFiles = length(A);

Gauges = [0.5, 2, 4, 6, 8, 10, 12, 14, 15];

Elev = zeros(nFiles, length(Gauges));
Time = zeros(nFiles, length(Gauges));
parfor i=1:length(Gauges)
    E = Elev(:,i);
    T = Time(:,i);
    ind = 1;
   
    A = dir([pathToDataFileStorage '*.txt']);
    nFiles = length(A);
    s = readInDataFile([pathToDataFileStorage A(ind).name]);
    I = find(s(:,4)==7);
    Freeparticles = s(I,1:2);
    dx = 0.05;
    
    location = Gauges(i);
    I1 = find(Freeparticles(:,1)>location-dx & Freeparticles(:,1)<= location+dx);
    GridFreeParticles = Freeparticles(I1,1:2);
    initial_height = max(GridFreeParticles(:,2));
    
    E(1) = 0;
    T(1) = 0;
    
    ind1=2;    
    while ind1 <= nFiles
        ind1;
        A = dir([pathToDataFileStorage '*.txt']);
        nFiles = length(A);
        
        s = readInDataFile([pathToDataFileStorage A(ind1).name]);
        
        tStep = (ind1-1)*storageStride;
        tTime = tStep*dt;
        
        I = find(s(:,4)==7);
        Freeparticles = s(I,1:2);

        I1 = find(Freeparticles(:,1)>location-dx & Freeparticles(:,1)<=location+dx);
        GridFreeParticles = Freeparticles(I1,1:2);
        surface_elevation = max(GridFreeParticles(:,2)) - initial_height;
        
        E(ind1) = surface_elevation;
        T(ind1) = tTime;
        
        ind1 = ind1+1;
    end
    Elev(:,i) = E;
    Time(:,i) = T;
end

writematrix(Elev, "Elevation_Leapfrog_shallowWater.txt");
writematrix(Time, "Time_Leapfrog_shallowWater.txt");

% f = fit(Time,Elevation,'smoothingspline','SmoothingParam',0.95);
% h = plot(f);
% xi = get(h,'XData');
% yi = get(h,'YData');
% 
% t = linspace(0,10,500);
% x  = 6.5;
% %ele = a * sin(k*x - w*t) - (k*a^2 /2)*((3 - (tanh(k * d))^2)/(2 * (tanh(k * d))^2)) * cos(2*k*x - 2 *w*t) - k*a^2/(2 *sinh(2*k*d));
% % ele = a * cos(w*t - k*x);
% % plot(t , ele ,'r');
% % hold on
% plot(xi,yi,'k') ;
% grid on;
% legend('Theoretical Stokes 1st order','Simulation results')
% xlabel('Time (seconds)')
% ylabel('\eta (metres)')
% title('Surface elevation comparison between theory and simulation results')
% f = gcf;
% saveas(f,'Elevation_latest.png')
   