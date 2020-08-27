%change in concentration per time step
%intialize: T = 700, 20000s. T = 800K, 5000s. 900K, 1000s
Tcount = 1;
for T = 900:100:1200
T
deltax = 5; %nm
timesteps = 1000000;
sites = 200;
index = 1:timesteps;
space = 1:sites;
k = 8.62*10^-5; %eV
epsilon = 10; %dielectric constant
e0 = 8.85*10^-12;
electroncharge = 1.6*10^-19;
bandgap = 3.5*electroncharge; %V
iconc = 2*10^16; %cm^-3
MgHconc = 2*10^19;
Mgconc = 2*10^16;
bulkH = 1.998*10^19;
Ea = 1; %eV
hiconc = bulkH*exp(-Ea/(k*T))*MgHconc/Mgconc;
hconc = bulkH*exp(-Ea/(k*T))*MgHconc/Mgconc*ones(sites,1);
global chargedepth;
chargedepth = round(10^9*sqrt(2*epsilon*e0*bandgap/(electroncharge^2*iconc*10^6)));
%chargedepth = 440; %nm

Ed = 0.7; %eV
D0 = 1.2*10^-3; %cm^2/s
D0 = 10^4; %nm^2/ 10^-5s 
Dh = D0*exp(-Ed/(k*T));
%Dh = 1;
%hconc = ones(sites,1);
hconcmat = ones(sites,timesteps);
Efmat = ones(sites,timesteps);
bulkHmat = ones(timesteps,1);
hconc(1) = 0;

%Efmax = 2/250000;
lattice = linspace(1,sites*deltax-(deltax-1),sites);
Efpinmax = electroncharge^2*Mgconc*10^6/(2*epsilon*e0);
Efpin = 1.7; %eV

count = 0;
for t = 1:timesteps
    Efmax = Efpin/(chargedepth+1)^2;
    Efpinmax = electroncharge^2*Mgconc*10^6/(2*epsilon*e0);
    Ef = Efpinmax*(lattice-chargedepth).^2; %eV
    %set right side of fermi level potential parabola  = 0
    Ef(floor(chargedepth/deltax):sites) = 0;

    href = hconc(sites-1);
    [hconc,forces] = newficks(hconc,Ef,sites,T,deltax,Dh);
    deltabulkH = 10^5*(-href + hconc(sites-1));
    if deltabulkH < 0;
        MgHconc = MgHconc + deltabulkH;
        Mgconc = Mgconc - deltabulkH;
        bulkH = bulkH + deltabulkH;
    end
    Efpinmax = electroncharge^2*Mgconc*10^6/(2*epsilon*e0);
    hconc(sites-1) = bulkH*exp(-Ea/(k*T))*MgHconc/Mgconc;
    chargedepth = round(10^9*sqrt(2*epsilon*e0*bandgap/(electroncharge^2*Mgconc*10^6)));
    hconcmat(:,t) = hconc;
    Efmat(:,t) = Ef;
    count = count + 1;
    bulkHmat(t) = bulkH;
    
 
end

hsum = 0;
for i  = 1:length(hconc)
    hsum = hsum + hconc(i);
end
%fractionretained = hsum/(sites)
Hconcentration(T/100-8) = hconcmat(sites,timesteps);
fractionremaining(T/100-8) = hconcmat(sites,timesteps)/hiconc;
hchange(T/100-8) = hiconc - hconcmat(sites,timesteps);
fractionhbulkchange(T/100-8) = bulkHmat(timesteps)/bulkHmat(1);

hconcmatT(:,:,Tcount) = hconcmat;
EfmatT(:,:,Tcount) = Efmat;

finall(T/100-8,:) = hconcmat(sites,:);
finallbulk(:,T/100-8) = bulkHmat;
Tcount = Tcount + 1;
end

%%
%plots
figure(1) %plot h concentration at desired time 'plottime'
plottime = 1000;
plot(lattice,hconcmatT(:,plottime,1),lattice,hconcmatT(:,plottime,2),lattice,hconcmatT(:,plottime,3),lattice,hconcmatT(:,plottime,4))
legend('900K','1000K','1100K','1200K')
title(['Hydrogen concentration after ' num2str(plottime) ' seconds'])
xlabel('Lattice Depth (nm)')
ylabel('Hydrogen Concentration (atoms/cm^3)')

sz = size(hconcmat);
%plot h concentration as it evolves with time
maxtime = 100000; %set length of anneal here
figure(2)
for k = 1:100:maxtime %sz(2)
    hold off
    plot(lattice,hconcmatT(:,k,1),lattice,hconcmatT(:,k,2),lattice,hconcmatT(:,k,3),lattice,hconcmatT(:,k,4))
    ylim([0, 2*10^18]); %comment out this line to have the axis limits rescale with concentration 
    legend('900K','1000K','1100K','1200K')
    title(['Anneal time = ' num2str(k-1) ' sec'])
    xlabel('Lattice Depth (nm)')
    ylabel('Hydrogen Concentration (atoms/cm^3)')
    drawnow
end