function [hconc,forces] = newficks(hconc,Ef,sites,T,deltax,Dh)
deltax = 5;
k = 8.62*10^-5; %eV
Dspec = 1;
forces = zeros(sites,3);
global chargedepth

for i = 2:sites-1
    
    hdif = (1/deltax)^2*(hconc(i-1)+hconc(i+1)-2*hconc(i));
    fermi = 1/(deltax)^2*((hconc(i)+hconc(i+1))/2*(Ef(i+1)-Ef(i))/(k*T)...
    - (hconc(i-1)+hconc(i))/2*(Ef(i)-Ef(i-1))/(k*T));
    forces(i,2) = hdif;
    forces(i,3) = fermi;
    
    if i == floor(chargedepth/deltax)
        fermi = 0;
    end
    deltaH = Dh*(hdif + fermi);
    
    hconc(i) = hconc(i) + deltaH;
    forces(i,1) = deltaH;
    
end

deltaH = Dh*(hdif);
hconc(sites) = hconc(sites-1);