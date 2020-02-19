function [gam, wgamma1] = gammamahdi(mapplan, mapactual, dta, dc,resolusiplan, resolusiactual,treshold)
% Muhammad Mahdi Ramadhan, Fisika 2015, FMIPA, Universitas Indonesia
% petunjuk penggunaan:
% fungsi untuk menghitung gamma index
% mapplan = citra dicom planning
% mapactual = citra dicom yang hasil pengukuran
% dta = distance to agreement (misal 3)
% dc =  dosedifference (misal 0.03)
%%
resolusiplan_x = resolusiplan(2);
resolusiplan_y = resolusiplan(1);
resolusiactual_x = resolusiactual(2);
resolusiactual_y = resolusiactual(1);
pixel = ceil(dta/resolusiactual_y);
[row,col] = find(mapactual(:,:));
x1 = min(col)+20;
x2 = max(col)-20;
y1 = min(row)+20;
y2 = max(row)-20;

[rowa,cola] = size(mapactual);
[rowp,colp] = size(mapplan);
 wgamma1 = zeros(rowa,cola);%sesuai pengukuran aktual
wadah = zeros(rowp,colp);%sesuai plan
wadah(:,:) = 100;
Npass = 0;
Nfailed = 0;
mapplan = mapplan - treshold;
mapactual = mapactual -treshold;
mapplan = max(mapplan,0.0001);
mapactual = max(mapactual,0.0001);
%% start counting
for i = y1 : y2
    for j = x1 : x2
        ya = resolusiactual_y*i-resolusiactual_y/2;
        xa = resolusiactual_x*j-resolusiactual_x/2;
        indexpy = round(ya/resolusiplan_y);
        indexpx = round(xa/resolusiplan_x);
        doseplan = mapplan(indexpy,indexpx);
            for k = indexpy-pixel-1:indexpy+pixel+1
                for l = indexpx-pixel-1:indexpx+pixel+1
                    deltadose = abs((mapplan(k,l) - mapactual(i,j))/doseplan);
                    yp = resolusiplan_y*k-resolusiplan_y/2;
                    xp = resolusiplan_x*l-resolusiplan_x/2;
                    deltad = hypot(yp-ya,xp-xa);
                    gvalue = sqrt((deltadose/dc)^2+(deltad/dta)^2);
                    wadah(k,l) = gvalue;
                end
            end
        if min(min(wadah)) < 1
            Npass = Npass +1;
            wgamma1(i,j) = 60;
        else
                Nfailed = Nfailed +1;
                wgamma1(i,j) = 100;
        end
        wadah(:,:) = 100;
    end
   sprintf('iterasi ke %i',i)
end
gam = 100*(Npass/(Npass+Nfailed));
end