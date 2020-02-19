function [gammaindeks, gammaNilai] = gammaindeks(roi,tps,epid,dta,dd,sptps,spepid)
%% petunjuk penggunaan
% roi = ukuran lapangan dalam mm
% tps = citra plan
% epid = citra aktual
% spepid = resolusi epid (matrix 2x1)
% sptps = resolusi tps (matrix 2x1)
% dta =  distance to agreement (3 untuk 3%)
% dd = dose difference (0.03 untuk 3 mm)

[rtps,ctps] = size(tps);
[repid,cepid] = size(epid);

tengahytps = round(rtps/2);
tengahxtps = round(ctps/2);

tengahyepid = round(repid/2);
tengahxepid = round(cepid/2);

%% roi 5 cm 
hroi = 0.5*roi;
%spepid = infoepid.ImagePlanePixelSpacing;
%sptps = infotps.PixelSpacing;%%
%spepid=(infoepid.RadiationMachineSAD/infoepid.RTImageSID)*infoepid.ImagePlanePixelSpacing;
kiriepid = tengahxepid - round(hroi/spepid(2));
kananepid = tengahxepid + round(hroi/spepid(2));
atasepid = tengahyepid - round(hroi/spepid(1));
bawahepid = tengahyepid + round(hroi/spepid(1));

kiritps = tengahxtps - round(hroi/sptps(2));
kanantps = tengahxtps + round(hroi/sptps(2));
atastps = tengahytps - round(hroi/sptps(1));
bawahtps = tengahytps + round(hroi/sptps(1));


%% image dalam roi

tpsbaru = zeros(round(roi/sptps(1)),round(roi/sptps(2)));
[ir,ic] = size(tpsbaru);
for irtps = 1:ir
    for ictps = 1:ic
        tpsbaru(irtps,ictps) = tps(atastps+irtps,kiritps+ictps);
    end
end

epidbaru = zeros(round(roi/spepid(1)),round(roi/spepid(2)));
[jr,jc] = size(epidbaru);
for jrepid = 1:jr
    for jcepid = 1:jc
        epidbaru(jrepid,jcepid) = epid(atasepid+jrepid,kiriepid+jcepid);
    end
end

%% === membuat gamma indeks sesuai mm ===================
barisDcm=round(roi/sptps(1));
kolomDcm=round(roi/sptps(2));
gammaNilai = zeros(barisDcm,kolomDcm);
Dta = round(dta/spepid(1));
awalB = 1;
awalK = 1;

for bk = 1:barisDcm
    for kk = 1:kolomDcm
        kriteria = zeros(2*Dta);
           
        %bagian kiri atas
        krib = 1; krik = 1;
        for ik = awalB:-1:(-Dta+awalB+1)
            for jk = awalK:-1:(-Dta+awalK+1)
                if ik < 1 || jk <1
                    break
                end
                rx = abs(awalK-jk);
                ry = abs(awalB-ik);
                r = sqrt(rx^2+ry^2);
                kriteria(krib,krik) = sqrt((r/Dta)^2+...
                    (tpsbaru(ik,jk)-epidbaru(bk,kk))^2/(dd*epidbaru(bk,kk))^2);% 0.03 adalah dc
                
                krik = krik+1;
            end
            krik = 1;
            krib = krib+1;
        end 
        
        %bagian kiri bawah
       for ikkb = awalB+1:(Dta+awalB)
            for jkkb = awalK:-1:(-Dta+awalK+1)
                if ikkb >= barisDcm || jkkb <1
                    break
                end
                rxkb = abs(awalK-jkkb);
                rykb = abs(awalB-ikkb);
                rkb = sqrt(rxkb^2+rykb^2);
                kriteria(krib,krik) = sqrt((rkb/Dta)^2+...
                    (tpsbaru(ikkb,jkkb)-epidbaru(bk,kk))^2/(0.03*epidbaru(bk,kk))^2);
                
                krik = krik+1;
            end
            krik = 1;
            krib = krib+1;
       end
        
        %bagian kanan atas
        krik = Dta+1; krib = 1;
       for ika = awalB:-1:(-Dta+awalB+1)
            for jka = awalK+1:(Dta+awalK)
                if ika <=1 || jka >=kolomDcm
                    break
                end
                rxa = abs(awalK-jka);
                rya = abs(awalB-ika);
                ra = sqrt(rxa^2+rya^2);
                kriteria(krib,krik) = sqrt((ra/Dta)^2+...
                    (tpsbaru(ika,jka)-epidbaru(bk,kk))^2/(0.03*epidbaru(bk,kk))^2);
                
                krik = krik+1;
            end
            krik = Dta+1;
            krib = krib+1;
       end
       
       %bagian kanan bawah
       for ikb = awalB+1:(Dta+awalB)
            for jkb = awalK+1:(Dta+awalK)
                if ikb >=barisDcm || jkb >=kolomDcm
                    break
                end
                rxb = abs(awalK-jkb);
                ryb = abs(awalB-ikb);
                rb = sqrt(rxb^2+ryb^2);
                kriteria(krib,krik) = sqrt((rb/Dta)^2+...
                    (tpsbaru(ikb,jkb)-epidbaru(bk,kk))^2/(0.03*epidbaru(bk,kk))^2);
                
                krik = krik+1;
            end
            krik = Dta+1;
            krib = krib+1;
       end
       
       
        awalK = awalK+1;
        if min(min(kriteria))==0
            [ba,ka,va] = find(kriteria);

%             raw = find(va==min(va));
           
            gammaNilai(bk,kk) = min(va);
        else
            gammaNilai(bk,kk) = min(min(kriteria));
        end
        
    end 
    awalK = 1; 
    awalB = awalB+1;
end

satu = find(gammaNilai<=1);
gammaindeks = length(satu)/(barisDcm*kolomDcm)*100