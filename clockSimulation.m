%%2-D simulation using a hexagonal structure%%
clear;clc;
%drugs=[0.05 0.15 0.45 1.35]; %for Ext. Data Fig. 5 parameter search
%for d=1:4 %for Ext. Data Fig. 5 parameter search
     %drug=1+drugs(d); %for Ext. Data Fig. 5 parameter search
%for d=1:40 %for Ext. Data Fig. 10d parameter search
     %drug=1+d*0.05; %for Ext. Data Fig. 10d parameter search
     drug=1; %1 for No Pulse Treatment %1.25 for Fig.2 10+30 simulation 
     pulse=4; % drug pulse freq in each cycle, for Fig.2 10+30 simulation
     
%for pulse=2:6 %for Ext. Data Fig. 5 parameter search

%with/without clock
clockamp=0; % clockamp=35 for clock, 0 for w/o clock
pcda = 0.4; % Dephosphorylation Rate due to Clock

%%space-time settings%%
L=44; % PSM size in cell #
prec=30; % posterior clock period in mins
perc=1.333; % drug pulse cycle / 30 min ratio

Vg = 0.055; % Tailbud Growth Speed
delt=prec/71;
tcell=delt/Vg; %shift clock matrix due to tail elongation
x=1:L;
R = 4; % PSM width in cell #
tailbud = 8; % # of tailbud cells
T = 2218; % simulation duration

randDirection = zeros(10^8, 1); % pre-allocates random direction matrix

% Clock's coupling to ppERK
%pcd=pcda*exp(-x./6+0.1); % use to simulate clock's effect on ERK
%exponentially decreasing along the PSM
pcd=pcda*ones(1,L);


%%Clock Matrix Settings%%
CRMatrix = zeros(L,T); % empty clock matrix
phi=zeros(L,T); %clock phase matrix
f=1./(1+exp(3*(x-35)/L));
f=2*pi/prec*f;
for i=1:T
    
    CRMatrix(:,i)=clockamp*(sin(phi(:,i))+1)/2; % clock signal
    phi(:,i+1)=phi(:,i)+delt*f(:);
    
    if mod(i,round(tcell/delt,0))==0
        phi=circshift(phi,1,1);
        phi(1,:)=phi(2,:);
    end
end

%%Simulation%%
tic

pb = .01; % Ligand Receptor Binding Rate
Diff = 2; % Ligand Diffusion Speed

pid  = 0.4; % Dephosphorylation Rate due to Inhibitor
ppp = 1.0; % Phosphorylation Rate

% description of phosphorylation rate with drug pulses
pppm=ppp*ones(1,T);
for t=15:T
    if mod(t-14,71*perc)<floor(71*perc/pulse)
        pppm(t)=pppm(t)/drug;
    end
end

pip = 0.5; % Inhibitor Protein Translation Rate
pir = -0.5; % Inhibitor Protein Decay Rate / Translation Rate

trlatedel = 14; % Trl+Secr Delay of FGF
actdel = 14; % Activation Delay from Complex to ERK
inhdel = 42; % Inhibitor Feedback (Trx+Trl) Delay

randDirection = randi(3, [10^8 1]); % direction for neighborFinder, only enough for 1 parameter set
dirCounter = 0; % Chooses which random num

pbd = -1; % Comp Decay Rate
rd = -0.025; % Ligand RNA Decay Rate
pfp = 1; % Free Ligand Protein Translation Rate
pfd = -0.25; % Free Ligand Protein Decay Rate

pd = -0.1; % Inactive Signaling Protein Decay Rate
pp = 12.5; % Inactive Signaling Protein Translation Rate
rp = 10; % Ligand RNA Transcription Rate

r0 = 0; % Initial Ligand RNA
pf0 = 0; % Initial Free Ligand Protein
pb0 = 0; % Initial Bound Ligand Protein
pi0 = 0; % Initial Inhibitor Protein
pp0 = 0; % Initial Active Signaling Protein
p0 = 0; % Initial Inactive Signaling Protein

RMatrix = r0 * ones(L, T, R); % Ligand RNA Matrix (mLIG)
PFMatrix = r0 * ones(L, T, R); % Ligand Protein Matrix (LIG)
PBMatrix = r0 * ones(L, T, R); % Bound Ligand Receptor Complex Matrix (COMP)
PIMatrix = r0 * ones(L, T, R); % Inhibitor Protein Matrix (INH)
PMatrix = r0 * ones(L, T, R); % Inactive Signaling Protein Matrix (SIG-)
PPMatrix = r0 * ones(L, T, R); % Active Signaling Protein Matrix (SIG+)
StatMatrix = PPMatrix; % Stationary (Tailbud Frame) Matrix

ry0matrix = r0 * ones(L, 1, R);
pfy0matrix = pf0 * ones(L, 1, R);
pby0matrix = pb0 * ones(L, 1, R);
piy0matrix = pi0 * ones(L, 1, R);
ppy0matrix = pp0 * ones(L, 1, R);
y0matrix = p0 * ones(L, 1, R);

for ti = 0:T-1
    %Euler's Method for dy/dt=f(y,t)
    N = 100;
    for j = 1:R
        for k = 1:tailbud
            
            tf = ti + 1; % initial and final ti
            h = (tf - ti) / N; % Time step
            ry = zeros(N + 1, 1);
            pfy = zeros(N + 1, 1);
            pby = zeros(N + 1, 1);
            piy = zeros(N + 1, 1);
            py = zeros(N + 1, 1);
            ppy = zeros(N + 1, 1);
            
            if ti == 0
                
                ry(1) = ry0matrix(k, 1, j);
                pfy(1) = pfy0matrix(k, 1, j);
                pby(1) = pby0matrix(k, 1, j);
                piy(1)=piy0matrix(k,1,j);
                py(1) = y0matrix(k, 1, j);
                ppy(1) = ppy0matrix(k, 1, j);
                
            else
                
                ry(1) = RMatrix(k, ti, j);
                pfy(1) = PFMatrix(k, ti, j);
                pby(1) = PBMatrix(k, ti, j);
                piy(1) = PIMatrix(k, ti, j);
                py(1) = PMatrix(k, ti, j);
                ppy(1) = PPMatrix(k, ti, j);
                
            end
            
            for n = 1:N
                
                ry(n + 1) = ry(n) + h * rp + h * rd * ry(n);
                
                if ti > actdel && ti > trlatedel
                    
                    pfy(n + 1) = pfy(n) + h * (pfp * RMatrix(k, ti - trlatedel, j) + pfd * pfy(n) - pb * pfy(n));
                    pby(n + 1) = pby(n) + h * (pb * pfy(n) + pbd * pby(n));
                    
                else
                    
                    pfy(n + 1) = pfy(n);
                    pby(n + 1) = pby(n);
                    
                end
                
                if ti > inhdel
                    
                    piy(n + 1) = piy(n) + h * (pip * PPMatrix(k, ti - inhdel, j) + pir * pip * piy(n));
                    
                    if py(n) < ppy(n)
                        
                        ppy(n + 1) = ppy(n) + h * (pppm(ti) * PBMatrix(k, ti - actdel, j) * py(n) + pd * ppy(n) - pid * piy(n) * ppy(n) - pcd(k) * CRMatrix(k, ti) * ppy(n));
                        py(n + 1) = py(n) + ppy(n) + h * (pp + pd * (py(n) + ppy(n))) - ppy(n + 1);
                        
                    else
                        
                        py(n + 1) = py(n) + h * (pp + pcd(k) * CRMatrix(k, ti) * ppy(n) + pid * piy(n) * ppy(n) + pd * py(n) - pppm(ti) * PBMatrix(k, ti - actdel, j) * py(n));
                        ppy(n + 1) = ppy(n) + py(n) + h * (pp + pd * (ppy(n) + py(n))) - py(n + 1);
                        
                    end
                    
                    
                    if py(n + 1) < 0 || py(n + 1) > 50000 || ppy(n + 1) < 0 || ppy(n + 1) > 50000
                        
                        py(n + 1) = 0;
                        ppy(n + 1) = 0;
                        
                    end
                    
                else
                    
                    piy(n + 1) = piy(n);
                    ppy(n + 1) = ppy(n);
                    py(n + 1) = py(n);
                    
                end
                
            end
            
            ry0 = ry(N + 1);
            pfy0 = pfy(N + 1);
            pby0 = pby(N + 1);
            piy0 = piy(N + 1);
            py0 = py(N + 1);
            ppy0 = ppy(N + 1);
            
            ry0matrix(k, 1, j) = ry0;
            pfy0matrix(k, 1, j) = pfy0;
            pby0matrix(k, 1, j) = pby0;
            piy0matrix(k, 1, j) = piy0;
            y0matrix(k, 1, j) = py0;
            ppy0matrix(k, 1, j) = ppy0;
            
            RMatrix(k, tf, j) = ry0;
            PFMatrix(k, tf, j) = pfy0;
            PBMatrix(k, tf, j) = pby0;
            PIMatrix(k, tf, j) = piy0;
            PMatrix(k, tf, j) = py0;
            PPMatrix(k, tf, j)=ppy0;
        end
    end
end


r0 = ry0matrix(tailbud, 1, j);
pf0 = pfy0matrix(tailbud, 1, j);
pb0 = pby0matrix(tailbud, 1, j);
pi0 = piy0matrix(tailbud, 1, j);
pp0 = ppy0matrix(tailbud, 1, j);
p0 = y0matrix(tailbud, 1, j);


for j = 1:R 
    for k = tailbud + 1:L
        
        N = 100; % Number of ti steps
        ry0 = r0;
        pfy0 = pf0;
        pby0 = pb0;
        piy0 = pi0;
        py0 = p0;
        ppy0 = pp0;
        
        past = floor((k - tailbud) / Vg);
        for ti = 1:past
            
            %Euler's Method for dy/dt=f(y,t)
            tf = ti + 1; % initial and final ti
            h = (tf - ti) / N; % Time step
            ry = zeros(N + 1, 1);
            pfy = zeros(N + 1, 1);
            pby = zeros(N + 1, 1);
            piy = zeros(N + 1, 1);
            py = zeros(N + 1, 1);
            ppy = zeros(N + 1, 1);
            
            ry(1) = ry0;
            pfy(1) = pfy0;
            pby(1) = pby0;
            piy(1) = piy0;
            py(1) = py0;
            ppy(1) = ppy0;
            
            for n = 1:N
                
                ry(n + 1) = ry(n) + h * rd * ry(n);
                
                if ti > actdel && ti > trlatedel
                    
                    pfy(n + 1) = pfy(n) + h * (pfp * RMatrix(k, ti - trlatedel, j) + pfd * pfy(n) - pb * pfy(n));
                    pby(n + 1) = pby(n) + h * (pb * pfy(n) + pbd * pby(n));
                    
                else
                    
                    pfy(n + 1) = pfy(n);
                    pby(n + 1) = pby(n);
                    
                end
                
                if ti > inhdel
                    
                    piy(n + 1) = piy(n) + h * (pip * PPMatrix(k, ti - inhdel, j) + pir * pip * piy(n));
                    
                    if py(n) < ppy(n)
                        
                        ppy(n + 1) = ppy(n) + h * (pppm(ti) * PBMatrix(k, ti - actdel, j) * py(n) + pd * ppy(n) - pid * piy(n) * ppy(n) - pcd(k) * CRMatrix(k, ti) * ppy(n));
                        py(n + 1) = py(n) + ppy(n) + h * (pp + pd * (py(n) + ppy(n))) - ppy(n + 1);
                        
                    else
                        
                        py(n + 1) = py(n) + h * (pp + pcd(k) * CRMatrix(k, ti) * ppy(n) + pid * piy(n) * ppy(n) + pd * py(n) - pppm(ti) * PBMatrix(k, ti - actdel, j) * py(n));
                        ppy(n + 1) = ppy(n) + py(n) + h * (pp + pd * (ppy(n) + py(n))) - py(n + 1);
                        
                    end
                else
                    
                    piy(n + 1) = piy(n);
                    ppy(n + 1) = ppy(n);
                    py(n + 1) = py(n);
                    
                end
                
            end
            
            ry0 = ry(N + 1);
            pfy0 = pfy(N + 1);
            pby0 = pby(N + 1);
            piy0 = piy(N + 1);
            py0 = py(N + 1);
            ppy0 = ppy(N + 1);
            
            ry0matrix(k, 1, j) = ry0;
            pfy0matrix(k, 1, j) = pfy0;
            pby0matrix(k, 1, j) = pby0;
            piy0matrix(k, 1, j) = piy0;
            y0matrix(k, 1, j) = py0;
            ppy0matrix(k, 1, j) = ppy0;
            
            RMatrix(k, tf, j) = ry0;
            PFMatrix(k, tf, j) = pfy0;
            PBMatrix(k, tf, j) = pby0;
            PIMatrix(k, tf, j) = piy0;
            PMatrix(k, tf, j) = py0;
            PPMatrix(k, tf, j) = ppy0;
            
        end
    end
end

for ti = 0:T - 1
    %Euler's Method for dy/dt=f(y,t)
    
    for j = 1:R
        for k = 1:L - 1
            
            tf = ti + 1; % initial and final ti
            h = (tf - ti) / N; % Time step
            ry = zeros(N + 1, 1);
            pfy = zeros(N + 1, 1);
            pby = zeros(N + 1, 1);
            piy = zeros(N + 1, 1);
            py = zeros(N + 1, 1);
            ppy = zeros(N + 1, 1);
            
            if ti == 0
                
                ry(1) = ry0matrix(k, 1, j);
                pfy(1) = pfy0matrix(k, 1, j);
                pby(1) = pby0matrix(k, 1, j);
                piy(1) = piy0matrix(k, 1, j);
                py(1) = y0matrix(k, 1, j);
                ppy(1) = ppy0matrix(k, 1, j);
                
            else
                
                ry(1) = RMatrix(k, ti, j);
                pfy(1) = PFMatrix(k, ti, j);
                pby(1) = PBMatrix(k, ti, j);
                piy(1) = PIMatrix(k,ti,j);
                py(1) = PMatrix(k, ti, j);
                ppy(1) = PPMatrix(k, ti, j);
            end
            
            for n = 1:N
                if k <= tailbud
                    ry(n + 1) = ry(n) + h * (rp + rd * ry(n));
                else
                    ry(n + 1) = ry(n) + h * rd * ry(n);
                end
                
                if ti > trlatedel && ti > actdel
                    value = pfy(n);
                    
                    dirCounter = dirCounter+1;
                    direction = randDirection(dirCounter);
                    c = ceil((k - 1) / L);
                    
                    switch direction % determines adjacent cells for hexagonal matrix
                        case 1 % 12 o'clock/6 o'clock
                            
                            if k == 1
                                neighbor1 = 0;
                            else
                                neighbor1 = PFMatrix(k - 1, ti, j); % 12 o'clock
                            end
                            
                            if k >= L
                                neighbor2 = 0;
                            else
                                neighbor2 = PFMatrix(k + 1, ti, j); % 6 o'clock
                            end
                            
                        case 2 % 2 o'clock/8 o'clock
                            
                            if (mod(j, 2) == 0 && k == 1) || j == R
                                neighbor1 = 0; % right (2 o'clock) not possible
                            elseif mod(j,2) == 0
                                neighbor1 = PFMatrix(k - 1, ti, j + 1); % 2 o'clock (even j)
                            else
                                neighbor1 = PFMatrix(k, ti, j + 1); % 2 o'clock (odd j)
                            end
                            
                            if j == 1 || (mod(j, 2) ~= 0 && k >= L)
                                neighbor2 = 0; % left (8 o'clock) not possible
                            elseif mod(j, 2) == 0
                                neighbor2 = PFMatrix(k, ti, j - 1); % 8 o'clock (even j)
                            else
                                neighbor2 = PFMatrix(k + 1, ti, j-1); % 8 o'clock (odd j)
                            end
                            
                        case 3 % 4 o'clock/10 o'clock
                            
                            if (mod(j, 2) ~= 0 && k >= L) || j == R
                                neighbor1 = 0; % right (4 o'clock) not possible
                            elseif mod(j, 2) == 0
                                neighbor1 = PFMatrix(k, ti, j + 1); % 4 o'clock (even j)
                            else
                                neighbor1 = PFMatrix(k + 1, ti, j + 1); % 4 o'clock (odd j)
                            end
                            
                            if j == 1 || (mod(j, 2) == 0 && k == 1)
                                neighbor2 = 0; % left (10 o'clock) not possible
                            elseif mod(j, 2) == 0
                                neighbor2 = PFMatrix(k - 1, ti, j - 1); % 10 o'clock (even j)
                            else
                                neighbor2 = PFMatrix(k, ti, j - 1); % 10 o'clock (odd j)
                            end
                    end
                    
                    switch c % k==1
                        case 0
                            if j == 1
                                switch direction
                                    case 1
                                        x = ((2 * neighbor2) - 2 * value);
                                    otherwise
                                        x = ((2 * neighbor1) - 2 * value);
                                end
                            elseif j == R
                                if mod(j, 2) == 0 && direction == 3
                                    x = 0;
                                else
                                    x = ((2 * neighbor2) - 2 * value);
                                end
                            else
                                if mod(j, 2) == 0
                                    switch direction
                                        case 3
                                            x = ((2 * neighbor1) - 2 *value);
                                        otherwise
                                            x = ((2 * neighbor2) - 2 * value);
                                    end
                                else
                                    switch direction
                                        case 1
                                            x = ((2 * neighbor2) - 2 * value);
                                        otherwise
                                            x = (neighbor1 + neighbor2 - 2 * value);
                                    end
                                end
                            end
                            
                        case 2 % k > L (k>M)
                            if j == 1
                                switch direction
                                    case 3
                                        x = (-2 * value);
                                    otherwise
                                        x = (neighbor1 - 2 * value);
                                end
                            elseif j == R
                                if mod(j, 2) == 0
                                    switch direction
                                        case 1
                                            x = (neighbor1 - 2 * value);
                                        otherwise
                                            x = (neighbor2 - 2 * value);
                                    end
                                else
                                    switch direction
                                        case 1
                                            x = (neighbor1 - 2 * value);
                                        case 2
                                            x = (-2 * value);
                                        case 3
                                            x = (neighbor2 - 2 * value);
                                    end
                                end
                            else
                                if mod(j, 2) == 0
                                    switch direction
                                        case 1
                                            x = (neighbor1 - 2 * value);
                                        otherwise
                                            x = (neighbor1 + neighbor2 - 2 * value);
                                    end
                                else
                                    switch direction
                                        case 3
                                            x = (neighbor2 - 2 * value);
                                        otherwise
                                            x = (neighbor1 - 2 * value);
                                    end
                                end
                            end
                            
                        otherwise % k is between 1 and L-1
                            if j == 1
                                switch direction
                                    case 1
                                        x = (neighbor1 + neighbor2 - 2 * value);
                                    otherwise
                                        x = (2 * neighbor1 - 2 * value);
                                end
                            elseif j == R
                                switch direction
                                    case 1
                                        x = (neighbor1 + neighbor2 - 2 * value);
                                    otherwise
                                        x = (2 * neighbor2 - 2 * value);
                                end
                            else
                                x = (neighbor1 + neighbor2 - 2 * value);
                            end
                            
                    end
                    pfy(n + 1) = pfy(n) + h * (pfp * RMatrix(k, ti - trlatedel, j) + pfd * pfy(n) - pb * pfy(n) + Diff * x);
                    pby(n + 1) = pby(n) + h * (pb * pfy(n) + pbd * pby(n));
                    
                else
                    pfy(n + 1) = pfy(n);
                    pby(n + 1) = pby(n);
                end
                
                    
               
                if ti > inhdel
                    piy(n + 1) = piy(n) + h * (pip * PPMatrix(k, ti - inhdel, j) + pir*pip * piy(n));
                    
                    if py(n) < ppy(n)
                        ppy(n+1)=ppy(n)+h*(pppm(ti)*PBMatrix(k, ti - actdel, j)*py(n)+pd*ppy(n)-pid*piy(n)*ppy(n) - pcd(k) * CRMatrix(k,ti) * ppy(n));
                        py(n+1)=py(n) + ppy(n)+h*(pp+pd*(py(n)+ppy(n))) - ppy(n+1);
                        
                    else
                        py(n+1)=py(n)+h*(pp + pcd(k) * CRMatrix(k,ti) * ppy(n) + pid*piy(n)*ppy(n)+pd*py(n)-pppm(ti)*PBMatrix(k, ti - actdel, j)*py(n));
                        ppy(n+1)=ppy(n)+py(n)+h*(pp + pd*(ppy(n)+py(n))) - py(n+1);
                    end
                    
                else
                    piy(n + 1) = piy(n);
                    ppy(n + 1) = ppy(n);
                    py(n + 1) = py(n);
                end
                
            end
            
            
            
            ry0 = ry(N + 1);
            pfy0 = pfy(N + 1);
            pby0 = pby(N + 1);
            piy0 = piy(N + 1);
            py0 = py(N + 1);
            ppy0 = ppy(N + 1);
            
            ry0matrix(k, 1, j) = ry0;
            pfy0matrix(k, 1, j) = pfy0;
            pby0matrix(k, 1, j) = pby0;
            piy0matrix(k, 1, j) = piy0;
            y0matrix(k, 1, j) = py0;
            ppy0matrix(k, 1, j) = ppy0;
            
            RMatrix(k, tf, j) = ry0;
            PFMatrix(k, tf, j) = pfy0;
            PBMatrix(k, tf, j) = pby0;
            PIMatrix(k, tf, j) = piy0;
            PMatrix(k, tf, j) = py0;
            PPMatrix(k, tf, j) = ppy0;
            
        end
        k = L;
        if ti > 0
            tf = ti + 1; % initial and final ti
            h = (tf - ti) / N; % Time step
            ry = zeros(N + 1, 1);
            pfy = zeros(N + 1, 1);
            pby = zeros(N + 1, 1);
            piy = zeros(N + 1, 1);
            py = zeros(N + 1, 1);
            ppy = zeros(N + 1, 1);
            
            
            ry(1) = RMatrix(k, ti, j);
            pfy(1) = PFMatrix(k, ti, j);
            pby(1) = PBMatrix(k, ti, j);
            piy(1) = PIMatrix(k, ti, j);
            py(1) = PMatrix(k, ti, j);
            ppy(1) = PPMatrix(k, ti, j);
            
            
            for n = 1:N
                
                ry(n + 1) = ry(n) + h * rd * ry(n);
                
                if ti > trlatedel && ti > actdel
                    value = pfy(n);
                    dirCounter = dirCounter + 1;
                    direction = randDirection(dirCounter);
                    c = ceil((k - 1) / L);
                    
                    switch direction % determines adjacent cells for hexagonal matrix
                        case 1 % 12 o'clock/6 o'clock
                            
                            if k == 1
                                neighbor1 = 0;
                            else
                                neighbor1 = PFMatrix(k - 1, ti, j); % 12 o'clock
                            end
                            
                            if k >= L
                                neighbor2 = 0;
                            else
                                neighbor2 = PFMatrix(k + 1, ti, j); % 6 o'clock
                            end
                            
                        case 2 % 2 o'clock/8 o'clock
                            
                            if (mod(j, 2) == 0 && k == 1) || j == R
                                neighbor1 = 0; % right (2 o'clock) not possible
                            elseif mod(j, 2) == 0
                                neighbor1 = PFMatrix(k - 1, ti, j + 1); % 2 o'clock (even j)
                            else
                                neighbor1 = PFMatrix(k, ti, j + 1); % 2 o'clock (odd j)
                            end
                            
                            if j == 1 || (mod(j, 2) ~= 0 && k >= L)
                                neighbor2 = 0; % left (8 o'clock) not possible
                            elseif mod(j, 2) == 0
                                neighbor2 = PFMatrix(k, ti, j - 1); % 8 o'clock (even j)
                            else
                                neighbor2 = PFMatrix(k + 1, ti, j-1); % 8 o'clock (odd j)
                            end
                            
                        case 3 % 4 o'clock/10 o'clock
                            
                            if (mod(j, 2) ~= 0 && k >= L) || j == R
                                neighbor1 = 0; % right (4 o'clock) not possible
                            elseif mod(j, 2) == 0
                                neighbor1 = PFMatrix(k, ti, j + 1); % 4 o'clock (even j)
                            else
                                neighbor1 = PFMatrix(k + 1, ti, j + 1); % 4 o'clock (odd j)
                            end
                            
                            if j == 1 || (mod(j, 2) == 0 && k == 1)
                                neighbor2 = 0; % left (10 o'clock) not possible
                            elseif mod(j, 2) == 0
                                neighbor2 = PFMatrix(k - 1, ti, j - 1); % 10 o'clock (even j)
                            else
                                neighbor2 = PFMatrix(k, ti, j - 1); % 10 o'clock (odd j)
                            end
                    end
                    
                    switch c % k==1
                        case 0
                            if j == 1
                                switch direction
                                    case 1
                                        x = ((2 * neighbor2) - 2 * value);
                                    otherwise
                                        x = ((2 * neighbor1) - 2 * value);
                                end
                            elseif j == R
                                if mod(j, 2) == 0 && direction == 3
                                    x = 0;
                                else
                                    x = ((2 * neighbor2) - 2 * value);
                                end
                            else
                                if mod(j, 2) == 0
                                    switch direction
                                        case 3
                                            x = ((2 * neighbor1) - 2 * value);
                                        otherwise
                                            x = ((2 * neighbor2) - 2 * value);
                                    end
                                else
                                    switch direction
                                        case 1
                                            x = ((2 * neighbor2) - 2 * value);
                                        otherwise
                                            x = (neighbor1 + neighbor2 - 2 * value);
                                    end
                                end
                            end
                            
                        case 2 % k > L (k>M)
                            if j == 1
                                switch direction
                                    case 3
                                        x = (-2 * value);
                                    otherwise
                                        x = (neighbor1 - 2 * value);
                                end
                            elseif j == R
                                if mod(j, 2) == 0
                                    switch direction
                                        case 1
                                            x = (neighbor1 - 2 * value);
                                        otherwise
                                            x = (neighbor2 - 2 * value);
                                    end
                                else
                                    switch direction
                                        case 1
                                            x = (neighbor1 - 2 * value);
                                        case 2
                                            x = (-2 * value);
                                        case 3
                                            x = (neighbor2 - 2 * value);
                                    end
                                end
                            else
                                if mod(j, 2) == 0
                                    switch direction
                                        case 1
                                            x = (neighbor1 - 2 * value);
                                        otherwise
                                            x = (neighbor1 + neighbor2 - 2 * value);
                                    end
                                else
                                    switch direction
                                        case 3
                                            x = (neighbor2 - 2 * value);
                                        otherwise
                                            x = (neighbor1 - 2 * value);
                                    end
                                end
                            end
                            
                        otherwise % k is between 1 and L-1
                            if j == 1
                                switch direction
                                    case 1
                                        x = (neighbor1 + neighbor2 - 2 * value);
                                    otherwise
                                        x = (2 * neighbor1 - 2 * value);
                                end
                            elseif j == R
                                switch direction
                                    case 1
                                        x = (neighbor1 + neighbor2 - 2 * value);
                                    otherwise
                                        x = (2 * neighbor2 - 2 * value);
                                end
                            else
                                x = (neighbor1 + neighbor2 - 2 * value);
                            end
                            
                    end
                    pfy(n + 1) = pfy(n) + h * (pfp * RMatrix(k, ti - trlatedel, j) + pfd * pfy(n) - pb * pfy(n) + Diff * x);
                    pby(n +1 ) = pby(n) + h * (pb * pfy(n) + pbd * pby(n));
                    
                else
                    pfy(n + 1) = pfy(n);
                    pby(n + 1) = pby(n);
                end
                
                
                if ti > inhdel
                    piy(n + 1) = piy(n) + h * (pip * PPMatrix(k, ti - inhdel, j) + pir*pip * piy(n));
                    
                    if py(n) < ppy(n)
                        ppy(n+1)=ppy(n)+h*(pppm(ti)*PBMatrix(k, ti - actdel, j)*py(n)+pd*ppy(n)-pid*piy(n)*ppy(n) - pcd(k) * CRMatrix(k,ti) * ppy(n));
                        py(n+1)=py(n) + ppy(n)+h*(pp+pd*(py(n)+ppy(n))) - ppy(n+1);
                        
                    else
                        py(n+1)=py(n)+h*(pp + pcd(k) * CRMatrix(k,ti) * ppy(n) + pid*piy(n)*ppy(n)+pd*py(n)-pppm(ti)*PBMatrix(k, ti - actdel, j)*py(n));
                        ppy(n+1)=ppy(n)+py(n)+h*(pp + pd*(ppy(n)+py(n))) - py(n+1);
                    end
                    
                else
                    piy(n + 1) = piy(n);
                    ppy(n + 1) = ppy(n);
                    py(n + 1) = py(n);
                end
                
            end
            
            
            
            ry0 = ry(N + 1);
            pfy0 = pfy(N + 1);
            pby0 = pby(N + 1);
            piy0 = piy(N + 1);
            py0 = py(N + 1);
            ppy0 = ppy(N + 1);
            
            ry0matrix(k, 1, j) = ry0;
            pfy0matrix(k, 1, j) = pfy0;
            pby0matrix(k, 1, j) = pby0;
            piy0matrix(k, 1, j) = piy0;
            y0matrix(k, 1, j) = py0;
            ppy0matrix(k, 1, j) = ppy0;
            
            RMatrix(k, tf, j) = ry0;
            PFMatrix(k, tf, j) = pfy0;
            PBMatrix(k, tf, j) = pby0;
            PIMatrix(k, tf, j) = piy0;
            PMatrix(k, tf, j) = py0;
            PPMatrix(k, tf, j) = ppy0;
            
        end
    end
    if mod(tf, floor(1 / Vg)) == 0
        StatMatrix(:, tf + 1 - floor(1 / Vg):tf, :) = PPMatrix(:, tf + 1 - floor(1 / Vg):tf, :);
        RMatrix = circshift(RMatrix, [1 0 0]); % shifts matrix 1 position along 1st dimension
        PFMatrix = circshift(PFMatrix, [1 0 0]);
        PBMatrix = circshift(PBMatrix, [1 0 0]);
        PIMatrix = circshift(PIMatrix, [1 0 0]);
        PMatrix = circshift(PMatrix, [1 0 0]);
        PPMatrix = circshift(PPMatrix, [1 0 0]);
        
        PFMatrix(1, :, :) = PFMatrix(2,:,:); % accounts for matrix shift
        RMatrix(1, :, :)=RMatrix(2,:,:);
        PIMatrix(1, :, :)=PIMatrix(2,:,:);
        PBMatrix(1, :, :)=PBMatrix(2,:,:);
        PMatrix(1, :, :)=PMatrix(2,:,:);
        PPMatrix(1, :, :)=PPMatrix(2,:,:);
    end
end

    zmatrix = zeros(L - tailbud, T);
    for i = 1:L - tailbud
        for j = 1:T
            zmatrix(i, j ) = PPMatrix(tailbud + i - 1, j)- 0.02 * PPMatrix(tailbud, j);
        end
    end

toc

%PPMatrix=zmatrix;

%%
%{
Keep Here for Normalization of ppERK Matrix
for i=1:size(PPMatrix,2)
    PPMatrix(:,i)=PPMatrix(:,i)/max(PPMatrix(1,:));
end
%}
%%
for j=T+1:-1:T-400
    PPMatrix(1:floor((T+1-j)*Vg),j)=NaN;
end

SFC=zmatrix;
SFC=SFC*0;
for i=1:size(zmatrix,2)
    j=1;
    while SFC(j,i)<=1 && j<size(zmatrix,1)
        j=j+1;
        SFC(j,i)=max(zmatrix(j-1,i),zmatrix(j,i))/min(zmatrix(j-1,i),zmatrix(j,i))-1;
    end
    SFC(1,i)=SFC(2,i);
    siz=size(SFC(SFC(:,i)<=0,i));
    x=NaN(siz);
    SFC(SFC(:,i)<=0,i)=x;
    
    SFC(isnan(SFC(:,i)),i)=1;
    
    siz2=size(SFC(SFC(:,i)>=1,i));
    z=ones(siz2);
    SFC(SFC(:,i)>=1,i)=z;
    
end

for j=T+1:-1:T-400
    SFC(1:floor((T+1-j)*Vg),j)=NaN;
    for tb=1:5
        if isnan(SFC(floor((T+1-j)*Vg)+tb,j))
            SFC(floor((T+1-j)*Vg)+tb,j)=0;
        end
    end
end


SFC0=SFC/0.44;
SFC0(SFC0>1)=1;
SFC1=SFC0;
SFC1(isnan(SFC1))=0;
H = fspecial('disk',3);
SFC1 = imfilter(SFC1,H,'replicate'); 
SFC1(isnan(SFC0))=NaN;

osc=ceil((1-min(smooth(StatMatrix(1,T-300:T-100),20))/max(smooth(StatMatrix(1,T-300:T-100),20)))*100);

abcds=SFC1(:,T-400:T-10)*100;
[X Y] = meshgrid(T-400:T-10,1:L-tailbud);
[xx yy] = meshgrid(T-400:71/30:T-10,1:1/7:L-tailbud);
abcs=interp2(X,Y,abcds,xx,yy);
figure
colormap(jet(64));
[Ms,cs]=contourf(abcs,[50 50]);

ylabel('PSM (um)');

% Create xlabel
xlabel('Time (min)');
%% Create title
%title(strcat('Clock A=',num2str(clockamp),' str=',num2str(pcda),' Osc=',num2str(ceil((1-min(StatMatrix(1,T-300:T-100))/max(StatMatrix(1,T-300:T-100)))*100)),'%'),'FontName','Arial', 'FontSize', 28);
%title(strcat('SFC clock=',num2str(clockamp),' pcd=',num2str(pcd)),'FontName','Arial', 'FontSize', 22);
%title(strcat('contour SFC drug=',num2str(drug),' pulse=',num2str(round(perc*30/pulse,0))),'FontName','Arial', 'FontSize', 22);
%%

%% Pick one from below to save SFC threshold border into folder
%savet=strcat('contour_SFC_clock_',num2str(clockamp),'_osc_',num2str(osc),'_',num2str(pcd),'.tif');
%savet=strcat('contour_SFC_clock_constant_',num2str(clockamp),'_',num2str(pcda),'_osc_',num2str(osc),'_Diff_',num2str(Diff),'.tif');
%savet=strcat('contour_SFC_drug_',num2str(drug),'_pulse_',num2str(round(perc*30/pulse,0)),'.tif');

%saveas(gcf,savet);
%%

%open the lut for SFC colors
load('matlab_SFCLUT.mat');

figure
mesh(SFC1(:,T-400:T-10))

%use below to convert any mesh plot into 2-D figures
fig = gcf;
axObjs = gca;
dataObjs = axObjs.Children;
x = dataObjs(1).XData;
y = dataObjs(1).YData;
zdata1 = dataObjs(1).ZData;

% Create figure
figure('InvertHardcopy','off','Color',[1 1 1],...
    'OuterPosition',[276 287 576 513]);
colormap(LUTSFC); %import the lut first!

% Create axes
axes1=gca;

% Create surf
surf1 = surf(zdata1,zdata1,'FaceLighting','none','EdgeLighting','flat',...
    'FaceColor',[1 1 1],...
    'EdgeColor','interp');


% Create ylabel
ylabel('PSM (cells)');

% Create xlabel
xlabel('Time (min)');

%% Create title
%%use for different inhibitory drug pulses
%strdrug=strcat('Pulse d=1/',num2str(pulse),' drug=',num2str((drug-1)/0.5));
%title(strdrug,'FontName','Arial', 'FontSize', 28); 
%title(strcat('SFC clock=',num2str(clockamp),' Pcd=',num2str(pcda)),'FontName','Arial', 'FontSize', 28);
%title(strcat('SFC drug=',num2str(drug),' pulse=',num2str(round(perc*30/pulse,0))),'FontName','Arial', 'FontSize', 22);

%%use for different clock parameters
%strclock=strcat('with Clock A=',num2str(clockamp),' Pcd=',num2str(pcda)); 
%title(strclock,'FontName','Arial', 'FontSize', 28);

%%use for raw ppERK plots
%title('Raw ppERK','FontName','Arial', 'FontSize', 28); 
%%

axes1=gca;
xlim(axes1,[0 400]);
ylim(axes1,[0 40]);
zlim(axes1,[-1 1]);
view(90,90);
% Set the remaining axes properties
set(gca,'BoxStyle','full','FontName','Arial','FontSize',20,'LineWidth',...
    1.5,'XGrid','on','XTick',[0 71 142 213 284 355],'YGrid','on','YTick',...
    [0 5 10 15 20 25 30 35 40]);
xticklabels({'0','30','60','90','120','150'});
caxis([0 1]);
% Create colorbar
colorbar(gca,'FontSize',16);
colorbar('Ticks',[0,0.25,0.5,0.75,1],'TickLabels',{'0','0.11','0.22','0.33','0.44'}); %use for rescaled SFC plots

%% Pick one from below to save SFC kymograph into folder
%savet2=strcat('SFC_clock_',num2str(clockamp),'_osc_',num2str(osc),'_',num2str(pcd),'.tif');
%savet2=strcat('SFC_clock_constant_',num2str(clockamp),'_',num2str(pcda),'_osc_',num2str(osc),'_Diff_',num2str(Diff),'.tif');
%savet2=strcat('SFC_drug_',num2str(drug),'_pulse_',num2str(round(perc*30/pulse,0)),'.tif');

%saveas(gcf,savet2);
%%


figure
mesh(PPMatrix(:,T-400:T-10))

%use below to convert any mesh plot into 2-D figures
fig = gcf;
axObjs = gca;
dataObjs = axObjs.Children;
x = dataObjs(1).XData;
y = dataObjs(1).YData;
zdata1 = dataObjs(1).ZData;

% Create figure
figure('InvertHardcopy','off','Color',[1 1 1],...
    'OuterPosition',[276 287 576 513]);
colormap(jet(64)); %import the lut first!

% Create axes
axes1=gca;

% Create surf
surf1 = surf(zdata1,zdata1,'FaceLighting','none','EdgeLighting','flat',...
    'FaceColor',[1 1 1],...
    'EdgeColor','interp');


% Create ylabel
ylabel('PSM (cells)');

% Create xlabel
xlabel('Time (min)');

%% Create title
%%use for different inhibitory drug pulses
%strdrug=strcat('Pulse d=1/',num2str(pulse),' drug=',num2str((drug-1)/0.5));
%title(strdrug,'FontName','Arial', 'FontSize', 28); 
%title(strcat('ppERK inh=',num2str(round(drug-1,2)),' pulse=',num2str(round(30*perc/pulse,1)),"'-",num2str(round(30*perc*(1-1/pulse),1)),"'"),'FontName','Arial', 'FontSize', 28);
%title(strcat('ppERK drug=',num2str(drug),' pulse=',num2str(round(perc*30/pulse,0))),'FontName','Arial', 'FontSize', 22);

%%use for different clock parameters
%strclock=strcat('with Clock A=',num2str(clockamp),' Pcd=',num2str(pcda)); 
%title(strclock,'FontName','Arial', 'FontSize', 28);

%%use for raw ppERK plots
%title('Raw ppERK','FontName','Arial', 'FontSize', 28);
%%
axes1=gca;
xlim(axes1,[0 400]);
ylim(axes1,[0 40]);
view(90,90);
% Set the remaining axes properties
set(gca,'BoxStyle','full','FontName','Arial','FontSize',20,'LineWidth',...
    1.5,'XGrid','on','XTick',[0 71 142 213 284 355],'YGrid','on','YTick',...
    [0 5 10 15 20 25 30 35 40]);
xticklabels({'0','30','60','90','120','150'});
caxis([0 1]); %60 for not normalized ppERK, 1 for normalized ppERK
% Create colorbar
colorbar(gca,'FontSize',16);
colorbar('Ticks',[0,20,40,60],'TickLabels',{'0','20','40','60'}); %use for rescaled SFC plots

%% Pick one from below to save ppERK kymograph into folder
%savet3=strcat('ppERK_clock_',num2str(clockamp),'_pcd_',num2str(pcd),'.tif');
%savet3=strcat('ppERK_clock_constant_',num2str(clockamp),'_',num2str(pcda),'_osc_',num2str(osc),'_Diff_',num2str(Diff),'.tif');
%savet3=strcat('NppERK_clock_constant_',num2str(clockamp),'_',num2str(pcda),'_osc_',num2str(osc),'_Diff_',num2str(Diff),'.tif');
%savet3=strcat('ppERK_drug_',num2str(drug),'_pulse_',num2str(round(perc*30/pulse,0)),'.tif');
%savet3=strcat('NppERK_drug_',num2str(drug),'_pulse_',num2str(round(perc*30/pulse,0)),'.tif');

%%saveas(gcf,savet3);
%%

close all
   % end %for Ext. Data Fig. 5 parameter search
%end %for Ext. Data Fig. 5 parameter search
%end %for Ext. Data Fig. 10d parameter search
