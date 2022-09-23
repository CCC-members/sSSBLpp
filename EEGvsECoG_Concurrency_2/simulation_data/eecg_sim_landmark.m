function [sim_data] = eecg_sim_landmark(sim_data,bands)
%% Parameters
F_data1        = sim_data.functional.data1.F;
Svvdata1       = sim_data.functional.data1.Svv;
Lvj            = sim_data.structural.L_data13D;
No             = length(Lvj);
W              = sim_data.structural.DNinv;
[Svvdata1,Lvj] = applying_reference(Svvdata1,Lvj);
ind            = sim_data.structural.surface.indL;
ind3D          = [3*ind-2 3*ind-1 3*ind];
ind3D          = ind3D';
ind3D          = ind3D(:);
Ngen           = length(ind);
Lvj            = Lvj(:,ind3D);
W              = W(ind3D,ind);
LvjW           = Lvj*W;
p1             = size(LvjW,1);
q              = size(LvjW,2);
Ip             = eye(p1);
Iq             = eye(q);
scaleL         = sqrt(sum(abs(diag(LvjW*LvjW')))/p1);
LvjW           = LvjW/scaleL;
WLjv           = LvjW';
LvjW2Ljv       = LvjW*WLjv;
gamma_grid     = 0:0.1:2;
gcv            = zeros(length(gamma_grid),1);
%%
for band = 1:size(bands,1)
    F1                     = bands(band,1); %frequency band lower limit
    F2                     = bands(band,2); %frequency band upper limit
    %% data1 Inverse solution filter
    [val,idf1]             = min(abs(F_data1-F1));
    [val,idf2]             = min(abs(F_data1-F2));
    %
    Svvdata1band           = mean(Svvdata1(:,:,idf1:idf2),3);
    scaleV                 = (sum(abs(diag(Svvdata1band)))/p1);
    Svvdata1band           = Svvdata1band/scaleV;
    %
    count                  = 1;
    for gamma = gamma_grid
        Txiv               = Ip - (10^gamma)*LvjW2Ljv + (10^gamma)*LvjW2Ljv*(Ip/((10^gamma)*LvjW2Ljv + Ip))*((10^gamma)*LvjW2Ljv);
        gcv(count)         = (1/p1)*sum(abs(diag(Txiv*Svvdata1band*Txiv')))/((1/p1)*sum(abs(diag(Txiv))))^2;
        disp(['eloreta gcv param_' , num2str(gamma)])
        count              = count + 1;
    end
    %
    [gcv_opt,idx_gamma]    = min(gcv);
    gamma                  = gamma_grid(idx_gamma);
    sigma2jW               = (10^gamma)*W';
    sigma2jWLjv            = (10^gamma)*WLjv;
    LvjWsigma2j            = (10^gamma)*LvjW;
    LvjWsigma2jW           = LvjWsigma2j*W';
    sigma2j_post0          = (W*sigma2jWLjv)/(LvjW*sigma2jWLjv+Ip);
    % Only save the diagonals of the Posterior Covariance
    for count_gen = 1:size(W,1)
        sigma2j_post(count_gen) = W(count_gen,:)*sigma2jW(:,count_gen) - sigma2j_post0(count_gen,:)*LvjWsigma2jW(:,count_gen);
    end
    Tjv                        = (W*sigma2jWLjv-sigma2j_post0*(LvjWsigma2j*WLjv));
    SvvTvj                     = Svvdata1band*Tjv';
    for count_gen = 1:size(Tjv,1)
        s2j(count_gen)         = abs(Tjv(count_gen,:)*SvvTvj(:,count_gen));
    end
    s2j                    = sum(reshape(abs(s2j),3,length(s2j)/3),1)';
    sigma2j_post           = sum(reshape(abs(sigma2j_post),3,length(sigma2j_post)/3),1)';
    stat                   = zeros(No/3,1);
    stat(ind)              = sqrt(s2j./sigma2j_post);
    indms                  = find(stat > 0);
    J                      = zeros(No/3,1);
    s2jfull                = zeros(No/3,1);
    s2jfull(ind)           = s2j;
    J(indms)               = s2jfull(indms);
    T                      = zeros(No,p1);
    T(ind3D,:)             = Tjv;
    %
    Svv1{band}             = Svvdata1band;
    Tdata1{band}           = T;
    Jdata1{band}           = J;
    statdata1{band}        = stat;
    indmsdata1{band}       = indms;
end
sim_data.functional.data1.T         = Tdata1;
sim_data.functional.data1.J         = Jdata1;
sim_data.functional.data1.stat      = statdata1;
sim_data.functional.data1.indms     = indmsdata1;
sim_data.functional.data1.Svvdata1  = Svv1;
end