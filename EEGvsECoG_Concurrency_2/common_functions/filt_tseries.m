function [sim_data] = filt_tseries(sim_data,bands)
Data1    = sim_data.functional.data1.Data;
Data1    = Data1';
Fs1      = sim_data.functional.data1.Fs;
Data2    = sim_data.functional.data2.Data;
Data2    = Data2';
Fs2      = sim_data.functional.data2.Fs;
Data_TS1 = cell(1,length(bands));
Data_TS2 = cell(1,length(bands));
for band = 1:length(bands)
    disp(strcat('-->> processing band: ', num2str(band)));
    f1             = bands(band,1); %frequency band lower limit
    f2             = bands(band,2); %frequency band upper limit  
    tmp            = bandpass(Data1,[f1 f2],Fs1);
    Data_TS1{band} = tmp';
    tmp            = bandpass(Data2,[f1 f2],Fs2);  
    Data_TS2{band} = tmp';
end
sim_data.functional.data1.Data_TS = Data_TS1';
sim_data.functional.data2.Data_TS = Data_TS2';
end