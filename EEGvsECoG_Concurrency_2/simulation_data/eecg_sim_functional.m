function [sim_data] = eecg_sim_functional(sim_data,preproced_data_path)
%% Loading data1
preproced_data = load(fullfile(preproced_data_path,'ECoG02notch.mat'));
Data           = preproced_data.ECoGnotch;
%% data1
disp('-->> Creating data1 cross-spectrum');
deltaf      = 0.25;
Fs          = 1000;
Nw          = 3;
rej_chan    = sim_data.structural.rej_chan;
Data(rej_chan,:) = [];
[F,Svv,PSD,Data_FC,Ntpoints,Nsegments] = cross_spectrum(Data,deltaf,Fs,Nw);

%% Saving data1
sim_data.functional.data1.F            = F;
sim_data.functional.data1.Fs           = Fs;
sim_data.functional.data1.Svv          = Svv;
sim_data.functional.data1.PSD          = PSD;
sim_data.functional.data1.Data_FC      = Data_FC;
sim_data.functional.data1.Data         = Data;
sim_data.functional.data1.deltaf       = deltaf;
sim_data.functional.data1.Ntpoints     = Ntpoints;
sim_data.functional.data1.Nsegments    = Nsegments;

%% Loading data2
preproced_data = load(fullfile(preproced_data_path,'EEG02notch.mat'));
Data           = preproced_data.EEGnotch;
%% data2
disp('-->> Creating data2 cross-spectrum');
deltaf      = 0.25;
Fs          = 1000;
Nw          = 6;
[F,Svv,PSD,Data_FC,Ntpoints,Nsegments] = cross_spectrum(Data,deltaf,Fs,Nw);

%% Saving data2
sim_data.functional.data2.F            = F;
sim_data.functional.data2.Fs           = Fs;
sim_data.functional.data2.Svv          = Svv;
sim_data.functional.data2.PSD          = PSD;
sim_data.functional.data2.Data_FC      = Data_FC;
sim_data.functional.data2.Data         = Data;
sim_data.functional.data2.deltaf       = deltaf;
sim_data.functional.data2.Ntpoints     = Ntpoints;
sim_data.functional.data2.Nsegments    = Nsegments;
end