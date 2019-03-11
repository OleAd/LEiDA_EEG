function [eigenvectors, design]=LEiDA_EEG_eigenvectors(data, group, FS, N_areas, window_size, frequency)

if contains(frequency,'alpha')
    flp=8;
    fhi=12;
elseif contains(frequency,'beta')
    flp=12;
    fhi=30;
elseif contains(frequency,'gamma')
    flp=30;
    fhi=80;
else
    error('Choose frequency correctly.')   
end

if mod(size(data,2), window_size)~=0
    warning('Window size does not fit data size, will discard last incomplete window.')
end




disp('Filtering data.');

TR=1/FS;
Tmax=size(data, 2);
eigenvectors=zeros((Tmax-2), N_areas);  

fnq=1/(2*TR);                 % Nyquist frequency
Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
k=6;                          % 6th order butterworth filter
[bfilt,afilt]=butter(k,Wn);   % construct the filter
clear fnq flp fhi Wn k

design=[];
t_all=0;


Phase_Data=zeros(N_areas, size(data,2));

for seed=1:N_areas
    data(seed,:)=data(seed,:)-mean(data(seed,:));
    signal_filt=filtfilt(bfilt,afilt,data(seed,:));
    Phase_Data(seed,:)=angle(hilbert(signal_filt));
end

repetitions=floor(size(data,2)/window_size);
repArray=1:window_size:size(data,2);

disp('Calculating eigenvectors');
tic
for t=2:repetitions-1

    iFC=zeros(N_areas);
    for n=1:N_areas
        for p=1:N_areas
          iFC(n,p)=circ_mean((cos(Phase_Data(n,repArray(t-1):repArray(t))-Phase_Data(p,repArray(t-1):repArray(t))))');
        end
    end
    
    [V1,~]=eigs(iFC,1);
    t_all=t_all+1;
    eigenvectors(t_all,:)=V1;
    
end
elapsedTime=toc;
msg=sprintf('This subject took %d seconds', round(elapsedTime));
disp(msg)
design=repelem(group, size(eigenvectors,1));

end