function Main_test(n)
delete(gcp('nocreate'));
if isempty(gcp('nocreate'))
    c=parcluster;
    Corenum=c.NumWorkers;
    threads=c.NumThreads;
    c.NumThreads = 1;
    parpool(Corenum);%36
else
    delete(gcp('nocreate'));
    c=parcluster;
    Corenum=c.NumWorkers;
    threads=c.NumThreads;
    c.NumThreads=1;
    %     c.NumWorkers=Corenum;%
    parpool(44);%36
end

% Random stream
stream=RandStream('mrg32k3a','Seed','shuffle');%'mlfg6331_64','mrg32k3a','philox4x32_10','threefry4x64_20',mldfg6331
File_name=['N=',num2str(n)];

A=zeros(1,1000);
parfor i=1:1000
    A(i)=log(i);
end
end