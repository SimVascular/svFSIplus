clear all; close all; clc; %#ok<CLALL>

N     = 64;
fhdr  = sprintf('../N%03d', N);
ntime = 2;
nsd   = 2;
tmin  = 0.0;
tmax  = 100.0;
nNo   = (N+1)^nsd;

u = zeros(nNo,nsd);
fname = sprintf('%s/csv/bforce.csv',fhdr);
V = csvread(fname);
u(:,1) = reshape(V, nNo, 1);

t = linspace(tmin, tmax, ntime);
fname = [fhdr '/bforce.dat'];
fid = fopen(fname,'w');
fprintf(fid,'%d   %d   %d\n', nsd, ntime, nNo);
for i=1:ntime
    fprintf(fid,'%.9f\n', t(i));
end
for a=1:nNo
    fprintf(fid,'%d\n',a);
    for i=1:ntime
        for j=1:nsd
            fprintf(fid,'%.18e ', u(a,j));
        end
        fprintf(fid,'\n');
    end
end
fclose(fid);


