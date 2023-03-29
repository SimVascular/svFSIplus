clear all; close all; clc; %#ok<CLALL>

N     = 64;
fldr  = sprintf('../N%03d', N);
ntime = 2;
nsd   = 2;
tmin  = 0.0;
tmax  = 100.0;

nFa   = 4;
flist = {'bottom', 'top', 'left', 'right'};

for fa=1:nFa
    fhdr   = flist{fa};
    fname  = sprintf('%s/csv/bc_%s_nodeid.csv', ...
        fldr, fhdr);
    nodeId = csvread(fname);
    nNo = size(nodeId,1);

    v = zeros(nNo,nsd);
    for j=1:nsd
        fname = sprintf('%s/csv/bc_%s_v%d.csv',...
            fldr, fhdr,j);
        if ~exist(fname,'file')
            continue;
        end
        C = csvread(fname)';
        v(:,j) = reshape(C, nNo, 1);
    end

    t = linspace(tmin, tmax, ntime);
    fname = sprintf('%s/%s_vbc.dat', fldr, fhdr);
    fid = fopen(fname,'w');
    fprintf(fid,'%d   %d   %d\n', nsd, ntime, nNo);
    for i=1:ntime
        fprintf(fid,'%.9f\n', t(i));
    end
    for a=1:nNo
        fprintf(fid,'%d\n', nodeId(a));
        for i=1:ntime
            for j=1:nsd
                fprintf(fid,'%.18e ', v(a,j));
            end
            fprintf(fid,'\n');
        end
    end
    fclose(fid);
end

