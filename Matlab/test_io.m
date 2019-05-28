%IO
rng(1)

dat_x = randn(3, 4)
dat_y = randn(3, 4)
dat_t = randn(3, 4)

dat_x(3, 2) = 1.0
dat_y(1, 2) = -2.0
dat_t(3, 4) = 3.0

dat_tt = -1.2e-1;

%%
filename = 'test';

N = floor(log10(dat_tt));
assert(abs(N) < 100)
ttmp = sprintf('%23.16E', dat_tt);

assert(ttmp(end-3) == 'E');
dat_tstr = [ttmp(1:end-2), '0', ttmp(end-1:end)];

fileID = fopen(['../unittest/',filename,'.hdr'],'w');
fprintf(fileID, '%s\n', 'serial');
fprintf(fileID, '%s\n', 'traj');
fprintf(fileID, '%d %d\n',  size(dat_x, 1), size(dat_y, 2));
fprintf(fileID, '%d\n',  2);
fprintf(fileID, '%s\n', dat_tstr);
fclose(fileID);

fileID = fopen(['../unittest/',filename,'_x.dat'],'w');
fwrite(fileID,dat_x,'double');
fclose(fileID);

fileID = fopen(['../unittest/',filename,'_y.dat'],'w');
fwrite(fileID,dat_y,'double');
fclose(fileID);

fileID = fopen(['../unittest/',filename,'_t.dat'],'w');
fwrite(fileID,dat_t,'double');
fclose(fileID);