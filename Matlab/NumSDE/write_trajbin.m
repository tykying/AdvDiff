function [] = write_trajbin(filepath, dat_x, dat_y, dat_ts, dat_t1) 

% Ensure the time text is printed with 3 digits
N = floor(log10(dat_t1));
assert(abs(N) < 100)
ttmp = sprintf('%23.16E', dat_t1);

assert(ttmp(end-3) == 'E');
dat_tstr = [ttmp(1:end-2), '0', ttmp(end-1:end)];

fileID = fopen([filepath,'.hdr'],'w');
fprintf(fileID, '%s\n', 'serial');
fprintf(fileID, '%s\n', 'traj');
fprintf(fileID, '%d %d\n',  size(dat_x, 1), size(dat_y, 2));
fprintf(fileID, '%d\n',  2);
fprintf(fileID, '%s\n', dat_tstr);
fclose(fileID);

fileID = fopen([filepath,'_x.dat'],'w');
fwrite(fileID,dat_x,'double');
fclose(fileID);

fileID = fopen([filepath,'_y.dat'],'w');
fwrite(fileID,dat_y,'double');
fclose(fileID);

fileID = fopen([filepath,'_t.dat'],'w');
fwrite(fileID,dat_ts,'double');
fclose(fileID);

end