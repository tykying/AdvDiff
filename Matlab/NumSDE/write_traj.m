%% Convert the trajectories and write binary file
function [] = write_traj(x, y, ts_list, t_1, filename)
dat_x = x';
dat_y = y';
dat_ts = zeros(size(dat_x));
for part = 1:size(dat_ts, 1)
  dat_ts(part, :) = ts_list';
end
dat_t1 = t_1;

write_trajbin(filename, dat_x, dat_y, dat_ts, dat_t1);

disp(['Finish writing to ', filename]);
end