Qkdat = dlmread('../../unittest/Qk16/Qk_r1.txt',',');
%  Qkdat = dlmread('../../unittest/Qk16/Qk.txt',',');
Qkdat(:,end) = [];

assert(any(Qkdat(:) ~= 0));
nullbasis = null(Qkdat);

for cmp = 1:3
Qk11dat = Qkdat(:, cmp:3:end);
nullbasis = null(Qk11dat);

% Qk11dat = Qkdat;
% nullbasis = null(Qk11dat);
% nullbasis = nullbasis(cmp:3:end, :);


%% Find out zero part of the null space
q = sum((abs(nullbasis) < 100*eps), 2);
nullbasisS = sum(nullbasis, 2);
%nullbasisS(abs(nullbasisS) > 100*eps) = 1;
q = nullbasisS;
figure(cmp)
ax = gca;
q2d = reshape(q, sqrt(size(q,1)), sqrt(size(q,1)));
imagesc(gca, q2d)

title_string = {"Ker K11", "Ker K22", "Ker K12"}
title(title_string{cmp})

%%
% figure(1)
for i = 1:size(nullbasis, 2)
  ax = gca;
  q = nullbasis(:, i);
  q2d = reshape(q, sqrt(size(q,1)), sqrt(size(q,1)));
%   imagesc(gca, q2d)
%   drawnow
%   pause
end

end