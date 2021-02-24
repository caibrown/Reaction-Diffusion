%%%%%%%%%%%%%%%%%%%% Pattern Formation in Nature %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Two component Reaction Diffusion system %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Gray-Scott Model %%%%%%%%%%%%%%%%%%%%%%%%%%

% "Feed" and "kill" rates of chemicals (morphogens)
%f=.015; k=.052;   % Example 1 "Wandering bubbles"
%f=.055; k=.062;   % Example 2 "Pufferfish skin"
%f=.025; k=.06;    % Example 3 "Polka dots"
f=.022; k=.05;    % Example 4 "Leopard spots"
%f=.027; k=.055;    % Example 5 psychedeliccc

da = 1; db = .5;   % Diffusion rates
width = 384;       % Size of grid

% Steps per second / Simulation seconds
dt = .25; tf = 10000;

% A and B are the concentrations of our two chemicals..
[t, A, B] = initial_conditions(width);

% Make figure
axes('Position',[0 0 1 1])
axis off
colormap hsv
hi = image(B); hi.CDataMapping = 'scaled'; % Scaled-color image

% Text object to show the current time
ht = text(3,width-8,'Time = 0','fontsize',24); ht.Color = [.95 .2 .8];

while t < tf
    anew = A + (da*laplac(A) - A.*B.^2 + f*(1-A))*dt;
    bnew = B + (db*laplac(B) + A.*B.^2 - (k+f)*B)*dt;
    A = anew;
    B = bnew;
    hi.CData = B;
    t = t+dt;
    ht.String = ['Time = ' num2str(t)];
    drawnow limitrate % makes MATLAB draw faster
end


% Local functions
function out = laplac(in)
  % Use circshift to impose "toroidal wrap", i.e. cells on the
  % top/bottom and right/left are treated as if they're adjacent
  out = -in ...
      + .20*(circshift(in,[ 1, 0]) + circshift(in,[-1, 0])  ...
      +      circshift(in,[ 0, 1]) + circshift(in,[ 0,-1])) ...
      + .05*(circshift(in,[ 1, 1]) + circshift(in,[-1, 1])  ...
      +      circshift(in,[-1,-1]) + circshift(in,[ 1,-1]));
end

function [t, A, B] = initial_conditions(n)
  t = 0;
  % Initialize A to one
  A = ones(n);
  % Initialize B to zero with central clump of ones
  B = zeros(n);
  B(187:192,187:192) = 1;
  B(181:130,191:130) = 1;
end
