function main_fmm

clear
rng('default')
rng(0)

%%%%%%%%%%%%%%% Set the geometry.
flag_geo = 4;  % Set this parameter to 1, 2, or 3 to choose a geometry.
if (flag_geo == 1)
  nmax = 50;
  ntot = 3000;
  box_geom_root = [1,0.5,0.5];
  xx = [box_geom_root(2) + box_geom_root(1)*(rand(1,ntot)-0.5);...
        box_geom_root(3) + box_geom_root(1)*(rand(1,ntot)-0.5)];
elseif (flag_geo == 2)
  nmax = 50;
  ntot = 900;
  tt = 2*pi*rand(1,ntot);
  rr = 1 + 0.025*randn(1,ntot);
  xx = [rr.*cos(tt);...
        rr.*sin(tt)];
  x1min = min(xx(1,:));
  x1max = max(xx(1,:));
  x2min = min(xx(2,:));
  x2max = max(xx(2,:));
  len   = (1 + 1e-10)*max(x1max - x1min,x2max - x2min);
  box_geom_root = [len,0.5*(x1min+x1max),0.5*(x2min+x2max)];
elseif (flag_geo == 3)
  nside     = 720;
  nmax      = 25;
  ntot      = nside*nside;
  h         = 1/nside;
  [xx1,xx2] = meshgrid(linspace(0.5*h, 1 - 0.5*h, nside));
  xx        = [reshape(xx1,1,ntot);...
               reshape(xx2,1,ntot)];
  box_geom_root = [1,0.5,0.5];
elseif (flag_geo == 4)
  nside     = 720;
  nmax      = 200;
  ntot      = nside*nside;
  ntot4     = round(ntot/4);
  ntot2     = ntot - ntot4 * 2;
  h         = 1/nside;
  r_inner   = 0.5;
  r_outer   = 1;

  rr        = sqrt(rand(1, ntot2) * (r_outer^2-r_inner^2) + r_inner^2);
  theta     = rand(1, ntot2) * pi;
  xxl1      = rand(1, ntot4) * (r_outer - r_inner) - 1;
  xxr1      = - rand(1, ntot4) * (r_outer - r_inner) + 1;
  xxl2      = rand(1, ntot4) - 1;
  xxr2      = rand(1, ntot4) - 1;
  xx        = [rr .* cos(theta), xxl1, xxr1;...
               rr .* sin(theta), xxl2, xxr2];
  box_geom_root = [2,0,0];
end

%%%%%%%%%%%%%%%% Set various parameters
flag_precomp = 0;    % Determines whether to precompute the operators.
                     % (This applies _only_ to operators that depend on xx.)
flag_mode    = '11'; % This parameter specifies the type of the 
                     % charges and the sources.
                     %   flag_mode = '00'   monopoles given, potentials sought
                     %   flag_mode = '01'   monopoles given, fields sought
                     %   flag_mode = '10'   dipoles given,   potentials sought
                     %   flag_mode = '11'   dipoles given,   fields sought
p_fmm        = 30;   % The order of the multipole expansions used.
params       = ntot; % The total number of points.
fprintf(1,'======= FMM test code ==============\n')
fprintf(1,'ntot      = %d\n',ntot)
fprintf(1,'p_fmm     = %d\n',p_fmm)
fprintf(1,'nbox_max  = %d\n',nmax)
fprintf(1,'flag_mode = %s\n',flag_mode)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Perform the actual compression %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
[NODES,T_OPS,indvec] = ...
     LOCAL_FMM01_init(xx,nmax,p_fmm,flag_precomp,flag_mode,params);
t_init = toc;
nlevels = NODES{2,size(NODES,2)};                               
fprintf(1,'t_init    = %0.2e (sec)\n',t_init)
mem1 = whos('NODES');
mem2 = whos('T_OPS');
mem  = mem1.bytes + mem2.bytes;
fprintf(1,'memory    = %0.2e (MB)\n',mem/(2^20))
fprintf(1,'nlevels   = %d\n',nlevels)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Set up and solve a test problem %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (flag_mode(2) == '0') % Monopole charges.
  qq = randn(ntot,3);
elseif (flag_mode(2) == '1') % Dipole charges.
  qq = randn(ntot,3) + 1i*randn(ntot,3);
end
t_ = tic;
[uu, t_ofs, t_ofo, t_ifo, t_ifi, t_tfi, t_close] = LOCAL_FMM01_apply(xx,NODES,T_OPS,indvec,qq,flag_precomp,flag_mode,params);
t_apply = toc(t_);
fprintf(1,'t_apply   = %0.2e    (nrhs = %d)\n',t_apply,size(qq,2))
t_tot = t_init + t_apply;

% tic
% if (ntot <= 2000)
%   A        = LOCAL_A_diag(xx,flag_mode,params);
%   uu_exact = A*qq;
%   if (flag_mode(1) == '0')
%     uu_exact = real(uu_exact);
%   end
%   fprintf(1,'Error     = %0.7e (relative error, full comparison)\n',norm(uu_exact-uu)/norm(uu_exact))
% else
%   n_samp = 20;
%   [~,ind_samp] = sort(rand(1,ntot));
%   ind_samp = ind_samp(1:n_samp);
%   uu_exact = LOCAL_applyrows(xx,flag_mode,params,ind_samp,qq);
%   fprintf(1,'Error     = %0.7e  (relative error, sampled)\n',norm(uu(ind_samp,:) - uu_exact)/norm(uu_exact))
% end
% t_direct = toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Display data in various ways %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Graphically give the locations and numbers of various boxes.
% LOCAL_draw_boxes(NODES,xx,box_geom_root);
% saveas(gcf, sprintf('plot/grid-ntot-%d-nmax-%d-flag_geo-%d.epsc', ntot, nmax, flag_geo))
% return

% Print various lists.
%lists_to_be_printed = 3:4;
%LOCAL_print_lists(NODES,lists_to_be_printed)

save(sprintf('data/ntot-%d-nmax-%d-flag_geo-%d.mat', ntot, nmax, flag_geo), ...
    't_tot', 't_init', 't_ofs', 't_ofo', 't_ifo', 't_ifi', 't_tfi', 't_close')

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [NODES,T_OPS,indvec] = LOCAL_FMM01_init(xx_orig,nmax,p_fmm,flag_precomp,flag_mode,params)

%%% Construct a root box.
x1max = max(max(xx_orig(1,:)));
x1min = min(min(xx_orig(1,:)));
x2max = max(max(xx_orig(2,:)));
x2min = min(min(xx_orig(2,:)));
len   = (1 + 1e-10)*max(x1max - x1min, x2max - x2min);
box_geom_root    = zeros(1,3);
box_geom_root(1) = len;
box_geom_root(2) = 0.5*(x1min + x1max);
box_geom_root(3) = 0.5*(x2min + x2max);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1: Set up a tree structure.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[NODES,nlevels,indvec] = LOCAL_get_tree(xx_orig,box_geom_root,nmax);
nboxes                 = size(NODES,2);
xx                     = xx_orig(:,indvec);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 2: Create the interaction lists.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Construct a list "leaf_list" of all childless boxes.
leaf_list  = zeros(1,nboxes);
ileaf      = 0;
for ibox = 1:nboxes
  if (length(NODES{04,ibox}) == 0)
    ileaf = ileaf+1;
    leaf_list(ileaf) = ibox;
  end
end
nleaf = ileaf;
leaf_list = leaf_list(1:nleaf);

% Initialize all lists.
for ibox = 1:nboxes
  for ilist = 10:2:18
    NODES{ilist,ibox} = [];
  end
end

% Create lists 2 and 5
for ibox = 2:nboxes
  lside          = NODES{1,ibox}(1);
  x1c            = NODES{1,ibox}(2);
  x2c            = NODES{1,ibox}(3);
  ifath          = NODES{3,ibox};
  fath_coll_list = [NODES{18,ifath},ifath];
  for ifathcoll = fath_coll_list
    for jbox = NODES{04,ifathcoll}
      if ~(ibox == jbox)
        y1c  = NODES{1,jbox}(2);
        y2c  = NODES{1,jbox}(3);
        dist = sqrt((y1c-x1c)*(y1c-x1c) + (y2c-x2c)*(y2c-x2c));
        if (dist < (1.9*lside))
          % jbox is a colleague of ibox
          NODES{18,ibox} = [NODES{18,ibox},jbox];
        else
          NODES{12,ibox} = [NODES{12,ibox},jbox];
        end
      end
    end
  end
end

% Create lists 1 and 3
for ibox = 2:nboxes
  ibox_geom = NODES{1,ibox};
  list1   = [];
  list3   = [];
  if (length(NODES{04,ibox}) == 0)
    for jbox = NODES{18,ibox}
      [list1_fromjbox,list3_fromjbox] = LOCAL_get_list1_and_list3(NODES,ibox_geom,jbox,xx);
      list1 = [list1,list1_fromjbox];
      list3 = [list3,list3_fromjbox];
    end
  end
  NODES{10,ibox} = list1;
  NODES{14,ibox} = list3;
end

% At this point, NODES{10,ibox} lists only touching boxes at the same
% or finer levels as ibox. We next add the missing boxes.
for ibox = 2:nboxes
  for jbox = NODES{10,ibox}
    lside_ibox = NODES{1,ibox}(1);
    lside_jbox = NODES{1,jbox}(1);
    if ((lside_ibox / lside_jbox) > 1.0001)
      NODES{10,jbox} = [NODES{10,jbox},ibox];
    end
  end
end

% Create list 4.
% It is the dual of list 3.
for ibox = 2:nboxes
  for jbox = NODES{14,ibox}
    NODES{16,jbox} = [NODES{16,jbox},ibox];
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 3: Construct the translation operators
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T_OPS = LOCAL_construct_translationops(box_geom_root,nlevels,p_fmm);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 4: Construct the "FMM data" (if requested)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (flag_precomp == 1)
    
  % Construct all matrices required for leaf interactions:
  %    - The self-interaction matrix.
  %    - T_ofs
  %    - T_tfi
  for ibox = 1:nboxes
    if (length(NODES{04,ibox}) == 0)
      ind            = NODES{06,ibox} - 1 + (1:NODES{07,ibox});
      box_geom       = NODES{1,ibox};
      xx_loc         = xx(:,ind);
      NODES{40,ibox} = LOCAL_A_diag(xx(:,ind),flag_mode,params);
      NODES{46,ibox} = LOCAL_T_ofs(box_geom,p_fmm,xx_loc,flag_mode);
      NODES{47,ibox} = LOCAL_T_tfi(xx_loc,box_geom,p_fmm,flag_mode);
    end
  end

  % Create the matrices for list 1.
  for ibox = 1:nboxes
    if (length(NODES{10,ibox}) > 0)
      J_trg = NODES{06,ibox} - 1 + (1:NODES{07,ibox});
      for j = 1:length(NODES{10,ibox})
        jbox              = NODES{10,ibox}(j);
        J_src             = NODES{06,jbox} - 1 + (1:NODES{07,jbox});
        NODES{41,ibox}{j} = LOCAL_A_offd(xx(:,J_trg),xx(:,J_src),flag_mode,params);
      end
    end
  end
  
  % Create the matrices for list 3.
  for ibox = 1:nboxes
    if (length(NODES{14,ibox}) > 0)
      J_trg  = NODES{06,ibox} - 1 + (1:NODES{07,ibox});
      xx_trg = xx(:,J_trg);
      for j = 1:length(NODES{14,ibox})
        jbox              = NODES{14,ibox}(j);
        box_geom_src      = NODES{1,jbox};
        NODES{43,ibox}{j} = LOCAL_T_tfo(xx_trg,box_geom_src,p_fmm,flag_mode);
      end
    end
  end

  % Create the matrices for list 4.
  for ibox = 1:nboxes
    if (length(NODES{16,ibox}) > 0)
      box_geom_trg = NODES{1,ibox};
      for j = 1:length(NODES{16,ibox})
        jbox   = NODES{16,ibox}(j);
        J_src  = NODES{06,jbox} - 1 + (1:NODES{07,jbox});
        xx_src = xx(:,J_src);
        NODES{44,ibox}{j} = LOCAL_T_ifs(box_geom_trg,p_fmm,xx_src,flag_mode);
      end
    end
  end

end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Auxiliary functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [list1,list3] = LOCAL_get_list1_and_list3(NODES,ibox_geom,jbox,xx)

lside_i = ibox_geom(1);
lside_j = NODES{1,jbox}(1);

lsep    = max(abs(NODES{1,jbox}(2)-ibox_geom(2)),abs(NODES{1,jbox}(3)-ibox_geom(3)));
          
is_touching  = (lsep < 0.50001*(lside_i + lside_j));
has_children = (length(NODES{04,jbox})>0);

if ~is_touching
  list1 = [];
  list3 = [jbox];
else
  if ~has_children
    list1 = [jbox];
    list3 = [];
  else
    list1 = [];
    list3 = [];
    for kbox = NODES{04,jbox}
      [list1_fromk,list3_fromk] = LOCAL_get_list1_and_list3(NODES,ibox_geom,kbox,xx);
      list1 = [list1,list1_fromk];
      list3 = [list3,list3_fromk];
    end
  end
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [NODES,nlevels,indvec] = LOCAL_get_tree(xx,box_geom,nmax)

ntot = size(xx,2);
len  = box_geom(1);
x1c  = box_geom(2);
x2c  = box_geom(3);


% BOXES is a temporary work array for setting up the tree-structure
BOXES       = zeros(18,100);
lenBOXES    = 100;
BOXES( 1,1) = NaN;
BOXES( 2,1) = 0;
BOXES( 3,1) = NaN;
BOXES( 6,1) = 1;
BOXES( 7,1) = ntot;
BOXES(10,1) = NaN;
BOXES(11,1) = x1c - 0.5*len;
BOXES(12,1) = x1c + 0.5*len;
BOXES(13,1) = x2c - 0.5*len;
BOXES(14,1) = x2c + 0.5*len;
BOXES(15,1) = -1;
BOXES(16,1) = -1;
BOXES(17,1) = -1;
BOXES(18,1) = -1;

indvec = 1:ntot;

% Create the tree structure by quartering any
% box that holds more than nmax nodes.
%
% We create the boxes one level at a time.
% The following loop is over LEVELS.
ibox_last  = 0;
ibox_new   = 1;
ilevel     = 0;
while (ibox_new > ibox_last)

  ibox_first = ibox_last+1;
  ibox_last  = ibox_new;
  ilevel     = ilevel+1;
  
  % Loop over all boxes on the level that was last created.
  % All newly created boxes are temporarily stored in the array TMPBOXES.
  for ibox = ibox_first:ibox_last
    
    % If ibox holds more than nmax nodes, it will be partitioned.
    if (BOXES(07,ibox) > nmax) 

      qfirst = BOXES(06,ibox);
      x1min  = BOXES(11,ibox);
      x1max  = BOXES(12,ibox);
      x2min  = BOXES(13,ibox);
      x2max  = BOXES(14,ibox);
      x1half = (x1min + x1max)/2;
      x2half = (x2min + x2max)/2;
      J      = qfirst - 1 + (1:BOXES(07,ibox));
      indloc = indvec(J);
      J_sw   = find( (xx(1,indloc) <= x1half) .* (xx(2,indloc) <= x2half) );
      J_nw   = find( (xx(1,indloc) <= x1half) .* (xx(2,indloc) >  x2half) );
      J_se   = find( (xx(1,indloc) >  x1half) .* (xx(2,indloc) <= x2half) );
      J_ne   = find( (xx(1,indloc) >  x1half) .* (xx(2,indloc) >  x2half) );
           
      indvec(J) = [indloc(J_sw),indloc(J_nw),indloc(J_se),indloc(J_ne)];
      
      npart_added = 0;
      % If there is not enough space to save the 4 children in
      % the array BOXES, then double the size of BOXES.
      if ((ibox_new + 4) > lenBOXES)
        BOXES  = [BOXES,zeros(size(BOXES,1),4+size(BOXES,2))];
        lenBOXES = size(BOXES,2);
      end
      if ~isempty(J_sw)
        ibox_new = ibox_new + 1;
        n_sw = length(J_sw);
        BOXES(:,ibox_new) = [NaN,ilevel,ibox,NaN,NaN,qfirst+npart_added,...
                             n_sw,NaN,NaN,NaN,x1min,x1half,x2min,x2half,...
                             -1,-1,-1,-1]';
        BOXES(15,ibox) = ibox_new;
        npart_added = npart_added + n_sw;
      end
      if ~isempty(J_nw)
        ibox_new = ibox_new + 1;
        n_nw = length(J_nw);
        BOXES(:,ibox_new) = [NaN,ilevel,ibox,NaN,NaN,qfirst+npart_added,...
                             n_nw,NaN,NaN,NaN,x1min,x1half,x2half,x2max,...
                             -1,-1,-1,-1]';
        BOXES(16,ibox) = ibox_new;
        npart_added = npart_added + n_nw;
      end
      if ~isempty(J_se)
        ibox_new = ibox_new + 1;
        n_se = length(J_se);
        BOXES(:,ibox_new) = [NaN,ilevel,ibox,NaN,NaN,qfirst+npart_added,...
                             n_se,NaN,NaN,NaN,x1half,x1max,x2min,x2half,...
                             -1,-1,-1,-1]';

        BOXES(17,ibox) = ibox_new;
        npart_added = npart_added + n_se;
      end
      if ~isempty(J_ne)
        ibox_new = ibox_new + 1;
        n_ne = length(J_ne);
        BOXES(:,ibox_new) = [NaN,ilevel,ibox,NaN,NaN,qfirst+npart_added,...
                             n_ne,NaN,NaN,NaN,x1half,x1max,x2half,x2max,...
                             -1,-1,-1,-1]';
        BOXES(18,ibox) = ibox_new;
        npart_added = npart_added + n_ne;
      end
    
    end
    
  end
  
end

nlevels = ilevel - 1;

% Let nboxes denote the number of boxes created,
% create the object NODES which will hold all
% relevant information, and transfer the 
% information in BOXES to NODES.
% We also delete the object BOXES.
nboxes  = ibox_new;
NODES = cell(45,nboxes);
for ibox = 1:nboxes
  NODES{01,ibox} = [     BOXES(12,ibox)-BOXES(11,ibox),...
                    0.5*(BOXES(11,ibox)+BOXES(12,ibox)),...
                    0.5*(BOXES(13,ibox)+BOXES(14,ibox))];
  NODES{02,ibox} = BOXES(02,ibox);
  NODES{03,ibox} = BOXES(03,ibox);
  NODES{04,ibox} = [];
  for j = 15:18
    if (BOXES(j,ibox) > 0)
      NODES{04,ibox} = [NODES{04,ibox},BOXES(j,ibox)];
    end
  end
  NODES{06,ibox} = BOXES(06,ibox);
  NODES{07,ibox} = BOXES(07,ibox);
end


return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function T_OPS = LOCAL_construct_translationops(box_geom_root,nlevels,p_fmm)

T_OPS = cell(3,nlevels);

for ilevel = 1:nlevels
    
  a_loc = box_geom_root(1)/(2^ilevel);
  
  % Construct the out-from-out and in-from-in operators
  T_OPS{1,ilevel}    = cell(1,4);
  T_OPS{3,ilevel}    = cell(1,4);
  box_geom_dad       = [    a_loc,           0,           0];
  box_geom_son1      = [0.5*a_loc, -0.25*a_loc, -0.25*a_loc];
  box_geom_son2      = [0.5*a_loc, -0.25*a_loc,  0.25*a_loc];
  box_geom_son3      = [0.5*a_loc,  0.25*a_loc, -0.25*a_loc];
  box_geom_son4      = [0.5*a_loc,  0.25*a_loc,  0.25*a_loc];
  T_OPS{1,ilevel}{1} = LOCAL_T_ofo(box_geom_dad,box_geom_son1,p_fmm);
  T_OPS{1,ilevel}{2} = LOCAL_T_ofo(box_geom_dad,box_geom_son2,p_fmm);
  T_OPS{1,ilevel}{3} = LOCAL_T_ofo(box_geom_dad,box_geom_son3,p_fmm);
  T_OPS{1,ilevel}{4} = LOCAL_T_ofo(box_geom_dad,box_geom_son4,p_fmm);
  T_OPS{3,ilevel}{1} = LOCAL_T_ifi(box_geom_son1,box_geom_dad,p_fmm);
  T_OPS{3,ilevel}{2} = LOCAL_T_ifi(box_geom_son2,box_geom_dad,p_fmm);
  T_OPS{3,ilevel}{3} = LOCAL_T_ifi(box_geom_son3,box_geom_dad,p_fmm);
  T_OPS{3,ilevel}{4} = LOCAL_T_ifi(box_geom_son4,box_geom_dad,p_fmm);

  % Construct the translation operators (the ones associated with "list 2")
  T_OPS{2,ilevel} = cell(1,49);
  box_geom_trg    = a_loc*[1,0,0];
  for i1 = (-3):3
    for i2 = (-3):3
      if (max(abs(i1),abs(i2)) >= 2)
        box_geom_src          = a_loc*[1,i1,i2];
        i_nr                  = (i1+3)*7 + i2 + 4; 
        T_OPS{2,ilevel}{i_nr} = LOCAL_T_ifo(box_geom_trg,box_geom_src,p_fmm);
      end
    end
  end
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function T = LOCAL_T_ofs(box_geom_src,p_fmm,xxsrc,flag_mode)

if (nargin == 3)
  flag_mode = '00';
end

nsrc = size(xxsrc,2);
dd   = ones(p_fmm,1)*((xxsrc(1,:)-box_geom_src(2)) + 1i*(xxsrc(2,:) - box_geom_src(3)));
PP   = ((1:p_fmm)')*ones(1,nsrc);

if (flag_mode(2) == '0')      % The sources are monopoles.

  T = [ ones(1,nsrc); -(1./PP).*(dd.^PP)];
  
elseif (flag_mode(2) == '1')  % The sources are dipoles.
   
  T = [ zeros(1,nsrc); dd.^(PP-1)];
  
else
    
  fprintf(1,'Error in LOCAL_T_ofs: This mode not implemented.\n')
  keyboard
  
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function T = LOCAL_T_tfi(xxtrg,box_geom,p_fmm,flag_mode)

if (nargin == 3)
  flag_mode = '00';
end

ntrg = size(xxtrg,2);

if (flag_mode(1) == '0')

  dd = ((xxtrg(1,:)-box_geom(2)) + 1i*(xxtrg(2,:) - box_geom(3))).' * ones(1,p_fmm);
  PP = ones(ntrg,1)*(1:p_fmm);
  T  = [ones(ntrg,1),dd.^PP];

elseif (flag_mode(1) == '1')
    
  dd = ((xxtrg(1,:)-box_geom(2)) + 1i*(xxtrg(2,:) - box_geom(3))).' * ones(1,p_fmm-1);
  PP = ones(ntrg,1)*(1:(p_fmm-1));
  T  = [zeros(ntrg,1),ones(ntrg,1),(dd.^PP).*(ones(ntrg,1)*(2:p_fmm))];
    
else
    
  fprintf(1,'Error in LOCAL_T_ofs: This mode not implemented.\n')
  keyboard
  
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function T = LOCAL_T_tfo(xxtrg,box_geom,p_fmm,flag_mode)

if (nargin == 3)
  flag_mode = '00';
end

ntrg = size(xxtrg,2);
PP   = ones(ntrg,1)*(1:p_fmm);
dd   = ((xxtrg(1,:)-box_geom(2)) + 1i*(xxtrg(2,:)-box_geom(3))).';

if (flag_mode(1) == '0')

  T  = [log(dd), 1./((dd*ones(1,p_fmm)).^PP)];
  
elseif (flag_mode(1) == '1')
    
  T  = [1./dd, -PP./((dd*ones(1,p_fmm)).^(PP+1))];
    
else
    
  fprintf(1,'Error in LOCAL_T_ofs: This mode not implemented.\n')
  keyboard
  
end

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function T = LOCAL_T_ifs(box_geom_trg,p_fmm,xxsrc,flag_mode)

if (nargin == 3)
  flag_mode = '00';
end

nsrc = size(xxsrc,2);

if (flag_mode(2) == '0')

  T = zeros(p_fmm+1,nsrc);
  dd            = (xxsrc(1,:) - box_geom_trg(2)) + 1i*(xxsrc(2,:)-box_geom_trg(3));
  PP               = ((1:p_fmm)')*ones(1,nsrc);
  T = [log(-dd); -1./(PP.*((ones(p_fmm,1)*dd).^PP))];
  
elseif (flag_mode(2) == '1')
    
  dd = ones(p_fmm+1,1) * ((xxsrc(1,:) - box_geom_trg(2)) + 1i*(xxsrc(2,:)-box_geom_trg(3)));
  PP = ((0:p_fmm)')*ones(1,nsrc);
  T  = -1./(dd.^(PP+1));
  
else
    
  fprintf(1,'Error in LOCAL_T_ofs: This mode not implemented.\n')
  keyboard

end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function T = LOCAL_T_ofo(box_geom_dad,box_geom_son,p_fmm)

cdad1  = box_geom_dad(2);
cdad2  = box_geom_dad(3);

cson1  = box_geom_son(2);
cson2  = box_geom_son(3);

dd     = (cson1-cdad1) + 1i*(cson2-cdad2);

[jj,ii] = meshgrid(0:(p_fmm-1));
X       = tril(dd.^max(ii-jj,0));
Y       = LOCAL_binomial(ii,min(ii,jj));
i_vec   = (1:p_fmm)';
t_vec   = -(1./i_vec).*(dd.^i_vec);
T       = [1,     zeros(1,p_fmm);...
           t_vec, X.*Y];

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function T = LOCAL_T_ifi(box_geom_son,box_geom_dad,p_fmm)

cdad1  = box_geom_dad(2);
cdad2  = box_geom_dad(3);

cson1  = box_geom_son(2);
cson2  = box_geom_son(3);

dd      = (cson1-cdad1) + 1i*(cson2-cdad2);
[jj,ii] = meshgrid(0:p_fmm);
X       = triu(dd.^max(jj-ii,0));
Y       = LOCAL_binomial(jj,min(ii,jj));
T       = X.*Y;

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function T = LOCAL_T_ifo(box_geom_trg,box_geom_src,p_fmm)

ct1  = box_geom_trg(2);
ct2  = box_geom_trg(3);

cs1  = box_geom_src(2);
cs2  = box_geom_src(3);

T    = ones(p_fmm+1,p_fmm+1);
dd   = (cs1-ct1) + 1i*(cs2-ct2);

% Construct the 1-1-element.
T(1,1) = log(-dd);

% Construct the top row;
jj = 1:p_fmm;
T(1,2:(p_fmm+1)) = ((-1).^jj)./((dd*ones(1,p_fmm)).^jj);

% Construct the left-most diagonal. 
ii = (1:p_fmm)';
T(2:(p_fmm+1),1) = -1./(ii.*((dd*ones(p_fmm,1)).^ii));

% Construct the "bulk" of the matrix
[JJ,II] = meshgrid(1:p_fmm);
T(2:(p_fmm+1),2:(p_fmm+1)) = ((-1).^JJ) .* LOCAL_binomial(II+JJ-1,JJ-1) .* (1./((dd*ones(p_fmm,p_fmm)).^(II+JJ)));

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function c = LOCAL_binomial(n, k)

nk = [n(:); k(:)];
c = exp(gammaln(n+1) - gammaln(k+1) - gammaln(n-k+1));  % binomial coefficient
i = n == floor(n + .5) & k == floor(k + .5);
c(i) = floor(c(i) + .5);                                % number of combinations

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function A = LOCAL_A_offd(xxtrg,xxsrc,flag_mode,params)

ntrg   = size(xxtrg,2);
nsrc   = size(xxsrc,2);
zztrg  = xxtrg(1,:) + 1i*xxtrg(2,:);
zzsrc  = xxsrc(1,:) + 1i*xxsrc(2,:);
dd     = zztrg.' * ones(1,nsrc) - ones(ntrg,1)*zzsrc;

if (flag_mode == '00')

  ddsq = dd.*conj(dd);
  A    = 0.5*log(ddsq);
  
elseif (flag_mode == '01')
    
  A    = 1./dd;

elseif (flag_mode == '10')
    
  A    = 1./dd;

elseif (flag_mode == '11')
    
  A    = -1./(dd.*dd);

else
    
  fprintf(1,'ERROR in LOCAL_A_offd: This choice of "flag_mode" is not implemented.\n')
  keyboard
  
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function A = LOCAL_A_diag(xx,flag_mode,params)

nloc = size(xx,2);
xxs  = ones(nloc,1)*(xx(1,:) + 1i * xx(2,:));
xxt  = (xx(1,:) + 1i * xx(2,:)).' * ones(1,nloc);
dd   = xxt - xxs + eye(nloc);

if (flag_mode == '00')

  ddsq = dd.*conj(dd);
  A    = 0.5*log(ddsq);
  for j = 1:nloc
    A(j,j) = 0;
  end

elseif (flag_mode == '01')
    
  A = 1./dd;
  for j = 1:nloc
    A(j,j) = 0;
  end
  
elseif (flag_mode == '10')
    
  A = 1./dd;
  for j = 1:nloc
    A(j,j) = 0;
  end
  
elseif (flag_mode == '11')
    
  A = -1./(dd.*dd);
  for j = 1:nloc
    A(j,j) = 0;
  end
  
else
    
  fprintf(1,'Error in LOCAL_A_diag - option not implemented.\n')
  keyboard
  
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [uu_orig, t_ofs, t_ofo, t_ifo, t_ifi, t_tfi, t_close] = LOCAL_FMM01_apply(xx_orig,NODES,T_OPS,indvec,qq_orig,flag_precomp,flag_mode,params)

p_fmm = size(T_OPS{2,1}{1},1)-1;

xx = xx_orig(:,indvec);
qq = qq_orig(indvec,:);

nboxes = size(NODES,2);
ntot   = size(qq,1);
nrhs   = size(qq,2);
uu     = zeros(ntot,nrhs);

% Initialize the temporary fields.
FIELDS = cell(4,nboxes);
for ibox = 1:nboxes
  FIELDS{2,ibox} = zeros(p_fmm+1,nrhs);
  FIELDS{4,ibox} = zeros(p_fmm+1,nrhs);
end

% Apply all T_ofs
tic
for ibox = 1:nboxes
  if (length(NODES{04,ibox}) == 0)
    ind      = NODES{06,ibox} - 1 + (1:NODES{07,ibox});
    if (flag_precomp == 0)
      box_geom = NODES{1,ibox};
      xx_src   = xx(:,ind);
      T_ofs    = LOCAL_T_ofs(box_geom,p_fmm,xx_src,flag_mode);
      FIELDS{2,ibox} = T_ofs*qq(ind,:);
    else
      FIELDS{2,ibox} = NODES{46,ibox}*qq(ind,:);
    end
  end
end
t_ofs = toc;

% Upwards pass.
toc
for ibox = nboxes:(-1):2
  if (length(NODES{04,ibox}) > 0)
    for ison = NODES{04,ibox}
      ilevel         = NODES{02,ibox};
      i_location     = LOCAL_which_son(NODES,ibox,ison); % i_location codifies the relative locations of ibox and ison.
      FIELDS{2,ibox} = FIELDS{2,ibox} + T_OPS{1,ilevel}{i_location}*FIELDS{2,ison};
    end
  end
end
t_ofo = toc;

% Process contributions from list 1.
% (Leaves that are touching - direct evaluation.)
tic
for ibox = 1:nboxes
  if (length(NODES{10,ibox}) > 0)
    J_trg = NODES{06,ibox} - 1 + (1:NODES{07,ibox});
    for j = 1:length(NODES{10,ibox})
      jbox  = NODES{10,ibox}(j);
      J_src = NODES{06,jbox} - 1 + (1:NODES{07,jbox});
      if (flag_precomp == 0)
        B           = LOCAL_A_offd(xx(:,J_trg),xx(:,J_src),flag_mode,params);
        uu(J_trg,:) = uu(J_trg,:) + B*qq(J_src,:);
      else
        uu(J_trg,:) = uu(J_trg,:) + NODES{41,ibox}{j}*qq(J_src,:);
      end
    end
  end
end
t_ofs = toc;

% Process contributions from list 2 --- the "interaction lists"
tic
for ibox = 1:nboxes
  if (length(NODES{12,ibox}) > 0)
    ilevel       = NODES{2,ibox};
    for j = 1:length(NODES{12,ibox})
      jbox   = NODES{12,ibox}(j);
      i1     = round((NODES{1,jbox}(2) - NODES{1,ibox}(2))/NODES{1,ibox}(1));
      i2     = round((NODES{1,jbox}(3) - NODES{1,ibox}(3))/NODES{1,ibox}(1));
      i_nr   = (i1+3)*7 + i2 + 4; % i_nr codifies the relative locations of ibox and jbox.
      FIELDS{4,ibox} = FIELDS{4,ibox} + T_OPS{2,ilevel}{i_nr}*FIELDS{2,jbox};
    end
  end
end
t_ifo = toc;

% Process contributions from list 3.
% For the leaf ibox, list 3 indicates the boxes that broadcast their
% outgoing expansions directly to the potential on ibox.
tic
for ibox = 1:nboxes
  if (length(NODES{14,ibox}) > 0)
    J_trg  = NODES{06,ibox} - 1 + (1:NODES{07,ibox});
    for j = 1:length(NODES{14,ibox})
      jbox         = NODES{14,ibox}(j);
      if (flag_precomp == 0)
        xx_trg       = xx(:,J_trg);
        box_geom_src = NODES{1,jbox};
        T_tfo        = LOCAL_T_tfo(xx_trg,box_geom_src,p_fmm,flag_mode);
        uu(J_trg,:)  = uu(J_trg,:) + T_tfo*FIELDS{2,jbox}; 
      else
        uu(J_trg,:)  = uu(J_trg,:) + NODES{43,ibox}{j}*FIELDS{2,jbox}; 
      end
    end
  end
end

% Process contributions from list 4.
% For box ibox, list 4 indicates the boxes jbox whose charges directly
% contribute to the incoming expansion on ibox.
for ibox = 1:nboxes
  if (length(NODES{16,ibox}) > 0)
    for j = 1:length(NODES{16,ibox})
      jbox  = NODES{16,ibox}(j);
      J_src = NODES{06,jbox} - 1 + (1:NODES{07,jbox});
      if (flag_precomp == 0)
        box_geom_trg   = NODES{1,ibox};
        xx_src         = xx(:,J_src);
        T_ifs          = LOCAL_T_ifs(box_geom_trg,p_fmm,xx_src,flag_mode);
        FIELDS{4,ibox} = FIELDS{4,ibox} + T_ifs*qq(J_src,:);
      else
        FIELDS{4,ibox} = FIELDS{4,ibox} + NODES{44,ibox}{j}*qq(J_src,:);
      end
    end
  end
end
t_close = toc;

% Downwards pass - apply all in-to-in operators.
tic
for ibox = 2:nboxes
  if (length(NODES{04,ibox}) > 0)
    for ison = NODES{04,ibox}
      ilevel         = NODES{02,ibox};
      i_location     = LOCAL_which_son(NODES,ibox,ison); % i_location codifies the relative locations of ibox and ison.
      FIELDS{4,ison} = FIELDS{4,ison} + T_OPS{3,ilevel}{i_location}*FIELDS{4,ibox};
    end
  end
end
t_ifi = toc;

% Expand all incoming expansions on the leaves to potentials.
tic
for ibox = 2:nboxes
  if (length(NODES{04,ibox}) == 0)
    ind       = NODES{06,ibox} - 1 + (1:NODES{07,ibox});
    if (flag_precomp == 0)
      xx_trg    = xx(:,ind);
      box_geom  = NODES{01,ibox};
      T_tfi     = LOCAL_T_tfi(xx_trg,box_geom,p_fmm,flag_mode);
      uu(ind,:) = uu(ind,:) + T_tfi*FIELDS{4,ibox};
    else
      uu(ind,:) = uu(ind,:) + NODES{47,ibox}*FIELDS{4,ibox};
    end
  end
end
t_tfi = toc;

% Add all self-interactions.
tic
for ibox = 1:nboxes
  if (length(NODES{04,ibox}) == 0)
    ind       = NODES{06,ibox} - 1 + (1:NODES{07,ibox});
    if (flag_precomp == 0)
      D         = LOCAL_A_diag(xx(:,ind),flag_mode,params);
      uu(ind,:) = uu(ind,:) + D*qq(ind,:);
    else
      uu(ind,:) = uu(ind,:) + NODES{40,ibox}*qq(ind,:);
    end
  end
end
t_close = t_close + toc;

if (flag_mode(1) == '0')
  uu = real(uu);
end

uu_orig = zeros(size(uu));
uu_orig(indvec,:) = uu;

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The four sons of a box are labelled using the integer i_location
% as follows
%
%    son 2     son 4
%    son 1     son 3
%
% This function simply retrieves this number.
function i_location = LOCAL_which_son(NODES,ibox,ison)

if (sum(NODES{4,ibox} == ison) <= 0)
  disp('ERROR in "LOCAL_which_son": The given boxes do not form a son/father pair.')
  keyboard
end

xc1_son =  NODES{1,ison}(2);
xc2_son  = NODES{1,ison}(3);
xc1_fath = NODES{1,ibox}(2);
xc2_fath = NODES{1,ibox}(3);

if (xc1_son < xc1_fath)
    if (xc2_son < xc2_fath)
        i_location = 1;
    else
        i_location = 2;
    end
else
    if (xc2_son < xc2_fath)
        i_location = 3;
    else
        i_location = 4;
    end
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function uu = LOCAL_applyrows(xx,flag_mode,params,ind_samp,qq)

if (length(xx)*length(ind_samp) > 4000000)
  disp('in LOCAL_applyrows')
  disp('The problem size seems a bit large ... breaking.')
  keyboard
end

ntot   = length(xx);
indtmp = 1:ntot;
indtmp(ind_samp) = 2*ntot;
indtmp = sort(indtmp);
ind_offd = indtmp(1:(ntot-length(ind_samp)));

A_offd = LOCAL_A_offd(xx(:,ind_samp),xx(:,ind_offd),flag_mode,params);
A_diag = LOCAL_A_diag(xx(:,ind_samp),flag_mode,params);

uu = A_offd*qq(ind_offd,:) + A_diag*qq(ind_samp,:);

if (flag_mode(1) == '0')
  uu = real(uu);
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function LOCAL_draw_boxes(NODES,xx,box_geom)

len = box_geom(1);

box = [box_geom(2) - 0.5*len,...
       box_geom(2) + 0.5*len,...
       box_geom(3) - 0.5*len,...
       box_geom(3) + 0.5*len];

nlevels = 1;
for ibox = 2:size(NODES,2)
  nlevels = max(NODES{2,ibox},nlevels);
end
for ilevel = 1:4
  subplot(2,2,ilevel)
  hold off
  plot(box([1,2,2,1,1]),box([3,3,4,4,3]),'k')
  hold on
  axis equal
  axis off
  plot(xx(1,:),xx(2,:),'g.');
  title(sprintf('Level %d',ilevel))
end

flag_warning = 0;

for ibox = 2:size(NODES,2)
  ilevel = NODES{2,ibox};
  if (ilevel <= 4)
    subplot(2,2,ilevel)
    x1min = NODES{1,ibox}(2) - 0.5*NODES{1,ibox}(1);
    x1max = NODES{1,ibox}(2) + 0.5*NODES{1,ibox}(1);
    x2min = NODES{1,ibox}(3) - 0.5*NODES{1,ibox}(1);
    x2max = NODES{1,ibox}(3) + 0.5*NODES{1,ibox}(1);
    plot([x1min,x1max,x1max,x1min,x1min],...
         [x2min,x2min,x2max,x2max,x2min],'r')
    hold on
    text(0.5*(x1min + x1max),...
         0.5*(x2min + x2max),...
         sprintf('%d',ibox))
  else
    flag_warning = 1;
  end
end

for ilevel = 1:min(nlevels,4)
  subplot(2,2,ilevel)
  hold off
end

if (flag_warning == 1)
  fprintf(1,'WARNING in LOCAL_draw_boxes: More than 4 levels!!')
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function LOCAL_print_lists(NODES,listvec)

nboxes = size(NODES,2);

for ilist = listvec

  fprintf(1,'=== Printing all list %d\n',ilist)
  for ibox = 2:nboxes
    if (length(NODES{10+2*(ilist-1),ibox})>0)
      list_string = [];
      for j = 1:length(NODES{10+2*(ilist-1),ibox})
        list_string = [list_string,sprintf(' %4d',NODES{10+2*(ilist-1),ibox}(j))];
      end
      disp([sprintf(' Box %4d has list-%1d:',ibox,ilist),list_string])
    end
  end
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
