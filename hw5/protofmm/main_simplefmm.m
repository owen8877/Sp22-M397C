function main_simplefmm

clear
rng('default')
rng(0)

%%%%%%%%%%%%%%% Set the geometry.
flag_geo = 3;  % Set this parameter to 1, 2, or 3 to choose a geometry.
if (flag_geo == 1)
  nmax = 50;
  ntot = 1000;
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
  nside     = 40;
  nmax      = 25;
  ntot      = nside*nside;
  h         = 1/nside;
  [xx1,xx2] = meshgrid(linspace(0.5*h, 1 - 0.5*h, nside));
  xx        = [reshape(xx1,1,ntot);...
               reshape(xx2,1,ntot)];
  box_geom_root = [1,0.5,0.5];
end

%%%%%%%%%%%%%%%% Set various parameters
p_fmm  = 20;  % The order of the multipole expansions.
fprintf(1,'======= FMM test code ==============\n')
fprintf(1,'ntot      = %d\n',ntot)
fprintf(1,'p_fmm     = %d\n',p_fmm)
fprintf(1,'nbox_max  = %d\n',nmax)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Apply the FMM
%%% 
%%% You can call the function to let it specify "box_geom_root" 
%%% by itself, or you can force it to use a given box.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
qq = randn(ntot,3);
uu = LOCAL_SIMPLEFMM(xx,qq,nmax,p_fmm,box_geom_root);
%uu = LOCAL_SIMPLEFMM(xx,qq,nmax,p_fmm);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Verify the solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A        = LOCAL_A_diag(xx);
uu_exact = A * qq;
fprintf(1,'||uu - uu_exact|| / ||uu_exact|| = %18.10e\n',...
        norm(uu - uu_exact)/norm(uu_exact))

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a very basic 2D FMM for the log-kernel.
%
% NOTE: It is VERY slow and should be used only for learning purposes.
%
% It evaluates a matrix-vector product 
%
%   uu_orig = A * qq_orig 
%
% where A is a matrix of size ntot x ntot whose ij entry is
%
%   A(i,j) = log(norm(xx_orig(:,i) - xx_orig(:,j))))
%
% The object xx_orig specifies the particle locations, its size is 2 x ntot.
%
% (The moniker "_orig" is used since the vectors will be reordered inside
%  the function to conform with the ordering of the tree.
%
% INPUTS:  xx_orig        - particle locations, a matrix of size 2 x ntot
%          qq_orig        - charges, a vector of size ntot x nrhs
%          nmax           - maximum number of particles in a leaf
%          p_fmm          - order of multipole expansions
%          box_geom_root  - the geometry of the box holding the particles xx
%                           (this input is optional)
%
% OUTPUT:  uu_orig        - the potentials
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function uu_orig = LOCAL_SIMPLEFMM(xx_orig,qq_orig,nmax,p_fmm,box_geom_root)

%%% If not "box_geom_root" object is specified, then create it:
if (nargin == 4)
  x1min = min(xx_orig(1,:));
  x1max = max(xx_orig(1,:));
  x2min = min(xx_orig(2,:));
  x2max = max(xx_orig(2,:));
  len   = (1 + 1e-10)*max(x1max-x1min,x2max-x2min);
  box_geom_root = [len,...
                   0.5 * (x1min + x1max),...
                   0.5 * (x2min + x2max)];
end

%%% Extract some parameters
ntot = size(xx_orig,2);
nrhs = size(qq_orig,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1: Set up a tree structure.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[NODES,indvec] = LOCAL_get_tree(xx_orig,box_geom_root,nmax);
nboxes         = size(NODES,2);
xx             = xx_orig(:,indvec);
figure(1)
LOCAL_draw_boxes(NODES,xx,box_geom_root);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 2: Construct the various lists needed.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NODES = LOCAL_get_lists(NODES,xx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize fields.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
qq = qq_orig(indvec,:); % we permute the indices so that the potentials associated with a box are contiguous. 
uu = zeros(ntot,nrhs); 
OUTFIELD = cell(1,nboxes);
INFIELD  = cell(1,nboxes);
for ibox = 1:nboxes
  OUTFIELD{ibox} = zeros(p_fmm+1,nrhs);
  INFIELD{ ibox} = zeros(p_fmm+1,nrhs);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Apply all T_ofs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ibox = 1:nboxes
  if (length(NODES{04,ibox}) == 0) % ibox has no children
    ind      = NODES{06,ibox} - 1 + (1:NODES{07,ibox});
    box_geom = NODES{1,ibox};
    xx_src   = xx(:,ind);
    T_ofs    = LOCAL_T_ofs(box_geom,p_fmm,xx_src);
    OUTFIELD{ibox} = T_ofs*qq(ind,:);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Apply all T_ofo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ibox = nboxes:(-1):2
  if ((NODES{3,ibox} > 1) && (length(NODES{12,ibox}) > 0))
    box_geom_dad = NODES{01,ibox}; 
    for ison = NODES{04,ibox}
      box_geom_son   = NODES{01,ison};
      T_ofo          = LOCAL_T_ofo(box_geom_dad,box_geom_son,p_fmm);
      OUTFIELD{ibox} = OUTFIELD{ibox} + T_ofo*OUTFIELD{ison};
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Apply all T_ifo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ibox = 1:nboxes
  if ((NODES{3,ibox} > 1) && (length(NODES{12,ibox}) > 0))
    box_geom_trg = NODES{1,ibox};
    for j = 1:length(NODES{12,ibox})
      jbox          = NODES{12,ibox}(j);
      box_geom_src  = NODES{1,jbox};
      T_ifo         = LOCAL_T_ifo(box_geom_trg,box_geom_src,p_fmm);
      INFIELD{ibox} = INFIELD{ibox} + T_ifo*OUTFIELD{jbox};
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Apply all T_tfo    ("list 3")
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ibox = 1:nboxes
  if (length(NODES{14,ibox}) > 0)
    J_trg  = NODES{06,ibox} - 1 + (1:NODES{07,ibox});
    for j = 1:length(NODES{14,ibox})
      jbox         = NODES{14,ibox}(j);
      xx_trg       = xx(:,J_trg);
      box_geom_src = NODES{1,jbox};
      T_tfo        = LOCAL_T_tfo(xx_trg,box_geom_src,p_fmm);
      uu(J_trg,:)  = uu(J_trg,:) + T_tfo*OUTFIELD{jbox}; 
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Apply all T_ifs    ("list 4")
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ibox = 1:nboxes
  if (length(NODES{16,ibox}) > 0)
    for j = 1:length(NODES{16,ibox})
      jbox          = NODES{16,ibox}(j);
      J_src         = NODES{06,jbox} - 1 + (1:NODES{07,jbox});
      box_geom_trg  = NODES{1,ibox};
      xx_src        = xx(:,J_src);
      T_ifs         = LOCAL_T_ifs(box_geom_trg,p_fmm,xx_src);
      INFIELD{ibox} = INFIELD{ibox} + T_ifs*qq(J_src,:);
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Apply all T_ifi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ibox = 2:nboxes
  if ((NODES{3,ibox} > 1) && (length(NODES{04,ibox}) > 0))
    box_geom_dad = NODES{01,ibox};
    for ison = NODES{04,ibox}
      box_geom_son  = NODES{01,ison};
      T_ifi         = LOCAL_T_ifi(box_geom_son,box_geom_dad,p_fmm);
      INFIELD{ison} = INFIELD{ison} + T_ifi*INFIELD{ibox};
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Apply all T_tfi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ibox = 2:nboxes
  if (length(NODES{04,ibox}) == 0) % ibox is a leaf
    ind       = NODES{06,ibox} - 1 + (1:NODES{07,ibox});
    xx_trg    = xx(:,ind);
    box_geom  = NODES{01,ibox};
    T_tfi     = LOCAL_T_tfi(xx_trg,box_geom,p_fmm);
    uu(ind,:) = uu(ind,:) + T_tfi*INFIELD{ibox};
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluate interactions with direct neighbors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ibox = 1:nboxes
  if (length(NODES{10,ibox}) > 0)
    J_trg = NODES{06,ibox} - 1 + (1:NODES{07,ibox});
    for j = 1:length(NODES{10,ibox})
      jbox        = NODES{10,ibox}(j);
      J_src       = NODES{06,jbox} - 1 + (1:NODES{07,jbox});
      B           = LOCAL_A_offd(xx(:,J_trg),xx(:,J_src));
      uu(J_trg,:) = uu(J_trg,:) + B*qq(J_src,:);
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluate leaf self-interactions 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ibox = 1:nboxes
  if (length(NODES{04,ibox}) == 0)
    ind       = NODES{06,ibox} - 1 + (1:NODES{07,ibox});
    B         = LOCAL_A_diag(xx(:,ind));
    uu(ind,:) = uu(ind,:) + B*qq(ind,:);
  end
end

% Undo the permutation of the indices.
uu_orig = zeros(size(uu));
uu_orig(indvec,:) = real(uu);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function NODES = LOCAL_get_lists(NODES,xx)

nboxes = size(NODES,2);

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

function [NODES,indvec] = LOCAL_get_tree(xx,box_geom,nmax)

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

% Let nboxes denote the number of boxes created,
% create the object NODES which will hold all
% relevant information, and transfer the 
% information in BOXES to NODES.
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

function T = LOCAL_T_ofs(box_geom_src,p_fmm,xxsrc)

nsrc = size(xxsrc,2);
dd   = ones(p_fmm,1)*((xxsrc(1,:)-box_geom_src(2)) + 1i*(xxsrc(2,:) - box_geom_src(3)));
PP   = ((1:p_fmm)')*ones(1,nsrc);
T    = [ ones(1,nsrc); -(1./PP).*(dd.^PP)];
  
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function T = LOCAL_T_tfi(xxtrg,box_geom,p_fmm)

ntrg = size(xxtrg,2);
dd = ((xxtrg(1,:)-box_geom(2)) + 1i*(xxtrg(2,:) - box_geom(3))).' * ones(1,p_fmm);
PP = ones(ntrg,1)*(1:p_fmm);
T  = [ones(ntrg,1),dd.^PP];

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function T = LOCAL_T_tfo(xxtrg,box_geom,p_fmm)

ntrg = size(xxtrg,2);
PP   = ones(ntrg,1)*(1:p_fmm);
dd   = ((xxtrg(1,:)-box_geom(2)) + 1i*(xxtrg(2,:)-box_geom(3))).';
T    = [log(dd), 1./((dd*ones(1,p_fmm)).^PP)];
  
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function T = LOCAL_T_ifs(box_geom_trg,p_fmm,xxsrc)

nsrc = size(xxsrc,2);
dd   = (xxsrc(1,:) - box_geom_trg(2)) + 1i*(xxsrc(2,:)-box_geom_trg(3));
PP   = ((1:p_fmm)')*ones(1,nsrc);
T    = [log(-dd); -1./(PP.*((ones(p_fmm,1)*dd).^PP))];

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function A = LOCAL_A_offd(xxtrg,xxsrc)

ntrg   = size(xxtrg,2);
nsrc   = size(xxsrc,2);
zztrg  = xxtrg(1,:) + 1i*xxtrg(2,:);
zzsrc  = xxsrc(1,:) + 1i*xxsrc(2,:);
dd     = zztrg.' * ones(1,nsrc) - ones(ntrg,1)*zzsrc;

ddsq = dd.*conj(dd);
A    = 0.5*log(ddsq);
  
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function A = LOCAL_A_diag(xx)

nloc = size(xx,2);
xxs  = ones(nloc,1)*(xx(1,:) + 1i * xx(2,:));
xxt  = (xx(1,:) + 1i * xx(2,:)).' * ones(1,nloc);
dd   = xxt - xxs + eye(nloc);

ddsq = dd.*conj(dd);
A    = 0.5*log(ddsq);
for j = 1:nloc
  A(j,j) = 0;
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
  h = plot(xx(1,:),xx(2,:),'g.');
%  set(h,'MarkerSize',20)
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
  fprintf(1,'WARNING in FMM01_draw_boxes: More than 4 levels!!\n')
end

return