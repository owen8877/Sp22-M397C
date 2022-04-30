function main_T_ops_are_fun

rng('default')
rng(0)

nsrc  = 30;   % Number of sources
ntrg  = 50;   % Number of targets
p_fmm = 20;   % Order of multipole expansion

% Set the geometry of the source and target boxes.
box_geom_src = [1,0.5,0.5];
box_geom_trg = [1,2.5,0.5];

% Set the source and target locations by randomly dumping points in the boxes.
xxsrc = [box_geom_src(2) + box_geom_src(1)*(rand(1,nsrc)-0.5);...
         box_geom_src(3) + box_geom_src(1)*(rand(1,nsrc)-0.5)];
xxtrg = [box_geom_trg(2) + box_geom_trg(1)*(rand(1,ntrg)-0.5);...
         box_geom_trg(3) + box_geom_trg(1)*(rand(1,ntrg)-0.5)];

% Compute the various operators.
A     = LOCAL_A_offd(xxtrg,xxsrc);
T_ofs = LOCAL_T_ofs(box_geom_src,p_fmm,xxsrc);
T_ifs = LOCAL_T_ifs(box_geom_trg,p_fmm,xxsrc);
T_ifo = LOCAL_T_ifo(box_geom_trg,box_geom_src,p_fmm);
T_tfi = LOCAL_T_tfi(xxtrg,box_geom_trg,p_fmm);
T_tfo = LOCAL_T_tfo(xxtrg,box_geom_src,p_fmm);

% Evaluate and print the errors.
fprintf(1,'p_fmm     = %3d                      ||ERR||           ||ERR||/||A||     max(max(abs(ERR)))\n',p_fmm)
ERR = A - real(T_tfo * T_ofs);
fprintf(1,'ERR = A - T_tfo *         T_ofs :    %12.5e      %12.5e      %12.5e\n',norm(ERR),norm(ERR)/norm(A),max(max(abs(ERR))))
ERR = A - real(T_tfi * T_ifo * T_ofs);
fprintf(1,'ERR = A - T_tfo *         T_ofs :    %12.5e      %12.5e      %12.5e\n',norm(ERR),norm(ERR)/norm(A),max(max(abs(ERR))))
ERR = A - real(T_tfi * T_ifs);
fprintf(1,'ERR = A - T_tfo *         T_ofs :    %12.5e      %12.5e      %12.5e\n',norm(ERR),norm(ERR)/norm(A),max(max(abs(ERR))))

% Set up a geometry for testing T_ofo and T_ifi.
box_geom_child = [1,0.8,0.6];
xx_child       = [box_geom_child(2) + box_geom_child(1)*(rand(1,200)-0.5);...
                box_geom_child(3) + box_geom_child(1)*(rand(1,200)-0.5)];
box_geom_parent = [2,  1,  1];
xx_outer_n   = [box_geom_parent(2) + 3.0*box_geom_parent(1)*(rand(1,100)-0.5);...
                box_geom_parent(3) + 1.5*box_geom_parent(1)*ones(1,100)];
xx_outer_s   = [box_geom_parent(2) + 3.0*box_geom_parent(1)*(rand(1,100)-0.5);...
                box_geom_parent(3) - 1.5*box_geom_parent(1)*ones(1,100)];
xx_outer_w   = [box_geom_parent(2) - 1.5*box_geom_parent(1)*ones(1,100);...
                box_geom_parent(3) + 3.0*box_geom_parent(1)*(rand(1,100)-0.5)];
xx_outer_e   = [box_geom_parent(2) + 1.5*box_geom_parent(1)*ones(1,100);...
                box_geom_parent(3) + 3.0*box_geom_parent(1)*(rand(1,100)-0.5)];
xx_outer     = [xx_outer_n,xx_outer_w,xx_outer_s,xx_outer_e];

fprintf(1,'Testing out-to-out and in-to-in.\n')
A     = LOCAL_A_offd(xx_outer,xx_child);
T_ofs = LOCAL_T_ofs(box_geom_child,p_fmm,xx_child);
T_ofo = LOCAL_T_ofo(box_geom_parent,box_geom_child,p_fmm);
T_tfo = LOCAL_T_tfo(xx_outer,box_geom_parent,p_fmm);
ERR = A - real(T_tfo * T_ofo * T_ofs);
fprintf(1,'ERR = A - T_tfo * T_ofo * T_ofs :    %12.5e      %12.5e      %12.5e\n',norm(ERR),norm(ERR)/norm(A),max(max(abs(ERR))))
A     = LOCAL_A_offd(xx_child,xx_outer);
T_ifs = LOCAL_T_ifs(box_geom_parent,p_fmm,xx_outer);
T_ifi = LOCAL_T_ifi(box_geom_child,box_geom_parent,p_fmm);
T_tfi = LOCAL_T_tfi(xx_child,box_geom_child,p_fmm);
ERR = A - real(T_tfi * T_ifi * T_ifs);
fprintf(1,'ERR = A - T_tfi * T_ifi * T_ifs :    %12.5e      %12.5e      %12.5e\n',norm(ERR),norm(ERR)/norm(A),max(max(abs(ERR))))

figure(1)
plot(xxsrc(1,:),xxsrc(2,:),'rx',...
     xxtrg(1,:),xxtrg(2,:),'bx',...
     box_geom_src(2) + box_geom_src(1)*[-0.5,0.5,0.5,-0.5,-0.5],...
     box_geom_src(3) + box_geom_src(1)*[-0.5,-0.5,0.5,0.5,-0.5],'r',...
     box_geom_trg(2) + box_geom_trg(1)*[-0.5,0.5,0.5,-0.5,-0.5],...
     box_geom_trg(3) + box_geom_trg(1)*[-0.5,-0.5,0.5,0.5,-0.5],'b')
legend('Sources','Targets')
axis equal
figure(2)
plot(box_geom_child(2) + box_geom_child(1)*[-0.5,0.5,0.5,-0.5,-0.5],...
     box_geom_child(3) + box_geom_child(1)*[-0.5,-0.5,0.5,0.5,-0.5],'-r',...
     xx_child(1,:),xx_child(2,:),'r.',...
     box_geom_parent(2) + box_geom_parent(1)*[-0.5,0.5,0.5,-0.5,-0.5],...
     box_geom_parent(3) + box_geom_parent(1)*[-0.5,-0.5,0.5,0.5,-0.5],'-b',...
     xx_outer(1,:),xx_outer(2,:),'k.')
legend('child','xxchild','parent.','xxouter')
axis equal

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

function T = LOCAL_T_ofo(box_geom_parent,box_geom_child,p_fmm)

cparent1  = box_geom_parent(2);
cparent2  = box_geom_parent(3);

cchild1  = box_geom_child(2);
cchild2  = box_geom_child(3);

dd     = (cchild1-cparent1) + 1i*(cchild2-cparent2);

[jj,ii] = meshgrid(0:(p_fmm-1));
X       = tril(dd.^max(ii-jj,0));
Y       = LOCAL_binomial(ii,min(ii,jj));
i_vec   = (1:p_fmm)';
t_vec   = -(1./i_vec).*(dd.^i_vec);
T       = [1,     zeros(1,p_fmm);...
           t_vec, X.*Y];

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function T = LOCAL_T_ifi(box_geom_child,box_geom_parent,p_fmm)

cparent1  = box_geom_parent(2);
cparent2  = box_geom_parent(3);

cchild1  = box_geom_child(2);
cchild2  = box_geom_child(3);

dd      = (cchild1-cparent1) + 1i*(cchild2-cparent2);
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

