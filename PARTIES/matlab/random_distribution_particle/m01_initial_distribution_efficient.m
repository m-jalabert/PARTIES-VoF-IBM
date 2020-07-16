%creates initial random particle distribution
clear all;
close all;

%% simulation setup for N Ardekani 2018
phi=0.01;%% set phi
L_x      = 1.0;
L_y      = 1.0;
L_z      = 1.0;
R50      = 1/50;
Vol_mean= 4.0 * pi * R50 * R50 * R50 / 3.0;
N_p      = int32(L_x*L_y*L_z*phi/Vol_mean);
variance = 0.;
stdeviat = sqrt(variance);
mu       = log(R50^2/sqrt(variance+R50^2));
sigma    = sqrt(log(variance/R50^2 + 1));
h        = R50/10;
%% walls 1 || periodic 0
south = 1; % has to be 1 by default
north = 1; % has to be 1 by default
east  = 1;
west  = 1;
front = 1;
back  = 1;


for p =1:N_p
    %% normal distribution
    % rp(p) = normrnd(R50,stdeviat);
    %% log-normal distribution
    %   rp(p) = lognrnd(mu,sigma);
    %% Constant
    rp(p) = R50;
    V_p(p) = 4.0 * pi * rp(p) * rp(p) * rp(p) / 3.0;
end
[h,stats] = cdfplot(rp);

fprintf('min : %g max   : %g max / min: %g \n', stats.min,  stats.max, max(rp) / min(rp) )
fprintf('Mean: %g Median: %g Std      : %g \n', stats.mean, stats.median, stats.std)
fprintf('Volume fraction: %g               \n', sum(V_p) / (L_x * L_y * L_z))

%% Decomposition
PX   = 5;
PY   = 5;
PZ   = 5;
% PTOT = PX * PY * PZ;
for i = 1:PX
    xs(i) = (i-1) * L_x / PX;
    xe(i) =  i    * L_x / PX;
end
for j = 1:PY
    ys(j) = (j-1) * L_y / PY;
    ye(j) =  j    * L_y / PY;
end
for k = 1:PZ
    zs(k) = (k-1) * L_z / PZ;
    ze(k) =  k    * L_z / PZ;
end

%% compute random coordinates
pd = makedist('Uniform');
p    = 0;
ph   = 0;
n    = 0;
fail = 0;
np(1:PX,1:PY,1:PZ) = 0;
while (N_p > p)
    n    = n + 1;
    p    = p + 1;
    
    xp(p) = L_x * random(pd);
    yp(p) = L_y * random(pd);
    zp(p) = L_z * random(pd);
    %     xp(p) = 20.0;
    %     yp(p) = 19.9+(p-1)*2.0*rp(p);
    %     zp(p) = 20.;
    
    ii     = floor(xp(p) * PX / L_x)+1;
    jj     = floor(yp(p) * PY / L_y)+1;
    kk     = floor(zp(p) * PZ / L_z)+1;
    
    %% particle - wall
    fail = 0;
    %SOUTH
    if (south == 1 && yp(p) < rp(p))
        fail = 1;
    end
    %NORTH
    if (north == 1 && yp(p)+rp(p) > L_y)
        fail = 1;
    end
    %EAST
    if (east == 1 && xp(p) < rp(p))
        fail = 1;
    end
    %WEST
    if (west == 1 && xp(p)+rp(p) > L_x)
        fail = 1;
    end
    %FRONT
    if (front == 1 && zp(p) < rp(p))
        fail = 1;
    end
    %BACK
    if (back == 1 && zp(p)+rp(p) > L_z)
        fail = 1;
    end
    if fail == 1
        p = p-1;
        continue;
    end
    
    np(ii,jj,kk)    = np(ii,jj,kk) + 1;
    pp              = np(ii,jj,kk);
    x(ii,jj,kk,pp) = xp(p);
    y(ii,jj,kk,pp) = yp(p);
    z(ii,jj,kk,pp) = zp(p);
    r(ii,jj,kk,pp) = rp(p);
    [fail,np(ii,jj,kk)] = m01_check_overlap(...
        xp(p),         yp(p),         zp(p),        rp(p),         ...
        x(ii,jj,kk,:), y(ii,jj,kk,:), z(ii,jj,kk,:),r(ii,jj,kk,:), ...
        pp);
    if fail == 1
        p = p-1;
        continue;
    end
    %     fprintf('ii : %d jj : %d kk : %d %d\n', ii, jj, kk, fail)
    %     fprintf('x  : %g y  : %g z  : %g\n', x(ii,jj,kk,pp), y(ii,jj,kk,pp), z(ii,jj,kk,pp))
    %     fprintf('master particle\n\n')
    
    %     neighborhood      front:         center:       back:
    %     process map:
    %                       6 |  7|  8     15|16 |17     24| 25| 26
    %                       __|___|___     __|___|___    __|___|___
    %                       3 | 4 |  5     12|13 |14     21| 22| 23
    %                       __|___|___     __|___|___    __|___|___
    %                       0 | 1 |  2     9 |10 |11     18| 19| 20
    %                         |   |          |   |         |   |
    
    ext_wes = xp(p) - rp(p); %east
    ext_eas = xp(p) + rp(p); %west
    
    ext_sou = yp(p) - rp(p); %south
    ext_nor = yp(p) + rp(p); %north
    
    ext_fro = zp(p) - rp(p); %front
    ext_bac = zp(p) + rp(p); %back
    
    
    ig = 0;
    jg = 0;
    kg = 0;
    
    if (ext_wes < xs(ii) )
        ig = ig-1;
    end
    if (ext_eas > xe(ii) )
        ig = ig+1;
    end
    
    if (ext_sou < ys(jj) )
        jg = jg-1;
    end
    if (ext_nor > ye(jj) )
        jg = jg+1;
    end
    
    if (ext_fro < zs(kk) )
        kg = kg-1;
    end
    if (ext_bac > ze(kk) )
        kg = kg+1;
    end
    
    id  = (ig+1)  + 3 * (jg+1) + 9 * (kg+1);
    
    if id ~= 13;
        if (id ==  4 ||...                                     %front
                id == 10 || id == 12 || id == 14 || id == 16 ||... %center
                id == 22)                                          %back
            
            [i,iii,j,jjj,k,kkk] = m01_get_indices(PX,ii+ig,PY,jj+jg,PZ,kk+kg);
            np(iii,jjj,kkk) = np(iii,jjj,kkk) + 1;
            pp              = np(iii,jjj,kkk)    ;
            x(iii,jjj,kkk,pp) = xp(p) - i * L_x;
            y(iii,jjj,kkk,pp) = yp(p) - j * L_y;
            z(iii,jjj,kkk,pp) = zp(p) - k * L_z;
            r(iii,jjj,kkk,pp) = rp(p);
            
            [fail,np(iii,jjj,kkk)] = m01_check_overlap(...
                x(iii,jjj,kkk,pp),y(iii,jjj,kkk,pp),z(iii,jjj,kkk,pp),r(iii,jjj,kkk,pp), ...
                x(iii,jjj,kkk,:) ,y(iii,jjj,kkk,:), z(iii,jjj,kkk,:), r(iii,jjj,kkk,:) , ...
                pp);
            if fail == 1
                p = p-1;
                continue;
            end
            
            
            %              fprintf('iii: %d jjj: %d kkk: %d id: %d\n', iii, jjj, kkk, id)
            %              fprintf('x  : %g y  : %g z  : %g\n', x(iii,jjj,kkk,pp), y(iii,jjj,kkk,pp), z(iii,jjj,kkk,pp))
            %              fprintf('periodic particle\n\n')
            
            
        elseif(id ==  9 || id == 11 || id == 15 || id == 17 ) %xy
            [i,iii,j,jjj,k,kkk] = m01_get_indices(PX,ii+ig,PY,jj+jg,PZ,kk);
            np(iii,jjj,kkk) = np(iii,jjj,kkk) + 1;
            pp              = np(iii,jjj,kkk)    ;
            x(iii,jjj,kkk,pp) = xp(p) - i * L_x;
            y(iii,jjj,kkk,pp) = yp(p) - j * L_y;
            z(iii,jjj,kkk,pp) = zp(p) - k * L_z;
            r(iii,jjj,kkk,pp) = rp(p);
            [fail,np(iii,jjj,kkk)] = m01_check_overlap(...
                x(iii,jjj,kkk,pp),y(iii,jjj,kkk,pp),z(iii,jjj,kkk,pp),r(iii,jjj,kkk,pp), ...
                x(iii,jjj,kkk,:) ,y(iii,jjj,kkk,:), z(iii,jjj,kkk,:), r(iii,jjj,kkk,:) , ...
                pp);
            if fail == 1
                p = p-1;
                continue;
            end
            
            [i,iii,j,jjj,k,kkk] = m01_get_indices(PX,ii+ig,PY,jj,PZ,kk);
            np(iii,jjj,kkk) = np(iii,jjj,kkk) + 1;
            pp              = np(iii,jjj,kkk)    ;
            x(iii,jjj,kkk,pp) = xp(p) - i * L_x;
            y(iii,jjj,kkk,pp) = yp(p) - j * L_y;
            z(iii,jjj,kkk,pp) = zp(p) - k * L_z;
            r(iii,jjj,kkk,pp) = rp(p);
            [fail,np(iii,jjj,kkk)] = m01_check_overlap(...
                x(iii,jjj,kkk,pp),y(iii,jjj,kkk,pp),z(iii,jjj,kkk,pp),r(iii,jjj,kkk,pp), ...
                x(iii,jjj,kkk,:) ,y(iii,jjj,kkk,:), z(iii,jjj,kkk,:), r(iii,jjj,kkk,:) , ...
                pp);
            if fail == 1
                p = p-1;
                continue;
            end
            %              fprintf('iii: %d jjj: %d kkk: %d id: %d\n', iii, jjj, kkk, id)
            %              fprintf('x  : %g y  : %g z  : %g\n', x(iii,jjj,kkk,pp), y(iii,jjj,kkk,pp), z(iii,jjj,kkk,pp))
            %              fprintf('periodic particle\n\n')
            
            [i,iii,j,jjj,k,kkk] = m01_get_indices(PX,ii,PY,jj+jg,PZ,kk);
            np(iii,jjj,kkk) = np(iii,jjj,kkk) + 1;
            pp              = np(iii,jjj,kkk)    ;
            x(iii,jjj,kkk,pp) = xp(p) - i * L_x;
            y(iii,jjj,kkk,pp) = yp(p) - j * L_y;
            z(iii,jjj,kkk,pp) = zp(p) - k * L_z;
            r(iii,jjj,kkk,pp) = rp(p);
            [fail,np(iii,jjj,kkk)] = m01_check_overlap(...
                x(iii,jjj,kkk,pp),y(iii,jjj,kkk,pp),z(iii,jjj,kkk,pp),r(iii,jjj,kkk,pp), ...
                x(iii,jjj,kkk,:) ,y(iii,jjj,kkk,:), z(iii,jjj,kkk,:), r(iii,jjj,kkk,:) , ...
                pp);
            if fail == 1
                p = p-1;
                continue;
            end
            
        elseif(id ==  3 || id ==  5 || id == 21 || id == 23 ) %xz
            [i,iii,j,jjj,k,kkk] = m01_get_indices(PX,ii+ig,PY,jj,PZ,kk);
            np(iii,jjj,kkk) = np(iii,jjj,kkk) + 1;
            pp              = np(iii,jjj,kkk)    ;
            x(iii,jjj,kkk,pp) = xp(p) - i * L_x;
            y(iii,jjj,kkk,pp) = yp(p) - j * L_y;
            z(iii,jjj,kkk,pp) = zp(p) - k * L_z;
            r(iii,jjj,kkk,pp) = rp(p);
            [fail,np(iii,jjj,kkk)] = m01_check_overlap(...
                x(iii,jjj,kkk,pp),y(iii,jjj,kkk,pp),z(iii,jjj,kkk,pp),r(iii,jjj,kkk,pp), ...
                x(iii,jjj,kkk,:) ,y(iii,jjj,kkk,:), z(iii,jjj,kkk,:), r(iii,jjj,kkk,:) , ...
                pp);
            if fail == 1
                p = p-1;
                continue;
            end
            
            [i,iii,j,jjj,k,kkk] = m01_get_indices(PX,ii,PY,jj,PZ,kk+kg);
            np(iii,jjj,kkk) = np(iii,jjj,kkk) + 1;
            pp              = np(iii,jjj,kkk)    ;
            x(iii,jjj,kkk,pp) = xp(p) - i * L_x;
            y(iii,jjj,kkk,pp) = yp(p) - j * L_y;
            z(iii,jjj,kkk,pp) = zp(p) - k * L_z;
            r(iii,jjj,kkk,pp) = rp(p);
            [fail,np(iii,jjj,kkk)] = m01_check_overlap(...
                x(iii,jjj,kkk,pp),y(iii,jjj,kkk,pp),z(iii,jjj,kkk,pp),r(iii,jjj,kkk,pp), ...
                x(iii,jjj,kkk,:) ,y(iii,jjj,kkk,:), z(iii,jjj,kkk,:), r(iii,jjj,kkk,:) , ...
                pp);
            if fail == 1
                p = p-1;
                continue;
            end
            
            [i,iii,j,jjj,k,kkk] = m01_get_indices(PX,ii+ig,PY,jj,PZ,kk+kg);
            np(iii,jjj,kkk) = np(iii,jjj,kkk) + 1;
            pp              = np(iii,jjj,kkk)    ;
            x(iii,jjj,kkk,pp) = xp(p) - i * L_x;
            y(iii,jjj,kkk,pp) = yp(p) - j * L_y;
            z(iii,jjj,kkk,pp) = zp(p) - k * L_z;
            r(iii,jjj,kkk,pp) = rp(p);
            [fail,np(iii,jjj,kkk)] = m01_check_overlap(...
                x(iii,jjj,kkk,pp),y(iii,jjj,kkk,pp),z(iii,jjj,kkk,pp),r(iii,jjj,kkk,pp), ...
                x(iii,jjj,kkk,:) ,y(iii,jjj,kkk,:), z(iii,jjj,kkk,:), r(iii,jjj,kkk,:) , ...
                pp);
            if fail == 1
                p = p-1;
                continue;
            end
        elseif(id == 1 || id == 7 || id == 19 || id == 25  ) %yz
            [i,iii,j,jjj,k,kkk] = m01_get_indices(PX,ii,PY,jj+jg,PZ,kk);
            np(iii,jjj,kkk) = np(iii,jjj,kkk) + 1;
            pp              = np(iii,jjj,kkk)    ;
            x(iii,jjj,kkk,pp) = xp(p) - i * L_x;
            y(iii,jjj,kkk,pp) = yp(p) - j * L_y;
            z(iii,jjj,kkk,pp) = zp(p) - k * L_z;
            r(iii,jjj,kkk,pp) = rp(p);
            [fail,np(iii,jjj,kkk)] = m01_check_overlap(...
                x(iii,jjj,kkk,pp),y(iii,jjj,kkk,pp),z(iii,jjj,kkk,pp),r(iii,jjj,kkk,pp), ...
                x(iii,jjj,kkk,:) ,y(iii,jjj,kkk,:), z(iii,jjj,kkk,:), r(iii,jjj,kkk,:) , ...
                pp);
            if fail == 1
                p = p-1;
                continue;
            end
            
            [i,iii,j,jjj,k,kkk] = m01_get_indices(PX,ii,PY,jj,PZ,kk+kg);
            np(iii,jjj,kkk) = np(iii,jjj,kkk) + 1;
            pp              = np(iii,jjj,kkk)    ;
            x(iii,jjj,kkk,pp) = xp(p) - i * L_x;
            y(iii,jjj,kkk,pp) = yp(p) - j * L_y;
            z(iii,jjj,kkk,pp) = zp(p) - k * L_z;
            r(iii,jjj,kkk,pp) = rp(p);
            [fail,np(iii,jjj,kkk)] = m01_check_overlap(...
                x(iii,jjj,kkk,pp),y(iii,jjj,kkk,pp),z(iii,jjj,kkk,pp),r(iii,jjj,kkk,pp), ...
                x(iii,jjj,kkk,:) ,y(iii,jjj,kkk,:), z(iii,jjj,kkk,:), r(iii,jjj,kkk,:) , ...
                pp);
            if fail == 1
                p = p-1;
                continue;
            end
            
            [i,iii,j,jjj,k,kkk] = m01_get_indices(PX,ii,PY,jj+jg,PZ,kk+kg);
            np(iii,jjj,kkk) = np(iii,jjj,kkk) + 1;
            pp              = np(iii,jjj,kkk)    ;
            x(iii,jjj,kkk,pp) = xp(p) - i * L_x;
            y(iii,jjj,kkk,pp) = yp(p) - j * L_y;
            z(iii,jjj,kkk,pp) = zp(p) - k * L_z;
            r(iii,jjj,kkk,pp) = rp(p);
            [fail,np(iii,jjj,kkk)] = m01_check_overlap(...
                x(iii,jjj,kkk,pp),y(iii,jjj,kkk,pp),z(iii,jjj,kkk,pp),r(iii,jjj,kkk,pp), ...
                x(iii,jjj,kkk,:) ,y(iii,jjj,kkk,:), z(iii,jjj,kkk,:), r(iii,jjj,kkk,:) , ...
                pp);
            if fail == 1
                p = p-1;
                continue;
            end
        elseif (id ==  0 || id ==  2 || id ==  6 || id ==  8 || ...
                id == 18 || id == 20 || id == 24 || id == 26)
            [i,iii,j,jjj,k,kkk] = m01_get_indices(PX,ii+ig,PY,jj,PZ,kk);
            np(iii,jjj,kkk) = np(iii,jjj,kkk) + 1;
            pp              = np(iii,jjj,kkk)    ;
            x(iii,jjj,kkk,pp) = xp(p) - i * L_x;
            y(iii,jjj,kkk,pp) = yp(p) - j * L_y;
            z(iii,jjj,kkk,pp) = zp(p) - k * L_z;
            r(iii,jjj,kkk,pp) = rp(p);
            [fail,np(iii,jjj,kkk)] = m01_check_overlap(...
                x(iii,jjj,kkk,pp),y(iii,jjj,kkk,pp),z(iii,jjj,kkk,pp),r(iii,jjj,kkk,pp), ...
                x(iii,jjj,kkk,:) ,y(iii,jjj,kkk,:), z(iii,jjj,kkk,:), r(iii,jjj,kkk,:) , ...
                pp);
            if fail == 1
                p = p-1;
                continue;
            end
            
            [i,iii,j,jjj,k,kkk] = m01_get_indices(PX,ii,PY,jj+jg,PZ,kk);
            np(iii,jjj,kkk) = np(iii,jjj,kkk) + 1;
            pp              = np(iii,jjj,kkk)    ;
            x(iii,jjj,kkk,pp) = xp(p) - i * L_x;
            y(iii,jjj,kkk,pp) = yp(p) - j * L_y;
            z(iii,jjj,kkk,pp) = zp(p) - k * L_z;
            r(iii,jjj,kkk,pp) = rp(p);
            [fail,np(iii,jjj,kkk)] = m01_check_overlap(...
                x(iii,jjj,kkk,pp),y(iii,jjj,kkk,pp),z(iii,jjj,kkk,pp),r(iii,jjj,kkk,pp), ...
                x(iii,jjj,kkk,:) ,y(iii,jjj,kkk,:), z(iii,jjj,kkk,:), r(iii,jjj,kkk,:) , ...
                pp);
            if fail == 1
                p = p-1;
                continue;
            end
            
            [i,iii,j,jjj,k,kkk] = m01_get_indices(PX,ii,PY,jj,PZ,kk+kg);
            np(iii,jjj,kkk) = np(iii,jjj,kkk) + 1;
            pp              = np(iii,jjj,kkk)    ;
            x(iii,jjj,kkk,pp) = xp(p) - i * L_x;
            y(iii,jjj,kkk,pp) = yp(p) - j * L_y;
            z(iii,jjj,kkk,pp) = zp(p) - k * L_z;
            r(iii,jjj,kkk,pp) = rp(p);
            [fail,np(iii,jjj,kkk)] = m01_check_overlap(...
                x(iii,jjj,kkk,pp),y(iii,jjj,kkk,pp),z(iii,jjj,kkk,pp),r(iii,jjj,kkk,pp), ...
                x(iii,jjj,kkk,:) ,y(iii,jjj,kkk,:), z(iii,jjj,kkk,:), r(iii,jjj,kkk,:) , ...
                pp);
            if fail == 1
                p = p-1;
                continue;
            end
            
            [i,iii,j,jjj,k,kkk] = m01_get_indices(PX,ii+ig,PY,jj+jg,PZ,kk);
            np(iii,jjj,kkk) = np(iii,jjj,kkk) + 1;
            pp              = np(iii,jjj,kkk)    ;
            x(iii,jjj,kkk,pp) = xp(p) - i * L_x;
            y(iii,jjj,kkk,pp) = yp(p) - j * L_y;
            z(iii,jjj,kkk,pp) = zp(p) - k * L_z;
            r(iii,jjj,kkk,pp) = rp(p);
            [fail,np(iii,jjj,kkk)] = m01_check_overlap(...
                x(iii,jjj,kkk,pp),y(iii,jjj,kkk,pp),z(iii,jjj,kkk,pp),r(iii,jjj,kkk,pp), ...
                x(iii,jjj,kkk,:) ,y(iii,jjj,kkk,:), z(iii,jjj,kkk,:), r(iii,jjj,kkk,:) , ...
                pp);
            if fail == 1
                p = p-1;
                continue;
            end
            
            [i,iii,j,jjj,k,kkk] = m01_get_indices(PX,ii+ig,PY,jj,PZ,kk+kg);
            np(iii,jjj,kkk) = np(iii,jjj,kkk) + 1;
            pp              = np(iii,jjj,kkk)    ;
            x(iii,jjj,kkk,pp) = xp(p) - i * L_x;
            y(iii,jjj,kkk,pp) = yp(p) - j * L_y;
            z(iii,jjj,kkk,pp) = zp(p) - k * L_z;
            r(iii,jjj,kkk,pp) = rp(p);
            [fail,np(iii,jjj,kkk)] = m01_check_overlap(...
                x(iii,jjj,kkk,pp),y(iii,jjj,kkk,pp),z(iii,jjj,kkk,pp),r(iii,jjj,kkk,pp), ...
                x(iii,jjj,kkk,:) ,y(iii,jjj,kkk,:), z(iii,jjj,kkk,:), r(iii,jjj,kkk,:) , ...
                pp);
            if fail == 1
                p = p-1;
                continue;
            end
            
            [i,iii,j,jjj,k,kkk] = m01_get_indices(PX,ii,PY,jj+jg,PZ,kk+kg);
            np(iii,jjj,kkk) = np(iii,jjj,kkk) + 1;
            pp              = np(iii,jjj,kkk)    ;
            x(iii,jjj,kkk,pp) = xp(p) - i * L_x;
            y(iii,jjj,kkk,pp) = yp(p) - j * L_y;
            z(iii,jjj,kkk,pp) = zp(p) - k * L_z;
            r(iii,jjj,kkk,pp) = rp(p);
            [fail,np(iii,jjj,kkk)] = m01_check_overlap(...
                x(iii,jjj,kkk,pp),y(iii,jjj,kkk,pp),z(iii,jjj,kkk,pp),r(iii,jjj,kkk,pp), ...
                x(iii,jjj,kkk,:) ,y(iii,jjj,kkk,:), z(iii,jjj,kkk,:), r(iii,jjj,kkk,:) , ...
                pp);
            if fail == 1
                p = p-1;
                continue;
            end
            
            [i,iii,j,jjj,k,kkk] = m01_get_indices(PX,ii+ig,PY,jj+jg,PZ,kk+kg);
            np(iii,jjj,kkk) = np(iii,jjj,kkk) + 1;
            pp              = np(iii,jjj,kkk)    ;
            x(iii,jjj,kkk,pp) = xp(p) - i * L_x;
            y(iii,jjj,kkk,pp) = yp(p) - j * L_y;
            z(iii,jjj,kkk,pp) = zp(p) - k * L_z;
            r(iii,jjj,kkk,pp) = rp(p);
            [fail,np(iii,jjj,kkk)] = m01_check_overlap(...
                x(iii,jjj,kkk,pp),y(iii,jjj,kkk,pp),z(iii,jjj,kkk,pp),r(iii,jjj,kkk,pp), ...
                x(iii,jjj,kkk,:) ,y(iii,jjj,kkk,:), z(iii,jjj,kkk,:), r(iii,jjj,kkk,:) , ...
                pp);
            if fail == 1
                p = p-1;
                continue;
            end
        end
    end
    if (mod(p,100) == 0)
        fprintf('%d  %d %g  %g  %g %g\n', p, n, xp(p), yp(p), zp(p), rp(p))
    end
    
    
end % while loop

%% write data
filename = ['./p_mobile.inp'];
fileID = fopen(filename,'w');
fprintf(fileID,'  %d     \r\n', N_p);

for p =1:N_p
    fprintf(fileID,'%12.8f %12.8f %12.8f %12.8f \r\n',xp(p), yp(p), zp(p), rp(p));
end

fclose(fileID);

%% plot data
showplot  = 1;
if (showplot == 1)
    fprintf('hello');
    %[x y z] = sphere;
    %a = [];
    
    ncol = 5;
    col = lines(ncol);
    C = zeros(N_p,3);
    for i = 1:N_p
        C(i,:) =  col(randi([1 ncol],1),:)';
    end
    
    

    
   
    
    figure
    
    scatter3sph(xp,yp,zp,'size',rp(1),'transp',1,'color',C); hold on;
    
    [Y,Z] = meshgrid(linspace(0,1,80));
    X = zeros(size(Z));
    vortX = cos((X+Y)*2) .* cos((X+Z)/5);
    h = surf(X,Y,Z,vortX);
    set(h,'LineStyle','none');
    [X,Z] = meshgrid(linspace(0,1,80));
    Y = zeros(size(Z));
    vortY = cos((X+Y)*2) .* cos((X+Z)*3);
    h = surf(X,Y,Z,vortY);
    set(h,'LineStyle','none');
    
    [X,Y] = meshgrid(linspace(0,1,80));
    Z = zeros(size(Y));
    vortZ = cos((X+Y)*3) .* cos((X+Z)*5);
    h = surf(X,Y,Z,vortZ);
    set(h,'LineStyle','none');
    
    set(gcf,'position',[1000 254 1224 1072]);
    set(gcf,'color','w');
    xlim([0 1])
    ylim([0 1])
    zlim([0 1])
    box on;
    
    
    %{
for p = 1: N_p
    new_sphere = [xp(p) yp(p) zp(p) rp(p)];
    a = vertcat(a, new_sphere);
    s(p) = surf(x*a(p,4)+a(p,1),y*a(p,4)+a(p,2),z*a(p,4)+a(p,3));
    if (p==1)
        hold on
    end
   
end

daspect([1 1 1])
%view(30,10)
xlabel('x')
ylabel('y')
zlabel('z')
    %}
end %showplot




