% mu-diff - Copyright (C) 2014-2017 X. Antoine and B. Thierry, University of Lorraine, CNRS, France
%
% See the LICENSE.txt file for license information. Please report all
% bugs and problems to either (replace [at] and [dot] by arobase and dot)
%   Xavier.Antoine [at] univ-lorraine [dot] fr
%   BThierry [at] math.cnrs [dot] fr
%
% ~~~
% Create and randomly distribue N_scat disks, with random centers and random radii inside the bos [xmin, xmax]x[ymin, ymax].
% [O, a] = CreateRandomDisks(xmin, xmax, ymin, ymax, N_scat, ...)
%
% Input Arguments:
% ----------------
% - xmin   [1x1]      : minimal abscissa limit of the box
% - xmax   [1x1]      : maximal abscissa limit of the box
% - ymin   [1x1]      : minimal ordinate limit of the box
% - ymax   [1x1]      : maximal ordinate limit of the box
% - N_scat [1x1]      : number of desired circular scatterers
%
% Output Arguments:
% -----------------
% - O      [2xN_scat] : coordinates of the center of the circular obstacles
% - a      [1xN_scat] : radii of the obstacles
%
% OPTIONS :
% ---------
% [O, a] = CreateRandomDisks(xmin, xmax, ymin, ymax, N_scat, amin, amax, 
%                            dmin, dmax, O_avoid, a_void, dmin_avoid, dmax_avoid)
%  value name (default value) [size] : description
% - amin (=1)         [1x1]    : minimal radius of the disk
% - amax (=1)         [1x1]    : maximal radius of the disk
% - dmin (=realmin)   [1x1]    : minimal distance between two disks
% - dmax (=realmax)   [1x1]    : maximal distance between two disks (Carreful with this option!)
% - O_avoid           [2xN]    : center of N disks (possibly of zero-radius) that must stay out
%                                of every scatterers
% - a_avoid           [1xN]    : radii of the disks of center O_avoid
% - dmin              [1x1]    : minimal distance that disks must kept from the other disks (O_avoid, a_void)
% - dmin              [1x1]    : maximal distance that disks must kept from the other disks (O_avoid, a_void)
%
% Remarks:
% --------
% - if (dmax <= 0) then dmax = realmax (no limit)
% - if (dmin <= 0) then dmin = realmin (no limit)
% - Same for dmin_avoid and dmax_avoid
% - Be carreful with dmax. It is not the distance between the closest disk,
%   but with ALL OTHER disks.
% - 1000 trials are done to place a disk. If after 1000 trials, a new disk
%   cannot be inserted, the function assumes that it is impossible and an
%   error message is returned.

function [O, a] = CreateRandomDisks(xmin, xmax, ymin, ymax, N_scat, varargin)

    %Maximum number of trial (security)
    trial_max = 1000;

    %Initialization
    O = zeros(2,N_scat);
    a = zeros(1,N_scat);
    amin = 1;
    amax = 1;
    dmin = realmin;
    dmax = realmax;
    O_avoid = [];
    a_avoid = [];
    dmin_avoid = realmin;
    dmax_avoid = realmax;

    %Reading arguments
    opt_argin = length(varargin);
    if(opt_argin >= 1) % Setting a minimal radius
        amin = varargin{1};
    end
    if(opt_argin >= 2) % Setting a maximal radius
        amax = varargin{2};
    end
    if(opt_argin >= 3) % Setting a minimal distance between two disks
        dmin = varargin{3};
    end
    if(opt_argin >= 4) % Setting a maximal distance between two disks
        dmax = varargin{4};
    end
    if(opt_argin >= 5) % Points to avoided (point source,...) ...
        O_avoid = varargin{5}; % coordinates of the points to be avoided
    end
    if(opt_argin >= 6)% ... or not points but disks ! (for maybe old disks)
        a_avoid = varargin{6}; % radii of the disks to be avoided
    end
    if(opt_argin >= 7)% dmin for the avoiding point/disks
        dmin_avoid = varargin{7};
    end
    if(opt_argin >= 8)% dmax for the avoiding point/disks
        dmax_avoid = varargin{8};
    end
    
    %Verification    
    if(amin > amax)
       aux = amin;
       amin = amax;
       amax = aux;
       warning('amin was greater than amax, exchanging values...');
    end
    
    if(dmax <=0)
        dmax = realmax;
    end
    
    if(~isempty(O_avoid) && isempty(a_avoid))
       a_avoid=zeros(1,size(O_avoid,2)); 
    end
    % Placement
    for p = 1:N_scat
        trial = 1;
        [centre_new, radius_new] = BuildRandomDisk(xmin, xmax, ymin, ymax, amin, amax);
        isOK = testDistance(O(:,1:p-1), a(1:p-1), xmin, xmax, ymin, ymax, dmin, dmax, centre_new, radius_new, O_avoid, a_avoid, dmin_avoid, dmax_avoid);
        while ( ~isOK && trial < trial_max)
            [centre_new, radius_new] = BuildRandomDisk(xmin, xmax, ymin, ymax, amin, amax);
            isOK = testDistance(O(:,1:p-1), a(1:p-1), xmin, xmax, ymin, ymax, dmin, dmax, centre_new, radius_new, O_avoid, a_avoid, dmin_avoid, dmax_avoid);
            trial = trial+1;
        end     
        if (trial == trial_max)
            error(['Cannot place obstacle number ',num2str(p)]);
        end
        %Obstacle is ok to be placed !
        O(:,p) = centre_new;
        a(p) = radius_new;
    end
end

%%
%
function [centre_new, radius_new] = BuildRandomDisk(xmin, xmax, ymin, ymax, amin, amax)
    
    real_xmin = xmin + amin;
    real_xmax = xmax - amin;
    real_ymin = ymin + amin;
    real_ymax = ymax - amin;
    
    x = real_xmin + (real_xmax - real_xmin).*rand(1,1);
    y = real_ymin + (real_ymax - real_ymin).*rand(1,1);
    %A possibly new disk:
    centre_new = [x;y];
    radius_new = amin + (amax - amin).*rand(1,1);
end


%%
function isOK = testDistance(O, a, xmin, xmax, ymin, ymax, dmin, dmax, centre_new, radius_new, O_avoid, a_avoid, dmin_avoid, dmax_avoid)

    OO = [O, O_avoid];
    aa = [a, a_avoid];
    x = centre_new(1);
    y = centre_new(2);
    test_dist = CheckPlacement(OO, aa, dmin, dmax, centre_new, radius_new);
    test_xmin = (x - radius_new >= xmin);
    test_xmax = (x + radius_new <= xmax);
    test_ymin = (y - radius_new >= ymin);
    test_ymax = (y + radius_new <= ymax);
    isOK = test_xmin*test_xmax*test_ymin*test_ymax*test_dist;
end

