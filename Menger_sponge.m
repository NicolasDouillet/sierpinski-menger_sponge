function [V, T] = Menger_sponge(nb_it, printable_ready, option_display)
%% Menger_sponge : function to compute, display, and save
% the Sierpinski-Menger sponge at any iteration / depth level.
%
% Author : nicolas.douillet9 (at) gmail.com, 2019-2024.
%
%
% Syntax
%
% Menger_sponge(nb_it);
% Menger_sponge(nb_it, printable_ready);
% Menger_sponge(nb_it, printable_ready, option_display);
% [V, T] = Menger_sponge(nb_it, printable_ready, option_display);
%
%
% Description
%
% Menger_sponge(nb_it) computes and displays the nb_it
% Sierpinski-Menger sponge included in the unit sphere.
%
% Menger_sponge(nb_it, printable_ready) prevents from
% creating non manifold edges when printable_ready is set to true /
% logical 1, and remove duplicated vertices and faces when it is set to
% *false / logical 0. In this latter case, the model is lighter (less
% vertices, less faces), but at the cost of non manifoldness.
%
% Menger_sponge(nb_it, printable_ready, option_display)
% displays it when option_display is set to logical *true/1 (default),
% and doesn't when it is set to  logical false/0.
%
% [V,T] = Menger_sponge(nb_it, printable_ready, option_display) saves
% the resulting vertex coordinates in the array V, and the triangulation in the array T.
%
%
% Input arguments
%
% - nb_it : positive integer scalar double, the number of iterations / depth level.
%
% - printable_ready : either logical, true/*false or numeric 1/*0.
%
% - option_display : either logical, *true/false or numeric *1/0.
%
%
% Output arguments
%
%       [ |  |  |]
% - V = [Vx Vy Vz], real matrix double, the vertex coordinates. Size(V) = [nb_vertices,3].
%       [ |  |  |]
%
%       [ |  |  |]
% - T = [T1 T2 T3], positive integer matrix double, the triangulation. Size(T) = [nb_triangles,3].
%       [ |  |  |]
%
%
% Example #1
%
% Computes and displays the Menger sponge
% at iteration 2, with minimum vertex and face numbers
%
% Menger_sponge(2);
%
%
% Example #2
%
% Computes and saves the 3D printable ready
% Menger sponge at iteration 3
%
% [V,T] = Menger_sponge(3,true,false);


%% Input parsing
assert(nargin < 4,'Too many input arguments.');

if ~nargin
    nb_it = 3;
    printable_ready = false;
    option_display = true;
elseif nargin > 0
    assert(isnumeric(nb_it) && nb_it == floor(nb_it) && nb_it >= 0,'nb_it parameter value must be numeric positive or null integer.');
    if nargin > 1
        assert(islogical(printable_ready) || isnumeric(printable_ready),'printable_ready parameter type must be either logical or numeric.');
        if nargin > 2
            assert(islogical(option_display) || isnumeric(option_display),'option_display parameter type must be either logical or numeric.');
        else
            option_display = true;
        end
    else
        printable_ready = false;
        option_display = true;
    end
end

warning('on');
if option_display && nb_it > 3
    warning('%s triangles to display ! Make sure your graphic card has enough memory.',num2str(12*20^nb_it))    
end
warning('off');


%% Body
% Summits of original cube (living in the unit sphere R(O,1))
a = sqrt(3)/3;

V1 = [a a a];
V2 = [-a a a];
V3 = [-a -a a];
V4 = [a -a a];
V5 = -V1;
V6 = -V2;
V7 = -V3;
V8 = -V4;

C = cube(V1, V2, V3, V4, V5, V6, V7, V8);


% Loop on nb_it
p = 0;

while p ~= nb_it 
    
    new_C_array = repmat(C, [1 1 20]);
    
    for j = 1 : size(C,3)
        
        C_current = C(:,:,j);        
        [V_new, F_new] = split_cube(C_current);
        
        for m = 1:size(F_new,1)/6 %  20 = 120 / 6
            
            new_cube = cube(V_new(F_new(6*(m-1)+1,1),:),...
                            V_new(F_new(6*(m-1)+1,2),:),...
                            V_new(F_new(6*(m-1)+1,3),:),...
                            V_new(F_new(6*(m-1)+1,4),:),...
                            V_new(F_new(6*(m-1)+2,1),:),...
                            V_new(F_new(6*(m-1)+2,2),:),...
                            V_new(F_new(6*(m-1)+2,3),:),...
                            V_new(F_new(6*(m-1)+2,4),:)); % first two lines vertices in this order
                        
            new_C_array(:,:,20*(j-1) + m) = new_cube;
            
        end
        
    end
    
    C = new_C_array;        
    p = p+1;
    
end

% Squares to triangles conversion
[V,T] = squares2triangles(C);

if ~printable_ready
    
    % Remove duplicated vertices
    [V,T] = remove_duplicated_vertices(V,T);
        
    % Remove duplicated triangles
    T = unique(sort(T,2),'rows','stable');
    
end

%% Display
if option_display
    
    cmap = [1 1 0];
    disp_Menger_sponge(V,T,cmap);
    
end


end % Menger_sponge


%% Split cube subfunction
function [V_new, F_new] = split_cube(C)
% Input
%
% C : cube structure
%
%
% Outputs
%
% V_new : The 8 + 4*6 + 8 = 40 (8 old + 32 new) newly created vertices -coordinates-
% F_new : 6 x 20 = 120 newly created facets line indices


one_third = 1/3;
two_third = 2/3;

% Four possible values for each coordinate, X, Y, or Z
val1x = min(C.vertex(:,1)); % min value
val4x = max(C.vertex(:,1)); % max value
val2x = val1x + one_third * (val4x - val1x);
val3x = val1x + two_third * (val4x - val1x);

val1y = min(C.vertex(:,2)); % min value
val4y = max(C.vertex(:,2)); % max value
val2y = val1y + one_third * (val4y - val1y);
val3y = val1y + two_third * (val4y - val1y);

val1z = min(C.vertex(:,3)); % min value
val4z = max(C.vertex(:,3)); % max value
val2z = val1z + one_third * (val4z - val1z);
val3z = val1z + two_third * (val4z - val1z);


% Compute  (8) + 4*6 + 8 = 40 new vertices coordinates

V_new = [val4x val4y val4z; % Top layer
         val3x val4y val4z; 
         val2x val4y val4z;
         val1x val4y val4z;
         val4x val3y val4z;
         val3x val3y val4z;
         val2x val3y val4z;
         val1x val3y val4z;
         val4x val2y val4z;
         val3x val2y val4z;
         val2x val2y val4z;
         val1x val2y val4z;
         val4x val1y val4z;
         val3x val1y val4z;
         val2x val1y val4z;
         val1x val1y val4z;...
         
         val4x val4y val3z; % Bottom first layer
         val3x val4y val3z; 
         val2x val4y val3z;
         val1x val4y val3z;
         val4x val3y val3z;
         val3x val3y val3z;
         val2x val3y val3z;
         val1x val3y val3z;
         val4x val2y val3z;
         val3x val2y val3z;
         val2x val2y val3z;
         val1x val2y val3z;
         val4x val1y val3z;
         val3x val1y val3z;
         val2x val1y val3z;
         val1x val1y val3z;...
         
         val4x val4y val2z; % Bottom second layer
         val3x val4y val2z; 
         val2x val4y val2z;
         val1x val4y val2z;
         val4x val3y val2z;
         val3x val3y val2z;
         val2x val3y val2z;
         val1x val3y val2z;
         val4x val2y val2z;
         val3x val2y val2z;
         val2x val2y val2z;
         val1x val2y val2z;
         val4x val1y val2z;
         val3x val1y val2z;
         val2x val1y val2z;
         val1x val1y val2z;...
         
         val4x val4y val1z; % Bottom face
         val3x val4y val1z; 
         val2x val4y val1z;
         val1x val4y val1z;
         val4x val3y val1z;
         val3x val3y val1z;
         val2x val3y val1z;
         val1x val3y val1z;
         val4x val2y val1z;
         val3x val2y val1z;
         val2x val2y val1z;
         val1x val2y val1z;
         val4x val1y val1z;
         val3x val1y val1z;
         val2x val1y val1z;
         val1x val1y val1z]; 

% 6 x 20 = 120 new facets
% /_!_\ Counter clockwise sorted for squares /_!_\

% General model  for a (a b c d e f g h) cube
% a b c d
% e f g h
% a d h e
% b a e f
% c b f g
% d c g h

F_new = [1 2 6 5; % Top layer top right corner cube
         17 18 22 21;
         1 5 21 17;
         2 1 17 18;
         6 2 18 22;
         5 6 22 21;...
         
         2 3 7 6; % Top layer top cross cube
         18 19 23 22;
         2 6 22 18;
         3 2 18 19;
         7 3 19 23;
         6 7 23 22;...         
         
         3 4 8 7; % Top layer top left corner cube
         19 20 24 23;
         3 7 23 19;
         4 3 19 20;
         8 4 20 24;
         7 8 24 23;...
         
         5 6 10 9; % Top layer right cross cube
         21 22 26 25;
         5 9 25 21;
         6 5 21 22;
         10 6 22 26;
         9 10 26 25;...
         
         7 8 12 11; % Top layer left cross cube
         23 24 28 27;
         7 11 27 23;
         8 7 23 24;
         12 8 24 28;
         11 12 28 27;...
         
         9 10 14 13; % Top layer bottom right corner cube
         25 26 30 29;
         9 13 29 25;
         10 9 25 26;
         14 10 26 30;
         13 14 30 29;...
         
         10 11 15 14; % Top layer bottom cross cube
         26 27 31 30;
         10 14 30 26;
         11 10 26 27;
         15 11 27 31;
         14 15 31 30;...
         
         11 12 16 15; % Top layer bottom left corner cube
         27 28 32 31;
         11 15 31 27;
         12 11 27 28;
         16 12 28 32;
         15 16 32 31;...         
         
         17 18 22 21; % Middle layer top right corner cube
         33 34 38 37;
         17 21 37 33;
         18 17 33 34;
         22 18 34 38;
         21 22 38 37;...                      
         
         19 20 24 23; % Middle layer top left corner cube
         35 36 40 39;
         19 23 39 35;
         20 19 35 36;
         24 20 36 40;
         23 24 40 39;...                                    
         
         25 26 30 29; % Middle layer bottom right corner cube
         41 42 46 45;
         25 29 45 41;
         26 25 41 42;
         30 26 42 46;
         29 30 46 45;...                  
         
         27 28 32 31; % Middle layer bottom left corner cube
         43 44 48 47;
         27 31 47 43;
         28 27 43 44;
         32 28 44 48;
         31 32 48 47;...          
         
         33 34 38 37; % Bottom layer top right corner cube
         49 50 54 53;
         33 37 53 49;
         34 33 49 50;
         38 34 50 54;
         37 38 54 53;...
         
         34 35 39 38; % Bottom layer top cross cube
         50 51 55 54;
         34 38 54 50;
         35 34 50 51;
         39 35 51 55;
         38 39 55 54;...         
         
         35 36 40 39; % Bottom layer top left corner cube
         51 52 56 55;
         35 39 56 51;
         36 35 51 52;
         40 36 52 56;
         39 40 56 55;...
         
         37 38 42 41; % Bottom layer right cross cube
         53 54 58 57;
         37 41 57 53;
         38 37 53 54;
         42 38 54 58;
         41 42 58 57;...
         
         39 40 44 43; % Bottom layer left cross cube
         55 56 60 59;
         39 43 59 55;
         40 39 55 56;
         44 40 56 60;
         43 44 60 59;...
         
         41 42 46 45; % Bottom layer bottom right corner cube
         57 58 62 61;
         41 45 61 57;
         42 41 57 58;
         46 42 58 62;
         45 46 62 61;...
         
         42 43 47 46; % Bottom layer bottom cross cube
         58 59 63 62;
         42 46 62 58;
         43 42 58 59;
         47 43 59 63;
         46 47 63 62;...
         
         43 44 48 47; % Bottom layer bottom left corner cube
         59 60 64 63;
         43 47 63 59;
         44 43 59 60;
         48 44 60 64;
         47 48 64 63];

end % split_cube


%% Build cube structure subfunction
function [C] = cube(V1, V2, V3, V4, V5, V6, V7, V8)
%
% V1, V2, V3, V4, V5, V6, V7, V8 : line vectors of cube eight vertices coordinates
%
% Vn = [Vxn Vyn Vzn]

F1 = [1 2 3 4];
F2 = [5 6 7 8];
F3 = [1 4 8 5];
F4 = [2 1 5 6];
F5 = [2 3 7 6];
F6 = [3 4 8 7];

C = struct('vertex', [V1; V2; V3; V4; V5; V6; V7; V8], ...
           'facet', [F1; F2; F3; F4; F5; F6]);
       
end % cube


%% Squares to triangles conversion subfunction
function [V, T] = squares2triangles(C)
%
% Author : nicolas.douillet9 (at) gmail.com, 2017-2024.
%
% Split struct array into two arrays : vertices & facets

S = size(C,3);
V = zeros(8*S,3);
T = zeros(12*S,3);

for k = 1:S
    
    for i = 1:size(C(:,:,k).vertex,1)
        
        V(8*(k-1)+i,:) = C(:,:,k).vertex(i,:);
        
    end
    
    for j = 1:size(C(:,:,k).facet,1) % 6
        
        a = C(:,:,k).facet(j,1) + 8*(k-1);
        b = C(:,:,k).facet(j,2) + 8*(k-1);
        c = C(:,:,k).facet(j,3) + 8*(k-1);
        d = C(:,:,k).facet(j,4) + 8*(k-1);
        
        T1 = sort([a b c]);
        T2 = sort([a d c]);
        
        T(12*(k-1)+2*(j-1)+1,:) = T1;
        T(12*(k-1)+2*j,:) = T2;
        
    end
    
end

% Remove duplicated triangles
T = unique(sort(T,2),'rows','stable');

end % squares2triangles


%% Display subfunction
function [] = disp_Menger_sponge(V, T, cmap)
%
% Author : nicolas.douillet9 (at) gmail.com, 2017-2024.

figure;
set(gcf,'Color',[0 0 0]), set(gca,'Color',[0 0 0]);
trisurf(T,V(:,1),V(:,2),V(:,3),'EdgeColor',cmap), shading interp, hold on;
colormap(cmap);
axis square, axis equal, axis tight, axis off;
grid off;
ax = gca;
ax.Clipping = 'off';
camlight left;
view(-45,35);

end % disp_Menger_sponge


%% Remove duplicated vertices subfunction
function [V_out, T_out] = remove_duplicated_vertices(V_in, T_in)

tol = 1e4*eps;
[V_out,~,n] = uniquetol(V_in,tol,'ByRows',true);
T_out = n(T_in);

end % remove_duplicated_vertices