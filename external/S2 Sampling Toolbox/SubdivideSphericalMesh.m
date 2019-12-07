function [TR,G]=SubdivideSphericalMesh(TR,k,G,W)
% Subdivide triangular (or quadrilateral) surface mesh, representing 
% a zero-centered unit sphere, k times using triangular (or quadrilateral)
% quadrisection. See function 'TriQuad' (or 'QuadQuad') for more info.
%
% INPUT:
%   - TR   : surface mesh represented as an object of 'TriRep' class,
%            'triangulation' class, or a cell such that TR={Tri,X}, where
%            Tri is an M-by-3 array of faces and X is an N-by-3 array of 
%            vertex coordinates.
%   - k    : desired number of subdivisions. k=1 is default.
%   - G    : optional; scalar or vector field defined at the vertices of TR.
%   - W    : optional; positive weights associated with vertices of TR.
%            See function 'TriQuad' (or 'QuadQuad') for more info.
%
% OUTPUT:
%   - TR  : subdivided mesh. Same format as input mesh.
%
% AUTHOR: Anton Semechko (a.semechko@gmail.com)
%


if nargin<2 || isempty(k), k=1; end
if nargin<3, G=[]; end
if nargin<4, W=[]; end


% Get the data structure
[Tri,X,fmt]=GetMeshData(TR);

% Make sure that the input mesh is either triangular or quadrilateral 
if ~(size(Tri,2)==3 || size(Tri,2)==4) || ~isequal(round(Tri),Tri)
    error('''SubdivideSphericalMesh'' only works with triangular and quadrilateral meshes')
end

% Make sure vertices of the input mesh lie on a unit sphere
X=ProjectOnSn(X);

% Return mesh as is if k<1
k=round(k(1));
if k<1 
    TR=MeshOut(Tri,X,fmt);
    return
end

% Spherical subdivision
for i=1:k
    
    % Subdivide mesh
    if size(Tri,2)==3
        [TR,W,G]=TriQuad({Tri X},W,G);
    else
        [TR,~,G]=QuadQuad({Tri X},W,G);
    end
    N1=size(X,1);
    [Tri,X]=deal(TR{1},TR{2});
        
    % Project vertices onto unit sphere 
    N1=N1+1;
    N2=size(X,1);
    X(N1:N2,:)=ProjectOnSn(X(N1:N2,:));

end
TR=MeshOut(Tri,X,fmt);


function TR=MeshOut(Tri,X,fmt)

switch fmt
    case 1
        TR=triangulation(Tri,X);
    case 2
        TR=TriRep(Tri,X); %#ok<*DTRIREP>
    case 3
        TR={Tri X};
    case 4
        TR=struct('faces',Tri,'vertices',X);
end

