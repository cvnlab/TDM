function [TR,W,G,SM,idx_unq]=TriQuad(TR,W,G)
% Subdivide a triangular surface mesh using generalized triangular 
% quadrisection. Triangular quadrisection is a linear subdivision procedure
% which inserts new vertices at the edge midpoints of the input mesh, 
% thereby producing four new faces for every face of the original mesh.
% Illustration of this operation is provided below:
% 
%                     x3                        x3
%                    /  \      subdivision     /  \
%                   /    \        ====>       v3__v2
%                  /      \                  / \  / \
%                x1________x2              x1___v1___x2
%
%                      Original vertices : x1, x2, x3
%                      New vertices      : v1, v2, v3
%                      New faces         : [x1 v1 v3; x2 v2 v1; x3 v3 v2; v1 v2 v3] 
%
% In case of generalized triangular quadrisection, positions of the newly
% inserted vertices do not have to correspond to edge midpoints, and may be
% varied by assigning (positive) weights to the vertices of the original 
% mesh. For example, let xi and xj be two vertices connected by an edge, 
% and suppose that Wi and Wj are the corresponding vertex weights. Position
% of the new point on the edge (xi,xj) is defined as (Wi*xi+Wj*xj)/(Wi+Wj).
% Note that in order to avoid degeneracies and self-intersections, all 
% weights must be real numbers greater than zero.
%
% INPUT ARGUMENTS:
%   - TR   : surface mesh represented as an object of 'TriRep' class,
%            'triangulation' class, or a cell such that TR={Tri,X}, where
%            Tri is an M-by-3 array of faces and X is an N-by-3 array of 
%            vertex coordinates.
%   - W    : optional input argument. N-by-1 array of NON-ZERO, POSITIVE 
%            vertex weights used during interpolation of the new vertices, 
%            where N is the total number of the original mesh vertices. 
%   - G     : scalar or vector field defined on the mesh vertices.
%
% OUTPUT:
%   - TR  : subdivided mesh. Same format as the input.
%   - W   : new set of vertex weights.
%   - G   : interpolated scalar or vector field.
%   - SM  : subdivision matrix
%
% AUTHOR: Anton Semechko (a.semechko@gmail.com)
%


% Get the list of vertex co-ordinates and list of faces
[Tri,X,fmt]=GetMeshData(TR);

% Make sure that the mesh is composed entirely of triangles 
if size(Tri,2)~=3
    error('This function is meant for triangular surface meshes, not quad meshes')
end

% Check vertex weights
flag=false;
if nargin<2 || isempty(W), flag=true; end

if ~flag && (numel(W)~=size(X,1) || sum(W<=eps))
    error('W must be a N-by-1 array with non-negative entries, where N is the # of mesh vertices')
end
if ~flag, W=W(:)+eps; end

% Field?
if nargin==3 
    if ~isempty(G) && (size(G,1)~=size(X,1) || ~ismatrix(G))
        error('3-rd input argument must be a %u-by-d matrix where d>0',size(X,1))
    end
else
    G=[];
end


% Compute new vertex positions
if ~flag
    w=bsxfun(@rdivide,[W(Tri(:,1)),W(Tri(:,2))],W(Tri(:,1))+W(Tri(:,2)));
    V1=bsxfun(@times,X(Tri(:,1),:),w(:,1))+bsxfun(@times,X(Tri(:,2),:),w(:,2));
    if ~isempty(G)
        G1=bsxfun(@times,G(Tri(:,1),:),w(:,1))+bsxfun(@times,G(Tri(:,2),:),w(:,2));
    end
    w1=W(Tri(:,1)).*w(:,1)+W(Tri(:,2)).*w(:,2);
    
    w=bsxfun(@rdivide,[W(Tri(:,2)),W(Tri(:,3))],W(Tri(:,2))+W(Tri(:,3)));
    V2=bsxfun(@times,X(Tri(:,2),:),w(:,1))+bsxfun(@times,X(Tri(:,3),:),w(:,2));
    if ~isempty(G)
        G2=bsxfun(@times,G(Tri(:,2),:),w(:,1))+bsxfun(@times,G(Tri(:,3),:),w(:,2));
    end
    w2=W(Tri(:,2)).*w(:,1)+W(Tri(:,3)).*w(:,2);
    
    w=bsxfun(@rdivide,[W(Tri(:,3)),W(Tri(:,1))],W(Tri(:,3))+W(Tri(:,1)));
    V3=bsxfun(@times,X(Tri(:,3),:),w(:,1))+bsxfun(@times,X(Tri(:,1),:),w(:,2));
    if ~isempty(G)
        G3=bsxfun(@times,G(Tri(:,3),:),w(:,1))+bsxfun(@times,G(Tri(:,1),:),w(:,2));
    end
    w3=W(Tri(:,3)).*w(:,1)+W(Tri(:,1)).*w(:,2);
    
    W_new=[w1;w2;w3];
else
    V1=(X(Tri(:,1),:)+X(Tri(:,2),:))/2;
    V2=(X(Tri(:,2),:)+X(Tri(:,3),:))/2;
    V3=(X(Tri(:,3),:)+X(Tri(:,1),:))/2;
    
    if ~isempty(G)
        G1=(G(Tri(:,1),:)+G(Tri(:,2),:))/2;
        G2=(G(Tri(:,2),:)+G(Tri(:,3),:))/2;
        G3=(G(Tri(:,3),:)+G(Tri(:,1),:))/2;
    end
end
V=[V1;V2;V3];

% Remove repeating vertices 
E=[Tri(:,1) Tri(:,2);Tri(:,2) Tri(:,3);Tri(:,3) Tri(:,1)];
E=sort(E,2);
[~,idx_unq,idx]=unique(E,'rows','stable'); % setOrder='stable' ensures that identical results will be obtained for meshes with same connectivity
V=V(idx_unq,:);
if ~flag
    W=[W;W_new(idx_unq)]; 
end
if ~isempty(G)
    G_new=[G1;G2;G3];
    G=[G;G_new(idx_unq,:)];
end

% Generate a subdivision matrix if one is required
Nx=size(X,1);   % # of vertices
Nt=size(Tri,1); % # of faces
if nargout>3
    
    idx_1=[ones(Nt,1);2*ones(Nt,1);3*ones(Nt,1)];
    idx_1=idx_1(idx_unq);
    idx_2=idx_1+1;
    idx_2(idx_2>3)=1;
    
    tri=repmat(Tri,[3 1]);
    tri=tri(idx_unq,:);
    idx_1=sub2ind(size(tri),(1:size(tri,1))',idx_1);
    idx_2=sub2ind(size(tri),(1:size(tri,1))',idx_2);
    
    idx_1=tri(idx_1);
    idx_2=tri(idx_2);    
    
    idx_0=(1:size(V,1))';
    idx_0=[idx_0;idx_0];
    
    SM=sparse(idx_0,[idx_1;idx_2],0.5*ones(size(idx_0)));
    SM=cat(1,speye(Nx),SM);
    
end

% Assign indices to the new triangle vertices
V1= Nx + idx(1:Nt);
V2= Nx + idx((Nt+1):2*Nt);
V3= Nx + idx((2*Nt+1):3*Nt);
clear idx

% Define new faces
T1= [Tri(:,1) V1 V3];
T2= [Tri(:,2) V2 V1];
T3= [Tri(:,3) V3 V2];
T4= [V1       V2 V3];
clear V1 V2 V3

T1=permute(T1,[3 1 2]);
T2=permute(T2,[3 1 2]);
T3=permute(T3,[3 1 2]);
T4=permute(T4,[3 1 2]);

Tri=cat(1,T1,T2,T3,T4);
Tri=reshape(Tri,[],3,1);

% New mesh
X=[X;V]; 
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
if nargout>1 && flag, W=ones(size(X,1),1); end

