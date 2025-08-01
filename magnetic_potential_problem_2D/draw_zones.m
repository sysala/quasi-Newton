function draw_zones(coord,elem,heter)

% =========================================================================
%
%  This function draws elements with a plastic response
%
%  input data:
%    coord - coordinates of the nodes, size(coord)=(2,n_n) where n_n is a
%            number of nodes
%    elem  - array containing numbers of nodes defining each element, 
%            size(elem)=(n_p,n_e), n_e = number of elements
%    heter - a logical array indicating the heterogeneity
%
% ======================================================================
%

  figure
  hold on
  coord_aux = [coord;zeros(1,size(coord,2))];
  patch('Faces',elem(1:3,~heter)','Vertices',coord_aux','FaceVertexCData',...
           0*ones(size(coord,2),1),'FaceColor','red','EdgeColor','blue'); 
  patch('Faces',elem(1:3,heter)','Vertices',coord_aux','FaceVertexCData',...
           0*ones(size(coord,2),1),'FaceColor','white','EdgeColor','blue'); 

  axis equal;  % real ratios
  view(2);     % standard view in 2D
  hold off;
%   axis off;
end