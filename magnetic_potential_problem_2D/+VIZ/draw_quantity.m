function draw_quantity(coord,elem,U,Q_node,size_x,size_y)

% =========================================================================
%
%  This function depicts prescribed nodal quantity
%
%  input data:
%    coord - coordinates of the nodes, size(coord)=(2,n_n) where n_n is a
%            number of nodes
%    surf - n_p_s x n_s array containing numbers of nodes defining each
%           surface element, n_s = number of surface elements
%    U - nodal displacements, size(U)=(2,n_n) to catch deformed shape
%        if the deformed shape is not required then set 0*U
%    Q_node - prescribed nodal quantity, size(Q_node)=(1,n_n)
%    size_x - length of the rectangular body in x-direction
%    size_y - length of the rectangular body in y-direction
%
% ======================================================================
%

  figure;
  hold on;
  
  % visualization of the quantity
  s = patch('Faces',elem(1:3,:)','Vertices',coord',...
        'FaceVertexCData',Q_node','FaceColor','interp','EdgeColor','none'); 
  
%   alpha(s,.5);
  colorbar;
  
  % undeformed shape of the body
  plot([0,size_x],[0,0])
  plot([size_x,size_y],[0,size_y])
  plot([size_x,0],[size_y,size_y])
  plot([0,0],[size_x,0])
   
  %
  box on
  view(2);

  axis equal;
  hold off;
  axis off;
end