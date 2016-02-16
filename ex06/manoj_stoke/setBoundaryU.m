function [M,I] = setBoundaryU(M,I,nx,ny)

rightSide =  nx+1 : nx+1 : (nx+1) * (ny+2);
leftSide =   1:nx+1:(nx+1)*(ny+2);
topSide    = (nx+1)*(ny+2) -(nx+1)+1 : (nx+1)*(ny+2);
bottomSide = 1:(nx+1);

top_inside = topSide - (nx+1);
bottom_inside = bottomSide + (nx+1);

M(leftSide,: ) = I(leftSide,:);
M(:,leftSide) = I(:,leftSide);

M(rightSide,: ) = I(rightSide,:);
M(:,rightSide) = I(:,rightSide);

M(topSide ,:) = 0.5*(I(topSide,:) + I(top_inside,:));
M(bottomSide ,:) = 0.5*(I(bottomSide,:) + I(bottom_inside,:));


end