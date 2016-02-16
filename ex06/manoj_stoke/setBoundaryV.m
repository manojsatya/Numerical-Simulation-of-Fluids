function M = setBoundaryV(M,I,nx,ny)

rightSide =  nx+2 : nx+2 : (nx+2) * (ny+1);
leftSide =   1:nx+2:(nx+2)*(ny+1);
topSide    = (nx+2)*(ny+1) -(nx+2)+1 : (nx+2)*(ny+1);
bottomSide = 1:(nx+2);

rightSide = rightSide +(nx+1)*(ny+2);
leftSide = leftSide + (nx+1)*(ny+2);
topSide = topSide + (nx+1)*(ny+2);
bottomSide = bottomSide + (nx+1)*(ny+2);

left_inside = leftSide + 1 ;
right_inside = rightSide -1 ;

M(leftSide ,:) = 0.5*(I(leftSide,:) + I(left_inside,:));
M(rightSide ,:) = 0.5*(I(rightSide,:) + I(right_inside,:));

M(topSide,: ) = I(topSide,:);
M(:,topSide) = I(:,topSide);

M(bottomSide,: ) = I(bottomSide,:);
M(:,bottomSide) = I(:,bottomSide);

end
