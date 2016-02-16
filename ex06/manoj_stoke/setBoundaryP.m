function M = setBoundaryP(M,I,nx,ny)

rightSide =  nx+2 : nx+2 : (nx+2) * (ny+2);
leftSide =   1:nx+2:(nx+2)*(ny+2);
topSide    = (nx+2)*(ny+2) -(nx+2)+1 : (nx+2)*(ny+2);
bottomSide = 1:(nx+2);

%In Matrix M , adjusting to the array number with p

rightSide = rightSide + (nx+1)*(ny+2) + (nx+2)*(ny+1);
leftSide = leftSide + (nx+1)*(ny+2) + (nx+2)*(ny+1);
topSide = topSide + (nx+1)*(ny+2) + (nx+2)*(ny+1);
bottomSide = bottomSide + (nx+1)*(ny+2) + (nx+2)*(ny+1);

right_inside = rightSide - 1 ;
left_inside = leftSide + 1 ;
top_inside = topSide - (nx+2);
bottom_inside = bottomSide + (nx+2);


%%%%% Applying Neumann Boundary Conditions
M(leftSide,:)   = -I(leftSide,:) + I(left_inside,:);
M(rightSide,:)   = -I(rightSide,:) + I(right_inside,:);
M(topSide,:)   = -I(topSide,:) + I(top_inside,:);
M(bottomSide,:)   = -I(bottomSide,:) + I(bottom_inside,:);



end