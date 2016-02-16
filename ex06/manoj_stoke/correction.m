function M_corr = correction(M_temp,nx,ny)

u_row = (nx+2) :(nx+2) : (nx+2) * (ny+2); % deleting nx multiple for u
v_col = 2*(nx+2)*(ny+2) - (nx+1) : 2*(nx+2)*(ny+2); % deleting last columns in v

M_temp(u_row,:) = [];
M_temp(:,u_row) = [];

v_col = v_col - (nx+2);

M_temp(v_col,:) = [];
M_temp(:,v_col) = [];

M_corr = M_temp;


end