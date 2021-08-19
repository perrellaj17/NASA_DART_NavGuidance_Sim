

function [M, B_nom] = B_plane_targeting(r_vec, v_vec, mu, alpha)
v_vec
%alpha is size of perturbation 'i.e. the dV increment for checking'

N_hat = [0 0 1]; %vector normal to target planet equatorial plane DO NOT CONFUSE WITH LOWERCASE n HAT
%loop once for the unpurturbed state, once for perturbations in each direction
for i = 1:4
    
    switch i
        %skip case 1, v_vec is unperturbed and remains itself
        case 2
            v_vec = v_vec + [alpha 0 0]' %dvx
%             v_vec = [alpha 0 0]' %dvx
        case 3
            v_vec = v_vec + [0 alpha 0]' %dvy
%             v_vec =[0 alpha 0]' %dvy
        case 4
            v_vec = v_vec + [0 0 alpha]' %dvz
%             v_vec = [0 0 alpha]' %dvz
    end
    v_mag(i) = norm(v_vec)
    r_mag = norm(r_vec)
    % motion of angular momentum direction
    %same direction as normal plane of target body motion
    h_hat(i,:) = cross(r_vec, v_vec) / norm(cross(r_vec, v_vec))
    
    %eccentricity of target - points towards periapsis of hyperbola
    e_vec(i,:) = (1/mu)*(( v_mag(i)^2)*r_vec - dot(r_vec,v_vec)*v_vec)-r_vec/r_mag
    e_mag(i) = norm(e_vec(i,:))
    e_hat(i,:) = e_vec(i,:)/e_mag(i)
    
    %hyperbolic asymptote angle
    Beta(i) = acos(1/e_mag(i)) %radians
    
    %lecture 16 slide 7 uses this eq, 17 slide 4 uses h_hat cross e_hat?????
    
    %basis of the vectors for the B frame coordinate frame
    S_hat(i,:) = cos(Beta(i))*e_hat(i,:) + sin(Beta(i))*(cross(h_hat(i,:), e_hat(i,:))/norm(cross(h_hat(i,:),e_hat(i,:))))
%     S_hat(i,:) = cos(Beta(i))*e_hat(i,:) + sin(Beta(i))*(cross(h_hat(i,:), e_hat(i,:)))
    T_hat(i,:) = cross(S_hat(i,:), N_hat) / norm(cross(S_hat(i,:), N_hat))
%     R_hat(i,:) = cross(S_hat(i,:), T_hat(i,:)) / norm(cross(S_hat(i,:), T_hat(i,:)))
    R_hat(i,:) = cross(S_hat(i,:), T_hat(i,:))
    
    %semi major axis relative to the target body
    a(i) = -mu/(2*(((v_mag(i)^2)/2) - mu/r_mag))
    %semi-minor axis
    b(i) = -a(i)*sqrt(e_mag(i)^2 - 1)
    
    %calculate the true anomaly for position vector and beta vector
    theta(i)  = acos ( (-a(i)*(e_mag(i)^2 - 1))/(r_mag*e_mag(i)) - 1/e_mag(i))
    theta_b(i) = pi - (pi/2 + Beta(i))
    
    %calculate the hyperbolic anomaly for position vector and beta vector
    H_r(i)  = acosh((e_mag(i)+cos(theta(i)))/(1+e_mag(i)*cos(theta(i))))
    H_b(i) = acosh((e_mag(i)+cos(theta_b(i)))/(1+e_mag(i)*cos(theta_b(i))))
    
    %calculate time from periapsis for position vector and B_vector
    t_r(i) = e_mag(i)*sinh(H_r(i)) - H_r(i) / sqrt(mu/(-a(i)^3))
    t_p(i) = e_mag(i)*sinh(H_b(i)) - H_b(i) / sqrt(mu/(-a(i)^3))
    
    
    %(LTOF) is the time from the current epoch until periapsis on the approach hyperbola relative to the targetplanet
    LTOF(i) = t_r(i) - t_p(i)
    
    
    %compute sensitivity Matrix

    %B_vec = B_T*T_hat + B_R*R_hat
    B_vec(i,:) = b(i)*cross(S_hat(i,:), h_hat(i,:))/norm(cross(S_hat(i,:), N_hat(:)))
    
    
    B_T(i) = dot(B_vec(i,:), T_hat(i,:))
    B_R(i) = dot(B_vec(i,:), R_hat(i,:))
    
    %%% IGNORE LTOF %%%
%     B_plane_vec(i,:) = [B_T(i) B_R(i) 60]
    B_plane_vec(i,:) = [B_T(i) B_R(i) LTOF(i)]
    
    %sc position doenst correspond to B plane locations (sc trajectory is
    %curving)
    %the b plane coords talk about where hyperbolic ... the b plane

    % use unperturbed position (r)
    % compute B-vec associated with that, denote as B_1
    % partial B wrt v_x B_r (unpert) - B_1 (perturbetd) / alpha
    if i > 1
        %unperturbed spacecraft state, aka B_star, ignore LTOF
        B_nom =B_plane_vec(1,:)
        
        M(:,i-1) = (B_nom - B_plane_vec(i,:))/alpha      
        %1 is partialB/dx
        %2 is partialB/dy
        %3 is partialB/dz
    end
end

%sensitivity matrix


end
