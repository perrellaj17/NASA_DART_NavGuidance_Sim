clear all, clc, clf

G = 6.6743E-20; %km^3/(kg*s^2)
mu_S = 1.3271244004193938E11; %km^3/s^2
mu_E = 3.98600440E5; %km^3/s^2
AU = 149597870.7; %km
a_E = 0.9994295585*AU; %km
i_E = 0.0028065954; %deg earth inclination
r_E_eq = 6378.1363; %km mean equitorial radius of the Earth km
g = 9.80665; %m/s this is a constant for the ideal rocket equation

epsilon = 0.409092802283074;
T_HCI_ECI = [1, 0, 0;
    0 ,cos(epsilon), -sin(epsilon);
    0, sin(epsilon), cos(epsilon)];
%transform heliocentric inertial frame to earth centered intertial

%get ephemerides
Earth_ephem = xlsread('Earth_ephem.xlsx');
cb21_ephem = xlsread('cb21_ephem.xlsx');
didymos_ephem = xlsread('Didymos_ephem.xlsx');





%% Project Overview
%{


Overview 



%}
%%%%




%% Mission Parameters 

%{

Parameters

%}



%% Body Plane targeting

%{
B-plane (“Body plane”) targeting
Objectives: Calculate the mid-course corrections, particularly stochastic Trajectory
Correction Maneuvers (TCMs), in response to perturbations and
general errors

Exert control over the parameters (e.g., orientation, periapsis) of
capture orbit around destination planet

The SC is flying along trajectory w/ b plane target parameters we want to
achieve, before arrival, make estiamtes of how we are doing in achieving
the necessary b plane targets. If it looks like we are off, we use this
algorithm to determine what small change in velocity is needed in the
spacecraft to meet b plane targets. usually for small dv changes.



B-plane state at impact (nominally all 0 if we wish to target the exact center of the body)


%}




%B_plane_vec_star = B_plane_vec + M*dv
%now solve for difference in velocity needed to drive the B-plane
%coordinates to target values

%convergence is defined according to difference between specified B plane
%values and currently achieved values

%inertial positiona and velocity relative to the planet!!!!!
r_vec = [-64961.204483358
    -133190.207105847
    0];

v_init = [1.37629832487732;
    2.39127882923125;
0 ];

mu = 42828; %mars mu

alpha = 1E-4; %perturbation value

%note no LTOF requirement
%B vector star is the given targen values

B_targ = [6000 6000 500]; %km

epsilon = 1E-6; %convergance tolerance
dV = [0; 0; 0]; %initial perturbation value
dB = inf; %initialize for first loop

disp('THIS CURRENTLY DIVERGES VERY QUICKLY!')
%lecture 17, 22:12 he starts talking about this
for i = 1:2000
    %CHECK IF SET DV TO ZERO IT WOULD IMMEDIATELY EXIT
    %on first iteration dV is 0. Then it will grow if not converged
%     if 1 == 1
%         v = v_init;
%     end
%     
%     v = v + dV(:,i)

  
    v = v_init + dV(:,i)
    
    %calc Sensitivity matrix and  Bt Br and LTOF
    [M, B_nom] = B_plane_targeting(r_vec, v, mu, alpha)

    %B_nom is the unperturbed
    %B_nom = [B_T(i) B_R(i) LTOF(i)]
    
    %calculate the magnitude difference between target and B plane with
    %this current dV value and no perturbation
    dB = B_targ' - B_nom'
 
    
    %if convergence met, return successful dV value, and exit loop
    if norm(dB) < epsilon
        dV(:,i) %dV that triggered the convergence
        break
    else
        dV(:,i+1) =dV(:,i) + inv(M)*B_targ'        
        %if you havent converged yet, then take sensitivy matrix and iterate
%         dV(:,i+1) = dV(:,i) + inv(M)*dB
        %when we loop again, we check this new dV(i) to see if it results in the
        %b plane target.
    end
end



%% Gravity assist

%{
FLYBY SEQUENCE CONSTRUCTION
-Lambert targeting, or differential correction, to construct trajectories
between flybys
	-Target from initial location to planet
	-Can insert Deep Space Maneuver along the way and see if optimizer drives it to zero
	-Target from planet flyby to next location (final destination Didymos)
	-Targeter/optimizer drives velocity discontinuity to zero at flyby
-Powered vs. unpowered gravity assists
	-Unpowered gravity assist includes no deterministic spacecraft
maneuvers during the planetary flyby
	-Powered gravity assist does include a deterministic spacecraft
maneuver during the planetary flyby (usually at the flyby hyperbola
periapsis relative to the planet) 

WILL I INCORPORATE POWER INTO THE ASSIST?
        -Hard to know a-priori whether including the maneuver will serve to
reduce overall mission Dv; best to allow the maneuver and see if the
optimizer drives it to zero

The spacecraft can flyby behind the planet or in front of the planet,
relative to the Sun
    -Behind the planet: Increase spacecraft’s heliocentric speed
    -In front of the planet: Decrease spacecraft’s heliocentric speed

Energy is conserved on the planetary flyby hyperbola; if it is an
unpowered flyby, then the pre- and post-flyby spacecraft v1 (relative
to the planet) are the same
    -But the heliocentric inertial direction of the pre- and post-flyby ~v1 has
changed (due to the effect of the planet’s gravity) and this change in
~v1 direction on either side of the flyby is what provide the net
heliocentric dv for the spacecraft

Assumptions:
    -“Impulsive” gravity assist: The flyby planet’s heliocentric position and
velocity are constant during the flyby
    -Because of the “impulsive” assumption, and using the patched conics
framework, we can use the spacecraft’s known pre-flyby heliocentric
velocity, the known position of the planet at the flyby epoch, and
dvGA to calculate the spacecraft’s post-flyby heliocentric orbit


Maximizing dvGA is not always the objective; often the key to
optimzing a mission’s overall flight plan hinges on obtaining the
needed value of dvGA at the proper time
%}



%this should be known by state of s/c and planet at flyby epoch
v_inf_i = nan;%v inf before flyby periapsis, incoming v inf
v_inf_i_mag = norm(v_inf_i)
v_inf_f = nan; %v inf after flyby periapsis, outgoing v inf
v_inf_f_mag = norm(v_inf_f)
%NOTE THAT THE MAGNITUDE OF THESE VECTORS IS THE SAME, SO YOU CAN INPUT
%EITHER AS V INF MAG BECASUE ENERGY CONSERVATION
v_inf_mag  = v_inf_i

dV_GA = v_inf_f - v_inf_i


%eccentricity of flyby hyperbola
e_hyp = 1 + (r_p * v_inf_mag^2)/mu_p
% We can compute the hyperbolic turning angle, delta
delta = 2*asin(1/e)

%r_p is the periapsis radius of the flyby hyperbola
%MOST SELECTABLE PARAMETER HER IS RADIUS OF PERIAPSIS

%WHAT TO DO IF MU P IS UNKNOWN??

%can clculate either way
dV_GA_mag = 2*v_inf_mag/e_hyp

dV_GA_mag = 2*v_inf_mag*sin(delta/2)


% frame A, 1 x_hat points along v_inf_i, and z_hat is normal to plane
% spanned by incoming and outgoing v_inf
vec = [cos(-delta); 
    sin(-delta); 
    0];

v_inf_f_A = v_inf_i_mag*vec
%now we must transform that vector to HCI

%what are r and v vectors?
e_hat = ((( v_inf_i_mag^2)/mu_p - 1/r)*r_vec - dot(r_vec,v_vec)*v_vec/mu_p)/norm(((( v_inf_i_mag^2)/mu_p - 1/r)*r_vec - dot(r_vec,v_vec)*v_vec/mu_p))


x_hat = v_inf_i/v_inf_i_mag
z_hat = cross(x_hat, e_hat)/norm(cross(x_hat, e_hat))
y_hat = cross(z_hat, x_hat)/norm(cross(z_hat, x_hat))

T = [x_hat y_hat z_hat]

v_inf_f = T*v_inf_f_A

%patched conics construction to get post flyby heliocentric velocity
v_s_f = v_p + v_inf_f

%via impulsive flyby assumption, apply patched conics again
% r_s_f = r_s_i = r_p


%% Lambert

Earth_ephem(:,2) = []; %remove NaN column
cb21_ephem(:,2) = []; %remove NaN column
didymos_ephem(:,2) = []; %remove NaN column


Earth_ephem(:,1) = Earth_ephem(:,1)- 2400000.5;%convert JD to MJD
cb21_ephem(:,1) = cb21_ephem(:,1)- 2400000.5;%convert JD to MJD
didymos_ephem(:,1) = didymos_ephem(:,1) - 2400000.5; %convert JD to MJD
disp('----- Joe Perrella ENAE741 DART Analysis -----')
disp('     ')
disp('-- Mission Concept of Operations --')
disp('(1) Initial location - Earth.')
disp('(2) Intermediate gravity assist - Asteroid cb21.')
disp('(3) Kinetic impact - Didymos binary asteroid moon.')
disp('     ')
disp('-- Mission Requirements --')

%ephem(1) = JD
%ephem(2) = r_x
%ephem(3) = r_y
%ephem(4) = r_z
%ephem(5) = v_x
%ephem(6) = v_y
%ephem(7) = v_z

%organize data
r_E = Earth_ephem(:,2:4);
r_cb =cb21_ephem(:,2:4);
r_d = didymos_ephem(:,2:4);

v_E = Earth_ephem(:,5:7);
v_cb =cb21_ephem(:,5:7);
v_d = didymos_ephem(:,5:7);

MJD_E = Earth_ephem(:,1);
MJD_cb = cb21_ephem(:,1);
MJD_d = didymos_ephem(:,1);



plot(r_E(:,1),r_E(:,2),r_d(:,1),r_d(:,2),r_cb(:,1),r_cb(:,2), 0,0, 'o', 'MarkerSize',10,...
    'MarkerEdgeColor','y',...
    'MarkerFaceColor',[0.9290 0.6940 0.1250])
xlabel('X [km]')
ylabel('Y [km]')
legend('Earth','Didymos Barycenter','CB21');
grid on;
%axis equal;
axis square;

% 
% plot3(r_E(:,1),r_E(:,2),r_E(:,3),r_d(:,1),r_d(:,2),r_d(:,3),r_cb(:,1),r_cb(:,2),r_cb(:,3), ...
%     r_E(1,1),r_E(1,2),r_E(1,3), '>', r_d(1,1),r_d(1,2),r_d(1,3), '>', r_cb(1,1),r_cb(1,2),r_cb(1,3), '>', [0,0], [0,0] , [0,0], '*', ...
%     r_E(end,1),r_E(end,2),r_E(end,3), 's', r_d(end,1),r_d(end,2),r_d(end,3), 's', r_cb(end,1),r_cb(end,2),r_cb(end,3), 's')
% xlabel('X [km]')
% ylabel('Y [km]')
% zlabel('Z [km]')
% legend('Earth','Didymos Barycenter','CB21');
% grid on;
% %axis equal;
% axis tight;


%% NEEDS TO BE MODIFIED NOW WOOOOHOOO!!!

GM_central = mu_S;

% meas_flag = 0;% measure dV min
% meas_flag = 1;% measure mass payload max
meas_flag = 2; %measure mass payload max witout dV arrival calculation



% Minimum departure date (MJD).
min_dep_date = MJD_E(1); %Earlies dep date: 2020-Jan-01
% Step size for scanning departure dates (days).
dep_step = 1;
% Maximum departure date (MJD).
max_dep_date = MJD_E(1+365*2); %Latest dep date: 2022-Jan-01

%store initial and final possible values in mdys
[min_dep_ymd(1), min_dep_ymd(2),min_dep_ymd(3), hrs, mins, secs] = JD_to_mdys(min_dep_date+ 2400000.5);
[max_dep_ymd(1), max_dep_ymd(2),max_dep_ymd(3), hrs, mins, secs] = JD_to_mdys(max_dep_date + 2400000.5);

dep_array = min_dep_date:dep_step:max_dep_date;
%len_dep = abs(min_dep_date-max_dep_date)
n_dep = size(dep_array,2);



% Minimum Time Of Flight (TOF) (days).
min_tof = 100;
% Step size for scanning TOFs (days).
tof_step = 1;
% Maximum TOF (days).
max_tof = MJD_E(end) - MJD_E(1); %1461
%tof array
tof_array = min_tof:tof_step:max_tof;
n_tof = size(tof_array,2);

n_iterations = n_dep*n_tof;

fprintf('Exploring %d date/TOF pairs.\n',n_iterations);

% Maximum arrival date (MJD).
% max_arr_date =R_PDC(end,1)
max_arr_date = 59477; %max arrival date: 2021-09-20
[max_arr_ymd(1), max_arr_ymd(2),max_arr_ymd(3), hrs, mins, secs] = JD_to_mdys(max_arr_date + 2400000.5);

fprintf('Earliest Departure Date: %d-%d-%d [YYYY-MM-DD].\n',min_dep_ymd(1),min_dep_ymd(2),min_dep_ymd(3));
fprintf('Latest Departure Date: %d-%d-%d [YYYY-MM-DD].\n',max_dep_ymd(1),max_dep_ymd(2),max_dep_ymd(3));

fprintf('Latest Arrival Date: %d-%d-%d [YYYY-MM-DD].\n',max_arr_ymd(1), max_arr_ymd(2),max_arr_ymd(3));


nrevmax =2; %number of maximum revolutions
fprintf('%d revolutions allowed by lambert solver.\n',nrevmax);



if meas_flag == 0
    z_lim = 12; %km/s max dV allowed
    disp('Optimizing minimum dV solutions.');
elseif meas_flag ==1
    z_lim = inf; %kg max mass allowed
    disp('Optimizing maximum delivered mass solutions.');
else
    z_lim = inf; %kg max mass allowed
    disp('Optimizing maximum delivered mass solutions (for rendezvous, excluding dV_Arr).');
end


% max_arr_angle = 120; %[deg] arrival phase angle
% min_SES_angle = 3; %[deg] sun Earth Spacecraft angle at asteroid rondezvous
max_arr_angle = inf; %[deg] arrival phase angle
min_SES_angle = -inf; %[deg] sun Earth Spacecraft angle at asteroid rondezvous

fprintf('Maximum Arrival Angle: %f [deg].\n',max_arr_angle);
fprintf('Minimum Sun Earth Spacecraft Angle: %f [deg].\n',min_SES_angle);

% max_DLA = 28.5; %deg
% min_DLA = -28.5; %deg

max_DLA = inf; %deg
min_DLA = -inf; %deg
fprintf('Maximum DLA: %f [deg].\n',max_DLA);
fprintf('Minimum DLA: %f [deg].\n',min_DLA);


sign_tof = [1 -1];
m = [0 0];

%{
 Lambert_tgt: By default, the short-way solution is computed. The long way solution
 may be requested by giving a negative value to the corresponding
 time-of-flight [tf].

 For problems with |m| > 0, there are generally two solutions. By
 default, the right branch solution will be returned. The left branch
 may be requested by giving a negative value to the corresponding
 number of complete revolutions [m].
%}

sign_tof(1) = 1;
m(1) = 0;
sign_tof(2) = -1;
m(2) = 0;
f = 3;

if nrevmax > 0
    
    for inc = 1:nrevmax
        
        %increment through short way solution, long way solution, and right and left branch solutions
        sign_tof(f) = 1;
        m(f) = inc;
        f = f + 1;
        
        sign_tof(f) = -1;
        m(f) = inc;
        f = f + 1;
        
        sign_tof(f) = 1;
        m(f) = -inc;
        f = f + 1;
        
        sign_tof(f) = -1;
        m(f) = -inc;
        f = f + 1;
        
    end
    
end

len_sign = length(sign_tof);


%initialize dV_total DLA, dv_dep, dv_arr, c3, DLA, phi, SES, dV_t, TOF, m_p
dV_total = NaN(n_dep, n_tof);
DLA_dep_store = NaN(n_dep, n_tof);
dv_arr_mag_store = NaN(n_dep, n_tof);
C3_dep_store = NaN(n_dep, n_tof);
dV_dep_store = NaN(n_dep, n_tof);
dep_date_store = NaN(n_dep, n_tof);
DLA_dep = zeros(len_sign,1);
phi_a = zeros(len_sign,1);
SES = zeros(len_sign,1);
phi_a_store = NaN(n_dep, n_tof);
SES_store = NaN(n_dep, n_tof);
dV_t_store = NaN(n_dep, n_tof);
TOF_store = NaN(n_dep, n_tof);
m_p_store = NaN(n_dep, n_tof);

%determine spacecraft Initial conditions in CPO
h_cpo = 185; %km
r_cpo = r_E_eq + h_cpo;
ISP = 310; %s
%v_sc_E_CPO = sqrt(mu_E/r_cpo);

disp('     ')
disp('-- Interplanetary Exploration --')
disp('Beginning grid search. Please standby...')
%loop through date indexs
%do both targets need to be in the same loop?
%loop through date indexs
for i = 1:n_dep
    
    %     turn the MJD into an index
    %     grab 1st departure date
    
    [val,dep_index] = ismember(dep_array(i),MJD_cb(:,1));
    
    %Grab current departure date and corresponding state vectors
    dep_date = MJD_E(dep_index,1);
    r_E_dep = r_E(dep_index,:);
    v_E_dep = v_E(dep_index,:);
    
    %loop through the flight times to determine corresp. arrival date
    for j = 1:n_tof
        
        arr_date = dep_date + tof_array(j);
        
        %if within max arrival date, get PDC state at arrival
        %otherwise skip solution and move to next arr date
        %if the max arrival been exceeded, break TOF and move
        %to the next dep date
        
        if arr_date <= max_arr_date
            
            %find the index that corresponds to the arrival date
            
            [val,arr_index] = ismember(arr_date,MJD_cb(:,1));
            %val just tells 1, is member, 0 is not member
            
            %get state vectros at arrival
            r_cb_arr = r_cb(arr_index,:);
            v_cb_arr = v_cb(arr_index,:);
            
            
            %init min dV for this grid point
            dV_min = inf;
            m_max = -inf;
            %loop through short way and long way, left and right
            %solutions for each dep arrival pair'
            
            for k = 1:len_sign
                %loop right/left branch solutions
                num_orb = m(k);
                %loop through short/long solutions
                tf = sign_tof(k)*abs(arr_date-dep_date); %days
                
                %%%%% INPUTS %%%%%
                %   r1, r1       [km]     position vectors of the two terminal points.
                %     tf        [days]    time of flight to solve for
                %      m          [-]     specifies the number of complete orbits to complete
                %                         (should be an integer)
                % GM_central   [km3/s2]   std. grav. parameter (GM = mu) of the central body
                [vi, vf, extremal_distances, exitflag] = lambert_tgt(r_E_dep, r_cb_arr, tf, num_orb, GM_central);
                %%%%% OUTPUTS %%%%%
                %vi, vf              [km/s] terminal velocities at the end-points
                %extremal_distances  [km]   minimum(1) and maximum(2) distance of the
                %                             spacecraft to the central body.
                
                
                %heliocentric velocity of s/c at Earth SOI (dep Vinf HCI)
                v_inf_E_dep = vi - v_E_dep;
                
                %Convert HCI to ECI
                v_inf_ECI =  T_HCI_ECI * v_inf_E_dep';
                v_inf_ECI_mag = norm(v_inf_ECI);
                
                %calculate dep DLA and only accept trajectory solutions that are within
                DLA_dep(k) = asind(v_inf_ECI(3)/v_inf_ECI_mag); %deg
                
                if (DLA_dep(k) <= max_DLA) && (DLA_dep(k) >= min_DLA)
                    
                    v_inf_E_dep_mag = norm(v_inf_E_dep);
                    C3_dep = v_inf_E_dep_mag^2;
                    
                    %now we calculate the approach phase angle and see if we are
                    %within spec
                    %calculate approach phase angle of june
                    r_cb_arr_hat = r_cb_arr/norm(r_cb_arr);
                    
                    %get the state befor arrival of the spacecraft
                    %v_sc_arr = vf;
                    
                    v_rel = vf- v_cb_arr;
                    v_rel_mag = norm(v_rel);
                    v_rel_hat = v_rel/v_rel_mag;
                    
                    phi_a(k) = acosd(dot(v_rel_hat,r_cb_arr_hat)); %DEG
                    
                    if phi_a(k) <= max_arr_angle
                        
                        %sun to earth
                        r_S_E_arr = -r_E(arr_index,:);
                        %spacecraft to earth
                        r_SC_E = r_cb_arr + r_S_E_arr;
                        %sun earth spacecraft angle at arrival
                        SES(k) = acosd( dot(r_S_E_arr,r_SC_E)/(norm(r_S_E_arr)*norm(r_SC_E)));
                        
                        if SES(k) >= min_SES_angle
                            %now that you have met all the criteria, get
                            %dV & mass
                            %for this solution gridpoint
                            
                            dV_dep = sqrt(C3_dep + 2*mu_E/r_cpo) - sqrt(mu_E/r_cpo); %km/s
                            
                            if meas_flag ==2
                                dv_cb_arr_mag = 0;
                            else
                                dv_cb_arr_mag = norm( vf - v_PDC_arr);
                            end
                            
                            dV_t = dV_dep + dv_cb_arr_mag;
                            
                            
                            if meas_flag == 0 %search for dV min
                                
                                %check if this left/right short/long solution is
                                %the dV min at this gridpoint,
                                
                                if dV_min > dV_t
                                    %store new min dV
                                    dV_total(i,j) = dV_t;
                                    
                                    %now we override the lowest dV value for
                                    %this gridpoint sol for future comparisons
                                    dV_min = dV_t;
                                    dV_t_store(i,j) = dV_min;
                                    %store DLA, dv_dep, dv_arr, c3, phi, ses, dep date for this
                                    %gridpoint solution.
                                    dep_date_store(i,j) = dep_date;
                                    dv_arr_mag_store(i,j) = dv_cb_arr_mag;
                                    C3_dep_store(i,j) = C3_dep;
                                    dV_dep_store(i,j) = dV_dep;
                                    
                                    phi_a_store(i,j) = phi_a(k);
                                    SES_store(i,j) = SES(k);
                                    DLA_dep_store(i,j) = DLA_dep(k);
                                    TOF_store(i,j) = tof_array(j);
                                end
                                
                                %check if this left/right short/long solution is
                                %the max mass at this gridpoint,
                            elseif meas_flag ==1
                                
                                [m_C3] = CCAFS_Intermediate_LV_mass_C3_DLA (C3_dep, DLA_dep(k));
                                %if NaN, handled by PCC_plot


                                m_p = m_C3/exp(dv_cb_arr_mag/(g*ISP/1000));
                                
                                %if current mass solution is larger than
                                %existing max
                                if m_max < m_p
                                    %store new max mass
                                    m_p_store(i,j) = m_p;
                                    %now we override the largest mass value for
                                    %this gridpoint sol for future comparisons
                                    m_max = m_p;
                                    
                                    %store DLA, dv_dep, dv_arr, c3, phi, ses, dep date for this
                                    %gridpoint solution.
                                    dV_total(i,j) = dV_t;
                                    dep_date_store(i,j) = dep_date;
                                    dv_arr_mag_store(i,j) = dv_cb_arr_mag;
                                    C3_dep_store(i,j) = C3_dep;
                                    dV_dep_store(i,j) = dV_dep;
                                    phi_a_store(i,j) = phi_a(k);
                                    SES_store(i,j) = SES(k);
                                    DLA_dep_store(i,j) = DLA_dep(k);
                                    TOF_store(i,j) = tof_array(j);
                                end
                                
                            else
                                
                                [m_C3] = CCAFS_Intermediate_LV_mass_C3_DLA (C3_dep, DLA_dep(k));
                                %if NaN, handled by PCC_plot
                                m_p = m_C3;
                                
                                %if current mass solution is larger than
                                %existing max at this grid point
                                if m_max < m_p 
                                    %store new max mass
                                    m_p_store(i,j) = m_p;
                                    
                                    %now we write over the largest mass value for
                                    %this gridpoint sol for future comparisons
                                    m_max = m_p;
                                    
                                    %store DLA, dv_dep, dv_arr, c3, phi, ses, dep date for this
                                    %gridpoint solution.
                                    dV_total(i,j) = dV_t;
                                    dep_date_store(i,j) = dep_date;
                                    %store relative arrival speed instead
                                    %of dV arrival
                                    dv_arr_mag_store(i,j) = v_rel_mag;
                                    C3_dep_store(i,j) = C3_dep;
                                    dV_dep_store(i,j) = dV_dep;
                                    phi_a_store(i,j) = phi_a(k);
                                    SES_store(i,j) = SES(k);
                                    DLA_dep_store(i,j) = DLA_dep(k);
                                    TOF_store(i,j) = tof_array(j);
                                end
                            end %end setting dV/mass gridpoint solution % storing data
                        end %end SSE angle loop
                    end %end arrival phase loop
                end %end DLA loop
            end %end sign loop
        end %end arrival date check loop
    end %end of flight times loop
end%end of departure dates loop



if meas_flag == 0
    
    %Pork Chop Plot dV
    [dV_t_sol,i_ext,j_ext] = pcc_plot(dep_array, tof_array, dV_total, z_lim, meas_flag, TOF_store);
    disp(' ')
    disp('-- Minimum dV solution --')
else
    %Pork Chop Plot kg
    [mp_max,i_ext,j_ext] = pcc_plot(dep_array, tof_array, m_p_store, z_lim, meas_flag, TOF_store);
    m_p_sol = m_p_store(i_ext,j_ext);
    dV_t_sol = dV_total(i_ext,j_ext);
    disp(' ')
    disp('-- Maximum delivered mass solution --')
    fprintf('Delivered Spacecraft Mass: %f [kg].\n',m_p_sol);
end

% 
% [min_row, ind_row]   = min(TOF_store);
% [TOF_min, ind_col] = min(min_row)





dep_date_sol = dep_array(i_ext);
[dep_date_sol_mdys(1), dep_date_sol_mdys(2), dep_date_sol_mdys(3), ~, ~, ~] = ...
    JD_to_mdys(dep_date_sol + 2400000.5);

tof_sol = tof_array(j_ext);
arr_date_sol = dep_date_sol + tof_sol;

[arr_date_sol_mdys(1), arr_date_sol_mdys(2),arr_date_sol_mdys(3), ~, ~, ~] = ...
    JD_to_mdys(arr_date_sol + 2400000.5);

%store DLA, dv_dep, dv_arr, c3 for this solution
dV_dep_sol = dV_dep_store(i_ext,j_ext);
C3_dep_sol = C3_dep_store(i_ext,j_ext);
DLA_sol = DLA_dep_store(i_ext,j_ext);
dv_arr_mag_sol = dv_arr_mag_store(i_ext,j_ext);
phi_a_sol = phi_a_store(i_ext,j_ext);
SES_sol = SES_store(i_ext,j_ext);

%print results
fprintf('Earth Departure Date: %d-%d-%d [YYYY-MM-DD].\n',dep_date_sol_mdys(1), dep_date_sol_mdys(2),dep_date_sol_mdys(3));
fprintf('Time of flight:  %d [days].\n',tof_sol);
fprintf('PDC Arrival Date: %d-%d-%d [YYYY-MM-DD].\n',arr_date_sol_mdys(1), arr_date_sol_mdys(2),arr_date_sol_mdys(3));
fprintf('Earth departure dV: %f [km/s].\n',dV_dep_sol);
fprintf('Earth departure C3: %f [km^2/s^2].\n',C3_dep_sol);
fprintf('Earth departure DLA: %f [deg].\n',DLA_sol);
fprintf('PDC arrival dV: %f [km/s].\n',dv_arr_mag_sol);
fprintf('Total mission dV: %f [km/s].\n',dV_t_sol);
fprintf('Sun-Earth-Spacecraft angle at PDC arrival: %f [deg].\n',SES_sol);
fprintf('PDC arrival approach phase angle: %f [deg].\n',phi_a_sol);


% T = table(dV_total,duration_total,tf,24,t_days_arr, datetime(2024,12,07,0,00,00),C3_dep,v_inf_E_dep_mag,deltaV_dep,dv_PDC_arr_mag , dv_NEA_dep_mag, 0, v_entry,DLA_dep, DLA_arr)
%
% writetable(T ,'HW4.xlsx','Sheet',2)


% %THESE ARE IN THE HELIOCENTRIC INERTIAL FRAME
% v_inf_E_arr =  vf - v_E_arr  ;
% v_inf_E_arr_mag = norm(v_inf_E_arr);
%
% %CONVERT HCI TO ECI
% v_inf_ECI =  T_HCI_ECI * v_inf_E_arr;
% v_inf_ECI_mag = norm(v_inf_ECI);
%
% DLA_arr = asind(v_inf_ECI(3)/v_inf_ECI_mag); %deg
%
%
% v_entry = norm(sqrt(v_inf_ECI_mag^2 + (2*mu_E)/r_entry));
% %so long as v_entry is less than 12.9 km/s (for non human Earth atmospher
% %entry), there is no dV required.


