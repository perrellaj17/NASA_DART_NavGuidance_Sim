function [mass] = CCAFS_Intermediate_LV_mass_C3_DLA (C3, DLA)

C3_row = 10:10:100;
DLA_col = [28.5 38 48 57 90];
DLA = abs(DLA); %need absolute value of DLA



row_len = length(C3_row);
col_len = length(DLA_col);


m = [4230 4115 3975 3840 3010;
    3400 3305 3190 3080 2395;
    2735 2660 2570 2475 1920;
    2210 2150 2075 2000 1540;
    1790 1740 1680 1620 1245;
    1455 1415 1365 1315 1005;
    1190 1150 1110 1070 815;
    970 940 905 870 655;
    790 765 735 710 540;
    645 625 600 575 430];

%check if C3 and DLA are in range of data 

if ((C3 < 10 || C3 > 100) || (DLA < 0 || DLA > 90))
    mass = NaN;
else
    
    %if DLA>57 deg, then take DLA = 90 value
    %if DLA<=28.5, set dLA = 28.5 value

    if DLA <= 28.5
        DLA = 28.5;
    elseif DLA > 57
        DLA = 90;
    end
    
    
    %firct, cycle through the C3 rows
    for i = 1:(row_len-1)
        %account for inclusivity of values at interpolation points
        %check which row c3 is in
        if (C3 <= C3_row(i+1) &&  C3 > C3_row(i)) || (C3 < C3_row(i+1) &&  C3 >= C3_row(i))
            %we are between row i and i+1
            row_l =  C3_row(i);
            row_h =  C3_row(i+1);
            row_index_l = i;
            row_index_h = i+1;
            break;
        end
    end
    
    %now, cycle through the DLA columns
    for j = 1:(col_len-1) 
        %account for inclusivity of values at interpolation points
        %check which column the DLA is in
        if (DLA <= DLA_col(j+1) &&  DLA > DLA_col(j)) || (DLA < DLA_col(j+1) &&  DLA >= DLA_col(j))
            %we are between column j and j+1
            col_l =  DLA_col(j);
            col_h =  DLA_col(j+1);
            col_index_l = j;
            col_index_h = j+1;
            break;
        end
    end
    
    
    %interpolate using a linear system to solve for coefficients
    
    matrix = [ 1 row_l col_l row_l*col_l;
        1 row_l col_h row_l*col_h;
        1 row_h col_l row_h*col_l;
        1 row_h col_h row_h*col_h];
    
    mass_data = [ m(row_index_l, col_index_l);
        m(row_index_l, col_index_h);
        m(row_index_h, col_index_l);
        m(row_index_h, col_index_h)];
    
    coeff = inv(matrix)*mass_data;
    
    mass = coeff(1) + coeff(2)*C3 + coeff(3)*DLA + coeff(4)*C3*DLA;
    
    
end

end




