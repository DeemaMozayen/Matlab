
% Input weather data 
data = [0.25 20 3.1 6 17; 0.3 33 3.1 5.7 47; 0.32 62 5.2 6 75; 0.41 103 7.6 7 105; 0.41 134 10.6 9 135; 0.39 136 14 12 162; 0.4 137 15.8 14.5 198; 0.42 120 15.4 15 228; 0.39 80 13.2 14.2 258; 0.33 46 10 12 288; 0.32 15 6 9.2 318; 0.27 15 4.2 7 344];
% Days in each month
days = [31 28 31 30 31 30 31 31 30 31 30 31];
% Further parameters
rho = 0.2; phi = 52.3; beta = 45; eta = 0.8075;

% DeWinter Exchange Factor 
A_c = 6;            % Collector Area 
F_r = 0.9;          % Heat Removal Factor 
U_L = 3.5;          % Collector Heat Loss Coeff 
m_c = 280;          % Collector Side [FlowRate*SpecificHeat] 
epsilon = 0.8;      % Heat Exchanger Effectiveness  
m_min = 252;        % Store Side [FlowRate*SpecificHeat] 
Fr_prime = (F_r)/( 1 + ((F_r*U_L*A_c)/(m_c))*(((m_c)/((epsilon*m_min)))-1));

for n = 1
    % Days in each month
    N=days(n);
    % Number of seconds in a given month 
    delta_t = 3600*24*N; 
    % Data Column = Monthly Clearness Indexes
    K = data(n,1);
    % Data Column = Monthly Insolation
    H = data(n,2)*3.6*(10)^6/N;
    % Data Column = Ambient Temperature
    Ta = data(n,3);
    % Data Column = Mains Temperature
    Tm = data(n,4);
    % Data Column = Average Day Numbers    
    Avd = data(n,5);     
    
    % Monthly Heating Load Calculation 
    c = 4200;       % Specific Heat of Water raised
    m = 200*N;      % Monthly Mass of Water Heated 
    QLM = m*c*(60-Tm);
    
    % X calculation
    X1 = (A_c*Fr_prime*U_L*(100-Ta)*delta_t)/QLM;           % Original X
    XW = X1*(11.6 + 1.18*60 + 3.86*Tm -2.32*Ta)/(100-Ta);   % Correction 1
    X = XW*(((300)/(75*A_c))^(-0.25));                      % Correction 2
    
    % Declination Angle
    delta = 23.45*sind(360*((284+Avd)/365));
    % Sunset Angles
    omega_s = acosd(-tand(delta)*tand(phi));
    omega_ss = min (omega_s, acosd(-tand(phi-beta)*tand(delta)));
    % Diffuse Radiation 
    if omega_s <= 81.4
        Hd = (1.391 - 3.560*K + 4.189*(K^2) - 2.137*(K^3));
    else
        Hd = (1.311 - 3.022*K + 3.427*(K^2) - 1.821*(K^3));
    end
    % Ratio of Inclined to Horizontal Radiation
    Rbt= (((cosd(phi-beta)*cosd(delta)*sind(omega_ss))+((pi/180)*omega_ss*sind(phi-beta)*sind(delta)))/((cosd(phi)*cosd(delta)*sind(omega_s))+((pi/180)*omega_s*sind(phi)*sind(delta))));
    % Monthly Average Incident Radiation
    Ht = H*(((1 - Hd)*Rbt) + (Hd*((1+cosd(beta))/2)) + (rho*((1-cosd(beta))/2)));
    % Y calculation
    Y = (A_c*Fr_prime*eta*Ht*N)/QLM;
    
    % Fm calculation
    Fm = 1.029*Y - 0.065*X - 0.245*(Y)^2 + 0.0018*(X)^2 + 0.0215*(Y)^3;
   
    Fm_out(n)=Fm;
    QLM_out(n)=QLM/(3.6*10^6);
    Q_solar(n)=Fm*(QLM/(3.6*10^6));
    Results = [QLM_out;Fm_out;];
    
    Fm_mean = mean(Fm_out);
    Fm_min = min(Fm_out);

end

    % Cost Assessment
    
    % Yearly Heating Cost (Gas only)
    Cg = 0.06*(sum(QLM_out)/0.75);
    % Gas energy needed to supplement Collector
    Q_extragas = [343 245 205 98 66 55 50 58 105 199 330 352]/0.75;
    % Yearly cost of Gas needed (w/solar)
    Cgs = (sum(Q_extragas))*0.06;
    % Yearly bestowed RHI
    Cr = (sum(Q_solar))*0.2066;
    % Yearly net heating cost (w/solar)
    Cs = Cgs - Cr;
    % Yearly Saving (solar vs Gas only)
    S = Cg - Cs;
    % Payback time of system
    P = 4000/S;
    
