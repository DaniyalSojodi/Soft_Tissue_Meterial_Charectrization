function BMG5300_Project()

% The material constants are determined by nonlinear optimization. 
% Once the constants are found, the inner radius, the longitudinal
% stretch, the pressure and the longitudinal force are double-checked.

    clc;
    clear % clears all variables
    format long
    
    % Reading experimental data file
    experimentalDataFileName = fopen('referencerawdata2.txt', 'rt');
    fileContentByLine = textscan(experimentalDataFileName, '%f %f %f',1);
    Ro = fileContentByLine{1}(1);           % stress-free outer radius (mm)
    Ri = fileContentByLine{2}(1);           % stress-free outer radius (mm)
    theta0 = fileContentByLine{3}(1);       % opening angle (deg)
    theta0 = theta0*pi/180;                 % opening angle (rad)
    
    fileContentByLine = textscan(experimentalDataFileName, '%f %f %f %f');
    fclose(experimentalDataFileName);
    
    % loading experimental data into variables
    Pexp = fileContentByLine{1};            % luminal pressure (eg. MPa or kPa)
    riexp = fileContentByLine{2};           % inner radius under Pexp (mm)
    lambdaexp = fileContentByLine{3};       % longitudinal stretch ratio under Pexp (-)
    Fexp = fileContentByLine{4};            % longitudinal force under Pexp (e.g. N or mN)
    
    % options for Levenberg-Marquardt method of solution
    options = optimset('largescale','off','MaxFunEvals',1e100,...
                       'tolFun',1e-30,'TolX',1e-30,'MaxIter',5e3,...
                       'Algorithm','levenberg-marquardt','Display','on','Display','iter','PlotFcns', @optimplotfval);


%     ksi = [-sqrt(3/5),0,sqrt(3/5)];
%     w = [5/9,8/9,5/9];

    % abscissas for Gauss-Legendre integration
    ksi = [-0.238619186083197,-0.661209386466265,-0.932469514203152,...
            0.238619186083197,0.661209386466265,0.932469514203152];
    % weights for Gauss-Legendre integration
    w = [0.467913934572691,0.360761573048139,0.171324492379170,....
         0.467913934572691,0.360761573048139,0.171324492379170];
     
    model = 'fung';

    % initial guess for material constants
    iniconstants = [0.5, 0.15, 0.4, 0.9]; 
    RadLambexp = [riexp,lambdaexp]; 
    
    [Constants, fval] = lsqnonlin(@(x) errPF(x, RadLambexp, Pexp, Fexp, Ro, Ri, theta0, ksi, w, model), iniconstants, [], [], options);
    disp(['Final value of the objective function: ', num2str(fval)]);
    fprintf('\n');
    fprintf('FUNG MODEL CONSTANTS -- Black markers\n');
    fprintf('c1 = %f [same unit as Pexp]\n',Constants(1));
    fprintf('c2 = %f \n',Constants(2));
    fprintf('c3 = %f \n',Constants(3));
    fprintf('c4 = %f \n',Constants(4));

    % initial guesses for theoretical inner radius and longitudinal stretch
    iniguess = [riexp,lambdaexp]; 
    RadLamb = lsqnonlin(@(x) errPF(Constants,x,Pexp,Fexp,Ro,Ri,...
                        theta0,ksi,w,model),iniguess,[],[],options); 
    % use optimized values to get theoretical pressure and force
    [Pressure,Force] = PF(Constants,RadLamb,Ro,Ri,theta0,ksi,w,model); 

    RadLambFung = RadLamb;
    PressureFung = Pressure;
    ForceFung = Force;

    model = 'guccione';

    % initial guess for material constants
    iniconstants = [0.076,0.5,0.6]; 
    RadLambexp = [riexp,lambdaexp];     
    
    %Constants = lsqnonlin(@(x) errPF(x,RadLambexp,Pexp,Fexp,Ro,Ri,...
        %theta0,ksi,w,model),iniconstants,[],[],options);
    %fprintf('\n');
    %fprintf('GUCCIONE MODEL CONSTANTS -- Green markers\n');
    %fprintf('c1 = %f [same unit as Pexp]\n',Constants(1));
    %fprintf('c2 = %f \n',Constants(2));
    %fprintf('c3 = %f \n',Constants(3));

    % initial guess for theoretical inner radius and longitudinal stretch
    iniguess = [riexp,lambdaexp]; 
    RadLamb = lsqnonlin(@(x) errPF(iniconstants,x,Pexp,Fexp,Ro,Ri,...
                        theta0,ksi,w,model),iniguess,[],[],options);
    % use optimized values to get theoretical pressure and force
    [Pressure,Force] = PF(iniconstants,RadLamb,Ro,Ri,theta0,ksi,w,model);

    RadLambGuc = RadLamb;
    PressureGuc = Pressure;
    ForceGuc = Force;
    
    model = 'hgo';

    % initial guess for material constants
    iniconstants = [0.01,0.1,0.1,1*pi/180];
    
    Constants = lsqnonlin(@(x) errPF(x,RadLambexp,Pexp,Fexp,Ro,Ri,...
        theta0,ksi,w,model),iniconstants,[],[],options); 
    fprintf('\n');
    fprintf('HGO MODEL CONSTANTS -- Blue markers\n');
    fprintf('c1 = %f [same unit as Pexp]\n',Constants(1));
    fprintf('c2 = %f [same unit as Pexp]\n',Constants(2));
    fprintf('c3 = %f \n',Constants(3));
    fprintf('alpha (degrees) = %f \n',abs(mod((Constants(4)),pi))*180/pi);

    % initial guess for theoretical inner radius and longitudinal stretch
    iniguess = [riexp,lambdaexp];
    RadLamb = lsqnonlin(@(x) errPF(Constants,x,Pexp,Fexp,Ro,Ri,...
                        theta0,ksi,w,model),iniguess,[],[],options);
    % use optimized values to get theoretical pressure and force    
    [Pressure,Force] = PF(Constants,RadLamb,Ro,Ri,theta0,ksi,w,model); 

    RadLambHGO = RadLamb;
    PressureHGO = Pressure;
    ForceHGO = Force;

    % plotting figures
    figure(1)
    subplot(2,1,1)
    plot(Pexp*1000,riexp,'-r','linewidth',2); hold on;
    plot(PressureFung*1000,RadLambFung(:,1),'+k',...
        PressureGuc*1000,RadLambGuc(:,1),'+g',...
        PressureHGO*1000,RadLambHGO(:,1),'+b')
    ylabel('Inner Radius (mm)')

    subplot(2,1,2)
    plot(Pexp*1000,lambdaexp,'-r','linewidth',2); hold on;
    plot(PressureFung*1000,RadLambFung(:,2),'+k',...
        PressureGuc*1000,RadLambGuc(:,2),'+g',...
        PressureHGO*1000,RadLambHGO(:,2),'+b')
    ylabel('Longitudinal Stretch Ratio (-)')
    xlabel('Inner Pressure (kPa)')    
    
end


function error = errPF(C,RadLamb,Pexp,Fexp,Ro,Ri,theta0,ksi,w,model)
    [Pt,Ft] = PF(C,RadLamb,Ro,Ri,theta0,ksi,w,model);
    if strcmp(model,'fung')
        errp=(Pt-Pexp)*1; % use 1, 10, 100, 1000 as weight
        errf=(Ft-Fexp)*100;    
    elseif strcmp(model,'guccione')
        errp=(Pt-Pexp)*1; % use 1, 10, 100, 1000 as weight
        errf=(Ft-Fexp)*1;
    elseif strcmp(model,'hgo')
        errp=(Pt-Pexp)*100; % use 1, 10, 100, 1000 as weight
        errf=(Ft-Fexp)*1;
    end
    error=[errp;errf]; % combines errp and errf into one vector
end

function [Pt,Ft] = PF(C,RadLamb,Ro,Ri,theta0,ksi,w,model)
    ri=RadLamb(:,1);
    lambda=RadLamb(:,2);
    P=zeros(length(ri),1); % initial pressure
    F=zeros(length(ri),1); % initial force
    % outer radius in current state
    ro=sqrt(((Ro^2-Ri^2)./(pi*lambda))*theta0+ri.^2); 
    for i=1:length(ksi) % Gauss-Legendre points
            r=(ro+ri)./2+(ro-ri)./2.*ksi(i); % current radius  
            R=sqrt(Ri.^2+(((pi.*lambda)./theta0).*(r.^2-ri.^2)));     
            Ezz=0.5*(lambda.^2-1);      
            Ett=0.5*(((pi.*r)./(theta0.*R)).^2-1);       
            calcPartials = str2func(model);
            % calculate partial derivatives for current model            
            [dWtt,dWzz] = calcPartials(C,Ett,Ezz); 
    P=P+w(i).*(((pi.*r)./(theta0.*R)).^2.*dWtt)./r;   
    F=F+w(i).*(2.*lambda.^2.*dWzz-((pi.*r)./(theta0.*R)).^2.*dWtt).*r;
    end
    Pt=((ro-ri)./2).*P; % theoretical pressure
    Ft=pi.*((ro-ri)./2).*F; % theoretical force
end

function [dWtt,dWzz] = fung(C,Ett,Ezz)
% Computes the partial derivatives of W for the Fung model
    dWtt=C(1)*(C(2).*Ett+C(4).*Ezz).*...
         exp((C(2).*Ett.^2)+(C(3).*Ezz.^2)+(2*C(4).*Ett.*Ezz));
    dWzz=C(1)*(C(3).*Ezz+C(4).*Ett).*...
         exp((C(2).*Ett.^2)+(C(3).*Ezz.^2)+(2*C(4).*Ett.*Ezz));       
end

function [dWtt,dWzz] = guccione(C,Ett,Ezz)
% Computes the partial derivatives of W for the Guccione model
    delta=(1+2.*Ett).*(1+2.*Ezz);
    Q=C(2).*Ett.^2+C(3).*(Ezz.^2+(0.5.*(1./delta-1)).^2);  
    dWtt=2.*C(2).*Ett-C(3).*(1./delta-1)./(delta.*(1+2.*Ett));
    dWtt=0.5.*C(1).*dWtt.*exp(Q);
    dWzz=2.*C(3).*Ezz-C(3).*(1./delta-1)./(delta.*(1+2.*Ezz));
    dWzz=0.5.*C(1).*dWzz.*exp(Q); 
end

function [dWtt,dWzz] = hgo(C,Ett,Ezz)
% Computes the partial derivatives of W for the HGO model
    delta=(1+2.*Ett).*(1+2.*Ezz);
    alpha=C(4);
    I4=(cos(alpha))^2.*(2.*Ett+1)+(sin(alpha))^2.*(2.*Ezz+1);
    Q1=(I4-1).^2;
    dWtt=C(1)*(1-1./(delta.*(1+2.*Ett)))+...
        4*C(2)*(cos(alpha))^2.*((I4-1).*exp(C(3).*Q1));
    dWzz=C(1)*(1-1./(delta.*(1+2.*Ezz)))+...
        4*C(2)*(sin(alpha))^2.*((I4-1).*exp(C(3).*Q1));
end
