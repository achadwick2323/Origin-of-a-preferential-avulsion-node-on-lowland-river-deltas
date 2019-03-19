%% Dimensionless backwater morphodynamics over multiple avulsions
function [ Output ] = MultiDimless( varargin )
%  Austin Chadwick 3/18/19
%  Contact: achadwick@caltech.edu
%  This code runs the script "Dimless" for multiple avulsion cycles.

%% Settings

    % Set output filename.
    savebin = 1; 
    folder = '/Users/austinchadwick/Documents/Research/Deltas/BackwaterModel/Output_MultiDimless';
    filepath = fullfile(folder, ['output_' datestr(datetime('now'),'yymmdd-HHMMSS')]);      
        if savebin; disp(filepath); end
    
    % Set whether the code will plot (plotbin = 1), or not (plotbin = 0).
    plotbin = 0;
    
    % Set max allowable of computational threads. 
    Nthreads = 1;
    maxNumCompThreads(Nthreads)
    
    % Set subroutine that runs individual avulsion cycles
    Dimless_fun = @(parameters,InitialConditions) ...
                  Dimless(parameters,InitialConditions); % morphodynamic model

    % Set model parameters
    N = 4; % [-] number of lobes
    A = 13; % [-] number of phases (i.e. avulsion node steps)
    parameters = ParameterBuilder_fun; % input parameters for Dimless_fun subroutine
    InitialConditions = 'default'; 

    % Set the steepest slope allowed during bed relaxation
	Sthresh = parameters.SedimentTransport.Sa; 
        
    % Initialize the output matrix
    data = cell(A,N); % for recording morphodynamics output
    cycle = zeros(A,N); % for recording the order of avulsion cycles
    
    % the initial conditions for the first row. 
    tempparameters = parameters;
    tempparameters.hfill = 0.001; 
    tempparameters.nu_CFL = 0.1; 
    temp = Dimless_fun(tempparameters,InitialConditions); 
        temp.InitialConditions.xm = temp.InitialConditions.xs; 
    InitialConditionsForFirstAvulsion = temp.InitialConditions; 
        Parent = InitialConditionsForFirstAvulsion;
        Concat = InitialConditionsForFirstAvulsion;
           
%% Loop through avulsion cycles

    for i = 1:A
        
        % Display
        disp(['--> cycle ' num2str(i) ' of ' num2str(A)])
        
        % Get previous & current parent channel
        Parent_prev = Parent;
        Parent = Concat; 
        if i == 1
            ParentColumn = 1; 
        elseif i~=1;
            ParentColumn = DaughterColumn; % previous daughter is now the parent
        end
        
        % Get daughter channel
            % Evaluate the distance to each shoreline
            xm = nan(1,N);
            for j = 1:N;
                if j ~= ParentColumn
                    ind = find(~cellfun(@isempty,data(:,j)),1,'last'); % potential daughter profile 
                    if isempty(ind); temp = InitialConditionsForFirstAvulsion; % if we haven't avulsed to this lobe yet, use initial conditions
                    else temp = data{ind,j}.Avulsion; end % otherwise, its in the data final conditions
                    xm(j) = FindMouth_fun(temp.xs*parameters.sigma,temp.eta,temp.xs,temp.xb,temp.etab,temp.etaF,Parent_prev.xio,parameters); % notice we use the sea level form the previous parent to recompute the river mouth on a lobe that has been drowned since last active
                    %xm(j) = temp.xm;                   
                end
            end
            % Take the lobe with minimum shoreline distance as new daughter channel.
            [~,jmin] = nanmin(xm);
            DaughterColumn = jmin; 
            ind = find(~cellfun(@isempty,data(:,jmin)),1,'last');
            if isempty(ind); Daughter = InitialConditionsForFirstAvulsion; % if we haven't avulsed to this lobe yet, use initial conditions
            else Daughter = data{ind,jmin}.Avulsion; end % otherwise, its in the data final conditions         

        % Set initial conditions
        InitialConditions = Concat; % Initial condition for the current parent is the previous concatenated profile ...
            % ... except the reference elevations, which is for the upcoming daughter channel
            InitialConditions.etaFref = Daughter.etaFref; 
            InitialConditions.etaref  = Daughter.etaref; 
            InitialConditions.hcref = Daughter.hcref; 
            % ... and except sea-level, which is that of the end of the previous parent
            InitialConditions.xio = Parent_prev.xio(end); 

        % Run morphodynamics to determine parent of this cycle
        disp('     --> running morphodynamics ...')
        cycle(i,ParentColumn) = 1; 
            parameters0 = parameters;
        data{i,ParentColumn} = Dimless_fun(parameters0,InitialConditions);                      
            % Quit if we have not avulsed
            if ~isfield(data{i,ParentColumn},'Avulsion'); disp('     --> (!) NO AVULSION, quitting multi dimless'); return; 
            else                               disp('     --> parent has avulsed. '); end;
            % Define parent as the bed @ avulsion
            Parent = data{i,ParentColumn}.Avulsion; 
            
        % Concatenate parent & daughter
        disp('     --> concatenating daughter & parent ...')
        Concat = BedConcat_fun(Daughter,Parent,Sthresh,plotbin); 
            % Relax the slope of the step in the concatenated profile
            disp('     --> relaxing concatenated profile ...')
            Concat.eta = StepRelaxer_fun(Concat.xs,parameters.sigma,Concat.eta, Sthresh,plotbin);
            Concat.etaF= StepRelaxer_fun(Concat.xs,parameters.sigma,Concat.etaF,Sthresh,plotbin); if plotbin; subplot(2,1,1); title('levee relaxation'); end; % relax levees    

        % save as we go
        if savebin
        Output.data = data;
        Output.cycle = cycle;
        save(filepath,'Output','-v7.3');
        disp(['     --> saving ' filepath])
        end
                        
    end % end loop over avulsions

end

%% Model parameters
function [parameters] = ParameterBuilder_fun()

    % Discretization of time & space
    N  = 200;	% [-] number of topset nodes
    sigma = linspace(0,1,N)'; % [-] dimensionless downstream distance
    nu_CFL = 0.001; % [-] CFL parameter. If this is too big, then errors in mass balance go up and instabilities may occur. 
    au = 0.75; % [-] upwind coefficient (1 for pure upwind, 0 for forward, 0.5 for central)              
    dt_max = 1e-3; % [-] max allowable timestep. Please adjust manually such that this does not exceed the flow duration parameter Te. 
    Mprint = 1e4; % [-] number of timesteps to save (assuming no avulsion triggered)
    Mtoprint = 1000; % [-] number of timesteps between prints.
    
    % Bankfull flow properties
    Cf_bf  = 0.01; % [-] friction factor for bankfull condition
    Frn_bf = 0.17; %  [-] bankfull normal-flow Froude number
    taun_bf = 1; % [-] bankfull normal flow shields number, tau*_n,bf = (hn_bf*Sn_bf) / (R*D) = hn_bf/(R*D) * Cf_bf * Frn_bf

    % Widths
    Width.FlowWidthRelation = 'Flare'; % 'Constant' or 'Flare'
    Width.DepoWidthRelation = 'Step'; % % 'Constant', 'Flare', or 'Step
        % For flare case
        Width.Flare.WidthToDepthRatio = 20; % [-] Ratio of width to depth in normal flow reac
        Width.Flare.theta_Depo = 30; % [deg] spreading angle of sediment @ mouth with flare relation % (50)
        Width.Flare.theta_Flow = 15; % 0.1; % [deg] spreading angle of water @ mouth with flare relation % (3)
        % For "Step" case, there is a step-change in depositional width at the river mouth.
        Width.Step.WidthAtMouth = 1; % [-] Depositional width downstream of mouth normalized by depositional width upstream of mouth. If equals 1, it is same as constant width case. 

    % Sediment transport
    SedimentTransport.tauc = 0.0495; % [-] shields number for threshold of motion
    SedimentTransport.Sa = 5; % [-] avalanche slope for shock capturing
    SedimentTransport.Sthresh = 5; % 0.5 * SedimentTransport.Sa; 
    SedimentTransport.Relation = 'EngelundHansen'; % 'WongParker', 'EngelundHansen'
       
    % Flow Variability
    FlowRegime.Scheme = 'MultiPhase';
        % sampling timescale
        FlowRegime.MultiPhase.Te = 1e-3; % [-] dimensionless flow duration 
        FlowRegime.MultiPhase.numbins = 20; % [-] number of discharge bins
        FlowRegime.MultiPhase.Cfi = Cf_bf*ones([FlowRegime.MultiPhase.numbins,1]); % [-] assume friction factor is the same for all discharges
        % the distribution itself
        FlowRegime.MultiPhase.F_bf = 0.05; % [-] bankfull exceedance probability
        FlowRegime.MultiPhase.psi_std = 0.5; % [-] std of log of stage height, for log-normal distribution
            FlowRegime.MultiPhase.CV = sqrt(exp(FlowRegime.MultiPhase.psi_std.^2)-1); % [-] coefficient of variation of stage height
            FlowRegime.MultiPhase.psi_mean = psi_mean_calculator_fun(FlowRegime.MultiPhase.psi_std,FlowRegime.MultiPhase.F_bf,0); % [-] mean of log of stage height, for log-normal distribution
        FlowRegime.MultiPhase.hn_max = 100; % [-] max allowable normal-flow depth 
   
    % Avulsion parameters
    hfill = 0.5; % [-] dimensionless avulsion threshold  
        upbin = 1; % Set whether we allow (upbin=0) or do not allow (upbin=1) avulsions upstream of initial shoreline. We never observed avulsions downstream of initial shoreline, so setting to 1 allows code to assume this & speed up. 
    
    % sea-level parameters
    dxibdt = 0; % Relative sea-level rise rate
    s = 0; % [-] normalized subsidence rate    

    % initial conditions
    InitialConditions.xs  = 6; % [-] initial foreset 
    InitialConditions.hs  = 1; % [-] initial river-mouth depth
    InitialConditions.xio = 2; % [-] initial sea level
    InitialConditions.etab = 0; % [-] initial elevation of delta toe
    InitialConditions.Sff = 0; % [-] slope of the false floor
    InitialConditions.L = 20; % [-] length of the whole domain
    
    % pack up
    parameters = ws2struct;

end
    
%% Concatenate bed and levees of parent onto the daughter upstream of the avulsion node
function [Daughter] = BedConcat_fun(Daughter,Parent,Sthresh,plotbin)

    % grid is same for both data
    sigma = linspace(0,1,length(Daughter.eta))';

    % daughter path
    xs1  = Daughter.xs; %  shoreline
    xb1  = Daughter.xb; %  toe
    eta1 = Daughter.eta; %  bed elev
    etaF1= Daughter.etaF; % floodplain elev
    
    % parent path
    xs2  = Parent.xs;
    xb2  = Parent.xb;
    eta2 = Parent.eta;
    etaF2= Parent.etaF; 
    xA =  Parent.xA; % avulsion location
    
    % create new profile thru concatenation
    xs3 = xs1; % shoreline is the same as initial
    xb3 = xb1; % toe is the same as initial
        % riverbed
        eta3 = nan(size(sigma)); % the bed ...
            % upstream of avulsion node, use profile @ avulsion
            ind = xs3*sigma <= xA;
            eta3(ind) = interp1(xs2*sigma, eta2, xs3*sigma(ind));
            % downstream of avulsion node, use very first profile
            ind = xs3*sigma > xA; 
            eta3(ind) = interp1(xs1*sigma, eta1, xs3*sigma(ind));
        % floodplain
        etaF3 = nan(size(sigma));
            % upstream of avulsion node, use parent profile
            ind = xs3*sigma <= xA;
            etaF3(ind) = interp1(xs2*sigma, etaF2, xs3*sigma(ind));
            % downstream of avulsion node, use daughter profile
            ind = xs3*sigma > xA; 
            etaF3(ind) = interp1(xs1*sigma, etaF1, xs3*sigma(ind));      

    % output is the concatenated profile
    Daughter_old = Daughter; % remember old
    Daughter.xs = xs3; 
    Daughter.xb = xb3; 
    Daughter.eta = eta3;
    Daughter.etaF= etaF3;

    % plot if desired.
    if plotbin
        figure(); title('bed concatenation'); hold on
        plot(xs1*sigma,eta1,':','color',0.7*ones(1,3));
        plot(xs2*sigma,eta2,'--','color',0.5*ones(1,3));
        plot(xs3*sigma,eta3,'-k');
        plot(xs1*sigma,etaF1,':g');
        plot(xs2*sigma,etaF2,'g--');
        plot(xs3*sigma,etaF3,'g-');
            xlabel('x^*'); ylabel('\eta^*');
            legend('daughter path','parent path','concatenated path','initial levee','parent levee','concatenated levee')
    end % end if plotbin
    

end

%% Step relaxer
function [eta, Srelax] = StepRelaxer_fun(xs,sigma,eta,Sthresh,plotbin)

        % Unrelaxed bed
        eta_old = eta;

        % Calculate slope & avulsion location 
        [S,Sstep,Sshore] = SlopeCalc_fun(xs,sigma,eta);
            [~,indA] = nanmax(S); 
            xA = xs*sigma(indA); 
          
        % Relaxation slope is the input threshold slope or the shoreline slope: whichever one is steeper. 
        % this ensures that the shoreline elevation remains fixed.
        Srelax = max([Sshore Sthresh]); 
        if Srelax == Sshore; disp('           --> (!) Input threshold slope is too shallow. Using slope to shoreline instead.'); end
                
        % Relax profile downstream of the step
        x = xs*sigma;
        eta_relax = eta(indA) - Srelax.*(x - xA); % relaxation line 
        ind = (x >= xA) & (eta_relax >= eta); % indices to apply the line
        eta(ind) = eta_relax(ind); % apply relaxation line to the area downstream of the step
            % If necessary, relax profile upstream of the step too.
            if Srelax == Sshore;
                disp('           --> (!) Slope to the shoreline might be too steep. Relaxing bed upstream of the step using the threshold slope. ');
                eta_relax2 = eta(end) - Sthresh.*(x - x(end)); % relaxation line. note that it intersects the shoreline now, not hte avulsion node. 
                ind2 = (eta_relax2 <= eta); % indices to apply the line
                eta(ind2) = eta_relax2(ind2); % apply relax line to the area upstream of the step
            end
            
        % plot if desired
        if plotbin
            
            % recalculate slope, for plotting
            S_old = S;
            N = length(sigma);
            S = nan(size(eta));
            S(1)     = -1./xs .* (eta(2) - eta(1))./(sigma(2) - sigma(1)); % forward difference
            S(2:N-1) = -1./xs .* (eta(3:N) - eta(2:N-1))./(sigma(3:N) - sigma(2:N-1)); % centered difference
            S(N)     = -1./xs .* (eta(N) - eta(N-1))./(sigma(N) - sigma(N-1));

            % compare profiles
            figure();
            subplot(2,1,1); title('bed relaxation'); hold on
            plot(xs*sigma,eta_old,'ko-','markerfacecolor','w','markersize',6); hold on
            %plot(xs*sigma(ind),eta_relax(ind),'rx-','markersize',5)
            plot(xs*sigma,eta,'ko-','markerfacecolor','k','markersize',3.5); 
                xlabel('x^*');
                ylabel('\eta^*');
                legend('raw profile','relaxed profile')
            subplot(2,1,2);
            h(1) = plot(xs*sigma,S_old,'ko-','markerfacecolor','w'); hold on
            h(2) = plot(xs*sigma,S,'ko-','markerfacecolor','k');
            h(3) = plot(xlim,[Sstep Sstep],'k:');
            h(4) = plot(xlim,[Sthresh Sthresh],'r-','linewidth',1);
            h(5) = plot(xlim,[Sshore Sshore],'b-','linewidth',1);
                xlabel('x^*')
                ylabel('d\eta / dt')
                legend(h,'raw profile','relaxed profile','raw step','input threshold','to shoreline')
                
        end % end if plotbin
            
       
end

%% Minor functions
 
    % engelund & hansen (1967) total load formulation
    function [qs,qs_feed] = EngelundHansen_fun(tau,hn,Cf,parameters)
    % note that total load = bedload + suspended load, i.e. the transport regimes that make up a significant portion of the bed. does not include wash load. 

        % unpack
        taun_bf = parameters.taun_bf;
        tauc = parameters.SedimentTransport.tauc;
    
        % calculate
        qs = (0.05./Cf) .* tau.^(5/2);
            qs( tau <= tauc ) = 0;      % no sediment is transported if shields number is below the threshold of motion

        % ghost node for sediment feed
        taun = taun_bf * hn; % shields number @ upstream end
        qs_feed = (0.05./Cf) .* taun.^(5/2); % shields number @ upstream end
            %qs_feed( taun <= tauc) = 0; 


    end
    
    % Multi-phase flow regime
    function [hn_new, Cf_new, Tprev_flow_new, cond_newhn] = MultiPhaseFlow_fun(hn,Cf,Tprev_flow,parameters)
    
        % unpack
        Te =  parameters.FlowRegime.MultiPhase.Te; % event timescale
        Fi =  parameters.FlowRegime.MultiPhase.Fi; % fraction that each flow is run
        psii= parameters.FlowRegime.MultiPhase.psii; % nat log of the flow depths in bin i
        Cfi=  parameters.FlowRegime.MultiPhase.Cfi; % friction factor for each flow depth
        
        % the new flow is the same as the previous ...
        hn_new = hn; 
        Cf_new = Cf; 
        Tprev_flow_new = Tprev_flow; % time elapsed since previous flow ended
        cond_newhn = 0; % this is 1 if we have changed discharge, and zero if discharge is the same as prev time
        
        % ... unless we have reached the end of a flow, when we switch the
        % flow randomly, weighting by the fractions
        if  (Tprev_flow >= Te);
            cond_newhn = 1;
            Tprev_flow_new = Tprev_flow - Te; 
            ind = randp(Fi,1); % random bin index drawn from weights according to Fi
            hn_new = exp(psii(ind)); 
            Cf_new = Cfi(ind); 
            
        end % end if
        
    
    end
    
    % Constant width calculator    
    function [w,wx] = ConstantWidth_fun(x,xm,parameters)
    % calculate width as a function of downstream distance. Can be constant
    % width, variable width, etc. Make sure you comment out all the "scenarios"
    % except for the one you want. this is faster than using switch case :)

        % scenario 1) constant-width throughout
        w  = ones(size(x)); % [m] width is constant 
        wx = zeros(size(x));  % [-] dw/dx is zero
        
    end
    
    % Flaring width calculator
    function [w,wx] = FlareWidth_fun(x,xm,theta,parameters)
    
        % unpack
        phi = (parameters.Width.Flare.WidthToDepthRatio * parameters.Frn_bf.^2 * parameters.Cf_bf)^-1; % (Lb / w_n_bf) factor to relate the horizontal scales
                
        % initialize
        w  = nan(size(x));
        wx = nan(size(x)); 
                
        % for x <= xm, river width is constant
        ind = ( x <= xm );
        w(ind) = 1;
            wx(ind) = 0;
        
        % for x > xm, river width expands linearly above a depth H0 = 1,
        % and then width is weighted.
            % widths above and below submerged channel
            w0 = 1; 
            wtheta = 1 + 2*tand(theta).*(x - xm) * phi; 
            % degree of submergence
            H0 = 1; % thickness of submerged channel
            H = H0 + 1*(x-xm); % approximation of depth from sea-level to submerged channel bed
            dH = H-H0; % approximation of depth from sea-level to submerged floodplain
            % weighted width formulation
            wtemp = (H0./H)*w0 + (dH./H).*wtheta;
            wxtemp = 2*tand(theta)*phi + H0*w0./H.^2 - (H0./H).^2*2*tand(theta)*phi; 
            % plug it in only where it applies
            ind = (x > xm);
            w(ind) = wtemp(ind);
            wx(ind) = wxtemp(ind);
                        
    end
    
    % relationship between mean & std for a log-normal distribution
    function [mu] = psi_mean_calculator_fun(sigma,F_bf,plotbin)

        % assuming a lognormal distribution of flow depths and set exceedence
        % probability of bankfull
        hnstar_temp = logspace(-3,0,1000);
        F_bf_fun = @(mu,sigma) 1 - trapz(hnstar_temp,pdf(makedist('lognormal','mu',mu,'sigma',sigma),hnstar_temp)); % calculates exceedence prob of bankfull for a given lognormal distribution of hnstar
       
        % solve for relatinoship between geometric mean & standard deviation of flows
            % initialize
            mu = nan(size(sigma));        
            % solver setings
            options = optimset('display','off','tolx',1e-6); 
            mu_guess = -0.5; % first guess
            % solve for every value of gstd
            for i = fliplr(1:length(sigma))
               %i % display
               resid = @(mu) sqrt((F_bf - F_bf_fun(mu,sigma(i))).^2); 
               [mu_optimum,~,exitflag(i,1)] = fsolve(resid,mu_guess,options); 
                   if exitflag(i) == 1
                   mu(i,1) = mu_optimum; % write
                   mu_guess = mu_optimum; % redo guess
                   end
            end
            % use linear interpolation to get nan values
            ind = ~isnan(mu);
            p = polyfit(sigma(ind),mu(ind),1);
            mu_fit = polyval(p,sigma);
  
            % plot if you want
            if plotbin
                fig = figure();
                for i = 1:length(sigma)
                subplot(1,2,1); hold on; pbaspect([1 1 1]);
                    plot(mu(i),sigma(i),'.','markersize',75); 
                    plot(mu_fit(i),sigma(i),'kx','markersize',10);
                    plot(mu_fit,sigma,'k-');
                        xlabel('\psi_m_e_a_n = mean(log( h_n / h_n_,_b_f ))')
                        ylabel('\psi_s_t_d = std(log( h_n / h_n_,_b_f ))')
                        legend('solution found','solution interpolated');
                        title('relationship between geometric mean & std')
                subplot(1,2,2); hold on; pbaspect([1 1 1]);
                    if ~isnan(mu(i))
                    pd = makedist('lognormal','mu',mu(i),'sigma',sigma(i));
                    pdf_temp = pdf(pd,hnstar_temp);
                    plot(hnstar_temp,pdf_temp,'-'); 
                    else
                    pd = makedist('lognormal','mu',mu_fit(i),'sigma',sigma(i));
                    pdf_temp = pdf(pd,hnstar_temp);
                    plot(hnstar_temp,pdf_temp,'--k');                     
                    end
                        %set(gca,'xscale','log')
                        xlabel('h_n / h_n_,_b_f')
                        ylabel('probability density')
                        ylim([0 5]);
                        title('distributions')
                end
            end

            % replace nans with the fit
            ind = isnan(mu);
            mu(ind) = mu_fit(ind); 


    end
   
    % write workspace to structure
    function WStruct=ws2struct()
    
        WSVARS = evalin('caller', 'who');
        for wscon=1:size(WSVARS,1)
            thisvar=evalin('caller', WSVARS{wscon});
            THEWORKSPACE.(WSVARS{wscon})=thisvar;
        end

        WStruct=THEWORKSPACE;
        
    end

    % round to specified significant digits
    function y=roundsd(x,n,method)

if nargin < 2
	error('Not enough input arguments.')
end

if nargin > 3
	error('Too many input arguments.')
end

if ~isnumeric(x)
		error('X argument must be numeric.')
end

if ~isnumeric(n) || ~isscalar(n) || n < 0 || mod(n,1) ~= 0
	error('N argument must be a scalar positive integer.')
end

opt = {'round','floor','ceil','fix'};

if nargin < 3
	method = opt{1};
else
	if ~ischar(method) || ~ismember(opt,method)
		error('METHOD argument is invalid.')
	end
end

e = floor(log10(abs(x)) - n + 1);
og = 10.^abs(e);
y = feval(method,x./og).*og;
k = find(e<0);
if ~isempty(k)
	y(k) = feval(method,x(k).*og(k))./og(k);
end	
y(x==0) = 0;
    end
    
    % non-exceedence to fractions
    function [ Fi, psii ] = Nonexceed2Fraction_fun( Fnei, psibi )

    % indices
    N = length(psibi); 
    k = (1:(N-1))';
    
    % psi values are the mean of the bounds of the psi bins
    psii(k,1) = 0.5*( psibi(k) + psibi(k+1) ); 
    
    % fractions are the differences in the fractions finer
    Fi(k,1) = abs(Fnei(k) - Fnei(k+1));

    end

    % fractions to non-exceedence
    function [ Fnei ] = Fraction2Nonexceed_fun( Fi, psibi )

        % initialize
        Fnei = zeros(size(psibi)); 

        % indices
        N = length(psibi); 

        % fractions finer
        for k = 1:(N-1)
        Fnei(1) = 1; 
        Fnei(k+1) = Fnei(k) - Fi(k);
        end

    end
    
    % slope calculator
    function [S,Sstep,Sshore] = SlopeCalc_fun(xs,sigma,eta)

            % whole profile
            N = length(sigma);
            S = nan(size(eta));
            S(1)     = -1./xs .* (eta(2) - eta(1))./(sigma(2) - sigma(1)); % forward difference
            S(2:N-1) = -1./xs .* (eta(3:N) - eta(2:N-1))./(sigma(3:N) - sigma(2:N-1)); % centered difference
            S(N)     = -1./xs .* (eta(N) - eta(N-1))./(sigma(N) - sigma(N-1)); % backward difference

            % find maximum slope (should be just d/s of avulsion site)
            [~,indA] = nanmax(S);
            xA = xs*sigma(indA);
            Sstep = S(indA); 

            % find slope from avulsion to shoreline
            Sshore = -(eta(end) - eta(indA))/(xs*sigma(end) - xs*sigma(indA)); % [-] slope to shoreline

    end
